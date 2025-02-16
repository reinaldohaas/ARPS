!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE READNAMELIST                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE readnamelist(progopt,max_dom,input_from_file,                &
     hinfmt,adasbasfn,adashisfn,bdybasfn,hisfile,nhisfile,finexist,     &
     nx,ny,nz,nzsoil,nstyps,nprocx_in,nprocy_in,ncompressx,ncompressy,  &
     use_arps_grid,nx_wrf,ny_wrf,nz_wrf,zlevels_wrf,ptop,               &
     i_parent_start,j_parent_start,parent_id,                           &
     mapproj_wrf,sclfct_wrf,lattru_wrf,lontru_wrf,                      &
     ctrlat_wrf,ctrlon_wrf,dx_wrf,dy_wrf,dt,                            &
     base_pres,base_temp,base_lapse,iso_temp,                           &
     sfcinitopt,wrftrnopt,sfcdtfn,geogdir,start_date,silwt,wvln,        &
     create_bdy,mgrdbas,tintv_bdywrf,tintv_bdyin,spec_bdy_width,        &
     diff_opt,km_opt,khdif,kvdif,mp_physics,ra_lw_physics,              &
     ra_sw_physics,sf_sfclay_physics,sf_surface_physics,bl_pbl_physics, &
     cu_physics, nprocx_wrf,nprocy_wrf,frames_per_outfile,              &
     restart_interval,radt,cudt,ifsnow,w_damping,parent_time_step_ratio,&
     iorder,korder,io_form,qx_zero_out,create_namelist,wrfnamelist,     &
     io_form_history, io_form_restart, history_interval,                &
     indir,outdir,staticdir,output_interval_2d,output_dir_2d,           &
     moist_adv_opt, dampcoef,wrfversion, istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  To read and propagate namelist input.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  09/15/2005
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: progopt     ! 0 = ARPS2WRF
                                      ! 1 = WRFSTATIC
  INTEGER, INTENT(OUT) :: istatus

!----------------------------------------------------------------------
!
! ARPS grid variables
!
!---------------------------------------------------------------------

  INTEGER, PARAMETER :: nmax_domains = 100

  INTEGER, DIMENSION(nmax_domains), INTENT(OUT) :: nx,ny,nz     ! Grid dimensions for ARPS.
  INTEGER, DIMENSION(nmax_domains), INTENT(OUT) :: nzsoil       ! Soil levels
  INTEGER, DIMENSION(nmax_domains), INTENT(OUT) :: nstyps       ! Maximum number of soil types.

  INTEGER, PARAMETER :: nhisfile_max = 100
  INTEGER, PARAMETER :: max_vertical_levels = 100

  CHARACTER(LEN=*), INTENT(OUT) :: hisfile(nhisfile_max)
  CHARACTER(LEN=*), INTENT(OUT) :: bdybasfn(nhisfile_max)
  INTEGER,          INTENT(OUT) :: nhisfile

  CHARACTER(LEN=256) :: grdbasfn
  INTEGER            :: lengbf

  ! hisfile(1) and bdybasfn(1) are for WRF input file
  !      They can be any ARPS history dumps including ADAS output,
  !      ARPS output or EXT2ARPS output
  !
  ! hisfile(2:nhisfil), bdybasfn(2:nhisfile) are for
  !      WRF lateral bounday files. They can be either ARPS history
  !      dumps or EXT2ARPS outputs

  CHARACTER(LEN=4), PARAMETER :: finfmt(11) =                           &
                          (/'.bin','.asc','.hdf','.pak','.svi','.bn2',  &
                            '.net','.npk','.gad','.grb','.v5d'  /)

  INTEGER, INTENT(OUT) :: finexist(nhisfile_max)
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
! Variables for mpi jobs
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER   :: fzone_arps = 3, fzone_wrf = 1
  INTEGER, INTENT(OUT) :: ncompressx, ncompressy ! compression in x and y direction:
                                    ! ncompressx=nprocx_in/nproc_x
                                    ! ncompressy=nprocy_in/nproc_y

  INTEGER, INTENT(OUT) :: nprocx_in, nprocy_in

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
  INTEGER, DIMENSION(nmax_domains), INTENT(OUT) :: nx_wrf      ! = nx-2 if the same grid as ARPS
  INTEGER, DIMENSION(nmax_domains), INTENT(OUT) :: ny_wrf      ! = ny-2 if the same grid as ARPS
  INTEGER, INTENT(OUT) :: nz_wrf      ! = nz-2 if the same grid as ARPS
                                      ! All are staggered values

  REAL,    INTENT(OUT) :: lattru_wrf(2)  ! array of true latitude of WRF map projection
  REAL,    INTENT(OUT) :: lontru_wrf     ! true longitude of WRF map projection
                                      ! = trulon_wrf

! Namelist variable declaration

  CHARACTER(LEN=7),   DIMENSION(nmax_domains), INTENT(OUT) :: sfcinitopt
  CHARACTER(LEN=19),  DIMENSION(nmax_domains), INTENT(OUT) :: start_date
  CHARACTER(LEN=256),                          INTENT(OUT) :: geogdir
  INTEGER,   DIMENSION(nmax_domains), INTENT(OUT) :: wrftrnopt
  REAL                            :: silavwt_parm_wrf,toptwvl_parm_wrf
  REAL,               INTENT(OUT) :: silwt,wvln

  INTEGER, INTENT(OUT) :: create_bdy      ! Create WRF boundary file
  INTEGER, INTENT(OUT) :: create_namelist ! Dump WRF namelist.input

  INTEGER, DIMENSION(nmax_domains), INTENT(OUT) :: use_arps_grid
                                          ! Use ARPS horizontal grid as WRF grid

  INTEGER, INTENT(OUT) :: tintv_bdywrf    ! Desired WRF boundary file interval (in seconds)
  CHARACTER(LEN=256)   :: bdyfheader      ! ARPS boundary input file header
  INTEGER              :: tbgn_bdyin      ! ARPS boundary begin time (in seconds)
  INTEGER, INTENT(OUT) :: tintv_bdyin     ! ARPS boundary file interval (in seconds)
  INTEGER              :: tend_bdyin      ! Last ARPS boundary file time
  INTEGER, INTENT(OUT) :: mgrdbas         ! Options for grid base file
                                       ! = 0 share same grid base as initial state file
                                       ! = 1 All ARPS boundary files share one grd base
                                       !     file but it is difference from the inital
                                       !     base file as specified using grdbasfn
                                       ! = 2 Each file has its own grid base file

  CHARACTER(LEN=256), INTENT(OUT) :: wrfnamelist  ! file name for WRF namelist.input

  INTEGER, INTENT(OUT) :: mapproj_wrf     ! Type of map projection in WRF model grid
                             ! modproj = 1  Polar Stereographic
                             !              projection.
                             !         = 2  Mercator projection.
                             !         = 3  Lambert projection.

  REAL,    INTENT(OUT) :: sclfct_wrf      ! Map scale factor.
                             ! Distance on map, between two latitudes
                             ! trulat1 and trulat2,
                             ! is = (Distance on earth)*sclfct.
                             ! For ARPS model runs,
                             ! generally this is 1.0

  REAL              :: trulat1_wrf, trulat2_wrf, trulon_wrf
                             ! 1st, 2nd real true latitude and true longitude
                             ! of WRF map projection
  REAL, INTENT(OUT) :: ctrlat_wrf      ! Center latitude of WRF model domain (deg. N)
  REAL, INTENT(OUT) :: ctrlon_wrf      ! Center longitude of WRF model domain (deg. E)

  REAL, DIMENSION(nmax_domains), INTENT(OUT) :: dx_wrf      ! WRF Grid spacing in x-direction
  REAL, DIMENSION(nmax_domains), INTENT(OUT) :: dy_wrf      ! WRF Grid spacing in y-direction

  INTEGER           :: vertgrd_opt     ! WRF sigma level scheme

  REAL, DIMENSION(nmax_domains), INTENT(OUT) :: ptop        ! WRF atmosphere top pressure in Pascal

  REAL              :: pbot
  REAL, INTENT(OUT) :: zlevels_wrf(max_vertical_levels)
                             ! WRF mass levels from 1.0 at surfact to
                             ! 0.0 at atmosphere top

  INTEGER, DIMENSION(nmax_domains), INTENT(OUT) :: i_parent_start
  INTEGER, DIMENSION(nmax_domains), INTENT(OUT) :: j_parent_start

  INTEGER, INTENT(OUT) :: iorder          ! order of polynomial for horizontal
                             ! interpolation (1 or 2)
  INTEGER, INTENT(OUT) :: korder          ! vertical interpolation order (1 or 2)

  INTEGER :: dyn_opt         ! WRF dynamics option
                             ! only works for = 2 Eulerian mass coordinate
  INTEGER, INTENT(OUT) :: diff_opt        ! WRF diffusion option
  INTEGER, INTENT(OUT) :: km_opt          ! WRF eddy coefficient option
  REAL,    DIMENSION(nmax_domains), INTENT(OUT) :: khdif   ! Horizontal diffusion constant (m^2/s)
  REAL,    DIMENSION(nmax_domains), INTENT(OUT) :: kvdif   ! Vertical diffusion constant (m^2/s)
  INTEGER, DIMENSION(nmax_domains), INTENT(OUT) :: mp_physics      ! WRF microphysics options
                             != 2 Lin et al. scheme
                             !   (QVAPOR,QRAIN,QSNOW,QCLOUD,QICE,QGRAUP)
                             != 3 NCEP 3-class simple ice scheme
                             !   (QVAPOR,QCLOUD,QICE,QRAIN,QSNOW)
                             != 4 NCEP 5-class scheme
                             !   (QVAPOR,QCLOUD,QICE,QRAIN,QSNOW)
                             != 5 Ferrier (new Eta) microphysics
                             !   (QVAPOR,QCLOUD)
  INTEGER, DIMENSION(nmax_domains), INTENT(OUT) :: ra_lw_physics   ! Longwave radiaiton option
  INTEGER, DIMENSION(nmax_domains), INTENT(OUT) :: ra_sw_physics   ! Shortwave radiation option
  INTEGER, DIMENSION(nmax_domains), INTENT(OUT) :: sf_sfclay_physics  ! WRF surface-layer option
  INTEGER, DIMENSION(nmax_domains), INTENT(OUT) :: sf_surface_physics ! WRF land-surface option
                                ! = 0 no land-surface
                                !     (DO NOT use)
                                ! = 1 Thermal diffusion scheme
                                !     (nzsoil_wrf = 5)
                                ! = 2 OSU land-surface model
                                !     (nzsoil_wrf = 4)
                                ! = 3 Do not use
                                !     (nzsoil_wrf = 6)
  INTEGER, DIMENSION(nmax_domains), INTENT(OUT) :: bl_pbl_physics   ! boundary-layer option
  INTEGER, DIMENSION(nmax_domains), INTENT(OUT) :: cu_physics       ! cumulus option
  REAL,    INTENT(OUT) :: dt               ! time-step for advection
  INTEGER, INTENT(OUT) :: spec_bdy_width   ! number of rows for specified boundary values nudging
  REAL,    INTENT(OUT) :: base_pres, base_temp, base_lapse, iso_temp
  INTEGER, INTENT(OUT) :: nprocx_wrf      ! Number of X direction processors for WRF run
  INTEGER, INTENT(OUT) :: nprocy_wrf      ! Number of Y direction processors for WRF run
  INTEGER, DIMENSION(nmax_domains), INTENT(OUT) :: frames_per_outfile
  REAL,    DIMENSION(nmax_domains), INTENT(OUT) :: radt,cudt
  INTEGER, DIMENSION(nmax_domains), INTENT(OUT) :: parent_time_step_ratio
  INTEGER, INTENT(OUT) :: restart_interval
  INTEGER, INTENT(OUT) :: ifsnow, w_damping
  INTEGER, INTENT(OUT) :: io_form
  INTEGER, INTENT(OUT) :: qx_zero_out

  INTEGER,                                     INTENT(OUT) :: max_dom
  INTEGER,            DIMENSION(nmax_domains), INTENT(OUT) :: parent_id
  INTEGER,            DIMENSION(nmax_domains), INTENT(OUT) :: hinfmt
  CHARACTER(LEN=256), DIMENSION(nmax_domains), INTENT(OUT) :: adashisfn, adasbasfn
  CHARACTER(LEN=256), DIMENSION(nmax_domains), INTENT(OUT) :: sfcdtfn
  LOGICAL,            DIMENSION(nmax_domains), INTENT(OUT) :: input_from_file

  INTEGER :: readsplit_in

  NAMELIST /message_passing/ nproc_x, nproc_y, readsplit_in, nprocx_in, nprocy_in

  NAMELIST /domains/ max_dom,parent_id,input_from_file

  NAMELIST /history_data/ hinfmt, adasbasfn, adashisfn

  NAMELIST /sfcdt/ sfcinitopt,wrftrnopt,sfcfmt,sfcdtfn,                 &
                   geogdir,silavwt_parm_wrf,toptwvl_parm_wrf,start_date

  NAMELIST /bdyspc/ create_bdy,tintv_bdywrf,bdyfheader,tbgn_bdyin,      &
                    tintv_bdyin,tend_bdyin,mgrdbas

  REAL :: max_dz

  NAMELIST /wrf_grid/ use_arps_grid,nx_wrf,ny_wrf,                      &
                      mapproj_wrf, sclfct_wrf,                          &
                      trulat1_wrf, trulat2_wrf,trulon_wrf,              &
                      ctrlat_wrf, ctrlon_wrf,                           &
                      dx_wrf, dy_wrf, i_parent_start, j_parent_start,   &
                      ptop,vertgrd_opt,nz_wrf,pbot,zlevels_wrf,max_dz

  NAMELIST /interp_options/ iorder, korder

  INTEGER,            INTENT(OUT) :: io_form_history, io_form_restart
  CHARACTER(LEN=256), INTENT(OUT) :: staticdir,indir, outdir
  INTEGER,            INTENT(OUT) :: output_interval_2d
  CHARACTER(LEN=40),  INTENT(OUT) :: output_dir_2d
  !LOGICAL,            INTENT(OUT) :: pd_moist
  INTEGER,            INTENT(OUT) :: moist_adv_opt
  INTEGER, DIMENSION(nmax_domains), INTENT(OUT) :: history_interval
  REAL,               INTENT(OUT) :: wrfversion
  REAL,               INTENT(OUT) :: dampcoef

  NAMELIST /wrf_opts/ diff_opt, km_opt, khdif, kvdif,                   &
                      mp_physics, ra_lw_physics, ra_sw_physics,         &
                      sf_sfclay_physics, sf_surface_physics,            &
                      bl_pbl_physics, cu_physics,                       &
                      base_pres,base_temp, base_lapse,iso_temp,dt,      &
                      spec_bdy_width, nprocx_wrf, nprocy_wrf,           &
                      frames_per_outfile,restart_interval,radt,cudt,    &
                      ifsnow, w_damping, parent_time_step_ratio,        &
                      io_form_history,io_form_restart,history_interval, &
                      indir,outdir,staticdir,output_interval_2d,        &
                      output_dir_2d,moist_adv_opt,dampcoef

  NAMELIST /output/ wrfversion,dirname,readyfl,io_form,qx_zero_out,     &
                    create_namelist,wrfnamelist

!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,n,ifile
  INTEGER :: lenstr, ireturn

  INTEGER :: domid

  LOGICAL :: fexist
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
! Read mpi related block
!
!-----------------------------------------------------------------------

  nproc_x = 1
  nproc_y = 1
  readsplit_in = 1
  nprocx_in = 1
  nprocy_in = 1
  IF(myproc == 0) THEN
    READ(5,message_passing)
    WRITE(6,'(a)')'Namelist message_passing was successfully read.'
  END IF
  readsplit(:) = readsplit_in

  CALL mpupdatei(nproc_x,1)
  CALL mpupdatei(nproc_y,1)
  CALL mpupdatei(readsplit,FINDX_NUM)

  IF(readsplit_in > 0) THEN
    nprocx_in  = nproc_x
    nprocy_in  = nproc_y
  END IF
  CALL mpupdatei(nprocx_in,1)
  CALL mpupdatei(nprocy_in,1)

  ncompressx = nprocx_in/nproc_x
  ncompressy = nprocy_in/nproc_y

  IF (mp_opt > 0 .AND. (                              &
       (MOD(nprocx_in,nproc_x) /= 0 .OR. MOD(nprocy_in,nproc_y) /= 0)   &
       .OR. (ncompressx < 1 .OR. ncompressy < 1) )  ) THEN
    IF (myproc == 0) WRITE(6,'(3x,a/,2(3x,2(a,I2)/))')                  &
      'nprocx_in (nprocy_in) must be a multiplier of nproc_x(nproc_y).',&
      'nprocx_in = ',nprocx_in, ', nprocy_in = ',nprocy_in,             &
      'nproc_x   = ',nproc_x,   ', nproc_y   = ',nproc_y
    CALL arpsstop('unmatched dimension size.',1);
  END IF

  IF (mp_opt == 0) THEN
    ncompressx = 1
    ncompressy = 1
    nproc_x = 1
    nproc_y = 1
    nprocx_in = 1
    nprocy_in = 1
    myproc = 0
    loc_x = 1
    loc_y = 1
    readsplit(:) = 0
  ELSE
    CALL mpinit_var
  END IF

  max_dom   = 1
  parent_id = 1
  input_from_file = .TRUE.
  IF(myproc == 0) THEN
    READ(5,domains)
    WRITE(6,'(a)')'Namelist domains was successfully read.'
  END IF
  CALL mpupdatei(max_dom,        1)
  CALL mpupdatei(parent_id,      nmax_domains)
  CALL mpupdatel(input_from_file,nmax_domains)
  IF (.NOT. input_from_file(1)) THEN
    WRITE(6,'(1x,a)') 'ERROR: input_from_file for domain 1 must be .TRUE..'
    CALL arpsstop('Wrong input_from_file in domains.',1)
  END IF
!
!-----------------------------------------------------------------------
!
!  Get the names of the input data files.
!
!-----------------------------------------------------------------------
!
  IF(myproc == 0) THEN
    READ(5,history_data)
    WRITE(6,'(a)') 'Namelist history_data was successfully read.'

    WRITE(6,'(2x,a6,a36,a36)') 'Domain','  ADAS file   ','  Grid BASE file  '
    WRITE(6,'(2x,a6,a36,a36)') '======','==============','=================='
    DO domid = 1, max_dom
      IF (input_from_file(domid))  &
      WRITE(6,'(2x,I3,a3,a36,a36)') domid,':  ',TRIM(adashisfn(domid)),TRIM(adasbasfn(domid))
    END DO
    WRITE(6,*)

    hisfile(1)  = adashisfn(1)
    bdybasfn(1) = adasbasfn(1)

  END IF
  CALL mpupdatei(hinfmt,nmax_domains)
  CALL mpupdatec(adashisfn,256*nmax_domains)
  CALL mpupdatec(adasbasfn,256*nmax_domains)

  finexist(:) = 0
  nhisfile    = 1
  finexist(1) = 1

!-----------------------------------------------------------------------
!
! Now get ARPS dimensions for program ARPS2WRF
!
!-----------------------------------------------------------------------

  IF (progopt == 0) THEN

    DO domid = 1,max_dom
      IF (input_from_file(domid)) THEN
        IF(mp_opt > 0 .AND. readsplit(FINDX_H) <= 0) THEN
          CALL gtsplitfn(adasbasfn(domid),1,1,loc_x,loc_y,1,1, &
                         0,0,1,lvldbg,grdbasfn,istatus)
        ELSE
          WRITE(grdbasfn,'(a)') TRIM(adasbasfn(domid))
        END IF
        lengbf = len_trim(grdbasfn)

        IF (myproc == 0) THEN
          CALL get_dims_from_data(hinfmt(domid),grdbasfn(1:lengbf),       &
          nx(domid),ny(domid),nz(domid),nzsoil(domid),nstyps(domid),ireturn)
        END IF

        CALL mpupdatei(ireturn,1)
        IF( ireturn /= 0 ) THEN
          PRINT*,'Problem occured when trying to get dimensions from data.'
          PRINT*,'Program stopped.'
          CALL arpsstop('get_dims_from_data error.',1)
        END IF
      END IF
    END DO

    CALL mpupdatei(nx,    nmax_domains)
    CALL mpupdatei(ny,    nmax_domains)
    CALL mpupdatei(nz,    nmax_domains)
    CALL mpupdatei(nzsoil,nmax_domains)
    CALL mpupdatei(nstyps,nmax_domains)

  ELSE         ! wrfstatic does not use ARPS dimensions so far
    nx = 1
    ny = 1
    nz = 1
    nzsoil = 1
    nstyps = 1
  END IF


!----------------------------------------------------------------------
!
!  Get surface characteristics file options
!
!----------------------------------------------------------------------
!
  sfcinitopt(:)(1:7)= ' '
  sfcfmt            = 1
  sfcdtfn(:)(1:256) = ' '
  wrftrnopt(:)      = 0
  geogdir(1:256)    = ' '
  start_date(:)     = '1998-05-25_00:00:00'
  silavwt_parm_wrf  = 0.0
  toptwvl_parm_wrf  = 2.0

  IF (myproc == 0) THEN
    READ(5,sfcdt)
    WRITE(6,'(a)') 'Namelist sfcdt was successfully read.'

    IF (progopt == 0) THEN

      DO domid = 1,max_dom
        IF (input_from_file(domid)) THEN
          IF (sfcinitopt(domid) == 'ARPS' .AND. mp_opt > 0 .AND. readsplit(FINDX_T) == 0) THEN
            CALL gtsplitfn(sfcdtfn(domid),1,1,loc_x,loc_y,1,1,          &
                         0,0,1,lvldbg,grdbasfn,istatus)
          ELSE
            WRITE(grdbasfn,'(a)') TRIM(sfcdtfn(domid))
          END IF

          INQUIRE(FILE=trim(grdbasfn), EXIST = fexist )
          IF(.NOT. fexist) THEN
            WRITE(6,*) 'The file ',TRIM(grdbasfn),' you specified does not exist.'
            CALL arpsstop('File does not exist.',1)
          END IF
        END IF   ! input_from_file
      END DO
    END IF
  END IF

  lenstr = LEN_TRIM(geogdir)
  IF (geogdir(lenstr:lenstr) /= '/') THEN
    lenstr = lenstr + 1
    geogdir(lenstr:lenstr) = '/'
  END IF

  silwt = silavwt_parm_wrf
  wvln  = toptwvl_parm_wrf

  CALL mpupdatei(sfcfmt,1)
  CALL mpupdatec(sfcinitopt,  7*nmax_domains)
  CALL mpupdatec(sfcdtfn,   256*nmax_domains)

!----------------------------------------------------------------------
!
!  Get boundary file specifications
!
!----------------------------------------------------------------------
!
  create_bdy = 0
  tintv_bdywrf = 10800
  bdyfheader = './eta25may1998'
  tbgn_bdyin = 10800
  tintv_bdyin= 10800
  tend_bdyin = 21600
  mgrdbas    = 0

  n = 0
  IF (myproc == 0) THEN
    READ(5,bdyspc)
    WRITE(6,'(a)') 'Namelist bdyspc was successfully read.'
    IF(create_bdy >= 1 .AND. progopt == 0) THEN
      IF(tintv_bdyin < 1) THEN
        WRITE(6,'(/a,I6,a/)') 'ERROR: Boudnary interval (tintv_bdyin =',  &
                          tintv_bdyin,') is not correct. Terminating ...'
        CALL arpsstop('Wrong namelist input parameters.',1)
      END IF

      n = 1
      IF (create_bdy == 1) tintv_bdywrf = tintv_bdyin

      DO i = tbgn_bdyin,tend_bdyin,tintv_bdyin
        nhisfile = nhisfile + 1
        WRITE(hisfile(nhisfile),'(2a,I6.6)')                            &
                   TRIM(bdyfheader),finfmt(hinfmt(1)),i
        finexist(nhisfile) = 1

        IF (mp_opt > 0 .AND. readsplit(FINDX_H) <= 0 ) THEN
          CALL gtsplitfn(hisfile(nhisfile),1,1,loc_x,loc_y,1,1,         &
                         0,0,1,lvldbg,grdbasfn,istatus)
        ELSE
          WRITE(grdbasfn,'(a)') hisfile(nhisfile)   ! used as temporary string
        END IF
        INQUIRE(FILE=trim(grdbasfn), EXIST = fexist )
        IF(.NOT. fexist) THEN
          WRITE(6,'(/1x,3a,I6,a/)') 'WARNING: The file ',TRIM(grdbasfn),  &
                  ' does not exist. Boundary file at ',i,' was skipped.'
          finexist(nhisfile) = 0
          !nhisfile = nhisfile - 1
          !CYCLE
          !CALL arpsstop('File does not exist.',1)
        END IF

        IF(mgrdbas == 1 .AND. nhisfile == 2) THEN
          WRITE(bdybasfn(nhisfile),'(3a)')                              &
                   TRIM(bdyfheader),finfmt(hinfmt(1)),'grdbas'
          n = 2      ! maximum index when checking the existance of grid and base file
        ELSE IF(mgrdbas > 1) THEN
          WRITE(bdybasfn(nhisfile),'(3a,I2.2)')                         &
                   TRIM(bdyfheader),finfmt(hinfmt(1)),'grdbas.',nhisfile-1
          n = nhisfile
        END IF
      END DO

    END IF   ! create_bdy == 1

    !
    ! Check grid and base file availability
    ! only do check for the root processor because nproc_x may not equal nprocx_in
    !
    DO i = 1,n
      IF (mp_opt > 0 .AND. readsplit(FINDX_H) <= 0 ) THEN
        CALL gtsplitfn(bdybasfn(i),1,1,loc_x,loc_y,1,1,         &
                       0,0,1,lvldbg,grdbasfn,istatus)
      ELSE
        WRITE(grdbasfn,'(a)') TRIM(bdybasfn(i))   ! used as temporary string
      END IF
      INQUIRE(FILE=trim(grdbasfn), EXIST = fexist )
      IF(.NOT. fexist) THEN
        WRITE(6,'(1x,3a)') 'WARNING: The ARPS grid and base file ',TRIM(grdbasfn),' does not exist.'
        IF (i == 1) THEN          ! grid & base file must be exist for the first time level
          CALL arpsstop('File does not exist.',1)
        ELSE IF (finexist(i) == 1) THEN    ! the history data exsit
          finexist(i) = 0
          WRITE(6,'(10x,a,I2,a,/,10x,a,/,10x,a,I2,a)')                       &
                  'The boundary data at time level ',i,' is skipped.',  &
                  'Actually, for better results, you can avoid this by specifying mgrdbas = 0/1,', &
                  'since the time dependent data file at time level ',i,' is available.'
        !ELSE                              ! the history data is also not exist
        END IF

      END IF
    END DO

  END IF ! myproc == 0
  CALL mpupdatei(n,1)
  CALL mpupdatei(nhisfile,1)
  CALL mpupdatec(bdybasfn,nhisfile_max*256)    ! does not include processor app.
  CALL mpupdatec(hisfile, nhisfile_max*256)

  CALL mpupdatei(create_bdy,1)
  CALL mpupdatei(mgrdbas, 1)
  CALL mpupdatei(tintv_bdywrf,1)
  CALL mpupdatei(tintv_bdyin,1)
  CALL mpupdatei(finexist,nhisfile_max)

!-----------------------------------------------------------------------
!
!  Get WRF grid options
!
!-----------------------------------------------------------------------
!
  use_arps_grid = 0
  nx_wrf  = nx-2
  ny_wrf  = ny-2

  i_parent_start = 1
  j_parent_start = 1

  ptop    = 5000
  pbot    = 101300
  vertgrd_opt = 0
  nz_wrf  = 31

  zlevels_wrf = -1.0
  max_dz      = 1000.

  IF (myproc == 0) THEN
    READ(5,wrf_grid)
    WRITE(6,'(a)') 'Namelist wrf_grid was successfully read.'

    DO domid = 1, max_dom

      IF(use_arps_grid(domid) == 1 .AND. progopt == 0 .AND. input_from_file(domid)) THEN

      WRITE(6,'(/,3(1x,a/))')   &
  '***** You have chosen to use ARPS grid in the data as your   ******', &
  '***** WRF grid. Please note that the two fake points in each ******', &
  '***** direction of ARPS grid are not used.                   ******'

      IF (mp_opt > 0) THEN
        IF( readsplit(FINDX_H) > 0 ) THEN
          IF( MOD(nx(domid)-fzone_arps,nproc_x) /= 0 .OR.                      &
              MOD(ny(domid)-fzone_arps,nproc_y) /= 0) THEN
            WRITE(6,'(a/,a/,4(a,i5))')      &
             'The specification of nproc_x or nproc_y is not matched with nx or ny.',&
             'nx-3 and ny-3 must be multiples of nproc_x and nproc_y respectively.', &
             'nx = ', nx(domid), ' ny = ', ny(domid), ' nproc_x = ',nproc_x, ' nproc_y = ',nproc_y
            nx(domid) = 0
            ny(domid) = 0                 ! to be exit later
          ELSE
            nx(domid) = (nx(domid)-fzone_arps)/nproc_x + fzone_arps
            ny(domid) = (ny(domid)-fzone_arps)/nproc_y + fzone_arps
          END IF
        ELSE
          nx(domid) = (nx(domid)-fzone_arps)*ncompressx + fzone_arps
          ny(domid) = (ny(domid)-fzone_arps)*ncompressy + fzone_arps
        END IF
      END IF

      WRITE(6,'(5(a,i5))') '  nx = ',nx(domid),', ny = ',ny(domid),', nz = ',nz(domid),      &
                           ', nzsoil = ',nzsoil(domid),', nstyps = ',nstyps(domid)

      nx_wrf(domid)  = nx(domid) - 2
      ny_wrf(domid)  = ny(domid) - 2

    ELSE IF (mp_opt > 0) THEN

      WRITE(6,'(1x,2a,/,a,/,2a,/,a,/,a)') 'WARNING: At present, ',      &
                'arps2wrf_mpi only works when WRF horizontal grid and', &
                ' ARPS horizontal grid are the same.', 'Please set ',   &
                'use_arps_grid = 1 in arps2wrf.input.',                 &
                'Or use no-mpi version of arps2wrf.',                   &
                ' Program stopping ...'
      CALL arpsstop('MPI mode does not work.',1)

    END IF

    END DO  ! max_dom

  END IF  ! myproc == 0

  CALL mpupdatei(nx,    nmax_domains)
  CALL mpupdatei(ny,    nmax_domains)
  CALL mpupdatei(nz,    nmax_domains)
  CALL mpupdatei(nzsoil,nmax_domains)
  CALL mpupdatei(nstyps,nmax_domains)

  DO domid = 1,max_dom
    IF (input_from_file(domid)) THEN
      IF( nx(domid) <= 0 .OR. ny(domid) <= 0 ) CALL arpsstop('Wrong size of dimensions.',1);
    ELSE
      nx(domid) = -1
      ny(domid) = -1
    END IF
  END DO

  nstyp = nstyps(1)

  CALL mpupdatei(nx_wrf,        nmax_domains)
  CALL mpupdatei(ny_wrf,        nmax_domains)
  CALL mpupdatei(i_parent_start,nmax_domains)
  CALL mpupdatei(j_parent_start,nmax_domains)
  CALL mpupdatei(nz_wrf,1)
  CALL mpupdatei(use_arps_grid, nmax_domains)
  CALL mpupdater(ptop,          nmax_domains)
  CALL mpupdater(zlevels_wrf,max_vertical_levels)

  !
  ! Map projection parameters do not need to be updated because
  ! use_arps_grid must be 1 for mpi job
  !
  lattru_wrf(1) = trulat1_wrf
  lattru_wrf(2) = trulat2_wrf
  lontru_wrf    = trulon_wrf

!-----------------------------------------------------------------------
!
!  Get WRF options (NetCDF file needs them for global attributes)
!
!-----------------------------------------------------------------------

  dyn_opt            = 2

!  IF(dyn_opt /= 2) THEN
!    WRITE(6,*) 'NOTE: ARPS2WRF only works for WRF MASS core at present'
!    WRITE(6,*) '      dyn_opt has been reseted to 2.'
!  END IF

  diff_opt           = 0
  km_opt             = 1
  khdif              = 0.0
  kvdif              = 0.0
  mp_physics         = 1    ! must be specified, used for moist variable determination
  ra_lw_physics      = 1
  ra_sw_physics      = 1
  sf_sfclay_physics  = 1
  sf_surface_physics = 1    ! must be specified, used for soil layer determination
  bl_pbl_physics     = 1
  cu_physics         = 1
  base_pres          = 100000.
  base_temp          = 290.
  base_lapse         = 50.
  iso_temp           = 0.
  dt                 = 40
  spec_bdy_width     = 5
  nprocx_wrf         = -1           ! only root processor write namelist parameters
  nprocy_wrf         = -1
  frames_per_outfile = 1
  restart_interval   = 60
  radt               = 30.
  cudt               = 0.
  ifsnow             = 0
  w_damping          = 0
  parent_time_step_ratio = 1
  io_form_history    = 2
  io_form_restart    = 2
  history_interval(:)= 60
  indir              = './'
  outdir             = './'
  staticdir          = './'
  output_interval_2d = 0
  output_dir_2d      = './'

  !pd_moist           = .FALSE.
  moist_adv_opt      = 0
  dampcoef           = 0.0

  IF (myproc == 0) THEN
    READ(5,wrf_opts)
    WRITE(6,'(a)') 'Namelist wrf_opts was successfully read.'

    IF(progopt == 0 ) THEN
      DO domid = 1,max_dom
        IF (sf_surface_physics(domid) > 3 .OR. sf_surface_physics(domid) < 1 ) THEN
          WRITE(6,*) 'Not a valid sf_surface_physics option - ',sf_surface_physics
          WRITE(6,*) 'It must be either 1, 2 or 3.'
          CALL arpsstop('Bad WRF namelist parameters inside arps2wrf.input.',1)
        END IF

        IF ( (bl_pbl_physics(domid) == 2 .AND. sf_sfclay_physics(domid) /= 2) .OR. &
             (bl_pbl_physics(domid) /= 2 .AND. sf_sfclay_physics(domid) == 2) ) THEN
          WRITE(6,'(1x,2(a,/))')                                        &
                  'MYJ PBL scheme requires a matched sfclay scheme.',   &
                  'Please check sf_sfclay_physics in the namelist file.'
          CALL arpsstop('Not matched sfclay and pbl options.',1)
        END IF
      END DO
    END IF

    lenstr = LEN_TRIM(staticdir)
    IF(lenstr > 0) THEN
       IF(staticdir(lenstr:lenstr) /= '/') THEN
         staticdir(lenstr+1:lenstr+1) = '/'
         lenstr = lenstr + 1
       END IF
    ELSE
       staticdir = './'
    END IF

    lenstr = LEN_TRIM(indir)
    IF(lenstr > 0) THEN
       IF(indir(lenstr:lenstr) /= '/') THEN
         indir(lenstr+1:lenstr+1) = '/'
         lenstr = lenstr + 1
       END IF
    ELSE
       indir = './'
    END IF

    lenstr = LEN_TRIM(outdir)
    IF(lenstr > 0) THEN
       IF(outdir(lenstr:lenstr) /= '/') THEN
         outdir(lenstr+1:lenstr+1) = '/'
         lenstr = lenstr + 1
       END IF
    ELSE
       outdir = './'
    END IF

  END IF
  CALL mpupdatei(mp_physics,nmax_domains)
  CALL mpupdatei(sf_surface_physics,nmax_domains)
  CALL mpupdatei(spec_bdy_width,1)
  CALL mpupdater(base_temp,1)

!-----------------------------------------------------------------------
!
! Compute vertical levels based on user's choice
!
!-----------------------------------------------------------------------

  IF (progopt == 0) THEN   ! for ARPS2WRF only
    IF (myproc == 0) CALL compute_eta(vertgrd_opt,nz_wrf,max_dz,        &
                     pbot,ptop,max_vertical_levels,zlevels_wrf,istatus)

    CALL mpupdatei(nz_wrf,1)
    CALL mpupdater(zlevels_wrf,max_vertical_levels)
  END IF
!
!-----------------------------------------------------------------------
!
!  Get interpolation options
!
!-----------------------------------------------------------------------
!
   iorder = 3
   korder = 2

   IF (myproc == 0) THEN
     READ(5,interp_options)
     WRITE(6,'(a)') 'Namelist interp_options was successfully read.'
   END IF
   CALL mpupdatei(iorder,1)
   CALL mpupdatei(korder,1)
!
!-----------------------------------------------------------------------
!
!  Get output options
!
!-----------------------------------------------------------------------
!
   wrfversion      = 2.2
   dirname         = './'
   readyfl         = 0
   create_namelist = 0

   IF (myproc == 0) THEN
     READ(5,output)
     WRITE(6,'(a)') 'Namelist output was successfully read.'
     lenstr = LEN_TRIM(dirname)
     IF(lenstr > 0) THEN
       IF(dirname(lenstr:lenstr) /= '/') THEN
         dirname(lenstr+1:lenstr+1) = '/'
         lenstr = lenstr + 1
       END IF
     ELSE
       dirname = './'
     END IF
   END IF
   CALL mpupdater(wrfversion,1)
   CALL mpupdatei(io_form,1)
   CALL mpupdatei(qx_zero_out,1)

!-----------------------------------------------------------------------
!
! Successfully finished with namelist input
!
!-----------------------------------------------------------------------

  istatus = 0

  RETURN
END SUBROUTINE readnamelist

!#######################################################################
!#######################################################################
!####                                                               ####
!####   SUBROUTINE compute_eta                                      ####
!####                                                               ####
!#######################################################################
!#######################################################################

SUBROUTINE compute_eta(vertgrdopt,nzwrf,maxdz,pbot,ptop,num_max_levels, &
                       znw,istatus)

!-----------------------------------------------------------------------
!
! Purpose
!   Compute WRF eta levels
!
!-----------------------------------------------------------------------

  USE wrf_metadata

  IMPLICIT NONE

  INTEGER, INTENT(IN)    :: vertgrdopt
  INTEGER, INTENT(INOUT) :: nzwrf
  REAL,    INTENT(IN)    :: maxdz
  REAL,    INTENT(IN)    :: pbot, ptop
  INTEGER, INTENT(IN)    :: num_max_levels
  REAL,    INTENT(INOUT) :: znw(num_max_levels)
  INTEGER, INTENT(OUT)   :: istatus

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: nlevels
  INTEGER :: k
  REAL    :: sigma1, sigma2, plevel1, plevel2

  LOGICAL :: compute_const_dz = .FALSE.

  INTEGER , PARAMETER              :: prac_levels = 17
  REAL,     DIMENSION(prac_levels) :: znw_prac, znu_prac, dnw_prac

  REAL    :: p00, t00, a, pb, p_surf, temp, t_init
  REAL    :: alb(nzwrf), phb(nzwrf)
  REAL    :: mub, ztop, ztop_pbl, dz
  INTEGER :: loop, loop1


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begin of executable code below
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    IF(vertgrdopt == 0) THEN

!-----------------------------------------------------------------------
!
! Check the validation of vertical level specifications
!
!-----------------------------------------------------------------------

      ! Make sure first value is 1.0
      IF (znw(1) /= 1.0) THEN   ! User ask us to compute znw
                                ! nzwrf is used here.
!        WRITE(6,'(A,F6.3)') 'Bad first level: ',zlevels_wrf(1)
!        WRITE(6,'(A)')      'Mass coordinate must range from 1.0 at '     &
!                            //'surface to 0.0 at top.'
!        CALL arpsstop('Bad WRF vertical levels.',1)
        !
        !  Compute eta levels assuming a constant delta z above the PBL.
        !
        compute_const_dz = .TRUE.

         p00 = base_pres
         t00 = base_temp
         a   = base_lapse
         !  Compute top of the atmosphere with some silly levels.  We just want to
         !  integrate to get a reasonable value for ztop.  We use the planned PBL-esque
         !  levels, and then just coarse resolution above that.  We know p_top, and we
         !  have the base state vars.

         p_surf = p00

         znw_prac = (/ 1.000 , 0.993 , 0.983 , 0.970 , 0.954 , 0.934 , 0.909 , &
                       0.88 , 0.8 , 0.7 , 0.6 , 0.5 , 0.4 , 0.3 , 0.2 , 0.1 , 0.0 /)

         DO k = 1 , prac_levels - 1
            znu_prac(k) = ( znw_prac(k) + znw_prac(k+1) ) * 0.5
            dnw_prac(k) = znw_prac(k+1) - znw_prac(k)
         END DO

         DO k = 1, prac_levels-1
            pb = znu_prac(k)*(p_surf - ptop) + ptop
            temp =             t00 + A*LOG(pb/p00)
            t_init = temp*(p00/pb)**(r_d/cp_wrf) - t0
            alb(k) = (r_d/p1000mb)*(t_init+t0)*(pb/p1000mb)**cvpm
         END DO

         !  Base state mu is defined as base state surface pressure minus p_top

         mub = p_surf - ptop

         !  Integrate base geopotential, starting at terrain elevation.

         phb(1) = 0.
         DO k  = 2,prac_levels
               phb(k) = phb(k-1) - dnw_prac(k-1)*mub*alb(k-1)
         END DO

         !  So, now we know the model top in meters.  Get the average depth above the PBL
         !  of each of the remaining levels.  We are going for a constant delta z thickness.

         ztop     = phb(prac_levels) / g_wrf
         ztop_pbl = phb(8          ) / g_wrf
         dz = ( ztop - ztop_pbl ) / REAL ( nzwrf - 8 )

         !  Standard levels near the surface so no one gets in trouble.

         DO k = 1 , 8
            znw(k) = znw_prac(k)
         END DO

         !  Using d phb(k)/ d eta(k) = -mub * alb(k), eqn 2.9
         !  Skamarock et al, NCAR TN 468.  Use full levels, so
         !  use twice the thickness.

         DO k = 8, nzwrf-1
            pb = znw(k) * (p_surf - ptop) + ptop
            temp =             t00 + A*LOG(pb/p00)
            t_init = temp*(p00/pb)**(r_d/cp_wrf) - t0
            alb(k) = (r_d/p1000mb)*(t_init+t0)*(pb/p1000mb)**cvpm
            znw(k+1) = znw(k) - dz*g_wrf / ( mub*alb(k) )
         END DO
         znw(nzwrf) = 0.000

         !  There is some iteration.  We want the top level, ztop, to be
         !  consistent with the delta z, and we want the half level values
         !  to be consistent with the eta levels.  The inner loop to 10 gets
         !  the eta levels very accurately, but has a residual at the top, due
         !  to dz changing.  We reset dz five times, and then things seem OK.

         DO loop1 = 1 , 5
            DO loop = 1 , 10
               DO k = 8, nzwrf-1
                  pb = (znw(k)+znw(k+1))*0.5 * (p_surf - ptop) + ptop
                  temp =             t00 + A*LOG(pb/p00)
                  t_init = temp*(p00/pb)**(r_d/cp_wrf) - t0
                  alb(k) = (r_d/p1000mb)*(t_init+t0)*(pb/p1000mb)**cvpm
                  znw(k+1) = znw(k) - dz*g_wrf / ( mub*alb(k) )
               END DO
               IF ( ( loop1 .EQ. 5 ) .AND. ( loop .EQ. 10 ) ) THEN
                  print *,'Converged znw(kte) should be 0.0 = ',znw(nzwrf)
               END IF
               znw(nzwrf) = 0.000
            END DO

            !  Here is where we check the eta levels values we just computed.

            DO k = 1, nzwrf-1
               pb = (znw(k)+znw(k+1))*0.5 * (p_surf - ptop) + ptop
               temp =             t00 + A*LOG(pb/p00)
               t_init = temp*(p00/pb)**(r_d/cp_wrf) - t0
               alb(k) = (r_d/p1000mb)*(t_init+t0)*(pb/p1000mb)**cvpm
            END DO

            phb(1) = 0.
            DO k  = 2,nzwrf
                  phb(k) = phb(k-1) - (znw(k)-znw(k-1)) * mub*alb(k-1)
            END DO

            !  Reset the model top and the dz, and iterate.

            ztop = phb(nzwrf)/g_wrf
            ztop_pbl = phb(8)/g_wrf
            dz = ( ztop - ztop_pbl ) / REAL ( nzwrf - 8 )
         END DO

         IF ( dz .GT. maxdz ) THEN
           print *,'z (m)            = ',phb(1)/g_wrf
           do k = 2 ,nzwrf
             print *,'z (m) and dz (m) = ',phb(k)/g_wrf,(phb(k)-phb(k-1))/g_wrf
           end do
           print *,'dz (m) above fixed eta levels = ',dz
           print *,'namelist max_dz (m) = ',maxdz
           print *,'namelist p_top (Pa) = ',ptop
           WRITE( 0, '(2x,a)') 'You need one of three things:'
           WRITE( 0, '(4x,a)') '1) More eta levels to reduce the dz: e_vert'
           WRITE( 0, '(4x,a)') '2) A lower p_top so your total height is reduced: p_top_requested'
           WRITE( 0, '(4x,a)') '3) Increase the maximum allowable eta thickness: max_dz'
           WRITE( 0, '(7x,a)') 'All are namelist options'
           CALL arpsstop ( 'dz above fixed eta levels is too large in compute_eta',1)
         END IF

      ELSE         !  zlevels_wrf was specified explicitly

        nzwrf = 0

        nlevels = 1
        ! Make sure things are decreasing to 0.0 from bottom to top.
        level_check_loop:                                                   &
        DO k = 2, num_max_levels
          ! See if this is a valid level, if so increment nlevels
          ! and perform additional QC checks
          IF (znw(k) >= 0.) THEN
            ! Check for decreasing value of META (after 2nd value found)
            IF (znw(k) >= znw(k-1)) THEN
              PRINT '(A,I2,A,I2)', 'Level ',k, 'is >= level ' , k-1
              PRINT '(A)', 'Mass must be listed in descending order!'
              CALL arpsstop('Bad WRF vertical levels.',1)
            ELSE
              nlevels = nlevels + 1
            END IF
          ELSE
            IF (nlevels < 2) THEN
              PRINT '(A)', 'Not enough levels specified'
              CALL arpsstop('Bad WRF vertical levels.',1)
            ELSE
              nzwrf = nlevels
              ! Check to make sure top level is 0.0
              IF (znw(nzwrf) /= 0.) THEN
                PRINT '(A,F6.3)', 'Bad top level value: ',znw(nzwrf)
                PRINT '(A)', 'Top level must be 0.0 for Mass Eta.'
                CALL arpsstop('Bad WRF vertical levels.',1)
              END IF
              EXIT level_check_loop
            END IF
          END IF
        END DO level_check_loop
      END IF

    ELSE IF (vertgrdopt == 1) THEN
      IF(nzwrf < 15) nzwrf = 15
      DO k = 1,nzwrf
        znw(k) = (nzwrf-k) / (nzwrf - 1.)
      END DO
    ELSE IF (vertgrdopt == 2) THEN
      IF(nzwrf < 15) nzwrf = 15
      DO k = 1,nzwrf,1
        znw(k) = SQRT( (nzwrf-k)/(nzwrf-1.) )
      END DO
    ELSE IF (vertgrdopt == 3) THEN
      IF(nzwrf < 15) nzwrf = 15
      nlevels = 1

      DO k = 1, nzwrf
        sigma1 =  (k-1) / (nzwrf-1.)
        sigma2 =  SQRT( (k-1) /(nzwrf-1.) )

        plevel1 = sigma1 * (pbot - ptop) + ptop
        plevel2 = sigma2 * (pbot - ptop) + ptop

        IF (plevel1 < 33333.33) THEN
          znw(nlevels) = sigma1
          nlevels = nlevels + 1
        END IF
        IF (plevel2 > 33333.33) THEN
          znw(nlevels) = sigma2
          nlevels = nlevels + 1
        END IF
      END DO

      nzwrf = nlevels - 1
      WRITE(6,'(a,a,I3)') 'WARNING: the actual number of vertical',     &
                        ' levels has been changed to ',nzwrf
      !
      ! Sort the levels in decrease order (bubble sort for simplicity)
      !
      DO k = 1,nzwrf
        DO nlevels = 1, nzwrf-k
          IF (znw(nlevels+1) > znw(nlevels) ) THEN
            sigma1 = znw(nlevels)
            znw(nlevels) = znw(nlevels+1)
            znw(nlevels+1) = sigma1
          END IF
        END DO
      END DO

    END IF

  !
  ! Ouput vertical levels for testing
  !
  WRITE(6,'(/1x,2(a,I3),a,/)') 'Number of vertical levels is ',nzwrf,   &
                      ', vertical sigma scheme is ',vertgrdopt,'.'
  IF (compute_const_dz) THEN
    WRITE(6,'(1x,a,/,6(4x,a,F12.3,/))') 'zlevels was computed using:',  &
      'ptop   = ',ptop, 'p00    = ',p00,'t00    = ',t00,'a      = ',a,  &
      'max_dz = ',maxdz,'dz     = ',dz
  ELSE IF (vertgrdopt == 3) THEN
    WRITE(6,'(1x,a,/,2(4x,a,F12.0,/))') 'The parameters used are:',     &
            'ptop   = ',ptop, 'pbot   = ',pbot
  END IF

  WRITE(6,FMT='(1x,a)',ADVANCE='NO') 'levels = '
  DO nlevels = 1, nzwrf
    IF ( nlevels /= 1 .AND. MOD(nlevels,5) == 1 ) &
    WRITE(6,FMT='(10x)',ADVANCE='NO')
    WRITE(6,FMT='(F8.5,A1)',ADVANCE='NO') znw(nlevels),', '
    IF ( MOD(nlevels,5) == 0 ) WRITE(6,*)
  END DO
  WRITE(6,*)

  istatus = 0
  RETURN
END SUBROUTINE compute_eta

SUBROUTINE get_check_grid_ratio(domid,nx_wrf,ny_wrf,dx_wrf,dy_wrf,      &
                       dx_wrf_parent, dy_wrf_parent,                    &
                       a_small_number,parent_grid_ratio, istatus)

!-----------------------------------------------------------------------
!
! Compute parent_grid_ratio and check its validation
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: domid
  INTEGER, INTENT(IN)  :: nx_wrf, ny_wrf
  REAL,    INTENT(IN)  :: dx_wrf, dy_wrf
  REAL,    INTENT(IN)  :: dx_wrf_parent, dy_wrf_parent
  REAL,    INTENT(IN)  :: a_small_number
  INTEGER, INTENT(OUT) :: parent_grid_ratio
  INTEGER, INTENT(OUT) :: istatus

  INCLUDE 'mp.inc'   ! require variable myproc

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  REAL    :: gridratiox, gridratioy
  INTEGER :: modx, mody

!-----------------------------------------------------------------------
!
! Beginning of executable code ....
!
!-----------------------------------------------------------------------

  gridratiox = dx_wrf_parent / dx_wrf
  gridratioy = dy_wrf_parent / dy_wrf
  parent_grid_ratio = NINT(gridratiox)

  IF ( ABS(gridratiox - gridratioy) > a_small_number) THEN     ! ensure to have the same ratio
    WRITE(6,'(1x,2(a,F7.2),/,4x,a,2(F12.2,a),/,4x,a,2(F12.2,a))')              &
    'gridratiox = ',gridratiox,', gridratioy = ',gridratioy,                   &
    'Parent grid has dx_wrf/dy_wrf = ',dx_wrf_parent,' / ',dy_wrf_parent,'.',  &
    'dx_wrf/dy_wrf                 = ',dx_wrf,       ' / ',dy_wrf,'.'
    WRITE(6,'(4x,a,/)') 'parent_grid_ratio is different in x and y direction.'
    istatus = -1
    RETURN
  END IF

  IF (ABS(gridratiox-parent_grid_ratio) > a_small_number .OR.           &
      parent_grid_ratio < 1 .OR. MOD(parent_grid_ratio,2) == 0) THEN
                                    !  ensure parent_grid_ratio is odd
    WRITE(6,'(1x,2(a,I2),/,4x,a,2(F12.2,a),/,4x,a,2(F12.2,a)  )')       &
    'Parent_grid_ratio = ',parent_grid_ratio,' is not acceptable for domain ',domid, &
    'Parent grid has dx_wrf/dy_wrf = ',dx_wrf_parent,' / ',dy_wrf_parent,'.',        &
    'dx_wrf/dy_wrf                 = ',dx_wrf,       ' / ',dy_wrf,'.'
    WRITE(6,'(4x,a,/)') 'parent_grid_ratio must be odd for real-data cases.'
    istatus = -2
    RETURN
  END IF

  IF (myproc == 0) WRITE(6,'(/,1x,a,I2,a,I2)')                          &
    'Domain - ',domid,', parent_grid_ratio = ',parent_grid_ratio

  IF (domid > 1) THEN
    modx = MOD((nx_wrf-1),parent_grid_ratio)
    mody = MOD((ny_wrf-1),parent_grid_ratio)

    IF ( modx /= 0 .OR. mody /= 0 ) THEN
      WRITE(6,'(1x,a,I2,/,2(4x,a,I5,a,I2,a,I2,/))')                   &
         'Both MOD(nx_wrf-1,parent_grid_ratio) and '//                  &
         'MOD(ny_wrf-1,parent_grid_ratio) must be 0 for domain ',domid, &
         'MOD(nx_wrf-1,parent_grid_ratio) = MOD(',nx_wrf-1,',',         &
          parent_grid_ratio,') = ',modx,                                &
          'MOD(ny_wrf-1,parent_grid_ratio) = MOD(',ny_wrf-1,',',        &
          parent_grid_ratio,') = ',mody
      istatus = -3
      RETURN
    END IF
  END IF

  istatus = 0
  RETURN
END SUBROUTINE get_check_grid_ratio

SUBROUTINE set_mp_physics_variables(domid,mp_physics,istatus)

!-----------------------------------------------------------------------
!
! PURPOSE:
!   Set microphysics variable based on namelist parameter mp_physics
!
! NOTE:
!
!   Based on Register.EM in WRFV2.2
!
!-----------------------------------------------------------------------

  USE wrf_metadata

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: domid
  INTEGER, INTENT(IN)  :: mp_physics
  INTEGER, INTENT(OUT) :: istatus

  INCLUDE 'globcst.inc'

!-----------------------------------------------------------------------
!
! Begin of executable code below ....
!
!-----------------------------------------------------------------------

  num_moist = 0

  P_QT_wrf  = 0
  P_QNI_wrf = 0
  num_scalar = 0

  istatus = 0

  SELECT CASE (mp_physics)
  CASE (0)                           ! passiveqv
    P_QV_wrf  = 1
    num_moist = 1
  CASE (1,3,98)                      ! kesslerscheme, wsm3scheme, ncepcloud3
    P_QV_wrf  = 1
    P_QC_wrf  = 2
    P_QR_wrf  = 3
    num_moist = 3
  CASE (2,6)                         ! linscheme, wsm6scheme
    P_QV_wrf  = 1
    P_QC_wrf  = 2
    P_QR_wrf  = 3
    P_QI_wrf  = 4
    P_QS_wrf  = 5
    P_QG_wrf  = 6
    num_moist = 6
  CASE (4,14,99)                    ! wsm5scheme,ncepcloud5
    P_QV_wrf  = 1
    P_QC_wrf  = 2
    P_QR_wrf  = 3
    P_QI_wrf  = 4
    P_QS_wrf  = 5
    num_moist = 5
  CASE (5)                           ! etampnew
    P_QV_wrf  = 1
    P_QC_wrf  = 2
    P_QR_wrf  = 3
    P_QI_wrf  = 4
    P_QS_wrf  = 5
    P_QG_wrf  = 6
    num_moist = 6

    P_QT_wrf    = 1
    num_scalar  = 1
  CASE (8,10,16)                     ! thompson
    P_QV_wrf  = 1
    P_QC_wrf  = 2
    P_QR_wrf  = 3
    P_QI_wrf  = 4
    P_QS_wrf  = 5
    P_QG_wrf  = 6
    num_moist = 6

    P_QNI_wrf   = 1
    num_scalar  = 1

  CASE DEFAULT
    istatus = -1
    WRITE(6,'(/,1x,2(a,I2),/)') 'ERROR: Wrong parameter - mp_physics = ',&
                                mp_physics,' for dom - ',domid
  END SELECT

  nscalar_wrf = num_moist-1
  arps_Q_ptr(:)    = 0
  IF (P_QC_wrf > 0 ) THEN
    IF (P_QC > 0) arps_Q_ptr(P_QC_wrf-1) = P_QC
    qnames_wrf(P_QC_wrf-1) = 'QCLOUD'
    qdescp_wrf(P_QC_wrf-1) = 'Cloud water mixing ratio'
  END IF

  IF (P_QR_wrf > 0 ) THEN
    IF (P_QR > 0) arps_Q_ptr(P_QR_wrf-1) = P_QR
    qnames_wrf(P_QR_wrf-1) = 'QRAIN'
    qdescp_wrf(P_QR_wrf-1) = 'Rain water mixing ratio'
  END IF

  IF (P_QS_wrf > 0 ) THEN
    IF (P_QS > 0) arps_Q_ptr(P_QS_wrf-1) = P_QS
    qnames_wrf(P_QS_wrf-1) = 'QSNOW'
    qdescp_wrf(P_QS_wrf-1) = 'Snow mixing ratio'
  END IF

  IF (P_QI_wrf > 0 ) THEN
    IF (P_QI > 0) arps_Q_ptr(P_QI_wrf-1) = P_QI
    qnames_wrf(P_QI_wrf-1) = 'QICE'
    qdescp_wrf(P_QI_wrf-1) = 'Ice mixing ratio'
  END IF

  IF (P_QG_wrf > 0 ) THEN
    IF (P_QG > 0 .AND. P_QH > 0) THEN
      arps_Q_ptr(P_QG_wrf-1) = 22222  ! Special handling to join P_QG & P_QH
    ELSE IF (P_QG > 0) THEN
      arps_Q_ptr(P_QG_wrf-1) = P_QG
    ELSE IF (P_QH > 0) THEN
      arps_Q_ptr(P_QG_wrf-1) = P_QH
    END IF
    qnames_wrf(P_QG_wrf-1) = 'QGRAUP'
    qdescp_wrf(P_QG_wrf-1) = 'Graupel mixing ratio'
  END IF

  RETURN
END SUBROUTINE set_mp_physics_variables
