!
!##################################################################
!##################################################################
!######                                                      ######
!######                   PROGRAM WRF2ARPS                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
PROGRAM wrfextsnd
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Extracts soundings from the WRF files and variables
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  2007 Feb 15 Based on WRF2ARPS by Yunheng Wang, Dan Weber
!              and extsnd program from ARPS package.
!
!  MODIFICATION HISTORY:
!
!  Yunheng Wang (03/12/2007)
!  Merged getwrfds with getwrfd so that it has the same interface as WRF2ARPS.
!  Added a new parameter "grid_id" for WRF nesting outputs.
!  Stopped the program when a station location exceeds the WRF domain.
!
!  Yunheng Wang (03/26/2007)
!  Removed namelist parameter frames_per_outfile. It is not determined by
!  the program automatically.
!
!-----------------------------------------------------------------------
!
  USE module_wrf2arps_post

  IMPLICIT NONE

  INTEGER :: nx_wrf,ny_wrf,nz_wrf
  INTEGER :: nzsoil_wrf,nstyp_wrf
  INTEGER :: iproj_wrf              ! external data map projection
  REAL    :: scale_wrf              ! external data map scale factor
  REAL    :: trlon_wrf              ! external data true longitude
  REAL    :: trlat1_wrf,trlat2_wrf,latnot_wrf(2)
                                    ! external data true latitude(s)
  REAL    :: ctrlat_wrf, ctrlon_wrf

  REAL    :: x0_wrf,y0_wrf          ! external data origin
  REAL    :: dx_wrf,dy_wrf,dt_wrf

  INTEGER :: sf_surface_physics, sf_sfclay_physics, mp_physics
!
!-----------------------------------------------------------------------
!
!  WRF forecast variables
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: x_wrf(:)        ! external data x-coordinate
  REAL, ALLOCATABLE :: y_wrf(:)        ! external data y-coordinate
  REAL, ALLOCATABLE :: xu_wrf(:)       ! external data u x-coordinate
  REAL, ALLOCATABLE :: yu_wrf(:)       ! external data u y-coordinate
  REAL, ALLOCATABLE :: xv_wrf(:)       ! external data v x-coordinate
  REAL, ALLOCATABLE :: yv_wrf(:)       ! external data v y-coordinate

  REAL, ALLOCATABLE :: lat_wrf(:,:)    ! external data latidude
  REAL, ALLOCATABLE :: lon_wrf(:,:)    ! external data longitude
  REAL, ALLOCATABLE :: latu_wrf(:,:)   ! external data latidude (x-stag)
  REAL, ALLOCATABLE :: lonu_wrf(:,:)   ! external data longitude (x-stag)
  REAL, ALLOCATABLE :: latv_wrf(:,:)   ! external data latidude (y-stag)
  REAL, ALLOCATABLE :: lonv_wrf(:,:)   ! external data longitude (y-stag)

  REAL, ALLOCATABLE :: zp_wrf(:,:,:)   ! external data physical height (m)
  REAL, ALLOCATABLE :: hgt_wrf(:,:,:)  ! Height (m) of scalar points
  REAL, ALLOCATABLE :: zpsoil_wrf(:,:,:)! Height (m) of soil layers

  REAL, ALLOCATABLE :: p_wrf(:,:,:)    ! Pressure (Pascals)
  REAL, ALLOCATABLE :: pt_wrf(:,:,:)   ! Potential Temperature (K)
  REAL, ALLOCATABLE :: t_wrf(:,:,:)    ! Temperature (K)
  REAL, ALLOCATABLE :: u_wrf(:,:,:)    ! Eastward wind component
  REAL, ALLOCATABLE :: v_wrf(:,:,:)    ! Northward wind component
  REAL, ALLOCATABLE :: vatu_wrf(:,:,:) ! NOTE: only used when use_wrf_grid /= 1
  REAL, ALLOCATABLE :: uatv_wrf(:,:,:) !
  REAL, ALLOCATABLE :: w_wrf(:,:,:)    ! Vertical wind component
  REAL, ALLOCATABLE :: qv_wrf(:,:,:)   ! Specific humidity (kg/kg)

  REAL, ALLOCATABLE :: qscalar_wrf(:,:,:,:)

  REAL, ALLOCATABLE :: rainf(:,:)   ! Total accumulated rainfall
  REAL, ALLOCATABLE :: cumrain(:,:)   ! Total accumulated cumulus rainfall
  REAL, ALLOCATABLE :: snowf(:,:)   ! Total accumulated snow and ice
  REAL, ALLOCATABLE :: graupel(:,:)   ! Total accumulated graupel

  REAL, ALLOCATABLE :: tsoil_wrf  (:,:,:,:)   ! Soil temperature  (K)
  REAL, ALLOCATABLE :: qsoil_wrf  (:,:,:,:)   ! Soil moisture (m3/m3)
  REAL, ALLOCATABLE :: wetcanp_wrf(:,:,:)     ! Canopy water amount

  REAL, ALLOCATABLE :: snowdpth_wrf(:,:)      ! Snow depth (m)
  REAL, ALLOCATABLE :: trn_wrf     (:,:)      ! External terrain (m)

  INTEGER, ALLOCATABLE :: soiltyp_wrf (:,:,:) ! Soil type
  REAL,    ALLOCATABLE :: stypfrct_wrf(:,:,:) ! Soil type fraction
  INTEGER, ALLOCATABLE :: vegtyp_wrf  (:,:)   ! Vegetation type
  REAL,    ALLOCATABLE :: veg_wrf     (:,:)   ! Vegetation fraction

!-------------------------------------------------------------------
!
!  Working arrays
!
!-------------------------------------------------------------------

  REAL, ALLOCATABLE :: tc(:)
  REAL, ALLOCATABLE :: tdc(:)
  REAL, ALLOCATABLE :: usnd(:)
  REAL, ALLOCATABLE :: vsnd(:)
  REAL, ALLOCATABLE :: dir(:)
  REAL, ALLOCATABLE :: spd(:)

  REAL, ALLOCATABLE :: dxfld(:)         ! on WRF grid
  REAL, ALLOCATABLE :: dyfld(:)
  REAL, ALLOCATABLE :: rdxfld(:)
  REAL, ALLOCATABLE :: rdyfld(:)
  REAL, ALLOCATABLE :: dxfldu(:)
  REAL, ALLOCATABLE :: dyfldu(:)
  REAL, ALLOCATABLE :: rdxfldu(:)
  REAL, ALLOCATABLE :: rdyfldu(:)
  REAL, ALLOCATABLE :: dxfldv(:)
  REAL, ALLOCATABLE :: dyfldv(:)
  REAL, ALLOCATABLE :: rdxfldv(:)
  REAL, ALLOCATABLE :: rdyfldv(:)

  REAL, ALLOCATABLE :: tem1_wrf(:,:,:)    ! Temporary work array
  REAL, ALLOCATABLE :: tem2_wrf(:,:,:)    ! Temporary work array
  REAL, ALLOCATABLE :: tem3_wrf(:,:,:)    ! Temporary work array
  REAL, ALLOCATABLE :: tem4_wrf(:,:,:)    ! Temporary work array
  REAL, ALLOCATABLE :: tem5_wrf(:,:,:)    ! Temporary work array
  REAL, ALLOCATABLE :: tem1_2d(:,:)       ! Temporary work array

  REAL, ALLOCATABLE :: xa_wrf(:,:)        ! WRF x coordinate on ARPS grid
  REAL, ALLOCATABLE :: ya_wrf(:,:)        ! WRF y coordinate on ARPS grid

!
!-----------------------------------------------------------------------
!
!  include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'bndry.inc'
  INCLUDE 'adjust.inc'
  INCLUDE 'mp.inc'

  INTEGER, PARAMETER :: maxwrffil = 100
  INTEGER, PARAMETER :: maxsnd = 500
  INTEGER, PARAMETER :: maxtime = 1000
!
!-----------------------------------------------------------------------
!
!  NAMELIST parameters
!
!-----------------------------------------------------------------------
!
  INTEGER            :: nprocx_in, nprocy_in

  CHARACTER(LEN=256) :: dir_extd            ! directory of external data
  CHARACTER(LEN=256) :: dir_output
  INTEGER            :: io_form
  LOGICAL            :: multifile

  CHARACTER(LEN=19)  :: init_time_str,start_time_str,end_time_str
  CHARACTER(LEN=11)  :: history_interval

  INTEGER :: outsnd
  INTEGER :: outsfc

  REAL :: ugrid
  REAL :: vgrid
  INTEGER :: locopt
  INTEGER :: nsnd
  REAL :: xpt(maxsnd)
  REAL :: ypt(maxsnd)
  CHARACTER (LEN=512) :: uafile
  CHARACTER (LEN=512) :: sfcfile
  CHARACTER (LEN=8) :: stid(maxsnd)
  REAL :: slat(maxsnd)
  REAL :: slon(maxsnd)
  REAL :: selev(maxsnd)
  INTEGER :: ipt(maxsnd)
  INTEGER :: jpt(maxsnd)
  INTEGER :: istnm(maxsnd)
  LOGICAL :: ptingrid(maxsnd)

  REAL :: sfcpr(maxtime,maxsnd)
  REAL :: tsfcf(maxtime,maxsnd)
  REAL :: tdsfcf(maxtime,maxsnd)
  REAL :: wdir(maxtime,maxsnd)
  REAL :: wspd(maxtime,maxsnd)
  REAL :: prcgp(maxtime,maxsnd)
  REAL :: prccu(maxtime,maxsnd)
  REAL :: prctot(maxtime,maxsnd)

  REAL :: tmax(maxsnd)
  REAL :: tmin(maxsnd)
  REAL :: spdmax(maxsnd)
  REAL :: prc0(maxsnd)
  REAL :: prcsum(maxsnd)

  INTEGER :: dmp_out_joined

  INTEGER :: grid_id

  NAMELIST /message_passing/ nproc_x, nproc_y, max_fopen,               &
                             nprocx_in, nprocy_in

  NAMELIST /wrfdfile/ dir_extd,init_time_str,io_form,grid_id,           &
                      start_time_str,history_interval,end_time_str

  NAMELIST /sndloc/ ugrid,vgrid,locopt,nsnd,                            &
                    slat,slon,selev,xpt,ypt,stid,istnm

  NAMELIST /output/ dir_output,outsnd,outsfc

  INTEGER :: ncompressx, ncompressy ! compression in x and y direction:
                                    ! ncompressx=nprocx_in/nproc_x
                                    ! ncompressy=nprocy_in/nproc_y
!-----------------------------------------------------------------------
!
!  Non-dimensional smoothing coefficient
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: smfct1 = 0.5
  REAL, PARAMETER :: smfct2 = -0.5
  REAL, PARAMETER :: rhmax  = 1.0
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: ms2kts =1.9438
  REAL, PARAMETER :: mm2inch=(1.0/25.4)

  INTEGER, PARAMETER :: isec24h = 86400
  CHARACTER (LEN=300) :: systring

  INTEGER :: isnow,jsnow,ii,jj
  INTEGER :: is, idummy
  INTEGER :: kbot,ktop

  REAL :: amin,amax
  REAL :: xumin,xumax,yvmin,yvmax
  REAL :: qvmin,qvmax,qvval
  REAL :: scrhgt,annhgt,csconst,pconst,tv_wrf,tsfc,tf,tdf,tdsfc
  REAL :: tcsoil,delt,delz,spdkts,rain_g,rain_c,raint
  REAL :: pres,temp,qvsat,rh,tvbar,qvprt,qtot,dirsfc,spdsfc
  REAL :: tctof,tktof

  CHARACTER(LEN=256) :: extdname(MAXWRFFIL)
  CHARACTER(LEN=256) :: tmpstr

  INTEGER, ALLOCATABLE :: fHndl(:,:)
  CHARACTER(LEN=19)    :: timestr
  CHARACTER(LEN=15)    :: timefil(maxtime)
  CHARACTER(LEN=17)    :: timewrf(maxtime)

  REAL               :: latnot(2)
  REAL               :: deltaz,omegasnd

  INTEGER            :: i,j,k,nq,ksmth
  INTEGER            :: strlen
  INTEGER            :: istatus
  INTEGER            :: nextdfil
  INTEGER            :: ifile,itime,jtime,ntime,isnd
  INTEGER            :: iniotfu

  INTEGER            :: iextmn,iextmx,jextmn,jextmx

  INTEGER            :: iyr,imo,iday,ihr,imin,isec,jldy
  INTEGER            :: myr,modyr,initsec,iabssec,istverif,kftime
  LOGICAL            :: rewindyr

  LOGICAL            :: fexist
  LOGICAL            :: first_time
  LOGICAL            :: in_verif
  LOGICAL            :: wait_verif
  LOGICAL            :: exit_early

  CHARACTER(LEN=1)   :: ach

  INTEGER            :: idist
  DOUBLE PRECISION   :: ntmergeinv, mfac

  REAL    :: ppasc,pmb,theta,smix,e,bige,alge
  REAL    :: dd,dmin,latd,lond

  INTEGER :: abstimes,abstimee,abstimei

  INTEGER :: iloc, jloc
  INTEGER :: frames_per_outfile(MAXWRFFIL)

  INTEGER :: num_scalar
!
!-----------------------------------------------------------------------
!
! roufns & lai convert table (see src/arpssfc/sfc_winter.tbl)
!
!-----------------------------------------------------------------------

  REAL, PARAMETER :: rfnstbl(14) = (/0.002, 0.02, 0.01, 0.03, 0.10,     &
                                     0.20,  0.40, 2.00, 0.005,0.01,     &
                                     0.02,  0.06, 0.04, 0.002 /)

  REAL, PARAMETER :: laitbl(14)  = (/0.5,   0.5,  0.02, 0.02, 0.02,     &
                                     0.05,  0.5,  0.5,  0.0,  0.02,     &
                                     0.0,   0.5,  0.5,  0.0   /)
!
!-----------------------------------------------------------------------
!
!  Function f_qvsat and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_qvsat  ! compute saturation specific humidity defined in
                   ! thermolib3d.f90

!fpp$ expand (f_qvsat)
!!dir$ inline always f_qvsat
!*$*  inline routine (f_qvsat)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!-----------------------------------------------------------------------
!
!  MIsc Initializations
!
!-----------------------------------------------------------------------
!
  jtime=0
  ntime=0
  tmax=-999.0
  tmin=999.0
  spdmax=-999.0
  prc0=-99.0
  prcsum=0.0
  mp_physics=2
!
!-----------------------------------------------------------------------
!
!  Initialize message passing processors.
!
!-----------------------------------------------------------------------
!
  ! Non-MPI defaults: All others initialized in mpinit_var
  mp_opt = 0
  myproc = 0
  nprocx_in  = 1
  nprocy_in  = 1
  dumpstride = 1
  readstride = 1

  CALL mpinit_proc(0)

  IF(myproc == 0) WRITE(6,'(10(/5x,a),/)')                                  &
      '###################################################################',&
      '###################################################################',&
      '####                                                           ####',&
      '####                Welcome to WRFEXTSND                       ####',&
      '####                                                           ####',&
      '####   Program that reads in output from the WRF model         ####',&
      '####             and produces text soundings.                  ####',&
      '####                                                           ####',&
      '###################################################################',&
      '###################################################################'


  mgrid = 1
  nestgrd = 0
!
!-----------------------------------------------------------------------
!
!  Read in message passing options.
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) THEN
    READ (5,message_passing)
    WRITE(6,'(a)')'  Namelist block message_passing sucessfully read.'
  END IF
  CALL mpupdatei(nproc_x,1)
  CALL mpupdatei(nproc_y,1)
  CALL mpupdatei(max_fopen,1)
  CALL mpupdatei(nprocx_in,1)
  CALL mpupdatei(nprocy_in,1)
!
!-----------------------------------------------------------------------
!
!  Initialize message passing variables.
!
!-----------------------------------------------------------------------
!
  CALL mpinit_var

!
!-----------------------------------------------------------------------
!
!  Read in namelist &wrfdfile
!
!-----------------------------------------------------------------------
!
  dir_extd = './'

  init_time_str         = '0000-00-00_00:00:00'
  start_time_str        = '0000-00-00_00:00:00'
  history_interval      = '00_00:00:00'
  end_time_str          = '0000-00-00_00:00:00'

  io_form               = 7
  grid_id               = 1

  multifile             = .FALSE.          ! not in namelist
  IF (myproc == 0) THEN
    READ(5,wrfdfile)
    WRITE(6,'(a)') '  Namelist wrfdfile read in successfully.'

    strlen = LEN_TRIM(dir_extd)
    IF(strlen > 0) THEN
      IF(dir_extd(strlen:strlen) /= '/') THEN
        dir_extd(strlen+1:strlen+1) = '/'
        strlen = strlen + 1
      END IF
    ELSE
      dir_extd = './'
    END IF

    IF (io_form > 100) THEN
      io_form = MOD(io_form,100)
      IF (mp_opt > 0) multifile = .TRUE.
    END IF

  END IF   ! myproc == 0
  CALL mpupdatec(dir_extd,256)
  CALL mpupdatei(io_form,1)
  CALL mpupdatel(multifile,1)
  CALL mpupdatec(init_time_str,19)
  CALL mpupdatec(start_time_str,19)
  CALL mpupdatec(end_time_str,19)
  CALL mpupdatec(history_interval,11)
  CALL mpupdatei(grid_id, 1)

!
!-----------------------------------------------------------------------
!
!  Read in namelist &wrfdfile
!
!-----------------------------------------------------------------------
!
  nsnd=1
  omegasnd=0.0
  xpt=-999.0
  ypt=-999.0
  slat=-999.0
  slon=-999.0
  ipt=-99
  jpt=-99

  IF (myproc == 0 ) THEN
    READ(5,sndloc)
    IF( nsnd > maxsnd )  then
      WRITE(6,'(a,/a,i5)')                                                &
      'The number of sounding locations to be extracted exceeded maximum ',&
      'allowed. nsnd is reset to ', maxsnd
      nsnd = maxsnd
    ENDIF
    DO isnd=1,nsnd
      xpt(isnd)=1000.*xpt(isnd)
      ypt(isnd)=1000.*ypt(isnd)
    END DO
  ENDIF

  CALL mpupdater(ugrid,1)
  CALL mpupdater(vgrid,1)
  CALL mpupdatei(locopt,1)
  CALL mpupdatei(nsnd,1)
  CALL mpupdater(xpt,maxsnd)
  CALL mpupdater(ypt,maxsnd)
  CALL mpupdater(slat,maxsnd)
  CALL mpupdater(slon,maxsnd)
  CALL mpupdatec(stid,8*maxsnd)
  CALL mpupdatei(istnm,maxsnd)

  IF (myproc == 0 ) THEN
    READ(5,output)
  ENDIF
  CALL mpupdatec(dir_output,256)
  CALL mpupdatei(outsnd,1)
  CALL mpupdatei(outsfc,1)

  IF (myproc == 0 ) THEN
    WRITE(6,'(4x,a,a)')       '  Output directory: ',TRIM(dir_output)
    WRITE(6,'(2(4x,a,i3),/)') '  outsnd: ',outsnd,'  outsfc: ',outsfc
  END IF

!=======================================================================
!
! NAMELIST readings are done
!
!=======================================================================

  rewindyr = .FALSE.

  READ(start_time_str,'(I4.4,5(a,I2.2))')      &
                  year,ach,month,ach,day,ach,hour,ach,minute,ach,second
  IF (year < 1960) THEN   ! maybe ideal case
    myr  =  year
    year =  1960
    rewindyr = .TRUE.
  END IF
  CALL ctim2abss(year,month,day,hour,minute,second,abstimes)

  READ(end_time_str,'(I4.4,a,I2.2,a,I2.2,a,I2.2,a,I2.2,a,I2.2)')        &
                  year,ach,month,ach,day,ach,hour,ach,minute,ach,second
  IF (rewindyr)  year = 1960
  CALL ctim2abss(year,month,day,hour,minute,second,abstimee)

  READ(history_interval,'(I2.2,a,I2.2,a,I2.2,a,I2.2)')                  &
                                     day,ach,hour,ach,minute,ach,second
  abstimei = day*24*3600+hour*3600+minute*60+second

  IF (multifile) THEN
    IF (MOD(nprocx_in,nproc_x) /= 0 .OR. MOD(nprocy_in,nproc_y) /= 0) THEN
      WRITE(6,*) 'nprocx_in (nprocy_in) must be dividable by nproc_x (nproc_y).'
      CALL arpsstop('WRONG message passing parameter.',1)
    END IF

    ncompressx = nprocx_in/nproc_x
    ncompressy = nprocy_in/nproc_y

  ELSE     ! non-mpi or mpi with one file

    ncompressx = 1
    ncompressy = 1

  END IF

  ALLOCATE(fHndl(ncompressx,ncompressy), STAT = istatus)

!-----------------------------------------------------------------------
!
! Check the availability of files and get parameter frames_per_outfile
!
!-----------------------------------------------------------------------

  frames_per_outfile(:) = 1
  nextdfil = 0

  CALL check_wrf_files(multifile,MAXWRFFIL,grid_id,io_form,nprocx_in,   &
       ncompressx,ncompressy,abstimes,abstimei,abstimee,rewindyr,myr,   &
       dir_extd,extdname,nextdfil,frames_per_outfile,istatus)

  IF (istatus /= 0) CALL arpsstop('ERROR in check_wrf_files, See STDOUT for details',1)


  joindmp(:) = 1
  IF (mp_opt > 0) THEN        ! should moved into mpinit_var later
    dumpstride = max_fopen
    readstride = nprocs
    IF (ANY(joindmp > 0)) dumpstride = nprocs   ! join and dump
    IF (multifile)        readstride = max_fopen
  END IF

  IF (ANY(frames_per_outfile > 1) .AND. readstride < nprocs) THEN
    WRITE(6,'(3(a,/))') 'WARNING: WRF2ARPS does not support multi-frame in ',&
               '         one file for split-form WRF history files.',        &
               '         The program is stopping ... ...'
    CALL arpsstop('frames_per_outfile should be 1.',1)
  END IF

!-----------------------------------------------------------------------
!
! Get dimension parameters from the first input file
!
!-----------------------------------------------------------------------

!  blocking inserted for ordering i/o for message passing
  DO i=0,nprocs-1,readstride
    IF(myproc >= i .AND. myproc <= i+readstride-1) THEN

      CALL open_wrf_file(TRIM(extdname(1)),io_form,multifile,.TRUE.,    &
                         ncompressx,ncompressy,fHndl)

      CALL get_wrf_metadata(fHndl,io_form,multifile,.TRUE.,1,1,         &
                         nx_wrf,ny_wrf,nz_wrf,nzsoil_wrf,               &
                         iproj_wrf,trlat1_wrf,trlat2_wrf,trlon_wrf,     &
                         ctrlat_wrf,ctrlon_wrf,dx_wrf,dy_wrf,dt_wrf,    &
                         sf_surface_physics,sf_sfclay_physics,          &
                         mp_physics,num_scalar,istatus)

      CALL close_wrf_file(fHndl,io_form,multifile,.TRUE.,ncompressx,ncompressy)

    END IF
    IF (mp_opt > 0) CALL mpbarrier
  END DO

  scale_wrf = 1.0
  nstyp_wrf = 1

  IF(myproc == 0) WRITE(6,'(/a,4(a,I4),/)') '  WRF grid dimensions:',   &
                       ' nx_wrf = ',nx_wrf, ' ny_wrf = ',ny_wrf,        &
                       ' nz_wrf = ',nz_wrf, ' nzsoil_wrf = ', nzsoil_wrf

  ! WRF forecast started time
  READ(init_time_str,'(I4.4,a,I2.2,a,I2.2,a,I2.2,a,I2.2,a,I2.2)')       &
                 year,ach,month,ach,day,ach,hour,ach,minute,ach,second
  CALL ctim2abss(year,month,day,hour,minute,second,initsec)

  thisdmp = abstimei
  tstart  = 0.0
  tstop   = abstimee-initsec
  latitud = ctrlat

  CALL init_data2d(0,0,0,abstimei,' ',' ',0,0,                          &
                   0,0,0,0,                                             &
                   0,0,0,0,                                             &
                   0,0,0,0,                                             &
                   .FALSE., istatus)

!-----------------------------------------------------------------------
!
!  Allocate and initalize arrays based on dimension parameters
!  read in from the input file
!
!-----------------------------------------------------------------------

  ! Note that for MP version nx & ny here are global values.  They will
  ! be reassigned to their per-processor value below.

  IF (mp_opt > 0) THEN
    nx_wrf = (nx_wrf - 1)/nproc_x + 1     ! fake zone for WRF is 1
    ny_wrf = (ny_wrf - 1)/nproc_y + 1
  ELSE
    nproc_x = 1
    nproc_y = 1
    nprocs  = 1
    joindmp(:) = 0
  END IF
!
!-----------------------------------------------------------------------
!
!  Allocate and initialize external grid variables
!
!-----------------------------------------------------------------------
!
  ALLOCATE(x_wrf(nx_wrf),stat=istatus)
  ALLOCATE(y_wrf(ny_wrf),stat=istatus)
  ALLOCATE(xu_wrf(nx_wrf),stat=istatus)
  ALLOCATE(yu_wrf(ny_wrf),stat=istatus)
  ALLOCATE(xv_wrf(nx_wrf),stat=istatus)
  ALLOCATE(yv_wrf(ny_wrf),stat=istatus)

  ALLOCATE(lat_wrf(nx_wrf,ny_wrf),stat=istatus)
  ALLOCATE(lon_wrf(nx_wrf,ny_wrf),stat=istatus)
  ALLOCATE(latu_wrf(nx_wrf,ny_wrf),stat=istatus)
  ALLOCATE(lonu_wrf(nx_wrf,ny_wrf),stat=istatus)
  ALLOCATE(latv_wrf(nx_wrf,ny_wrf),stat=istatus)
  ALLOCATE(lonv_wrf(nx_wrf,ny_wrf),stat=istatus)

  ALLOCATE(zp_wrf(nx_wrf,ny_wrf,nz_wrf),stat=istatus)
  ALLOCATE(hgt_wrf(nx_wrf,ny_wrf,nz_wrf),stat=istatus)
  ALLOCATE(zpsoil_wrf(nx_wrf,ny_wrf,nzsoil_wrf),stat=istatus)

  ALLOCATE(p_wrf(nx_wrf,ny_wrf,nz_wrf),stat=istatus)
  ALLOCATE(t_wrf(nx_wrf,ny_wrf,nz_wrf),stat=istatus)
  ALLOCATE(u_wrf(nx_wrf,ny_wrf,nz_wrf),stat=istatus)
  ALLOCATE(v_wrf(nx_wrf,ny_wrf,nz_wrf),stat=istatus)
  ALLOCATE(vatu_wrf(nx_wrf,ny_wrf,nz_wrf),stat=istatus)
  ALLOCATE(uatv_wrf(nx_wrf,ny_wrf,nz_wrf),stat=istatus)
  ALLOCATE(w_wrf(nx_wrf,ny_wrf,nz_wrf),stat=istatus)
  ALLOCATE(pt_wrf(nx_wrf,ny_wrf,nz_wrf),stat=istatus)
  ALLOCATE(qv_wrf(nx_wrf,ny_wrf,nz_wrf),stat=istatus)
  ALLOCATE(qscalar_wrf(nx_wrf,ny_wrf,nz_wrf,nscalar),stat=istatus)

  ALLOCATE(tsoil_wrf   (nx_wrf,ny_wrf,nzsoil_wrf,0:nstyp_wrf),stat=istatus)
  ALLOCATE(qsoil_wrf   (nx_wrf,ny_wrf,nzsoil_wrf,0:nstyp_wrf),stat=istatus)
  ALLOCATE(wetcanp_wrf (nx_wrf,ny_wrf,0:nstyp_wrf),stat=istatus)
  ALLOCATE(soiltyp_wrf (nx_wrf,ny_wrf,nstyp_wrf),stat=istatus)
  ALLOCATE(stypfrct_wrf(nx_wrf,ny_wrf,nstyp_wrf),stat=istatus)

  ALLOCATE(vegtyp_wrf  (nx_wrf,ny_wrf),stat=istatus)
  ALLOCATE(veg_wrf     (nx_wrf,ny_wrf),stat=istatus)
  ALLOCATE(snowdpth_wrf(nx_wrf,ny_wrf),stat=istatus)
  ALLOCATE(trn_wrf     (nx_wrf,ny_wrf),stat=istatus)

  ALLOCATE(tc(nz_wrf))
  ALLOCATE(tdc(nz_wrf))
  ALLOCATE(usnd(nz_wrf))
  ALLOCATE(vsnd(nz_wrf))
  ALLOCATE(dir(nz_wrf))
  ALLOCATE(spd(nz_wrf))

  x_wrf=0.0
  y_wrf=0.0
  xu_wrf=0.0
  yu_wrf=0.0
  xv_wrf=0.0
  yv_wrf=0.0

  lat_wrf=0.0
  lon_wrf=0.0
  latu_wrf=0.0
  lonu_wrf=0.0
  latv_wrf=0.0
  lonv_wrf=0.0

  zp_wrf=0.0
  hgt_wrf=0.0
  zpsoil_wrf=0.0
  p_wrf=0.0
  t_wrf=0.0
  u_wrf=0.0
  v_wrf=0.0
  vatu_wrf=0.0
  uatv_wrf=0.0
  w_wrf=0.0
  pt_wrf=0.0
  qv_wrf=0.0
  qscalar_wrf=0.0

  trn_wrf     =0.0

  tsoil_wrf   =-999.0
  qsoil_wrf   =-999.0
  wetcanp_wrf =-999.0
  snowdpth_wrf=-999.0

  soiltyp_wrf = -999
  stypfrct_wrf= -999.0
  vegtyp_wrf  = -999
  veg_wrf     = -999.0

  CALL allocate_data2d_ext( nx_wrf,ny_wrf, istatus )

!-----------------------------------------------------------------------
!
! Allocate working arrays
!
!-----------------------------------------------------------------------

  ALLOCATE(dxfld(nx_wrf),stat=istatus)
  ALLOCATE(dyfld(ny_wrf),stat=istatus)
  ALLOCATE(rdxfld(nx_wrf),stat=istatus)
  ALLOCATE(rdyfld(ny_wrf),stat=istatus)
  ALLOCATE(dxfldu(nx_wrf),stat=istatus)
  ALLOCATE(dyfldu(ny_wrf),stat=istatus)
  ALLOCATE(rdxfldu(nx_wrf),stat=istatus)
  ALLOCATE(rdyfldu(ny_wrf),stat=istatus)
  ALLOCATE(dxfldv(nx_wrf),stat=istatus)
  ALLOCATE(dyfldv(ny_wrf),stat=istatus)
  ALLOCATE(rdxfldv(nx_wrf),stat=istatus)
  ALLOCATE(rdyfldv(ny_wrf),stat=istatus)

  ALLOCATE(tem1_wrf(nx_wrf,ny_wrf,nz_wrf),stat=istatus)
  ALLOCATE(tem2_wrf(nx_wrf,ny_wrf,nz_wrf),stat=istatus)
  ALLOCATE(tem3_wrf(nx_wrf,ny_wrf,nz_wrf),stat=istatus)
  ALLOCATE(tem4_wrf(nx_wrf,ny_wrf,nz_wrf),stat=istatus)
  ALLOCATE(tem5_wrf(nx_wrf,ny_wrf,nz_wrf),stat=istatus)

  ALLOCATE(tem1_2d(nx_wrf,ny_wrf),        stat=istatus)

  ALLOCATE(xa_wrf(nx_wrf,ny_wrf),stat=istatus)
  ALLOCATE(ya_wrf(nx_wrf,ny_wrf),stat=istatus)

  dxfld=0.0
  dyfld=0.0
  rdxfld=0.0
  rdyfld=0.0
  dxfldu=0.0
  dyfldu=0.0
  rdxfldu=0.0
  rdyfldu=0.0
  dxfldv=0.0
  dyfldv=0.0
  rdxfldv=0.0
  rdyfldv=0.0

  tem1_wrf=0.0
  tem2_wrf=0.0
  tem3_wrf=0.0
  tem4_wrf=0.0
  tem5_wrf=0.0

  xa_wrf=0.0
  ya_wrf=0.0
!
!-----------------------------------------------------------------------
!
!  Loop through the data times provided via NAMELIST.
!
!-----------------------------------------------------------------------
!
  iniotfu = 21  ! FORTRAN unit number used for data output

  first_time = .TRUE.
  wait_verif = .TRUE.
  in_verif   = .FALSE.

  readstride = readstride/(ncompressx*ncompressy)
  IF (readstride < 1) THEN
    WRITE(6,*) 'ERROR: readstride < 1, please check max_fopen in namelist file.'
    WRITE(6,*) '       Please remember that readstride = max_fopen/(ncompressx*ncompressy)'
    CALL arpsstop('max_fopen too small',1)
  END IF

  exit_early = .FALSE.
  DO ifile = 1,nextdfil

    IF (myproc == 0) WRITE(6,'(2x,2a,/)') 'Reading file ',TRIM(extdname(ifile))

    IF (frames_per_outfile(ifile) == 1) THEN  ! Finish read in the following block
      !
      !  blocking inserted for ordering i/o for message passing
      !
      DO i=0,nprocs-1,readstride
        IF(myproc >= i .AND. myproc <= i+readstride-1) THEN

          CALL open_wrf_file(TRIM(extdname(ifile)),io_form,multifile,   &
                             .FALSE.,ncompressx,ncompressy,fHndl)

!-----------------------------------------------------------------------
!
!  Retrieve global attributes from the file
!
!-----------------------------------------------------------------------

          CALL get_wrf_metadata(fHndl,io_form,multifile,.FALSE.,        &
                          ncompressx,ncompressy,                        &
                          nx_wrf,ny_wrf,nz_wrf,nzsoil_wrf,              &
                          iproj_wrf,trlat1_wrf,trlat2_wrf,trlon_wrf,    &
                          ctrlat_wrf,ctrlon_wrf,dx_wrf,dy_wrf,dt_wrf,   &
                          sf_surface_physics,sf_sfclay_physics,         &
                          mp_physics,num_scalar,istatus)

          scale_wrf = 1.0
          latnot_wrf(1) = trlat1_wrf
          latnot_wrf(2) = trlat2_wrf

          IF (mp_opt > 0) THEN
            nx_wrf = (nx_wrf - 1)/nproc_x + 1
            ny_wrf = (ny_wrf - 1)/nproc_y + 1
          END IF

          itime = 1

          CALL get_wrf_Times(fHndl,io_form,multifile,ncompressx,ncompressy, &
                             itime,timestr)
!
!-----------------------------------------------------------------------
!
!  Call getwrfd to reads and converts data to ARPS units
!
!  NOTE: u_wrf, v_wrf are just values extracted from data files. It may
!        need to be rotated or extend to be MPI valid.
!
!-----------------------------------------------------------------------
!
          CALL getwrfd(fHndl,io_form,multifile,ncompressx,ncompressy,itime, &
                   timestr,nx_wrf,ny_wrf,nz_wrf,nzsoil_wrf,             &
                   iproj_wrf,scale_wrf,trlon_wrf,latnot_wrf,            &
                   ctrlat_wrf,ctrlon_wrf,dx_wrf,dy_wrf,x0_wrf,y0_wrf,   &
                   sf_surface_physics,sf_sfclay_physics,num_scalar,     &
                   lat_wrf,lon_wrf,latu_wrf,lonu_wrf,latv_wrf,lonv_wrf, &
                   zp_wrf,hgt_wrf,zpsoil_wrf, p_wrf,t_wrf,              &
                   u_wrf,v_wrf, w_wrf,                                  &
                   qv_wrf,qscalar_wrf,tem1_wrf,                         &
                   tsoil_wrf(:,:,:,0),qsoil_wrf(:,:,:,0),               &
                   wetcanp_wrf(:,:,0),snowdpth_wrf,trn_wrf,             &
                   soiltyp_wrf(:,:,1),vegtyp_wrf,veg_wrf,               &
                   tem1_wrf,tem2_wrf,tem3_wrf,tem4_wrf,istatus)

          CALL get_wrf_data2d(fHndl,io_form,                            &
                   multifile,ncompressx,ncompressy,itime,timestr,       &
                   tem1_wrf,istatus)

          CALL close_wrf_file(fHndl,io_form,multifile,.FALSE.,          &
                                ncompressx,ncompressy)

        END IF
        IF (mp_opt > 0) CALL mpbarrier
      END DO

    ELSE     ! frames_per_outfile > 1, only open the file and read meta

      CALL open_wrf_file(TRIM(extdname(ifile)),io_form,multifile,      &
                         .FALSE.,ncompressx,ncompressy,fHndl)

!-----------------------------------------------------------------------
!
!  Retrieve global attributes from the file
!
!-----------------------------------------------------------------------

      CALL get_wrf_metadata(fHndl,io_form,multifile,.FALSE.,            &
                          ncompressx,ncompressy,                        &
                          nx_wrf,ny_wrf,nz_wrf,nzsoil_wrf,              &
                          iproj_wrf,trlat1_wrf,trlat2_wrf,trlon_wrf,    &
                          ctrlat_wrf,ctrlon_wrf,dx_wrf,dy_wrf,dt_wrf,   &
                          sf_surface_physics,sf_sfclay_physics,         &
                          mp_physics,num_scalar,istatus)

      scale_wrf = 1.0
      latnot_wrf(1) = trlat1_wrf
      latnot_wrf(2) = trlat2_wrf

      IF (mp_opt > 0) THEN
        nx_wrf = (nx_wrf - 1)/nproc_x + 1
        ny_wrf = (ny_wrf - 1)/nproc_y + 1
      END IF

    END IF

    DO itime = 1, frames_per_outfile(ifile)

      IF (frames_per_outfile(ifile) > 1) THEN   ! read data only when frames_per_outfile >1
                                         ! otherwise, it has been read before.
        CALL get_wrf_Times(fHndl,io_form,multifile,ncompressx,ncompressy,&
                           itime,timestr)
!
!-----------------------------------------------------------------------
!
!  Call getwrfds to reads and converts data to ARPS units
!
!  NOTE: u_wrf, v_wrf are just values extracted from data files. It may
!        need to be rotated or extend to be MPI valid.
!
!-----------------------------------------------------------------------
!
        CALL getwrfd(fHndl,io_form,multifile,ncompressx,ncompressy,itime,&
                   timestr,nx_wrf,ny_wrf,nz_wrf,nzsoil_wrf,             &
                   iproj_wrf,scale_wrf,trlon_wrf,latnot_wrf,            &
                   ctrlat_wrf,ctrlon_wrf,dx_wrf,dy_wrf,x0_wrf,y0_wrf,   &
                   sf_surface_physics,sf_sfclay_physics,num_scalar,     &
                   lat_wrf,lon_wrf,latu_wrf,lonu_wrf,latv_wrf,lonv_wrf, &
                   zp_wrf,hgt_wrf,zpsoil_wrf, p_wrf,t_wrf,              &
                   u_wrf,v_wrf, w_wrf,                                  &
                   qv_wrf,qscalar_wrf,tem1_wrf,                         &
                   tsoil_wrf(:,:,:,0),qsoil_wrf(:,:,:,0),               &
                   wetcanp_wrf(:,:,0),snowdpth_wrf,trn_wrf,             &
                   soiltyp_wrf(:,:,1),vegtyp_wrf,veg_wrf,               &
                   tem1_wrf,tem2_wrf,tem3_wrf,tem4_wrf,istatus)


        CALL get_wrf_data2d(fHndl,io_form,                              &
                 multifile,ncompressx,ncompressy,itime,timestr,         &
                 tem1_wrf,istatus)

      END IF


!-----------------------------------------------------------------------
!
! Post-reading processing, most of this block is inside getwrfd before
! It must be moved here because of ordering I/O for message passing
!
!-----------------------------------------------------------------------

      CALL adj_wrfuv(multifile,1,nx_wrf,ny_wrf,nz_wrf,                  &
                iproj_wrf,scale_wrf,trlon_wrf,latnot_wrf,x0_wrf,y0_wrf, &
                lonu_wrf,lonv_wrf,u_wrf,vatu_wrf,uatv_wrf,v_wrf,        &
                tem1_wrf,tem2_wrf,istatus)

      pt_wrf(:,:,:) = t_wrf(:,:,:)*((p0/p_wrf(:,:,:))**rddcp)

      stypfrct_wrf(:,:,1) = 1.
      tsoil_wrf(:,:,:,1) = tsoil_wrf(:,:,:,0)
      qsoil_wrf(:,:,:,1) = qsoil_wrf(:,:,:,0)
      wetcanp_wrf(:,:,1) = wetcanp_wrf(:,:,0)
!
!-----------------------------------------------------------------------
!
!  Time conversions.
!  Formats:  timestr='1998-05-25_18:00:00
!
!-----------------------------------------------------------------------
!
      jtime=jtime+1
      READ(timestr,'(I4.4,a1,I2.2,a1,I2.2,a1,I2.2,a1,I2.2,a1,I2.2)')    &
                 iyr,ach,imo,ach,iday,ach,ihr,ach,imin,ach,isec
      WRITE(timefil(jtime),'(i4.4,i2.2,i2.2,a1,3(i2.2))') &
                 iyr,imo,iday,'-',ihr,imin,isec
      WRITE(timewrf(jtime),'(i4.4,i2.2,i2.2,1x,i2.2,2(a1,i2.2))') &
                 iyr,imo,iday,ihr,':',imin,':',isec
      CALL julday(iyr,imo,iday,jldy)
      CALL ctim2abss(iyr,imo,iday,ihr,imin,isec,iabssec)
      IF(wait_verif) THEN
        IF(ihr == 6) THEN
          in_verif=.true.
          wait_verif=.false.
          CALL ctim2abss(iyr,imo,iday,6,0,0,istverif)
        END IF
      ELSE IF (in_verif) THEN
        IF((iabssec - istverif) > isec24h) in_verif=.false.
      END IF

      kftime=iabssec - initsec
      curtim=FLOAT(kftime)

      IF(myproc == 0) WRITE(6,'(/,a,a19,/,a,i12,a,i12,a,/)')            &
          ' Subroutine getwrfd was called for data valid at ',timestr,  &
          ' Which is ',iabssec,' abs seconds or ',kftime,               &
          ' seconds from the WRF/ARPS initial time.'

!-----------------------------------------------------------------------

      CALL a3dmax0(lat_wrf,1,nx_wrf,1,nx_wrf,1,ny_wrf,1,ny_wrf,         &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'lat_wrf_min= ', amin,', lat_wrf_max=',amax
      CALL a3dmax0(lon_wrf,1,nx_wrf,1,nx_wrf,1,ny_wrf,1,ny_wrf,         &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'lon_wrf_min= ', amin,', lon_wrf_max=',amax

      CALL a3dmax0(p_wrf  ,1,nx_wrf,1,nx_wrf,1,ny_wrf,1,ny_wrf,         &
                   1,nz_wrf,1,nz_wrf,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'p_wrf_min  = ', amin,', p_wrf_max  =',amax
      CALL a3dmax0(hgt_wrf,1,nx_wrf,1,nx_wrf,1,ny_wrf,1,ny_wrf,         &
                   1,nz_wrf,1,nz_wrf,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'hgt_wrf_min= ', amin,', hgt_wrf_max=',amax
      CALL a3dmax0(t_wrf  ,1,nx_wrf,1,nx_wrf,1,ny_wrf,1,ny_wrf,         &
                   1,nz_wrf,1,nz_wrf,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   't_wrf_min  = ', amin,', t_wrf_max  =',amax
      CALL a3dmax0(u_wrf  ,1,nx_wrf,1,nx_wrf,1,ny_wrf,1,ny_wrf,         &
                   1,nz_wrf,1,nz_wrf,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'u_wrf_min  = ', amin,', u_wrf_max  =',amax
      CALL a3dmax0(v_wrf  ,1,nx_wrf,1,nx_wrf,1,ny_wrf,1,ny_wrf,         &
                   1,nz_wrf,1,nz_wrf,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'v_wrf_min  = ', amin,', v_wrf_max  =',amax
      CALL a3dmax0(w_wrf  ,1,nx_wrf,1,nx_wrf,1,ny_wrf,1,ny_wrf,         &
                   1,nz_wrf,1,nz_wrf,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'w_wrf_min  = ', amin,', w_wrf_max  =',amax
      CALL a3dmax0(qv_wrf ,1,nx_wrf,1,nx_wrf,1,ny_wrf,1,ny_wrf,         &
                   1,nz_wrf,1,nz_wrf,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'qv_wrf_min = ', amin,', qv_wrf_max =',amax
      DO nq=1,nscalar
        CALL a3dmax0(qscalar_wrf(1,1,1,nq),1,nx_wrf,1,nx_wrf,1,ny_wrf,1,ny_wrf,   &
                   1,nz_wrf,1,nz_wrf,amax,amin)
        IF (myproc == 0) WRITE(6,'(1x,a,e13.6,2a,e13.6)')               &
              TRIM(qnames(nq))//'_ext_min = ', amin,', ',TRIM(qnames(nq))//'_ext_max =',amax
      END DO


      CALL a3dmax0(tsoil_wrf(1,1,1,0),1,nx_wrf,1,nx_wrf,1,ny_wrf,1,ny_wrf, &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'tsoil_sfc_wrf_min= ', amin,', tsoil_sfc_wrf_max=',amax
      CALL a3dmax0(tsoil_wrf(1,1,2,0),1,nx_wrf,1,nx_wrf,1,ny_wrf,1,ny_wrf, &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'tsoil_1st_wrf_min= ', amin,', tsoil_1st_wrf_max=',amax
      CALL a3dmax0(qsoil_wrf(1,1,1,0),1,nx_wrf,1,nx_wrf,1,ny_wrf,1,ny_wrf, &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'qsoil_sfc_wrf_min= ', amin,', qsoil_sfc_wrf_max=',amax
      CALL a3dmax0(qsoil_wrf(1,1,2,0),1,nx_wrf,1,nx_wrf,1,ny_wrf,1,ny_wrf, &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'qsoil_1st_wrf_min= ', amin,', qsoil_1st_wrf_max=',amax
      CALL a3dmax0(wetcanp_wrf,1,nx_wrf,1,nx_wrf,1,ny_wrf,1,ny_wrf,     &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'wetcanp_wrf_min= ', amin,', wetcanp_wrf_max=',amax
      CALL a3dmax0(snowdpth_wrf,1,nx_wrf,1,nx_wrf,1,ny_wrf,1,ny_wrf,    &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'snowd_wrf_min= ', amin,', snow_wrf_max=',amax

      CALL a3dmax0(trn_wrf    ,1,nx_wrf,1,nx_wrf,1,ny_wrf,1,ny_wrf,     &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'trn_wrf_min= ', amin,', trn_wrf_max=',amax

      CALL a3dmax0(veg_wrf,1,nx_wrf,1,nx_wrf,1,ny_wrf,1,ny_wrf,         &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'veg_wrf_min= ', amin,', veg_wrf_max=',amax
      IF(myproc == 0) PRINT*,' '
!
!-----------------------------------------------------------------------
!
!    First time through the time loop, calculate grid
!    transformation info.
!
!-----------------------------------------------------------------------
!
      IF(first_time) THEN
!
!-----------------------------------------------------------------------
!
!  Find x,y locations of WRF grid.
!
!-----------------------------------------------------------------------
!
        CALL setmapr(iproj_wrf,scale_wrf,latnot_wrf,trlon_wrf)
        CALL setorig(1,x0_wrf,y0_wrf)
        DO j=1,ny_wrf
          CALL lltoxy(1,1,lat_wrf(1,j),lon_wrf(1,j),                    &
                      x_wrf(1),y_wrf(j))
        END DO
        DO i=1,nx_wrf
          CALL lltoxy(1,1,lat_wrf(i,1),lon_wrf(i,1),                    &
                      x_wrf(i),y_wrf(1))
        END DO
        DO j=1,ny_wrf
          CALL lltoxy(1,1,latu_wrf(1,j),lonu_wrf(1,j),                  &
                      xu_wrf(1),yu_wrf(j))
        END DO
        DO i=1,nx_wrf
          CALL lltoxy(1,1,latu_wrf(i,1),lonu_wrf(i,1),                  &
                          xu_wrf(i),yu_wrf(1))
        END DO
        DO j=1,ny_wrf
          CALL lltoxy(1,1,latv_wrf(1,j),lonv_wrf(1,j),                  &
                          xv_wrf(1),yv_wrf(j))
        END DO
        DO i=1,nx_wrf
          CALL lltoxy(1,1,latv_wrf(i,1),lonv_wrf(i,1),                  &
                          xv_wrf(i),yv_wrf(1))
        END DO

!
!-----------------------------------------------------------------------
!
!  Find x,y of sounding locations in external grid.
!
!-----------------------------------------------------------------------
!
        DO isnd=1,nsnd
          IF(locopt == 1) THEN
            CALL lltoxy(1,1,slat(isnd),slon(isnd),xpt(isnd),ypt(isnd))
          END IF

          IF (xpt(isnd) > x_wrf(1) .AND. xpt(isnd) < x_wrf(nx_wrf) .AND.    &
              ypt(isnd) > y_wrf(1) .AND. ypt(isnd) < y_wrf(ny_wrf) ) THEN

            ptingrid(isnd)=.TRUE.
            dmin=((xpt(isnd)-x_wrf(1))*(xpt(isnd)-x_wrf(1))+                &
                  (ypt(isnd)-y_wrf(1))*(ypt(isnd)-y_wrf(1)))

            ipt(isnd)=1
            jpt(isnd)=1

            DO j=1,ny_wrf-1
              DO i=1,nx_wrf-1
                dd=((xpt(isnd)-x_wrf(i))*(xpt(isnd)-x_wrf(i))+              &
                    (ypt(isnd)-y_wrf(j))*(ypt(isnd)-y_wrf(j)))
                IF(dd < dmin) THEN
                  dmin=dd
                  ipt(isnd)=i
                  jpt(isnd)=j
                END IF
              END DO
            END DO
            WRITE(6,'(a,f10.4,f10.4,/a,i5,i5,a,f10.4,f10.4)')             &
              ' Nearest WRF pt to diagnostic lat,lon: ',                  &
              slat(isnd),slon(isnd),                                      &
              ' Diagnostic i,j: ', ipt(isnd),jpt(isnd),                   &
              ' lat,lon= ', lat_wrf(ipt(isnd),jpt(isnd)),                 &
                            lon_wrf(ipt(isnd),jpt(isnd))
            WRITE(6,'(1x,a,f10.2,a)') 'Terrain at WRF pt: ',              &
                              trn_wrf(ipt(isnd),jpt(isnd)),' meters.'

          ELSE
            ptingrid(isnd)=.FALSE.
            ipt(isnd)=1
            jpt(isnd)=1
!            WRITE(6,'(1x,a,I2,a)') 'ERROR: Station ',isnd,              &
!                    ' location is outside of the WRF domain.'
!            CALL arpsstop('Station location exceeds the WRF domain',1)
          END IF
        END DO ! isnd


        first_time = .FALSE.

      END IF ! first_time

      DO isnd=1,nsnd
      IF(ptingrid(isnd)) THEN

        iloc=ipt(isnd)
        jloc=jpt(isnd)

        WRITE(6,'(///2(a,I4),/2x,a)')                                    &
          ' External sounding at iloc,jloc = ',iloc,', ',jloc,           &
          'k   pres    hgt   temp   theta   dewp     u      v     dir    spd'
!
!-----------------------------------------------------------------------
!
!  Output sounding to look like GEMPAK SNLIST.FIL
!  to use the Skew-T and Hodograph programs
!  e.g., those by Rich Carpenter and NSHARP
!
!  Example of what the header looks like:
!
!23456789012345678901234567890123456789012345678901234567890
!SNPARM = PRES;HGHT;TMPC;DWPC;DRCT;SPED;OMEG
!
!STID = SEP        STNM =    72260   TIME = 920308/1500
!SLAT =    32.21   SLON =   -98.18   SELEV =   399.0
!
!  PRES     HGHT     TMPC     DWPC     DRCT     SPED     OMEG
!
!-----------------------------------------------------------------------
!
      IF(outsnd > 0) THEN
        WRITE(uafile,'(6a)') TRIM(dir_output),'/',TRIM(stid(isnd)),&
           '-',timefil(jtime),'.gemsnd'
        OPEN(31,FILE=TRIM(uafile),STATUS='unknown')
        modyr=mod(iyr,100)
        WRITE(31,'(/a/)') ' SNPARM = PRES;HGHT;TMPC;DWPC;DRCT;SPED;OMEG'
        WRITE(31,'(a,a8,3x,a,i8,3x,a,i2.2,i2.2,i2.2,a1,i2.2,i2.2)')     &
            ' STID = ',stid(isnd),                                      &
            'STNM = ',istnm(isnd),                                      &
            'TIME = ',modyr,imo,iday,'/',ihr,imin
        WRITE(31,'(a,f8.2,3x,a,f8.2,3x,a,f7.1/)')                       &
            ' SLAT = ',slat(isnd),'SLON = ',slon(isnd),                 &
            'SELV = ',selev(isnd)
        WRITE(31,'(6x,a)')                                              &
            'PRES     HGHT     TMPC     DWPC     DRCT     SPED     OMEG'
      END IF
!
!  Convert units of external data and write as a sounding.
!
        DO k=1,nz_wrf
          pmb=.01*p_wrf(iloc,jloc,k)
          tc(k)=t_wrf(iloc,jloc,k)-273.15
          theta=pt_wrf(iloc,jloc,k)
          IF( qv_wrf(iloc,jloc,k) > 0.) THEN
            smix=qv_wrf(iloc,jloc,k)/(1.-qv_wrf(iloc,jloc,k))
            e=(pmb*smix)/(0.62197 + smix)
            bige=e/( 1.001 + ( (pmb - 100.) / 900.) * 0.0034)
            alge = ALOG(bige/6.112)
            tdc(k) = (alge * 243.5) / (17.67 - alge)
          ELSE
            tdc(k) = tc(k)-30.
          END IF

          usnd=0.5*(u_wrf(iloc,jloc,k)+ &
                    u_wrf((iloc+1),jloc,k))
          vsnd=0.5*(v_wrf(iloc,jloc,k)+ &
                    v_wrf(iloc,(jloc+1),k))

          CALL uvrotdd(1,1,lon_wrf(iloc,jloc),                &
                       usnd(k),vsnd(k),dir(k),spd(k))

          WRITE(6,'(i4,f6.0,f9.0,f7.1,f7.1,f7.1,f7.1,f7.1,f7.1,f7.1)')  &
                        k,pmb,                                          &
                        hgt_wrf(iloc,jloc,k),                           &
                        tc(k),theta,tdc(k),                             &
                        usnd(k),vsnd(k),dir(k),spd(k)
          IF(outsnd > 0 ) WRITE(31,'(1x,7(F9.2))')                    &
                 pmb,hgt_wrf(iloc,jloc,k),                            &
                 tc(k),tdc(k),dir(k),spd(k),omegasnd
        END DO
        IF(outsnd > 0) CLOSE(31)
      ELSE  ! pt not in grid
!
!     Create files with an error message
!
        IF(outsnd > 0) THEN
          WRITE(uafile,'(6a)') TRIM(dir_output),'/',TRIM(stid(isnd)),&
           '-',timefil(jtime),'.gemsnd'
          OPEN(31,FILE=TRIM(uafile),STATUS='unknown')
          WRITE(31,'(3a)') ' Station ',stid(isnd),' not in model domain'
          CLOSE(31)
        END IF
      END IF
!
      scrhgt=selev(isnd)+2.0
      annhgt=selev(isnd)+10.0
      tcsoil=tsoil_wrf(iloc,jloc,1,0)-273.15

      IF(scrhgt > trn_wrf(iloc,jloc)) THEN
        DO k=1,nz_wrf
          IF(hgt_wrf(iloc,jloc,k) > scrhgt) THEN
            ktop=k
            EXIT
          END IF
        END DO
        IF(k > 1) THEN
          kbot=k-1
          delt=tc(ktop)-tc(kbot)
          delz=hgt_wrf(iloc,jloc,ktop)-hgt_wrf(iloc,jloc,kbot)
          tsfc=tc(kbot)+((scrhgt-hgt_wrf(iloc,jloc,kbot))*(delt/delz))
        ELSE
          delt=tc(1)-tcsoil
          delz=hgt_wrf(iloc,jloc,1)-trn_wrf(iloc,jloc)
          tsfc=tcsoil+((scrhgt-trn_wrf(iloc,jloc))*(delt/delz))
        END IF
      ELSE ! below terrain
        tsfc=tc(1)+0.0065*(scrhgt-hgt_wrf(iloc,jloc,1))
      END IF

      print *, ' terrain, scrhgt, hgt(1):', &
                 trn_wrf(iloc,jloc),scrhgt,hgt_wrf(iloc,jloc,1)

      print *, ' tcsoil,tsfc,tc(1):',tcsoil,tsfc,tc(1)

      tf=tctof(tsfc)
      tdf=tctof(tdc(1))
      spdkts=spd(1)*ms2kts
      rain_g=raing_ext(iloc,jloc)*mm2inch
      rain_c=rainc_ext(iloc,jloc)*mm2inch
      raint=rain_g+rain_c
!      WRITE(isfcunit(isnd),'(1x,a8,1x,i4.4,2i2.2,1x,i2.2,2(a,i2.2),4f6.0,3f8.3)') &
!            stid(isnd),iyr,imo,iday,ihr,':',imin,':',isec, &
!            tf,tdf,dir(1),spdkts,raing,rainc,raint

      pmb=0.01*psfc_ext(iloc,jloc)
      tf=tktof(t2m_ext(iloc,jloc))
      IF( qv2m_ext(iloc,jloc) > 0.) THEN
         smix=qv2m_ext(iloc,jloc)/(1.-qv2m_ext(iloc,jloc))
         e=(pmb*smix)/(0.62197 + smix)
         bige=e/( 1.001 + ( (pmb - 100.) / 900.) * 0.0034)
         alge = ALOG(bige/6.112)
         tdsfc = (alge * 243.5) / (17.67 - alge)
      ELSE
        tdsfc = t2m_ext(iloc,jloc)-303.15
      END IF
      tdf=tctof(tdsfc)
      CALL uvrotdd(1,1,lon_wrf(iloc,jloc),                &
               u10m_ext(iloc,jloc),v10m_ext(iloc,jloc),dirsfc,spdsfc)
      spdkts=spdsfc*ms2kts

      sfcpr(jtime,isnd)=pmb
      tsfcf(jtime,isnd)=tf
      tdsfcf(jtime,isnd)=tdf
      wdir(jtime,isnd)=dirsfc
      wspd(jtime,isnd)=spdkts
      prcgp(jtime,isnd)=rain_g
      prccu(jtime,isnd)=rain_c
      prctot(jtime,isnd)=raint

        IF(in_verif) THEN
          tmax(isnd)=max(tmax(isnd),tf)
          tmin(isnd)=min(tmin(isnd),tf)
          spdmax(isnd)=max(spdmax(isnd),spdkts)
          IF(prc0(isnd) < 0.) prc0(isnd)=raint
          prcsum(isnd)=raint-prc0(isnd)
        END IF

      END DO    ! isnd

      IF (iabssec >= abstimee) THEN       ! exit if time exceeds required in namelist
        exit_early = .TRUE.
        EXIT
      END IF

    END DO     ! itime

    IF (frames_per_outfile(ifile) > 1)  THEN
      CALL close_wrf_file(fHndl,io_form,multifile,.FALSE.,              &
                          ncompressx,ncompressy)
    END IF

    IF (exit_early) EXIT

  END DO  ! ifile
  ntime=jtime

  IF(outsfc > 0) THEN
    DO isnd=1, nsnd
      WRITE(sfcfile,'(6a)') TRIM(dir_output),'/',  &
          TRIM(stid(isnd)),'-',timefil(1),'.sfc'
      OPEN(32,FILE=TRIM(sfcfile),STATUS='unknown')
      IF(ptingrid(isnd)) THEN
      WRITE(32,'(2a,2(a,f7.0))')  &
          ' Station: ',stid(isnd),'   Elev(m): ',selev(isnd),  &
          '    WRF Terrain(m): ', trn_wrf(ipt(isnd),jpt(isnd))
      WRITE(32,'(a,a)') '   Date  Time(UTC) SfcPr  Tsfc Tdsfc', &
          '  WDir  WSpd  PrecGP  PrecCu PrecTot'
      DO jtime=1,ntime
        WRITE(32,'(1x,a,5f6.0,3f8.3)') timewrf(jtime), &
            sfcpr(jtime,isnd),tsfcf(jtime,isnd),tdsfcf(jtime,isnd), &
            wdir(jtime,isnd),wspd(jtime,isnd),prcgp(jtime,isnd), &
            prccu(jtime,isnd),prctot(jtime,isnd)
      END DO
      WRITE(32,'(a,f6.0,a,f6.0,a,f6.0,a,f8.3)') &
          ' Tmin=',tmin(isnd),'  Tmax=',tmax(isnd),&
          '  WSpd Max=',spdmax(isnd),' Precip=',prcsum(isnd)
      ELSE
      WRITE(6,'(3a)')  ' Station: ',stid(isnd),' is outside WRF domain'
      WRITE(32,'(3a)') ' Station: ',stid(isnd),' is outside WRF domain'
      END IF
      CLOSE(32)
    END DO
  END IF

  CALL io_shutdown(io_form)

!
!-----------------------------------------------------------------------
!
!  Friendly exit message.
!
!-----------------------------------------------------------------------
!
  IF(myproc == 0) WRITE(6,'(a/a,i4,a)')                                 &
          ' ==== Normal succesful completion of WRFEXTSND ====',        &
          '      Processed',nextdfil,' file(s)'

  IF (mp_opt > 0) CALL mpexit(0)

  STOP
!
!-----------------------------------------------------------------------
!
!  Error status returned from getwrfd
!
!-----------------------------------------------------------------------
!
  IF(myproc == 0) WRITE(6,'(a,i6)')                                     &
          ' Aborting, error reading external file. istatus=',istatus

  IF (mp_opt > 0) CALL mpexit(0)

  STOP
END PROGRAM wrfextsnd

FUNCTION tctof(tc)
  IMPLICIT NONE
  REAL :: tctof
  REAL :: tc
  REAL, PARAMETER :: fconst = 32.0
  REAL, PARAMETER :: fcratio = (9./5.)

  tctof = (tc*fcratio) + fconst

  RETURN
END

FUNCTION tktof(tk)
  IMPLICIT NONE
  REAL :: tktof
  REAL :: tk
  REAL, PARAMETER :: fconst = 32.0
  REAL, PARAMETER :: fcratio = (9./5.)

  tktof = ((tk-273.15)*fcratio) + fconst

  RETURN
END

