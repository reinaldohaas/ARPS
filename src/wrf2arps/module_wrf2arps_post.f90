MODULE module_wrf2arps_post
  IMPLICIT NONE

  !-------------- on ARPS grid -----------------------------------------
  REAL,    POINTER     :: t2m   (:,:)
  REAL,    POINTER     :: th2m  (:,:)
  REAL,    POINTER     :: qv2m  (:,:)
  REAL,    POINTER     :: u10m  (:,:)
  REAL,    POINTER     :: v10m  (:,:)

  REAL,    POINTER     :: raing(:,:)
  REAL,    POINTER     :: rainc(:,:)
  REAL,    POINTER     :: rainsnow(:,:)
  REAL,    POINTER     :: raingrpl(:,:)
  REAL,    POINTER     :: rainhail(:,:)

  REAL,    POINTER     :: raddn (:,:)     ! sfc total downward radiation flux (W/m*2)

  REAL,    POINTER     :: wspd10max  (:,:)
  REAL,    POINTER     :: w_up_max   (:,:)
  REAL,    POINTER     :: w_dn_max   (:,:)
  REAL,    POINTER     :: refd_max   (:,:)
  REAL,    POINTER     :: up_heli_max(:,:)
  REAL,    POINTER     :: grpl_max   (:,:)
  REAL,    POINTER     :: ltg1_max   (:,:)
  REAL,    POINTER     :: ltg2_max   (:,:)
  REAL,    POINTER     :: ltg3_max   (:,:)
  REAL,    POINTER     :: wspd1kmmax (:,:)
  REAL,    POINTER     :: crefd_max  (:,:)
  REAL,    POINTER     :: refdm10c_max(:,:)
  REAL,    POINTER     :: pblh       (:,:)
  REAL,    POINTER     :: up_heli16_max(:,:)

  REAL,    POINTER     :: snow_max   (:,:)
  REAL,    POINTER     :: grpl05_max (:,:)
  REAL,    POINTER     :: up_heli_maxe(:,:)
  REAL,    POINTER     :: up_heli_maxp(:,:)

  REAL, ALLOCATABLE, TARGET :: post2d(:,:)

  !-------------------- On External grid -------------------------------

  REAL,    POINTER     :: psfc_ext(:,:)
  REAL,    POINTER     :: t2m_ext (:,:)
  REAL,    POINTER     :: th2m_ext(:,:)
  REAL,    POINTER     :: qv2m_ext(:,:)
  REAL,    POINTER     :: u10m_ext(:,:)
  REAL,    POINTER     :: v10m_ext(:,:)

  REAL,    POINTER     :: raing_ext   (:,:)
                          ! PBL time-step total precipitation (mm)
                          ! whether is it "raing" in ARPS????
  REAL,    POINTER     :: rainc_ext   (:,:)
                          ! time-step cumulus precipitation (mm)
  REAL,    POINTER     :: raingrpl_ext(:,:)
  REAL,    POINTER     :: rainsnow_ext(:,:)
  REAL,    POINTER     :: rainhail_ext(:,:)

  REAL,    POINTER     :: raddn_ext (:,:)

  REAL,    POINTER     :: wspd10max_ext  (:,:)
  REAL,    POINTER     :: w_up_max_ext   (:,:)
  REAL,    POINTER     :: w_dn_max_ext   (:,:)
  REAL,    POINTER     :: refd_max_ext   (:,:)
  REAL,    POINTER     :: up_heli_max_ext(:,:)
  REAL,    POINTER     :: grpl_max_ext   (:,:)
  REAL,    POINTER     :: ltg1_max_ext   (:,:)
  REAL,    POINTER     :: ltg2_max_ext   (:,:)
  REAL,    POINTER     :: ltg3_max_ext   (:,:)
  REAL,    POINTER     :: wspd1kmmax_ext (:,:)
  REAL,    POINTER     :: crefd_max_ext  (:,:)
  REAL,    POINTER     :: refdm10c_max_ext(:,:)
  REAL,    POINTER     :: pblh_ext       (:,:)
  REAL,    POINTER     :: up_heli16_max_ext(:,:)

  REAL,    POINTER     :: snow_max_ext   (:,:)
  REAL,    POINTER     :: grpl05_max_ext (:,:)
  REAL,    POINTER     :: up_heli_maxe_ext(:,:)
  REAL,    POINTER     :: up_heli_maxp_ext(:,:)

  REAL, ALLOCATABLE, TARGET :: post2d_ext(:,:)

  !---------------------- Flags and Private variables ------------------

  INTEGER, PRIVATE :: nx, ny, nx_ext, ny_ext

  LOGICAL, PRIVATE :: call_postcore, SPRING_EXPERIMENT
  INTEGER, PRIVATE :: iaccu, iabssec1
  INTEGER, PRIVATE :: icape, iascii, i2dfmt, igempak

  CHARACTER(LEN=256), PRIVATE :: outheader,gemoutheader

  INTEGER, PRIVATE :: ilite,iltgci,icrtm,isatid,chbgn,chend,icitm,user_emis
  INTEGER, PRIVATE :: ibeg_offset,iend_offset,jbeg_offset,jend_offset

  !---------------------- PRIVATE SUBPROCEDURE ------------------

  PRIVATE :: get_wrf_rain, print_maxmin_rain

CONTAINS

  !#####################################################################

  SUBROUTINE init_data2d(i2dfmtin,igempakin,iaccuin,abstimei,           &
                         header1,header2,icapein,iasciiin,              &
                         ilitein,iltgciin,icrtmin,isatidin,             &
                         chbgnin,chendin,icitmin,user_emisin,           &
                         ibeg_offsetin,iend_offsetin,                   &
                         jbeg_offsetin,jend_offsetin,                   &
                         spring,istatus)

  !---------------------------------------------------------------------

    INTEGER, INTENT(IN)  :: i2dfmtin,igempakin
    INTEGER, INTENT(IN)  :: iaccuin, abstimei
    CHARACTER(LEN=*), INTENT(IN) :: header1, header2
    INTEGER, INTENT(IN)  :: icapein, iasciiin
    INTEGER, INTENT(IN)  :: ilitein,iltgciin,icrtmin,isatidin,          &
                            chbgnin,chendin,icitmin,user_emisin,        &
                            ibeg_offsetin,iend_offsetin,                &
                            jbeg_offsetin,jend_offsetin

    LOGICAL, INTENT(IN)  :: spring

    INTEGER, INTENT(OUT) :: istatus

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    i2dfmt  = i2dfmtin
    igempak = igempakin

    call_postcore = .FALSE.
    IF (i2dfmt > 0 .OR. igempak == 1) THEN
        call_postcore = .TRUE.
    END IF

    iaccu    = iaccuin
    iabssec1 = abstimei
    icape    = icapein
    iascii   = iasciiin

    outheader    = header1
    gemoutheader = header2

    ilite        = ilitein
    iltgci       = iltgciin
    icrtm        = icrtmin
    isatid       = isatidin
    chbgn        = chbgnin
    chend        = chendin
    icitm        = icitmin
    user_emis    = user_emisin
    ibeg_offset  = ibeg_offsetin
    iend_offset  = iend_offsetin
    jbeg_offset  = jbeg_offsetin
    jend_offset  = jend_offsetin

    SPRING_EXPERIMENT = spring

    RETURN
  END SUBROUTINE init_data2d

  !#####################################################################

  SUBROUTINE allocate_data2d_arps(nxin,nyin,istatus)

    INTEGER, INTENT(IN)  :: nxin, nyin
    INTEGER, INTENT(OUT) :: istatus

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    nx = nxin
    ny = nyin

    ALLOCATE( t2m(nx,ny),           STAT = istatus )
    ALLOCATE( th2m(nx,ny),          STAT = istatus )
    ALLOCATE( qv2m(nx,ny),          STAT = istatus )
    ALLOCATE( u10m(nx,ny),          STAT = istatus )
    ALLOCATE( v10m(nx,ny),          STAT = istatus )
    ALLOCATE( raddn(nx,ny),         STAT = istatus )

    ALLOCATE( raing(nx,ny),         STAT = istatus )
    ALLOCATE( rainc(nx,ny),         STAT = istatus )
    ALLOCATE( rainsnow(nx,ny),         STAT = istatus )
    ALLOCATE( raingrpl(nx,ny),         STAT = istatus )
    ALLOCATE( rainhail(nx,ny),         STAT = istatus )

    t2m = 0.0
    th2m = 0.0
    qv2m = 0.0
    u10m = 0.0
    v10m = 0.0

    raddn = 0.0

    raing = 0.0
    rainc = 0.0
    rainsnow = 0.0
    raingrpl = 0.0
    rainhail = 0.0

    IF ( call_postcore ) THEN

      IF ( SPRING_EXPERIMENT ) THEN
        ALLOCATE( wspd10max(nx,ny),     STAT = istatus )
        ALLOCATE( w_up_max(nx,ny),      STAT = istatus )
        ALLOCATE( w_dn_max(nx,ny),      STAT = istatus )
        ALLOCATE( refd_max(nx,ny),      STAT = istatus )
        ALLOCATE( up_heli_max(nx,ny),   STAT = istatus )
        ALLOCATE( grpl_max(nx,ny),      STAT = istatus )
        ALLOCATE( ltg1_max(nx,ny),      STAT = istatus )
        ALLOCATE( ltg2_max(nx,ny),      STAT = istatus )
        ALLOCATE( ltg3_max(nx,ny),      STAT = istatus )
        ALLOCATE( wspd1kmmax(nx,ny),    STAT = istatus )
        ALLOCATE( crefd_max(nx,ny),     STAT = istatus )
        ALLOCATE( refdm10c_max(nx,ny),  STAT = istatus )
        ALLOCATE( pblh(nx,ny),          STAT = istatus )
        ALLOCATE( up_heli16_max(nx,ny), STAT = istatus )

        ALLOCATE( snow_max(nx,ny),      STAT = istatus )
        ALLOCATE( grpl05_max(nx,ny),    STAT = istatus )
        ALLOCATE( up_heli_maxe(nx,ny),  STAT = istatus )
        ALLOCATE( up_heli_maxp(nx,ny),  STAT = istatus )


        wspd10max = 0.0
        w_up_max = 0.0
        w_dn_max = 0.0
        refd_max = 0.0
        up_heli_max = 0.0
        grpl_max = 0.0
        ltg1_max = 0.0
        ltg2_max = 0.0
        ltg3_max = 0.0
        wspd1kmmax = 0.0
        crefd_max = 0.0
        refdm10c_max = 0.0
        pblh = 0.0
        up_heli16_max = 0.0

        snow_max = 0.0
        grpl05_max = 0.0
        up_heli_maxe = 0.0
        up_heli_maxp = 0.0

      ELSE

        ALLOCATE(post2d(nx,ny), STAT = istatus)
        post2d = -999.9

        wspd10max      => post2d
        w_up_max       => post2d
        w_dn_max       => post2d
        refd_max       => post2d
        up_heli_max    => post2d
        grpl_max       => post2d
        ltg1_max       => post2d
        ltg2_max       => post2d
        ltg3_max       => post2d
        wspd1kmmax     => post2d
        crefd_max      => post2d
        refdm10c_max   => post2d
        pblh           => post2d
        up_heli16_max  => post2d

        snow_max       => post2d
        grpl05_max     => post2d
        up_heli_maxe   => post2d
        up_heli_maxp   => post2d
      END IF

    END IF

    CALL check_alloc_status(istatus, "allocate_data2d_arps:post2d")

    RETURN
  END SUBROUTINE allocate_data2d_arps

  !#####################################################################

  SUBROUTINE allocate_data2d_ext(nxin,nyin, istatus)

    INTEGER, INTENT(IN)  :: nxin, nyin
    INTEGER, INTENT(OUT) :: istatus

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    nx_ext = nxin
    ny_ext = nyin

    ALLOCATE(psfc_ext     (nx_ext,ny_ext),stat=istatus)
    ALLOCATE(t2m_ext      (nx_ext,ny_ext),stat=istatus)
    ALLOCATE(th2m_ext     (nx_ext,ny_ext),stat=istatus)
    ALLOCATE(qv2m_ext     (nx_ext,ny_ext),stat=istatus)
    ALLOCATE(u10m_ext     (nx_ext,ny_ext),stat=istatus)
    ALLOCATE(v10m_ext     (nx_ext,ny_ext),stat=istatus)

    ALLOCATE(raing_ext    (nx_ext,ny_ext), STAT = istatus)
    ALLOCATE(rainc_ext    (nx_ext,ny_ext), STAT = istatus)
    ALLOCATE(rainsnow_ext (nx_ext,ny_ext), STAT = istatus)
    ALLOCATE(raingrpl_ext (nx_ext,ny_ext), STAT = istatus)
    ALLOCATE(rainhail_ext (nx_ext,ny_ext), STAT = istatus)

    ALLOCATE(raddn_ext    (nx_ext,ny_ext),stat=istatus)

    psfc_ext  = 0.0

    raddn_ext = 0.0
    t2m_ext   = 0.0
    th2m_ext  = 0.0
    qv2m_ext  = 0.0
    u10m_ext  = 0.0
    v10m_ext  = 0.0

    raing_ext = 0.0
    rainc_ext = 0.0
    rainsnow_ext = 0.0
    raingrpl_ext = 0.0
    rainhail_ext = 0.0

    IF (call_postcore) THEN

      IF ( SPRING_EXPERIMENT ) THEN
        ALLOCATE(wspd10max_ext(nx_ext,ny_ext), stat=istatus)
        ALLOCATE(w_up_max_ext (nx_ext,ny_ext), stat=istatus)
        ALLOCATE(w_dn_max_ext (nx_ext,ny_ext), stat=istatus)
        ALLOCATE(refd_max_ext (nx_ext,ny_ext), stat=istatus)
        ALLOCATE(up_heli_max_ext(nx_ext,ny_ext), stat=istatus)
        ALLOCATE(grpl_max_ext   (nx_ext,ny_ext), stat=istatus)
        ALLOCATE(ltg1_max_ext   (nx_ext,ny_ext), stat=istatus)
        ALLOCATE(ltg2_max_ext   (nx_ext,ny_ext), stat=istatus)
        ALLOCATE(ltg3_max_ext   (nx_ext,ny_ext), stat=istatus)
        ALLOCATE(wspd1kmmax_ext (nx_ext,ny_ext), stat=istatus)
        ALLOCATE(crefd_max_ext  (nx_ext,ny_ext), stat=istatus)
        ALLOCATE(refdm10c_max_ext (nx_ext,ny_ext), stat=istatus)
        ALLOCATE(pblh_ext         (nx_ext,ny_ext), stat=istatus)
        ALLOCATE(up_heli16_max_ext(nx_ext,ny_ext), stat=istatus)

        ALLOCATE(snow_max_ext    (nx_ext,ny_ext), stat=istatus)
        ALLOCATE(grpl05_max_ext  (nx_ext,ny_ext), stat=istatus)
        ALLOCATE(up_heli_maxe_ext(nx_ext,ny_ext), stat=istatus)
        ALLOCATE(up_heli_maxp_ext(nx_ext,ny_ext), stat=istatus)

        wspd10max_ext   = 0.0
        w_up_max_ext    = 0.0
        w_dn_max_ext    = 0.0
        refd_max_ext    = 0.0
        up_heli_max_ext = 0.0
        grpl_max_ext    = 0.0
        ltg1_max_ext    = 0.0
        ltg2_max_ext    = 0.0
        ltg3_max_ext    = 0.0
        wspd1kmmax_ext  = 0.0
        crefd_max_ext   = 0.0
        refdm10c_max_ext  = 0.0
        pblh_ext          = 0.0
        up_heli16_max_ext = 0.0

        snow_max_ext     = 0.0
        grpl05_max_ext   = 0.0
        up_heli_maxe_ext = 0.0
        up_heli_maxp_ext = 0.0
      ELSE
        ALLOCATE(post2d_ext(nx_ext,ny_ext), STAT = istatus)
        post2d_ext = -999.9

        psfc_ext          => post2d_ext
        t2m_ext           => post2d_ext
        th2m_ext          => post2d_ext
        qv2m_ext          => post2d_ext
        u10m_ext          => post2d_ext
        v10m_ext          => post2d_ext
        raddn_ext         => post2d_ext

        wspd10max_ext     => post2d_ext
        w_up_max_ext      => post2d_ext
        w_dn_max_ext      => post2d_ext
        refd_max_ext      => post2d_ext
        up_heli_max_ext   => post2d_ext
        grpl_max_ext      => post2d_ext
        ltg1_max_ext      => post2d_ext
        ltg2_max_ext      => post2d_ext
        ltg3_max_ext      => post2d_ext
        wspd1kmmax_ext    => post2d_ext
        crefd_max_ext     => post2d_ext
        refdm10c_max_ext  => post2d_ext
        pblh_ext          => post2d_ext
        up_heli16_max_ext => post2d_ext

        snow_max_ext      => post2d_ext
        grpl05_max_ext    => post2d_ext
        up_heli_maxe_ext  => post2d_ext
        up_heli_maxp_ext  => post2d_ext
      END IF

    END IF

    CALL check_alloc_status(istatus, "allocate_data2d_ext:post2d_ext")

    RETURN
  END SUBROUTINE allocate_data2d_ext

  !##################################################################
  !
  SUBROUTINE get_wrf_data2d(fHndl,io_form,multifile,ncmprx,ncmpry,itime,timestr, &
                            tem2d1,istatus)

  !------------------------------------------------------------------------
  !
  !  PURPOSE:
  !
  !  Read in WRF 2D data in one time level (itime) from the opened NetCDF
  !  file (fHndl) and convert those data to ARPS units if needed.
  !
  !------------------------------------------------------------------------
  !
  !  AUTHOR:
  !  Yunheng Wang (04/22/2013)
  !
  !  MODIFICATION HISTORY:
  !
  !------------------------------------------------------------------------
    IMPLICIT NONE

  !------------------------------------------------------------------------
  !
  !  Include files
  !
  !------------------------------------------------------------------

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

    REAL,         INTENT(INOUT) :: tem2d1 (nx_ext,  ny_ext)

    INTEGER,      INTENT(OUT) :: istatus

  !------------------------------------------------------------------------
  !
  !  Misc. local variables
  !
  !------------------------------------------------------------------

    INTEGER  :: nxlg_ext, nylg_ext
    INTEGER  :: nxd, nyd           ! dimensions in data file
    INTEGER  :: nxd_stag, nyd_stag

    REAL, ALLOCATABLE :: temlg1(:,:)
    REAL, ALLOCATABLE :: temlg2(:,:)

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !
  !  Beginning of executable code...
  !
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

  !-----------------------------------------------------------------------
  !
  ! Allocate working arrays
  !
  !-----------------------------------------------------------------------

    nxlg_ext = (nx_ext-1)*nproc_x + 1
    nylg_ext = (ny_ext-1)*nproc_y + 1

    IF (mp_opt > 0 .AND. multifile) THEN
      nxd = (nx_ext - 1)/ncmprx
      nyd = (ny_ext - 1)/ncmpry
      nxd_stag = nxd
      nyd_stag = nyd
      IF (loc_x == nproc_x) nxd_stag = nxd + 1
      IF (loc_y == nproc_y) nyd_stag = nyd + 1

      ALLOCATE(temlg2(1,1),   STAT = istatus)  ! do not use it
    ELSE
      nxd = nxlg_ext - 1
      nyd = nylg_ext - 1
      nxd_stag = nxlg_ext
      nyd_stag = nylg_ext

      ALLOCATE(temlg2(nxlg_ext, nylg_ext ),  STAT = istatus)
    END IF

    ALLOCATE(  temlg1(nxd_stag, nyd_stag ),  STAT = istatus)

  !---------------------------------------------------------------------
  !
  ! Get Near-Sfc variables
  !
  !---------------------------------------------------------------------


    CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                    timestr,itime,1,'Q2','',                              &
                    'west_east','south_north',nx_ext,ny_ext,qv2m_ext(:,:),&
                    nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

    CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                    timestr,itime,1,'T2','',                              &
                    'west_east','south_north',nx_ext,ny_ext,t2m_ext(:,:), &
                    nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

    CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                    timestr,itime,1,'TH2','',                             &
                    'west_east','south_north',nx_ext,ny_ext,th2m_ext(:,:),&
                    nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

    CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                    timestr,itime,1,'PSFC','',                            &
                    'west_east','south_north',nx_ext,ny_ext,psfc_ext(:,:),&
                    nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

    CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                    timestr,itime,1,'U10','',                             &
                    'west_east','south_north',nx_ext,ny_ext,u10m_ext(:,:),&
                    nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

    CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                    timestr,itime,1,'V10','',                             &
                    'west_east','south_north',nx_ext,ny_ext,v10m_ext(:,:),&
                    nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

    IF ( call_postcore .AND. SPRING_EXPERIMENT ) THEN

        ! SPRING 2010 SPECIFIC fields
        CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                        timestr,itime,1,'WSPD10MAX','',                       &
                        'west_east','south_north',nx_ext,ny_ext,wspd10max_ext(:,:),&
                        nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

        CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                        timestr,itime,1,'W_UP_MAX','',                        &
                        'west_east','south_north',nx_ext,ny_ext,w_up_max_ext(:,:),&
                        nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

        CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                        timestr,itime,1,'W_DN_MAX','',                        &
                        'west_east','south_north',nx_ext,ny_ext,w_dn_max_ext(:,:),&
                        nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

        CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                        timestr,itime,1,'REFD_MAX','',                        &
                        'west_east','south_north',nx_ext,ny_ext,refd_max_ext(:,:),&
                        nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

        CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                        timestr,itime,1,'UP_HELI_MAX','',                     &
                        'west_east','south_north',nx_ext,ny_ext,up_heli_max_ext(:,:),&
                        nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

        CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                        timestr,itime,1,'GRPL_MAX','',                        &
                        'west_east','south_north',nx_ext,ny_ext,grpl_max_ext(:,:),&
                        nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

        ! CAPS SPRING 2011 SPECIFIC

        CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                        timestr,itime,1,'LTG1_MAX','',                        &
                        'west_east','south_north',nx_ext,ny_ext,ltg1_max_ext(:,:),&
                        nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

        CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                        timestr,itime,1,'LTG2_MAX','',                        &
                        'west_east','south_north',nx_ext,ny_ext,ltg2_max_ext(:,:),&
                        nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

        CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                        timestr,itime,1,'LTG3_MAX','',                        &
                        'west_east','south_north',nx_ext,ny_ext,ltg3_max_ext(:,:),&
                        nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

        ! CAPS SPRING 2012 SPECIFIC

        CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                        timestr,itime,1,'WSPD1KMMAX','',                      &
                        'west_east','south_north',nx_ext,ny_ext,wspd1kmmax_ext(:,:),&
                        nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

        CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                        timestr,itime,1,'CREFD_MAX','',                       &
                        'west_east','south_north',nx_ext,ny_ext,crefd_max_ext(:,:),&
                        nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

        CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                        timestr,itime,1,'REFDM10C_MAX','',                    &
                        'west_east','south_north',nx_ext,ny_ext,refdm10c_max_ext(:,:),&
                        nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

        CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                        timestr,itime,1,'PBLH1','',                           &
                        'west_east','south_north',nx_ext,ny_ext,pblh_ext(:,:),&
                        nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

        CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                        timestr,itime,1,'UP_HELI_MAX16','',                   &
                        'west_east','south_north',nx_ext,ny_ext,up_heli16_max_ext(:,:),&
                        nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

       ! CAPS SPRING 2013 SPECIFIC

        CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                        timestr,itime,1,'SNOW_MAX','',                        &
                        'west_east','south_north',nx_ext,ny_ext,snow_max_ext(:,:),&
                        nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

        CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                        timestr,itime,1,'GRPL05_MAX','',                      &
                        'west_east','south_north',nx_ext,ny_ext,grpl05_max_ext(:,:),&
                        nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

        CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                        timestr,itime,1,'UP_HELI_MAXE','',                    &
                        'west_east','south_north',nx_ext,ny_ext,up_heli_maxe_ext(:,:),&
                        nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

        CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                        timestr,itime,1,'UP_HELI_MAXP','',                    &
                        'west_east','south_north',nx_ext,ny_ext,up_heli_maxp_ext(:,:),&
                        nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)
    END IF
  !
  !-----------------------------------------------------------------------
  !
  ! Downward radiation flux at surface
  !
  !-----------------------------------------------------------------------

    CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,              &
                    timestr,itime,1,'SWDOWN','',                        &
                    'west_east','south_north',nx_ext,ny_ext,tem2d1,     &
                    nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

    CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,              &
                    timestr,itime,1,'GLW','',                           &
                    'west_east','south_north',nx_ext,ny_ext,raddn_ext,  &
                    nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

    raddn_ext = raddn_ext + tem2d1   ! Total downward radiation flux


  !
  !WDT 2004-01-10 GMB: switch from time-step precipitation to accumulated total precip:
  !from ncdump -h:
  ! RAINC:description = "ACCUMULATED TOTAL CUMULUS PRECIPITATION" ;
  ! RAINNC:description = "ACCUMULATED TOTAL GRID SCALE PRECIPITATION" ;
  ! RAINBL:description = "PBL TIME-STEP TOTAL PRECIPITATION" ;
  ! RAINCV:description = "TIME-STEP CUMULUS PRECIPITATION" ;
  !
  !---------------------------------------------------------------------
  !
  ! Cumulus total precipitation
  !
  !---------------------------------------------------------------------

    CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,              &
                    timestr,itime,1,'RAINC','',                         &
                    'west_east','south_north',nx_ext,ny_ext,rainc_ext,  &
                    nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

  !---------------------------------------------------------------------
  !
  ! Grid scale total precipitation
  !
  !---------------------------------------------------------------------

    CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,              &
                    timestr,itime,1,'RAINNC','',                        &
                    'west_east','south_north',nx_ext,ny_ext,raing_ext,  &
                    nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

    CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,              &
                    timestr,itime,1,'SNOWNC','',                        &
                    'west_east','south_north',nx_ext,ny_ext,rainsnow_ext(:,:),&
                    nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)


    CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,              &
                    timestr,itime,1,'GRAUPELNC','',                     &
                    'west_east','south_north',nx_ext,ny_ext,raingrpl_ext(:,:),&
                    nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

    CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,              &
                    timestr,itime,1,'HAILNC','',                        &
                    'west_east','south_north',nx_ext,ny_ext,rainhail_ext(:,:),&
                    nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

  !-----------------------------------------------------------------------
  !
  ! Deallocate the working arrays
  !
  !-----------------------------------------------------------------------

    DEALLOCATE( temlg1, temlg2 )


    RETURN
  END SUBROUTINE get_wrf_data2d

  !##################################################################

  SUBROUTINE get_wrf_rain(fHndl,io_form,multifile,ncmprx,ncmpry,itime,timestr,  &
                          istatus)

  !------------------------------------------------------------------------
  !
  !  PURPOSE:
  !
  !  Read in precipitation related fields in one time level (itime) from
  !  the opened NetCDF file (fHndl)
  !
  !------------------------------------------------------------------------
  !
  !  AUTHOR:
  !  Fanyou Kong (04/07/2007)
  !
  !  MODIFICATION HISTORY:
  !
  !  04/12/2012 (Fanyou Kong)
  !  Add raingrpl_ext and rainhail_ext
  !
  !  04/18/2013 (Fanyou Kong)
  !  Add rainsnow_ext
  !
  !  04/22/2013 (Y. Wang)
  !  Moved to this module and renamed from getwrfrain.
  !
  !------------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER,      INTENT(IN)  :: ncmprx, ncmpry
    INTEGER,      INTENT(IN)  :: fHndl(ncmprx,ncmpry)
    INTEGER,      INTENT(IN)  :: io_form      ! 5 - PHDF5
                                              ! 7 - NetCDF
    LOGICAL,      INTENT(IN)  :: multifile
    INTEGER,      INTENT(IN)  :: itime   ! Recorder number
    CHARACTER(*), INTENT(IN)  :: timestr

    INTEGER,      INTENT(OUT) :: istatus

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

    INTEGER  :: nxlg_ext, nylg_ext
    INTEGER  :: nxd, nyd, nzd            ! dimensions in data file
    INTEGER  :: nxd_stag, nyd_stag, nzd_stag

    REAL, ALLOCATABLE :: temlg1(:,:)
    REAL, ALLOCATABLE :: temlg2(:,:)


  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !
  !  Beginning of executable code...
  !
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

  !-----------------------------------------------------------------------
  !
  ! Allocate working arrays
  !
  !-----------------------------------------------------------------------

    nxlg_ext = (nx_ext-1)*nproc_x + 1
    nylg_ext = (ny_ext-1)*nproc_y + 1

    IF (mp_opt > 0 .AND. multifile) THEN
      nxd = (nx_ext - 1)/ncmprx
      nyd = (ny_ext - 1)/ncmpry
      nxd_stag = nxd
      nyd_stag = nyd
      IF (loc_x == nproc_x) nxd_stag = nxd + 1
      IF (loc_y == nproc_y) nyd_stag = nyd + 1

      ALLOCATE(temlg2(1,1),   STAT = istatus)  ! do not use it
    ELSE
      nxd = nxlg_ext - 1
      nyd = nylg_ext - 1
      nxd_stag = nxlg_ext
      nyd_stag = nylg_ext

      ALLOCATE(temlg2(nxlg_ext,nylg_ext),   STAT = istatus)
    END IF

    ALLOCATE(temlg1(nxd_stag, nyd_stag), STAT = istatus)

  !
  !WDT 2004-01-10 GMB: switch from time-step precipitation to accumulated total precip:
  !from ncdump -h:
  ! RAINC:description = "ACCUMULATED TOTAL CUMULUS PRECIPITATION" ;
  ! RAINNC:description = "ACCUMULATED TOTAL GRID SCALE PRECIPITATION" ;
  ! RAINBL:description = "PBL TIME-STEP TOTAL PRECIPITATION" ;
  ! RAINCV:description = "TIME-STEP CUMULUS PRECIPITATION" ;
  !
  !-----------------------------------------------------------------------
  !
  ! Cumulus total precipitation
  !
  !-----------------------------------------------------------------------

    CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                    timestr,itime,1,'RAINC','',                           &
                    'west_east','south_north',nx_ext,ny_ext,rainc_ext,    &
                    nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

  !-----------------------------------------------------------------------
  !
  ! Grid scale total precipitation
  !
  !-----------------------------------------------------------------------

    CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                    timestr,itime,1,'RAINNC','',                          &
                    'west_east','south_north',nx_ext,ny_ext,raing_ext,    &
                    nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)
  !
  !-----------------------------------------------------------------------
  !
  ! Grid scale total snow precipitation
  !
  !-----------------------------------------------------------------------

    CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                    timestr,itime,1,'SNOWNC','',                          &
                    'west_east','south_north',nx_ext,ny_ext,rainsnow_ext, &
                    nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)
  !
  !-----------------------------------------------------------------------
  !
  ! Grid scale total graupel precipitation
  !
  !-----------------------------------------------------------------------

    CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                    timestr,itime,1,'GRAUPELNC','',                       &
                    'west_east','south_north',nx_ext,ny_ext,raingrpl_ext, &
                    nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)
  !
  !-----------------------------------------------------------------------
  !
  ! Grid scale total hail precipitation
  !
  !-----------------------------------------------------------------------

    CALL get_wrf_2d(fHndl,io_form,multifile,ncmprx,ncmpry,                &
                    timestr,itime,1,'HAILNC','',                          &
                    'west_east','south_north',nx_ext,ny_ext,rainhail_ext, &
                    nxd,nyd,temlg1,nxlg_ext,nylg_ext,temlg2,istatus)

  !-----------------------------------------------------------------------
  !
  ! Deallocate the working arrays
  !
  !-----------------------------------------------------------------------

    DEALLOCATE(temlg1,temlg2)

    RETURN
  END SUBROUTINE get_wrf_rain

  !#####################################################################

  SUBROUTINE make_data2d(use_wrf_grid,iorder,iscl,jscl,x_ext,y_ext,     &
                         xs2d,ys2d,dxfld,dyfld,rdxfld,rdyfld,           &
                         tem1_ext,istatus)

    INTEGER, INTENT(IN)    :: use_wrf_grid, iorder
    REAL,    INTENT(IN)    :: x_ext(nx_ext),y_ext(ny_ext)
    REAL,    INTENT(IN)    :: xs2d(nx,ny), ys2d(nx,ny)  ! arps point coordinate in external grid
    INTEGER, INTENT(IN)    :: iscl(nx,ny), jscl(nx,ny)  ! arps point indices in external array
    REAL,    INTENT(IN)    :: dxfld(nx_ext), dyfld(ny_ext)
    REAL,    INTENT(IN)    :: rdxfld(nx_ext),rdyfld(ny_ext)
    REAL,    INTENT(INOUT) :: tem1_ext(nx_ext,ny_ext,4)

    INTEGER, INTENT(OUT)   :: istatus

  !---------------------------------------------------------------------

    INTEGER :: i, j

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    IF(use_wrf_grid /= 1) THEN

      CALL mkarps2d(nx_ext,ny_ext,nx,ny,                              &
                    iorder,iscl,jscl,x_ext,y_ext,                     &
                    xs2d,ys2d,t2m_ext,t2m,                            &
                    dxfld,dyfld,rdxfld,rdyfld,                        &
                    tem1_ext(:,:,1),tem1_ext(:,:,2),                  &
                    tem1_ext(:,:,3),tem1_ext(:,:,4))
      CALL mkarps2d(nx_ext,ny_ext,nx,ny,                              &
                    iorder,iscl,jscl,x_ext,y_ext,                     &
                    xs2d,ys2d,th2m_ext,th2m,                          &
                    dxfld,dyfld,rdxfld,rdyfld,                        &
                    tem1_ext(:,:,1),tem1_ext(:,:,2),                  &
                    tem1_ext(:,:,3),tem1_ext(:,:,4))
      CALL mkarps2d(nx_ext,ny_ext,nx,ny,                              &
                    iorder,iscl,jscl,x_ext,y_ext,                     &
                    xs2d,ys2d,qv2m_ext,qv2m,                          &
                    dxfld,dyfld,rdxfld,rdyfld,                        &
                    tem1_ext(:,:,1),tem1_ext(:,:,2),                  &
                    tem1_ext(:,:,3),tem1_ext(:,:,4))
      CALL mkarps2d(nx_ext,ny_ext,nx,ny,                              &
                    iorder,iscl,jscl,x_ext,y_ext,                     &
                    xs2d,ys2d,u10m_ext,u10m,                          &
                    dxfld,dyfld,rdxfld,rdyfld,                        &
                    tem1_ext(:,:,1),tem1_ext(:,:,2),                  &
                    tem1_ext(:,:,3),tem1_ext(:,:,4))
      CALL mkarps2d(nx_ext,ny_ext,nx,ny,                              &
                    iorder,iscl,jscl,x_ext,y_ext,                     &
                    xs2d,ys2d,v10m_ext,v10m,                          &
                    dxfld,dyfld,rdxfld,rdyfld,                        &
                    tem1_ext(:,:,1),tem1_ext(:,:,2),                  &
                    tem1_ext(:,:,3),tem1_ext(:,:,4))
      CALL mkarps2d(nx_ext,ny_ext,nx,ny,                              &
                    iorder,iscl,jscl,x_ext,y_ext,                     &
                    xs2d,ys2d,raddn_ext,raddn,                        &
                    dxfld,dyfld,rdxfld,rdyfld,                        &
                    tem1_ext(:,:,1),tem1_ext(:,:,2),                  &
                    tem1_ext(:,:,3),tem1_ext(:,:,4))


      !  Processing PBL grid scale precipitation

      !IF(rainout == 1) THEN

        CALL mkarps2d(nx_ext,ny_ext,nx,ny,                              &
                       iorder,iscl,jscl,x_ext,y_ext,                    &
                       xs2d,ys2d,raing_ext,raing,                       &
                       dxfld,dyfld,rdxfld,rdyfld,                       &
                       tem1_ext(:,:,1),tem1_ext(:,:,2),                 &
                       tem1_ext(:,:,3),tem1_ext(:,:,4))

        CALL mkarps2d(nx_ext,ny_ext,nx,ny,                              &
                      iorder,iscl,jscl,x_ext,y_ext,                     &
                      xs2d,ys2d,rainc_ext,rainc,                        &
                      dxfld,dyfld,rdxfld,rdyfld,                        &
                      tem1_ext(:,:,1),tem1_ext(:,:,2),                  &
                      tem1_ext(:,:,3),tem1_ext(:,:,4))

        CALL mkarps2d(nx_ext,ny_ext,nx,ny,                              &
                      iorder,iscl,jscl,x_ext,y_ext,                     &
                      xs2d,ys2d,rainsnow_ext,rainsnow,                  &
                      dxfld,dyfld,rdxfld,rdyfld,                        &
                      tem1_ext(:,:,1),tem1_ext(:,:,2),                  &
                      tem1_ext(:,:,3),tem1_ext(:,:,4))

        CALL mkarps2d(nx_ext,ny_ext,nx,ny,                              &
                      iorder,iscl,jscl,x_ext,y_ext,                     &
                      xs2d,ys2d,raingrpl_ext,raingrpl,                  &
                      dxfld,dyfld,rdxfld,rdyfld,                        &
                      tem1_ext(:,:,1),tem1_ext(:,:,2),                  &
                      tem1_ext(:,:,3),tem1_ext(:,:,4))

        CALL mkarps2d(nx_ext,ny_ext,nx,ny,                              &
                      iorder,iscl,jscl,x_ext,y_ext,                     &
                      xs2d,ys2d,rainhail_ext,rainhail,                  &
                      dxfld,dyfld,rdxfld,rdyfld,                        &
                      tem1_ext(:,:,1),tem1_ext(:,:,2),                  &
                      tem1_ext(:,:,3),tem1_ext(:,:,4))
      !END IF

    ELSE
      DO j = 2,ny-2
        DO i = 2,nx-2
          t2m(i,j)      = t2m_ext(i-1,j-1)
          th2m(i,j)     = th2m_ext(i-1,j-1)
          qv2m(i,j)     = qv2m_ext(i-1,j-1)
          u10m(i,j)     = u10m_ext(i-1,j-1)
          v10m(i,j)     = v10m_ext(i-1,j-1)
          raddn(i,j)    = raddn_ext(i-1,j-1)

          raing(i,j)    = raing_ext(i-1,j-1)
          rainc(i,j)    = rainc_ext(i-1,j-1)
          rainsnow(i,j) = rainsnow_ext(i-1,j-1)
          raingrpl(i,j) = raingrpl_ext(i-1,j-1)
          rainhail(i,j) = rainhail_ext(i-1,j-1)

        END DO
      END DO
      CALL edgfill(t2m,      1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
      CALL edgfill(th2m,     1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
      CALL edgfill(qv2m,     1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
      CALL edgfill(u10m,     1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
      CALL edgfill(v10m,     1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
      CALL edgfill(raddn,    1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)

      CALL edgfill(raing,    1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
      CALL edgfill(rainc,    1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
      CALL edgfill(rainsnow, 1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
      CALL edgfill(raingrpl, 1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
      CALL edgfill(rainhail, 1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)

    END IF

    raing(:,:)    = max(raing(:,:),    0.0)
    rainc(:,:)    = max(rainc(:,:),    0.0)
    rainsnow(:,:) = max(rainsnow(:,:), 0.0)
    raingrpl(:,:) = max(raingrpl(:,:), 0.0)
    rainhail(:,:) = max(rainhail(:,:), 0.0)

    IF (call_postcore .AND. SPRING_EXPERIMENT ) THEN

      IF ( use_wrf_grid /= 1) THEN

        CALL mkarps2d(nx_ext,ny_ext,nx,ny,                              &
                      iorder,iscl,jscl,x_ext,y_ext,                     &
                      xs2d,ys2d,wspd10max_ext,wspd10max,                &
                      dxfld,dyfld,rdxfld,rdyfld,                        &
                      tem1_ext(:,:,1),tem1_ext(:,:,2),                  &
                      tem1_ext(:,:,3),tem1_ext(:,:,4))
        CALL mkarps2d(nx_ext,ny_ext,nx,ny,                              &
                      iorder,iscl,jscl,x_ext,y_ext,                     &
                      xs2d,ys2d,w_up_max_ext,w_up_max,                  &
                      dxfld,dyfld,rdxfld,rdyfld,                        &
                      tem1_ext(:,:,1),tem1_ext(:,:,2),                  &
                      tem1_ext(:,:,3),tem1_ext(:,:,4))
        CALL mkarps2d(nx_ext,ny_ext,nx,ny,                              &
                      iorder,iscl,jscl,x_ext,y_ext,                     &
                      xs2d,ys2d,w_dn_max_ext,w_dn_max,                  &
                      dxfld,dyfld,rdxfld,rdyfld,                        &
                      tem1_ext(:,:,1),tem1_ext(:,:,2),                  &
                      tem1_ext(:,:,3),tem1_ext(:,:,4))
        CALL mkarps2d(nx_ext,ny_ext,nx,ny,                              &
                      iorder,iscl,jscl,x_ext,y_ext,                     &
                      xs2d,ys2d,refd_max_ext,refd_max,                  &
                      dxfld,dyfld,rdxfld,rdyfld,                        &
                      tem1_ext(:,:,1),tem1_ext(:,:,2),                  &
                      tem1_ext(:,:,3),tem1_ext(:,:,4))
        CALL mkarps2d(nx_ext,ny_ext,nx,ny,                              &
                      iorder,iscl,jscl,x_ext,y_ext,                     &
                      xs2d,ys2d,up_heli_max_ext,up_heli_max,            &
                      dxfld,dyfld,rdxfld,rdyfld,                        &
                      tem1_ext(:,:,1),tem1_ext(:,:,2),                  &
                      tem1_ext(:,:,3),tem1_ext(:,:,4))
        CALL mkarps2d(nx_ext,ny_ext,nx,ny,                              &
                      iorder,iscl,jscl,x_ext,y_ext,                     &
                      xs2d,ys2d,grpl_max_ext,grpl_max,                  &
                      dxfld,dyfld,rdxfld,rdyfld,                        &
                      tem1_ext(:,:,1),tem1_ext(:,:,2),                  &
                      tem1_ext(:,:,3),tem1_ext(:,:,4))

        CALL mkarps2d(nx_ext,ny_ext,nx,ny,                              &
                      iorder,iscl,jscl,x_ext,y_ext,                     &
                      xs2d,ys2d,ltg1_max_ext,ltg1_max,                  &
                      dxfld,dyfld,rdxfld,rdyfld,                        &
                      tem1_ext(:,:,1),tem1_ext(:,:,2),                  &
                      tem1_ext(:,:,3),tem1_ext(:,:,4))
        CALL mkarps2d(nx_ext,ny_ext,nx,ny,                              &
                      iorder,iscl,jscl,x_ext,y_ext,                     &
                      xs2d,ys2d,ltg2_max_ext,ltg2_max,                  &
                      dxfld,dyfld,rdxfld,rdyfld,                        &
                      tem1_ext(:,:,1),tem1_ext(:,:,2),                  &
                      tem1_ext(:,:,3),tem1_ext(:,:,4))
        CALL mkarps2d(nx_ext,ny_ext,nx,ny,                              &
                      iorder,iscl,jscl,x_ext,y_ext,                     &
                      xs2d,ys2d,ltg3_max_ext,ltg3_max,                  &
                      dxfld,dyfld,rdxfld,rdyfld,                        &
                      tem1_ext(:,:,1),tem1_ext(:,:,2),                  &
                      tem1_ext(:,:,3),tem1_ext(:,:,4))
        CALL mkarps2d(nx_ext,ny_ext,nx,ny,                              &
                      iorder,iscl,jscl,x_ext,y_ext,                     &
                      xs2d,ys2d,wspd1kmmax_ext,wspd1kmmax,              &
                      dxfld,dyfld,rdxfld,rdyfld,                        &
                      tem1_ext(:,:,1),tem1_ext(:,:,2),                  &
                      tem1_ext(:,:,3),tem1_ext(:,:,4))
        CALL mkarps2d(nx_ext,ny_ext,nx,ny,                              &
                      iorder,iscl,jscl,x_ext,y_ext,                     &
                      xs2d,ys2d,crefd_max_ext,crefd_max,                &
                      dxfld,dyfld,rdxfld,rdyfld,                        &
                      tem1_ext(:,:,1),tem1_ext(:,:,2),                  &
                      tem1_ext(:,:,3),tem1_ext(:,:,4))
        CALL mkarps2d(nx_ext,ny_ext,nx,ny,                              &
                      iorder,iscl,jscl,x_ext,y_ext,                     &
                      xs2d,ys2d,refdm10c_max_ext,refdm10c_max,          &
                      dxfld,dyfld,rdxfld,rdyfld,                        &
                      tem1_ext(:,:,1),tem1_ext(:,:,2),                  &
                      tem1_ext(:,:,3),tem1_ext(:,:,4))
        CALL mkarps2d(nx_ext,ny_ext,nx,ny,                              &
                      iorder,iscl,jscl,x_ext,y_ext,                     &
                      xs2d,ys2d,pblh_ext,pblh,                          &
                      dxfld,dyfld,rdxfld,rdyfld,                        &
                      tem1_ext(:,:,1),tem1_ext(:,:,2),                  &
                      tem1_ext(:,:,3),tem1_ext(:,:,4))
        CALL mkarps2d(nx_ext,ny_ext,nx,ny,                              &
                      iorder,iscl,jscl,x_ext,y_ext,                     &
                      xs2d,ys2d,up_heli16_max_ext,up_heli16_max,        &
                      dxfld,dyfld,rdxfld,rdyfld,                        &
                      tem1_ext(:,:,1),tem1_ext(:,:,2),                  &
                      tem1_ext(:,:,3),tem1_ext(:,:,4))

        CALL mkarps2d(nx_ext,ny_ext,nx,ny,                              &
                      iorder,iscl,jscl,x_ext,y_ext,                     &
                      xs2d,ys2d,snow_max_ext,snow_max,                  &
                      dxfld,dyfld,rdxfld,rdyfld,                        &
                      tem1_ext(:,:,1),tem1_ext(:,:,2),                  &
                      tem1_ext(:,:,3),tem1_ext(:,:,4))
        CALL mkarps2d(nx_ext,ny_ext,nx,ny,                              &
                      iorder,iscl,jscl,x_ext,y_ext,                     &
                      xs2d,ys2d,grpl05_max_ext,grpl05_max,              &
                      dxfld,dyfld,rdxfld,rdyfld,                        &
                      tem1_ext(:,:,1),tem1_ext(:,:,2),                  &
                      tem1_ext(:,:,3),tem1_ext(:,:,4))
        CALL mkarps2d(nx_ext,ny_ext,nx,ny,                              &
                      iorder,iscl,jscl,x_ext,y_ext,                     &
                      xs2d,ys2d,up_heli_maxe_ext,up_heli_maxe,          &
                      dxfld,dyfld,rdxfld,rdyfld,                        &
                      tem1_ext(:,:,1),tem1_ext(:,:,2),                  &
                      tem1_ext(:,:,3),tem1_ext(:,:,4))
        CALL mkarps2d(nx_ext,ny_ext,nx,ny,                              &
                      iorder,iscl,jscl,x_ext,y_ext,                     &
                      xs2d,ys2d,up_heli_maxp_ext,up_heli_maxp,          &
                      dxfld,dyfld,rdxfld,rdyfld,                        &
                      tem1_ext(:,:,1),tem1_ext(:,:,2),                  &
                      tem1_ext(:,:,3),tem1_ext(:,:,4))

      ELSE

        DO j = 2,ny-2
          DO i = 2,nx-2
            wspd10max(i,j)   = wspd10max_ext(i-1,j-1)
            w_up_max(i,j)    = w_up_max_ext(i-1,j-1)
            w_dn_max(i,j)    = w_dn_max_ext(i-1,j-1)
            refd_max(i,j)    = refd_max_ext(i-1,j-1)
            up_heli_max(i,j) = up_heli_max_ext(i-1,j-1)
            grpl_max(i,j)    = grpl_max_ext(i-1,j-1)
            ltg1_max(i,j)    = ltg1_max_ext(i-1,j-1)
            ltg2_max(i,j)    = ltg2_max_ext(i-1,j-1)
            ltg3_max(i,j)    = ltg3_max_ext(i-1,j-1)
            wspd1kmmax(i,j)  = wspd1kmmax_ext(i-1,j-1)
            crefd_max(i,j)   = crefd_max_ext(i-1,j-1)
            refdm10c_max(i,j)  = refdm10c_max_ext(i-1,j-1)
            pblh(i,j)          = pblh_ext(i-1,j-1)
            up_heli16_max(i,j) = up_heli16_max_ext(i-1,j-1)

            snow_max(i,j)     = snow_max_ext(i-1,j-1)
            grpl05_max(i,j)   = grpl05_max_ext(i-1,j-1)
            up_heli_maxe(i,j) = up_heli_maxe_ext(i-1,j-1)
            up_heli_maxp(i,j) = up_heli_maxp_ext(i-1,j-1)
          END DO
        END DO

        CALL edgfill(wspd10max,     1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
        CALL edgfill(w_up_max,      1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
        CALL edgfill(w_dn_max,      1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
        CALL edgfill(refd_max,      1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
        CALL edgfill(up_heli_max,   1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
        CALL edgfill(grpl_max,      1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
        CALL edgfill(ltg1_max,      1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
        CALL edgfill(ltg2_max,      1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
        CALL edgfill(ltg3_max,      1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
        CALL edgfill(wspd1kmmax,    1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
        CALL edgfill(crefd_max,     1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
        CALL edgfill(refdm10c_max,  1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
        CALL edgfill(pblh,          1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
        CALL edgfill(up_heli16_max, 1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)

        CALL edgfill(snow_max,      1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
        CALL edgfill(grpl05_max,    1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
        CALL edgfill(up_heli_maxe,  1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
        CALL edgfill(up_heli_maxp,  1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
      END IF

    END IF

    RETURN
  END SUBROUTINE make_data2d

  !#####################################################################

  SUBROUTINE print_maxmin_data2d_arps( myproc, istatus )
    INTEGER, INTENT(IN)  :: myproc
    INTEGER, INTENT(OUT) :: istatus

  !---------------------------------------------------------------------
    REAL :: amax, amin

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    CALL a3dmax0(t2m  ,1,nx,1,nx-1,1,ny,1,ny-1,                       &
                 1,1,1,1,amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                 't2m_min = ', amin,', t2m_max =',amax
    CALL a3dmax0(th2m  ,1,nx,1,nx-1,1,ny,1,ny-1,                      &
                 1,1,1,1,amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                 'th2m_min= ', amin,', th2m_max=',amax
    CALL a3dmax0(qv2m  ,1,nx,1,nx-1,1,ny,1,ny-1,                      &
                 1,1,1,1,amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                 'qv2m_min= ', amin,', qv2m_max=',amax
    CALL a3dmax0(u10m  ,1,nx,1,nx-1,1,ny,1,ny-1,                      &
                 1,1,1,1,amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                 'u10m_min= ', amin,', u10m_max=',amax
    CALL a3dmax0(v10m  ,1,nx,1,nx-1,1,ny,1,ny-1,                      &
                 1,1,1,1,amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                 'v10m_min= ', amin,', v10m_max=',amax

    CALL a3dmax0(raing  ,1,nx,1,nx-1,1,ny,1,ny-1,                     &
                 1,1,1,1,amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                 'raing_min= ', amin,', raing_max=',amax
    CALL a3dmax0(rainc  ,1,nx,1,nx-1,1,ny,1,ny-1,                     &
                 1,1,1,1,amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                 'rainc_min= ', amin,', rainc_max=',amax
    CALL a3dmax0(rainsnow  ,1,nx,1,nx-1,1,ny,1,ny-1,                  &
                 1,1,1,1,amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                 'rainsnow_min= ', amin,', rainsnow_max=',amax
    CALL a3dmax0(raingrpl  ,1,nx,1,nx-1,1,ny,1,ny-1,                  &
                 1,1,1,1,amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                 'raingrpl_min= ', amin,', raingrpl_max=',amax
    CALL a3dmax0(rainhail  ,1,nx,1,nx-1,1,ny,1,ny-1,                  &
                 1,1,1,1,amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                 'rainhail_min= ', amin,', rainhail_max=',amax

    CALL a3dmax0(raddn  ,1,nx,1,nx-1,1,ny,1,ny-1,                     &
                 1,1,1,1,amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                 'raddn_min= ', amin,', raddn_max=',amax

    IF (call_postcore) THEN

      CALL a3dmax0(wspd1kmmax  ,1,nx,1,nx-1,1,ny,1,ny-1,                &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'wspd1kmmax_min= ', amin,', wspd1kmmax_max=',amax
      CALL a3dmax0(crefd_max  ,1,nx,1,nx-1,1,ny,1,ny-1,                 &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'crefd_max_min= ', amin,', crefd_max_max=',amax
      CALL a3dmax0(refdm10c_max  ,1,nx,1,nx-1,1,ny,1,ny-1,              &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'refdm10c_max_min= ', amin,', refdm10c_max_max=',amax
      CALL a3dmax0(pblh  ,1,nx,1,nx-1,1,ny,1,ny-1,                      &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'pblh_min= ', amin,', pblh_max=',amax
      CALL a3dmax0(up_heli16_max  ,1,nx,1,nx-1,1,ny,1,ny-1,             &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'up_heli16_max_min= ', amin,', up_heli16_max_max=',amax

      CALL a3dmax0(snow_max  ,1,nx,1,nx-1,1,ny,1,ny-1,                  &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'snow_max_min= ', amin,', snow_max_max=',amax
      CALL a3dmax0(grpl05_max  ,1,nx,1,nx-1,1,ny,1,ny-1,                &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'grpl05_max_min= ', amin,', grpl05_max_max=',amax
      CALL a3dmax0(up_heli_maxe  ,1,nx,1,nx-1,1,ny,1,ny-1,              &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'up_heli_maxe_min= ', amin,', up_heli_maxe_max=',amax
      CALL a3dmax0(up_heli_maxp  ,1,nx,1,nx-1,1,ny,1,ny-1,              &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'up_heli_maxp_min= ', amin,', up_heli_maxp_max=',amax
    END IF

    RETURN
  END SUBROUTINE print_maxmin_data2d_arps

  !#####################################################################

  SUBROUTINE print_maxmin_data2d_ext( myproc, istatus )
    INTEGER, INTENT(IN)  :: myproc
    INTEGER, INTENT(OUT) :: istatus

  !---------------------------------------------------------------------
    REAL :: amax, amin

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    CALL a3dmax0(raddn_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,       &
                 1,1,1,1,amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                          &
                 'raddn_ext_min   = ', amin,', raddn_ext_max   = ',amax
    CALL a3dmax0(t2m_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,         &
                 1,1,1,1,amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                          &
                 't2m_ext_min    = ', amin,', t2m_ext_max    = ',amax
    CALL a3dmax0(th2m_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,        &
                 1,1,1,1,amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                          &
                 'th2m_ext_min   = ', amin,', th2m_ext_max   = ',amax
    CALL a3dmax0(qv2m_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,        &
                 1,1,1,1,amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                          &
                 'qv2m_ext_min   = ', amin,', qv2m_ext_max   = ',amax
    CALL a3dmax0(u10m_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,        &
                 1,1,1,1,amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                          &
                 'u10m_ext_min   = ', amin,', u10m_ext_max   = ',amax
    CALL a3dmax0(v10m_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,        &
                 1,1,1,1,amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                          &
                 'v10m_ext_min   = ', amin,', v10m_ext_max   = ',amax

    CALL a3dmax0(raing_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,       &
                 1,1,1,1,amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                          &
                 'raing_ext_min   = ', amin,', raing_ext_max   = ',amax
    CALL a3dmax0(rainc_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,       &
                 1,1,1,1,amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                          &
                 'rainc_ext_min   = ', amin,', rainc_ext_max   = ',amax
    CALL a3dmax0(rainsnow_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,    &
                 1,1,1,1,amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                          &
                 'rainsnow_ext_min   = ', amin,', rainsnow_ext_max   = ',amax
    CALL a3dmax0(raingrpl_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,    &
                 1,1,1,1,amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                          &
                 'raingrpl_ext_min   = ', amin,', raingrpl_ext_max   = ',amax
    CALL a3dmax0(rainhail_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,    &
                 1,1,1,1,amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                          &
                 'rainhail_ext_min   = ', amin,', rainhail_ext_max   = ',amax

    IF (call_postcore) THEN

      CALL a3dmax0(wspd10max_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext, &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'wspd10max_ext_min   = ', amin,', wspd10max_ext_max   = ',amax
      CALL a3dmax0(w_up_max_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,  &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'w_up_max_ext_min   = ', amin,', w_up_max_ext_max   = ',amax
      CALL a3dmax0(w_dn_max_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,  &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'w_dn_max_ext_min   = ', amin,', w_dn_max_ext_max   = ',amax
      CALL a3dmax0(refd_max_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,  &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'refd_max_ext_min   = ', amin,', refd_max_ext_max   = ',amax
      CALL a3dmax0(up_heli_max_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext, &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'up_heli_max_ext_min   = ', amin,', up_heli_max_ext_max   = ',amax
      CALL a3dmax0(grpl_max_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,  &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'grpl_max_ext_min   = ', amin,', grpl_max_ext_max   = ',amax
      CALL a3dmax0(ltg1_max_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,  &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'ltg1_max_ext_min   = ', amin,', ltg1_max_ext_max   = ',amax
      CALL a3dmax0(ltg2_max_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,  &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'ltg2_max_ext_min   = ', amin,', ltg2_max_ext_max   = ',amax
      CALL a3dmax0(ltg3_max_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,  &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'ltg3_max_ext_min   = ', amin,', ltg3_max_ext_max   = ',amax
      CALL a3dmax0(wspd1kmmax_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext, &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'wspd1kmmax_ext_min   = ', amin,', wspd1kmmax_ext_max   = ',amax
      CALL a3dmax0(crefd_max_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext, &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'crefd_max_ext_min   = ', amin,', crefd_max_ext_max   = ',amax
      CALL a3dmax0(refdm10c_max_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext, &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'refdm10c_max_ext_min   = ', amin,', refdm10c_max_ext_max   = ',amax
      CALL a3dmax0(pblh_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,      &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'pblh_ext_min   = ', amin,', pblh_ext_max   = ',amax
      CALL a3dmax0(up_heli16_max_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext, &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'up_heli16_max_ext_min   = ', amin,', up_heli16_max_ext_max   = ',amax

      CALL a3dmax0(snow_max_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,  &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'snow_max_ext_min   = ', amin,', snow_max_ext_max   = ',amax
      CALL a3dmax0(grpl05_max_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,&
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'grpl05_max_ext_min   = ', amin,', grpl05_max_ext_max   = ',amax
      CALL a3dmax0(up_heli_maxe_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,  &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'up_heli_maxe_ext_min   = ', amin,', up_heli_maxe_ext_max   = ',amax
      CALL a3dmax0(up_heli_maxp_ext  ,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext,  &
                   1,1,1,1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                        &
                   'up_heli_maxp_ext_min   = ', amin,', up_heli_maxp_ext_max   = ',amax
    END IF

    RETURN
  END SUBROUTINE print_maxmin_data2d_ext

  !#####################################################################

  SUBROUTINE print_maxmin_rain( myproc, nxin, nyin, message, affix,     &
                            varg, varc, vars, varp, varh,               &
                            istatus )
  !---------------------------------------------------------------------

    INTEGER, INTENT(IN)  :: myproc
    INTEGER, INTENT(IN)  :: nxin, nyin
    CHARACTER(LEN=*), INTENT(IN) :: message ! message beforehand
    CHARACTER(LEN=*), INTENT(IN) :: affix   ! variable affix
    REAL,    INTENT(IN)  :: varg(nxin,nyin)
    REAL,    INTENT(IN)  :: varc(nxin,nyin)
    REAL,    INTENT(IN)  :: vars(nxin,nyin)
    REAL,    INTENT(IN)  :: varp(nxin,nyin)
    REAL,    INTENT(IN)  :: varh(nxin,nyin)

    INTEGER, INTENT(OUT) :: istatus

  !---------------------------------------------------------------------
    REAL :: amax, amin

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    IF (myproc == 0) WRITE(6,'(1x,a)') TRIM(message)

    CALL a3dmax0(varg,1,nxin,1,nxin,1,nyin,1,nyin,1,1,1,1,amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                          &
                 'raing_'//TRIM(affix)//'_min   = ', amin,              &
               ', raing_'//TRIM(affix)//'_max   = ', amax

    CALL a3dmax0(varc,1,nxin,1,nxin,1,nyin,1,nyin,1,1,1,1,amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                          &
                 'rainc_'//TRIM(affix)//'_min   = ', amin,              &
               ', rainc_'//TRIM(affix)//'_max   = ', amax

    CALL a3dmax0(vars,1,nxin,1,nxin,1,nyin,1,nyin,1,1,1,1,amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                          &
                 'rainsnow_'//TRIM(affix)//'_min   = ', amin,           &
               ', rainsnow_'//TRIM(affix)//'_max   = ', amax

    CALL a3dmax0(varp,1,nxin,1,nxin,1,nyin,1,nyin,1,1,1,1,amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                          &
                 'raingrpl_'//TRIM(affix)//'_min   = ', amin,           &
               ', raingrpl_'//TRIM(affix)//'_max   = ', amax

    CALL a3dmax0(varh,1,nxin,1,nxin,1,nyin,1,nyin,1,1,1,1,amax,amin)
    IF(myproc == 0) WRITE(6,'(1x,2(a,e13.6))')                          &
                 'rainhail_'//TRIM(affix)//'_min   = ', amin,           &
               ', rainhail_'//TRIM(affix)//'_max   = ', amax

    RETURN
  END SUBROUTINE print_maxmin_rain

  !#####################################################################

  SUBROUTINE process_post_arps( curtim,hdmpfn,ldmpf,mp_physics,         & 
                       nz,nzsoil,nstyp,nscalar,x,y,z,zp,                &
                       uprt,vprt,w ,ptprt, pprt, qvprt, qscalar,        &
                       ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,    &
                       soiltyp,stypfrct,vegtyp,lai,roufns,veg,          &
                       tsoil,qsoil,wetcanp,snowdpth,                    &
                       tem1,tem2,tem3,tem4,tem5,tem6,istatus )

  !---------------------------------------------------------------------

    REAL,    INTENT(IN)  :: curtim
    CHARACTER(LEN=*), INTENT(IN) :: hdmpfn
    INTEGER, INTENT(IN)  :: ldmpf

    INTEGER, INTENT(IN)  :: nz, nzsoil, nstyp, nscalar, mp_physics

    REAL,    INTENT(IN) :: x(nx)
    REAL,    INTENT(IN) :: y(ny)
    REAL,    INTENT(IN) :: z(nz)
    REAL,    INTENT(IN) :: zp(nx,ny,nz)
    REAL,    INTENT(IN) :: uprt(nx,ny,nz)
    REAL,    INTENT(IN) :: vprt(nx,ny,nz)
    REAL,    INTENT(IN) :: w(nx,ny,nz)
    REAL,    INTENT(IN) :: pprt(nx,ny,nz)
    REAL,    INTENT(IN) :: ptprt(nx,ny,nz)
    REAL,    INTENT(IN) :: qvprt(nx,ny,nz)
    REAL,    INTENT(IN) :: qscalar(nx,ny,nz,nscalar)
    REAL,    INTENT(IN) :: pbar(nx,ny,nz)
    REAL,    INTENT(IN) :: ptbar(nx,ny,nz)
    REAL,    INTENT(IN) :: qvbar(nx,ny,nz)
    REAL,    INTENT(IN) :: ubar(nx,ny,nz)
    REAL,    INTENT(IN) :: vbar(nx,ny,nz)
    REAL,    INTENT(IN) :: wbar(nx,ny,nz)
    REAL,    INTENT(IN) :: rhobar(nx,ny,nz)
    REAL,    INTENT(IN) :: tsoil(nx,ny,nzsoil,0:nstyp)
    REAL,    INTENT(IN) :: qsoil(nx,ny,nzsoil,0:nstyp)
    REAL,    INTENT(IN) :: wetcanp(nx,ny,0:nstyp)
    REAL,    INTENT(IN) :: snowdpth(nx,ny)

    INTEGER, INTENT(IN) :: soiltyp (nx,ny,nstyp)
    REAL,    INTENT(IN) :: stypfrct(nx,ny,nstyp)
    INTEGER, INTENT(IN) :: vegtyp(nx,ny)

    REAL,    INTENT(IN) :: lai   (nx,ny)
    REAL,    INTENT(IN) :: roufns(nx,ny)
    REAL,    INTENT(IN) :: veg   (nx,ny)

    REAL,    INTENT(INOUT) :: tem1(nx,ny,nz),tem2(nx,ny,nz),tem3(nx,ny,nz)
    REAL,    INTENT(INOUT) :: tem4(nx,ny,nz),tem5(nx,ny,nz),tem6(nx,ny,nz)

    INTEGER, INTENT(OUT) :: istatus

  !---------------------------------------------------------------------

    INTEGER :: lenbin,lengem

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    IF(call_postcore)  THEN

      lenbin = len_trim(outheader)
      lengem = len_trim(gemoutheader)

      CALL postcore(nx,ny,nz,nzsoil,nstyp,curtim,mp_physics+100,        &
              hdmpfn(1:ldmpf),ldmpf,outheader(1:lenbin),                &
              lenbin,gemoutheader(1:lengem),lengem,                     &
              icape,iaccu,iascii,i2dfmt,igempak,                        &
              ilite,iltgci,icrtm,isatid,chbgn,chend,icitm,user_emis,    &
              ibeg_offset,iend_offset,jbeg_offset,jend_offset,          &
              x,y,z,zp,uprt,vprt,w,ptprt, pprt,                         &
              qvprt, qscalar,                                           &
              ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,             &
              soiltyp,stypfrct,vegtyp,lai,roufns,veg,                   &
              tsoil,qsoil,wetcanp,snowdpth,                             &
              raing,rainc,rainsnow,raingrpl,rainhail,                   &
              t2m,th2m,qv2m,u10m,v10m,raddn,                            &
              wspd10max,w_up_max,w_dn_max,refd_max,up_heli_max,grpl_max,&
              ltg1_max,ltg2_max,ltg3_max,                               &
              wspd1kmmax,crefd_max,pblh,                                &
              up_heli16_max,refdm10c_max,                               &
              snow_max,grpl05_max,up_heli_maxe,up_heli_maxp,            &
              tem1,tem2,tem3,tem4,tem5,tem6)
    ELSE
      post2d_ext(:,:) = -9999.9
    END IF

    RETURN
  END SUBROUTINE process_post_arps

  !##################################################################
  !
  SUBROUTINE process_rain(curtim, iabssec, multifile, use_wrf_grid,             &
                  grid_id, io_form, nprocx_in, ncmprx,ncmpry,           &
                  numdigits, rewindyr,myr, dir_ext,                     &
                  nz_ext,iorder,iscl,jscl,x_ext,y_ext,xs2d,ys2d,        &
                  dxfld,dyfld,rdxfld,rdyfld,                            &
                  tem1_ext,istatus)

  !---------------------------------------------------------------------
  !
  ! Convert precipitation saved in raing, rainc etc. as accumulated variables
  ! to difference between two continuous output time levels
  !
  !---------------------------------------------------------------------
    IMPLICIT NONE

    REAL,    INTENT(IN) :: curtim
    INTEGER, INTENT(IN) :: iabssec
    INTEGER, INTENT(IN) :: io_form
    INTEGER, INTENT(IN) :: use_wrf_grid, grid_id, nprocx_in
    LOGICAL, INTENT(IN) :: multifile, rewindyr
    INTEGER, INTENT(IN) :: ncmprx,ncmpry, numdigits
    INTEGER, INTENT(IN) :: myr
    CHARACTER(LEN=*), INTENT(IN) :: dir_ext

    INTEGER, INTENT(IN) :: nz_ext

    REAL,    INTENT(IN) :: x_ext(nx_ext),y_ext(ny_ext)
    REAL,    INTENT(IN) :: dxfld(nx_ext),dyfld(ny_ext),rdxfld(nx_ext),rdyfld(ny_ext)

    INTEGER, INTENT(IN)  :: iorder
    INTEGER, INTENT(IN)  :: iscl(nx,ny),jscl(nx,ny)
    REAL,    INTENT(IN)  :: xs2d(nx,ny),ys2d(nx,ny)

    REAL,    INTENT(OUT) :: tem1_ext (nx_ext,ny_ext,4)

    INTEGER, INTENT(OUT) :: istatus

  !---------------------- Working arrays -------------------------------

    REAL, ALLOCATABLE :: raing0(:,:),   rainc0(:,:),   rainsnow0(:,:)
    REAL, ALLOCATABLE :: raingrpl0(:,:),rainhail0(:,:)

  !-----------------------------------------------------------------------
  !
  ! Misc. local
  !
  !-----------------------------------------------------------------------

    INTEGER :: abstimes0, abstimee0
    CHARACTER(LEN=256) :: extdname(1)
    INTEGER :: nextdfil, idummy(1)

    INTEGER           :: fHndl(ncmprx,ncmpry)
    CHARACTER(LEN=19) :: timestr = ' '

    REAL    :: amax, amin

    INTEGER :: i,j

    INCLUDE 'mp.inc'

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !
  ! Begining of executable code ...
  !
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    IF(call_postcore)  THEN

      IF(iaccu == 1 .AND. int(curtim) >= iabssec1 ) THEN ! piecewise accppt

        IF(myproc==0) PRINT *,'Read previous time level precipitation ...'
        abstimes0 = iabssec - iabssec1
        abstimee0 = iabssec - iabssec1

        IF (io_form /= 7) CALL arpsstop('ERROR: works with netCDF file only.',1)

        CALL check_wrf_files(multifile,1,grid_id,io_form,               &
           nprocx_in,ncmprx,ncmpry,numdigits,                           &
           abstimes0,iabssec1,abstimee0,rewindyr,myr,                   &
           dir_ext,extdname,nextdfil,idummy,istatus)

        IF (istatus /= 0)    CALL arpsstop('ERROR in check_wrf_files, See STDOUT for details',1)
        IF (nextdfil /= 1)   CALL arpsstop('ERROR: Should find 1 file only.',1)
        IF (idummy(1) /= 1)  CALL arpsstop('ERROR: Only works with 1 frame in one file.',1)

        IF (myproc == 0) WRITE(6,'(/,1x,2a,/)') 'Reading file ',TRIM(extdname(1))

        !
        !  blocking inserted for ordering i/o for message passing
        !
        DO i=0,nprocs-1,readstride
          IF(myproc >= i .AND. myproc <= i+readstride-1) THEN

            CALL open_wrf_file(TRIM(extdname(1)),io_form,multifile,.FALSE., &
                               ncmprx,ncmpry,numdigits,fHndl)

            CALL get_wrf_Times(fHndl,io_form,multifile,ncmprx,ncmpry,    &
                               1,timestr)

            CALL get_wrf_rain(fHndl,io_form,multifile,ncmprx,ncmpry,1, timestr,         &
                              istatus)

            CALL close_wrf_file(fHndl,io_form,multifile,.FALSE.,ncmprx,ncmpry)

          END IF
          IF (mp_opt > 0) CALL mpbarrier
        END DO

        IF (istatus /= 0 ) RETURN

        IF(myproc == 0) WRITE(6,'(1x,2a)') 'Read in precipitation from WRF history dumps at ',timestr

        CALL print_maxmin_rain( myproc, nx_ext, ny_ext,                 &
                'WRF precipitation at previous output time', 'ext',     &
                raing_ext, rainc_ext, rainsnow_ext, raingrpl_ext, rainhail_ext, &
                istatus )

        ALLOCATE(raing0(nx,ny),    STAT = istatus)
        ALLOCATE(rainc0(nx,ny),    STAT = istatus)
        ALLOCATE(rainsnow0(nx,ny), STAT = istatus)
        ALLOCATE(raingrpl0(nx,ny), STAT = istatus)
        ALLOCATE(rainhail0(nx,ny), STAT = istatus)
        raing0 = 0.0
        rainc0 = 0.0
        rainsnow0 = 0.0
        raingrpl0 = 0.0
        rainhail0 = 0.0

  !-----------------------------------------------------------------------
  !
  !  Processing PBL grid scale precipitation
  !
  !-----------------------------------------------------------------------

        IF(use_wrf_grid /= 1) THEN
           CALL mkarps2d(nx_ext,ny_ext,nx,ny,                           &
                         iorder,iscl,jscl,x_ext,y_ext,                  &
                         xs2d,ys2d,raing_ext,raing0,                    &
                         dxfld,dyfld,rdxfld,rdyfld,                     &
                         tem1_ext(1,1,1),tem1_ext(1,1,2),               &
                         tem1_ext(1,1,3),tem1_ext(1,1,4))
        ELSE
          DO j = 2,ny-2
            DO i = 2,nx-2
              raing0(i,j) = raing_ext(i-1,j-1)
            END DO
          END DO
          CALL edgfill(raing0, 1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
        END IF

  !-----------------------------------------------------------------------
  !
  !  Processing cumulus precipitation
  !
  !-----------------------------------------------------------------------

        IF(use_wrf_grid /= 1) THEN
          CALL mkarps2d(nx_ext,ny_ext,nx,ny,                            &
                        iorder,iscl,jscl,x_ext,y_ext,                   &
                        xs2d,ys2d,rainc_ext,rainc0,                     &
                        dxfld,dyfld,rdxfld,rdyfld,                      &
                        tem1_ext(1,1,1),tem1_ext(1,1,2),                &
                        tem1_ext(1,1,3),tem1_ext(1,1,4))
        ELSE
          DO j = 2,ny-2
            DO i = 2,nx-2
              rainc0(i,j) = rainc_ext(i-1,j-1)
            END DO
          END DO
          CALL edgfill(rainc0, 1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
        END IF

  !-----------------------------------------------------------------------
  !
  !  Processing PBL grid scale snow precipitation
  !
  !-----------------------------------------------------------------------

        IF(use_wrf_grid /= 1) THEN
          CALL mkarps2d(nx_ext,ny_ext,nx,ny,                            &
                        iorder,iscl,jscl,x_ext,y_ext,                   &
                        xs2d,ys2d,rainsnow_ext,rainsnow0,               &
                        dxfld,dyfld,rdxfld,rdyfld,                      &
                        tem1_ext(1,1,1),tem1_ext(1,1,2),                &
                        tem1_ext(1,1,3),tem1_ext(1,1,4))
        ELSE
          DO j = 2,ny-2
            DO i = 2,nx-2
              rainsnow0(i,j) = rainsnow_ext(i-1,j-1)
            END DO
          END DO
          CALL edgfill(rainsnow0, 1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
        END IF

  !-----------------------------------------------------------------------
  !
  !  Processing PBL grid scale graupel precipitation
  !
  !-----------------------------------------------------------------------

        IF(use_wrf_grid /= 1) THEN
          CALL mkarps2d(nx_ext,ny_ext,nx,ny,                            &
                        iorder,iscl,jscl,x_ext,y_ext,                   &
                        xs2d,ys2d,raingrpl_ext,raingrpl0,               &
                        dxfld,dyfld,rdxfld,rdyfld,                      &
                        tem1_ext(1,1,1),tem1_ext(1,1,2),                &
                        tem1_ext(1,1,3),tem1_ext(1,1,4))
        ELSE
          DO j = 2,ny-2
            DO i = 2,nx-2
              raingrpl0(i,j) = raingrpl_ext(i-1,j-1)
            END DO
          END DO
          CALL edgfill(raingrpl0, 1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
        END IF

  !-----------------------------------------------------------------------
  !
  !  Processing PBL grid scale hail precipitation
  !
  !-----------------------------------------------------------------------

        IF(use_wrf_grid /= 1) THEN
           CALL mkarps2d(nx_ext,ny_ext,nx,ny,                           &
                         iorder,iscl,jscl,x_ext,y_ext,                  &
                         xs2d,ys2d,rainhail_ext,rainhail0,              &
                         dxfld,dyfld,rdxfld,rdyfld,                     &
                         tem1_ext(1,1,1),tem1_ext(1,1,2),               &
                         tem1_ext(1,1,3),tem1_ext(1,1,4))
        ELSE
          DO j = 2,ny-2
            DO i = 2,nx-2
              rainhail0(i,j) = rainhail_ext(i-1,j-1)
            END DO
          END DO
          CALL edgfill(rainhail0, 1,nx,2,nx-2,1,ny,2,ny-2,1,1,1,1)
        END IF

        IF (istatus /= 0) THEN
          PRINT *, 'ERROR: returned from readrain.'
          RETURN
        END IF

        CALL print_maxmin_rain( myproc, nx, ny,                         &
                'Precipitation at ARPS grid at earlier output time', 'pre',     &
                raing0, rainc0, rainsnow0, raingrpl0, rainhail0,        &
                istatus )

  !-----------------------------------------------------------------------
  !
  !  Get precipitation difference
  !
  !-----------------------------------------------------------------------

        raing    = max(raing - raing0, 0.0)
        rainc    = max(rainc - rainc0, 0.0)
        rainsnow = max(rainsnow - rainsnow0, 0.0)
        raingrpl = max(raingrpl - raingrpl0, 0.0)
        rainhail = max(rainhail - rainhail0, 0.0)

        CALL print_maxmin_rain( myproc, nx, ny,                         &
                'Precipitation difference at ARPS grid ', 'out',        &
                raing, rainc, rainsnow, raingrpl, rainhail,             &
                istatus )

        DEALLOCATE(raing0,    STAT = istatus)
        DEALLOCATE(rainc0,    STAT = istatus)
        DEALLOCATE(rainsnow0, STAT = istatus)
        DEALLOCATE(raingrpl0, STAT = istatus)
        DEALLOCATE(rainhail0, STAT = istatus)

      END IF
    END IF

    RETURN
  END SUBROUTINE process_rain

  !#####################################################################

  SUBROUTINE deallocate_data2d_arps( istatus )
    INTEGER, INTENT(OUT) :: istatus
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    DEALLOCATE( t2m,           STAT = istatus )
    DEALLOCATE( th2m,          STAT = istatus )
    DEALLOCATE( qv2m,          STAT = istatus )
    DEALLOCATE( u10m,          STAT = istatus )
    DEALLOCATE( v10m,          STAT = istatus )
    DEALLOCATE( raddn,         STAT = istatus )

    DEALLOCATE( raing,         STAT = istatus )
    DEALLOCATE( rainc,         STAT = istatus )
    DEALLOCATE( rainsnow,      STAT = istatus )
    DEALLOCATE( raingrpl,      STAT = istatus )
    DEALLOCATE( rainhail,      STAT = istatus )

    IF (call_postcore) THEN

      IF ( SPRING_EXPERIMENT ) THEN
        DEALLOCATE( wspd10max,     STAT = istatus )
        DEALLOCATE( w_up_max,      STAT = istatus )
        DEALLOCATE( w_dn_max,      STAT = istatus )
        DEALLOCATE( refd_max,      STAT = istatus )
        DEALLOCATE( up_heli_max,   STAT = istatus )
        DEALLOCATE( grpl_max,      STAT = istatus )
        DEALLOCATE( ltg1_max,      STAT = istatus )
        DEALLOCATE( ltg2_max,      STAT = istatus )
        DEALLOCATE( ltg3_max,      STAT = istatus )
        DEALLOCATE( wspd1kmmax,    STAT = istatus )
        DEALLOCATE( crefd_max,     STAT = istatus )
        DEALLOCATE( refdm10c_max,  STAT = istatus )
        DEALLOCATE( pblh,          STAT = istatus )
        DEALLOCATE( up_heli16_max, STAT = istatus )

        DEALLOCATE( snow_max,     STAT = istatus )
        DEALLOCATE( grpl05_max,   STAT = istatus )
        DEALLOCATE( up_heli_maxe, STAT = istatus )
        DEALLOCATE( up_heli_maxp, STAT = istatus )
      ELSE

        DEALLOCATE(post2d, STAT = istatus)

      END IF

    END IF

    RETURN
  END SUBROUTINE deallocate_data2d_arps

  !#####################################################################

  SUBROUTINE deallocate_data2d_ext( istatus )
    INTEGER, INTENT(OUT) :: istatus
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    DEALLOCATE(psfc_ext,     STAT=istatus)
    DEALLOCATE(t2m_ext ,     STAT=istatus)
    DEALLOCATE(th2m_ext,     STAT=istatus)
    DEALLOCATE(qv2m_ext,     STAT=istatus)
    DEALLOCATE(u10m_ext,     STAT=istatus)
    DEALLOCATE(v10m_ext,     STAT=istatus)

    DEALLOCATE(raing_ext,         STAT=istatus)
    DEALLOCATE(rainc_ext,         STAT=istatus)
    DEALLOCATE(rainsnow_ext,      STAT=istatus)
    DEALLOCATE(raingrpl_ext,      STAT=istatus)
    DEALLOCATE(rainhail_ext,      STAT=istatus)

    DEALLOCATE(raddn_ext,         STAT=istatus)

    IF (call_postcore) THEN

      IF ( SPRING_EXPERIMENT ) THEN
        DEALLOCATE(wspd10max_ext,     STAT=istatus)
        DEALLOCATE(w_up_max_ext,      STAT=istatus)
        DEALLOCATE(w_dn_max_ext,      STAT=istatus)
        DEALLOCATE(refd_max_ext,      STAT=istatus)
        DEALLOCATE(up_heli_max_ext,   STAT=istatus)
        DEALLOCATE(grpl_max_ext,      STAT=istatus)
        DEALLOCATE(ltg1_max_ext,      STAT=istatus)
        DEALLOCATE(ltg2_max_ext,      STAT=istatus)
        DEALLOCATE(ltg3_max_ext,      STAT=istatus)
        DEALLOCATE(wspd1kmmax_ext,    STAT=istatus)
        DEALLOCATE(crefd_max_ext,     STAT=istatus)
        DEALLOCATE(refdm10c_max_ext,  STAT=istatus)
        DEALLOCATE(pblh_ext,          STAT=istatus)
        DEALLOCATE(up_heli16_max_ext, STAT=istatus)

        DEALLOCATE(snow_max_ext,      STAT=istatus)
        DEALLOCATE(grpl05_max_ext,    STAT=istatus)
        DEALLOCATE(up_heli_maxe_ext,  STAT=istatus)
        DEALLOCATE(up_heli_maxp_ext,  STAT=istatus)

      ELSE
        DEALLOCATE(post2d_ext, STAT = istatus)

      END IF
    END IF

    RETURN
  END SUBROUTINE deallocate_data2d_ext

END MODULE module_wrf2arps_post
