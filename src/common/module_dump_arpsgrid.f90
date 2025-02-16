MODULE module_output_arpsgrid
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! This module contains data for the ARPS grid and it also handles ARPS
! related IO (dumps ARPS history files).
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  USE module_arpsgrid_constants

  TYPE type_arpsgrid

    INTEGER :: nx, ny, nz, nzsoil, nstyps
    INTEGER :: nxlg, nylg

    REAL    :: dx, dy

    INTEGER :: mapproj
    REAL    :: ctrlat,  ctrlon
    REAL    :: trulat1, trulat2, trulon

    REAL, ALLOCATABLE :: xlat(:,:),ylat(:,:),slat(:,:)
    REAL, ALLOCATABLE :: xlon(:,:),ylon(:,:),slon(:,:)

    REAL, ALLOCATABLE :: x(:), xs(:)
    REAL, ALLOCATABLE :: y(:), ys(:)
    REAL, ALLOCATABLE :: z(:)
    REAL, ALLOCATABLE :: zp(:,:,:), zpsoil(:,:,:), zps(:,:,:)

    REAL, ALLOCATABLE :: pprt (:,:,:)
    REAL, ALLOCATABLE :: ptprt(:,:,:)
    REAL, ALLOCATABLE :: u (:,:,:)
    REAL, ALLOCATABLE :: v (:,:,:)
    REAL, ALLOCATABLE :: w (:,:,:)
    REAL, ALLOCATABLE :: qv(:,:,:)
    REAL, ALLOCATABLE :: qscalar(:,:,:,:)
    REAL, ALLOCATABLE :: tke(:,:,:)

    !nmm2arps specific
    REAL, ALLOCATABLE :: pbar(:,:,:)  ! Base state pressure (Pascal)
    REAL, ALLOCATABLE :: ptbar(:,:,:) ! Base state potential temperature (K)
    REAL, ALLOCATABLE :: qvbar(:,:,:) ! Base state water vapor specific humidity
                                    ! (kg/kg)
    REAL, ALLOCATABLE :: ubar(:,:,:)  ! Base state u-velocity (m/s)
    REAL, ALLOCATABLE :: vbar(:,:,:)  ! Base state v-velocity (m/s)
    REAL, ALLOCATABLE :: wbar(:,:,:)  ! Base state w-velocity (m/s)

    REAL, ALLOCATABLE :: rhobar(:,:,:)! Base state density (kg/m3).

    INTEGER, ALLOCATABLE :: soiltyp(:,:,:)
    REAL,    ALLOCATABLE :: stypfrct(:,:,:)
    INTEGER, ALLOCATABLE :: vegtyp(:,:)
    REAL,    ALLOCATABLE :: veg(:,:), lai(:,:), roufns(:,:)
    REAL,    ALLOCATABLE :: tsoil(:,:,:,:)
    REAL,    ALLOCATABLE :: qsoil(:,:,:,:)
    REAL,    ALLOCATABLE :: wetcanp(:,:,:)
    REAL,    ALLOCATABLE :: snowdpth(:,:)

    REAL,    ALLOCATABLE :: raing (:,:)
    REAL,    ALLOCATABLE :: rainc (:,:)

!-----------------------------------------------------------------------
    REAL,    ALLOCATABLE :: zpu(:,:,:), zpv(:,:,:)

    REAL,    ALLOCATABLE :: hterain(:,:)
    REAL,    ALLOCATABLE :: mapfct(:,:,:)
    REAL,    ALLOCATABLE :: rhostr(:,:,:)! Base state density rhobar times j3.
    REAL,    ALLOCATABLE :: wcont(:,:,:) ! Contravariant vertical velocity (m/s)
    REAL,    ALLOCATABLE :: j1(:,:,:)
    REAL,    ALLOCATABLE :: j2(:,:,:)
    REAL,    ALLOCATABLE :: j3(:,:,:)
    REAL,    ALLOCATABLE :: j3soil(:,:,:)
    REAL,    ALLOCATABLE :: j3soilinv(:,:,:)
    REAL,    ALLOCATABLE :: aj3z(:,:,:)

    REAL,    ALLOCATABLE :: t2m(:,:),th2m(:,:),qv2m(:,:)
    REAL,    ALLOCATABLE :: u10m(:,:),v10m(:,:),raddn(:,:)
    REAL,    ALLOCATABLE :: wspd10max(:,:),w_up_max(:,:),w_dn_max(:,:)
    REAL,    ALLOCATABLE :: refd_max(:,:),up_heli_max(:,:),grpl_max(:,:)

    LOGICAL :: first_time

  END TYPE type_arpsgrid

  TYPE (type_arpsgrid) :: arpsgrid

  LOGICAL :: arps_allocated = .FALSE.

  CONTAINS

  SUBROUTINE arpsgrid_meta(iproj,truelat1,truelat2,truelong,sclf,       &
                 nx,ny,nz,nstyps,delx,dely,delz,cenlat,cenlong,         &
                 strhoptin,dzminin,zrefsfcin,dlayer1in,dlayer2in,       &
                 strhtunein,zflatin,nzsoil,dzsoilin,soilmodelopt,       &
                 terndtain,ternfmtin,runnamein,dirnamein,               &
                 nocmntin,cmntin,num_moist,num_qq,                      &
        qc_exist,qr_exist,qi_exist,qs_exist,qg_exist,qh_exist,          &
        nc_exist,nr_exist,ni_exist,ns_exist,ng_exist,nh_exist,          &
                 zr_exist,zi_exist,zs_exist,zg_exist,zh_exist,          &
        nn_exist,rim_exist,                                             &
                 dmp_out_joined,terndmpin,soildmpin,hdmpfmtin,          &
                 basoutin, grdoutin, varoutin, iceoutin, sfcoutin,      &
                 landoutin, radoutin, tkeoutin, flxoutin,               &
                 mstoutin, rainoutin, prcoutin, trboutin,               &
                 exbcdmpin, qcexoutin, qrexoutin, qiexoutin,            &
                 qsexoutin, qhexoutin,                                  &
                 hdfcomprin, filcmprsin, readyflin,                     &
                 timestr,initsec,abstimei,abstimee, dbglvl, istatus)

!#######################################################################
!
!  PURPOSE:
!
!    Set meta data in the arpsgrid and the including files from namelist
!    parameters
!
!#######################################################################

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: iproj
    REAL,    INTENT(IN)  :: truelat1,truelat2,truelong,sclf
    INTEGER, INTENT(IN)  :: nx,ny,nz,nstyps,nzsoil
    REAL,    INTENT(IN)  :: delx,dely,delz
    REAL,    INTENT(IN)  :: cenlat,cenlong
    INTEGER, INTENT(IN)  :: strhoptin
    REAL,    INTENT(IN)  :: dzminin,zrefsfcin,dlayer1in,dlayer2in
    REAL,    INTENT(IN)  :: strhtunein,zflatin,dzsoilin
    INTEGER, INTENT(IN)  :: soilmodelopt
    CHARACTER(LEN=MAXFILELEN), INTENT(IN) :: terndtain, dirnamein
    CHARACTER(LEN=80),         INTENT(IN) :: runnamein
    INTEGER, INTENT(IN)  :: ternfmtin
    INTEGER, INTENT(IN)  :: nocmntin
    INTEGER, INTENT(IN)  :: num_moist, num_qq
    INTEGER, INTENT(IN)  :: qc_exist,qr_exist,qi_exist,qs_exist,qg_exist,qh_exist
    INTEGER, INTENT(IN)  :: nc_exist,nr_exist,ni_exist,ns_exist,ng_exist,nh_exist,nn_exist
    INTEGER, INTENT(IN)  :: zr_exist,zi_exist,zs_exist,zg_exist,zh_exist
    INTEGER, INTENT(IN)  :: rim_exist
    CHARACTER(LEN=80), INTENT(IN) :: cmntin(50)
    INTEGER, INTENT(IN)  :: dmp_out_joined
    INTEGER, INTENT(IN)  :: terndmpin, soildmpin
    INTEGER, INTENT(IN)  :: hdmpfmtin
    INTEGER, INTENT(IN)  :: basoutin, grdoutin, varoutin, iceoutin, sfcoutin
    INTEGER, INTENT(IN)  :: landoutin, radoutin, tkeoutin, flxoutin
    INTEGER, INTENT(IN)  :: mstoutin, rainoutin, prcoutin, trboutin
    INTEGER, INTENT(IN)  :: exbcdmpin, qcexoutin, qrexoutin, qiexoutin, &
                            qsexoutin, qhexoutin
    INTEGER, INTENT(IN)  :: hdfcomprin, filcmprsin, readyflin
    CHARACTER(LEN=19), INTENT(IN) :: timestr
    INTEGER, INTENT(IN)  :: initsec, abstimei, abstimee

    INTEGER, INTENT(IN)  :: dbglvl
    INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

    INCLUDE 'mp.inc'
    INCLUDE 'globcst.inc'
    INCLUDE 'grid.inc'
    INCLUDE 'bndry.inc'

    INTEGER :: n,i
    CHARACTER(LEN=1) :: ach

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    lvldbg = dbglvl

!-----------------------------------------------------------------------

    mgrid = 1
    nestgrd = 0

!-----------------------------------------------------------------------

    nstyp   = nstyps
    strhopt = strhoptin
    dzmin   = dzminin
    zrefsfc = zrefsfcin
    dlayer1 = dlayer1in
    dlayer2 = dlayer2in
    strhtune = strhtunein
    zflat    = zflatin
    dx      = delx
    dy      = dely
    dz      = delz
    dzsoil  = dzsoilin
    soilmodel_option = soilmodelopt
    mapproj = iproj
    trulat1 = truelat1
    trulat2 = truelat2
    trulon  = truelong
    sclfct  = sclf
    ctrlat  = cenlat
    ctrlon  = cenlong

    terndmp  = terndmpin
    soildmp  = soildmpin
    hdmpfmt  = hdmpfmtin
    hdfcompr = hdfcomprin

    mstout  = mstoutin
    rainout = rainoutin
    prcout  = prcoutin
    trbout  = trboutin

    basout  = basoutin
    grdout  = grdoutin
    varout  = varoutin
    iceout  = iceoutin
    tkeout  = tkeoutin
    sfcout  = sfcoutin
    landout = landoutin
    radout  = radoutin
    flxout  = flxoutin

    exbcdmp = exbcdmpin
    qcexout = qcexoutin
    qrexout = qrexoutin
    qiexout = qiexoutin
    qsexout = qsexoutin
    qhexout = qhexoutin

    filcmprs = filcmprsin
    readyfl  = readyflin

    nocmnt  = nocmntin
    DO n = 1, nocmntin
      cmnt(n) = cmntin(n)
    END DO
!-----------------------------------------------------------------------
!  Order is important below
!
    nscalar = 0
    IF (qc_exist > 0) THEN
      nscalar = nscalar + 1
      P_QC = nscalar
      nscalarq = nscalar
      qnames(P_QC) = 'qc'; qdescp(P_QC) = 'Cloud water mixing ratio (kg/kg)'
    END IF

    IF (qr_exist > 0) THEN
      nscalar = nscalar + 1
      P_QR = nscalar
      nscalarq = nscalar
      qnames(P_QR) = 'qr'; qdescp(P_QR) = 'Rain  water mixing ratio (kg/kg)'
    END IF

    IF (qi_exist > 0) THEN
      nscalar = nscalar + 1
      P_QI = nscalar
      nscalarq = nscalar
      qnames(P_QI) = 'qi'; qdescp(P_QI) = 'Cloud ice   mixing ratio (kg/kg)'
    END IF

    IF (qs_exist > 0) THEN
      nscalar = nscalar + 1
      P_QS = nscalar
      nscalarq = nscalar
      qnames(P_QS) = 'qs'; qdescp(P_QS) = 'Snow mixing ratio (kg/kg)'
    END IF

    IF (qg_exist > 0) THEN
      nscalar = nscalar + 1
      P_QG = nscalar
      nscalarq = nscalar
      qnames(P_QG) = 'qg'; qdescp(P_QG) = 'Graupel mixing ratio (kg/kg)'
    END IF

    IF (qh_exist > 0) THEN
      nscalar = nscalar + 1
      P_QH = nscalar
      nscalarq = nscalar
      qnames(P_QH) = 'qh'; qdescp(P_QH) = 'Hail mixing ratio (kg/kg)'
    END IF

!-----------------------------------------------------------------------

    IF (nc_exist > 0) THEN
      nscalar = nscalar + 1
      P_NC = nscalar
      qnames(P_NC) = 'nc'; qdescp(P_NC) = 'Cloud water number concentration (#/kg)'
    END IF

    IF (nr_exist > 0) THEN
      nscalar = nscalar + 1
      P_NR = nscalar
      qnames(P_NR) = 'nr'; qdescp(P_NR) = 'Rain number concentration (#/kg)'
    END IF

    IF (ni_exist > 0) THEN
      nscalar = nscalar + 1
      P_NI = nscalar
      qnames(P_NI) = 'ni'; qdescp(P_NI) = 'Ice number concentration (#/kg)'
    END IF

    IF (ns_exist > 0) THEN
      nscalar = nscalar + 1
      P_NS = nscalar
      qnames(P_NS) = 'ns'; qdescp(P_NS) = 'Snow number concentration (#/kg)'
    END IF

    IF (ng_exist > 0) THEN
      nscalar = nscalar + 1
      P_NG = nscalar
      qnames(P_NG) = 'ng'; qdescp(P_NG) = 'Graupel number concentration (#/kg)'
    END IF

    IF (nh_exist > 0) THEN
      nscalar = nscalar + 1
      P_NH = nscalar
      qnames(P_NH) = 'nh'; qdescp(P_NH) = 'Hail number concentration (#/kg)'
    END IF

    IF (nn_exist > 0) THEN    ! CCN Number concentration
      nscalar = nscalar + 1
      P_NI = nscalar
      qnames(P_NI) = 'ni'; qdescp(P_NI) = 'CCN number concentration (#/kg)'
    END IF

    IF (rim_exist > 0) THEN    ! F_RIMEF
      nscalar = nscalar + 1
      qnames(nscalar) = 'F_RIMEF_PHY'; qdescp(nscalar) = 'F_RIMEF from NMM core'
    END IF

!-----------------------------------------------------------------------

    IF (zr_exist > 0) THEN
      nscalar = nscalar + 1
      P_ZR = nscalar
      qnames(P_ZR) = 'zr'; qdescp(P_ZR) = 'Rain reflectivity (m6/kg)'
    END IF

    IF (zi_exist > 0) THEN
      nscalar = nscalar + 1
      P_ZI = nscalar
      qnames(P_ZI) = 'zi'; qdescp(P_ZI) = 'Ice reflectivity (m6/kg)'
    END IF

    IF (zs_exist > 0) THEN
      nscalar = nscalar + 1
      P_ZS = nscalar
      qnames(P_ZS) = 'zs'; qdescp(P_ZS) = 'Snow reflectivity (m6/kg)'
    END IF

    IF (zg_exist > 0) THEN
      nscalar = nscalar + 1
      P_ZG = nscalar
      qnames(P_ZG) = 'zg'; qdescp(P_ZG) = 'Graupel reflectivity (m6/kg)'
    END IF

    IF (zh_exist > 0) THEN
      nscalar = nscalar + 1
      P_ZH = nscalar
      qnames(P_ZH) = 'zh'; qdescp(P_ZH) = 'Hail reflectivity (m6/kg)'
    END IF

!-----------------------------------------------------------------------
    arpsgrid%nz = nz
    arpsgrid%nzsoil = nzsoil
    arpsgrid%nstyps = nstyps

    arpsgrid%dx = dx
    arpsgrid%dy = dy

    arpsgrid%mapproj = iproj
    arpsgrid%trulat1 = truelat1
    arpsgrid%trulat2 = truelat2
    arpsgrid%trulon  = truelong
    arpsgrid%ctrlat  = cenlat
    arpsgrid%ctrlon  = cenlong

    arpsgrid%nx = nx
    arpsgrid%ny = ny

    arpsgrid%nxlg = (nx-3)*nproc_x+3
    arpsgrid%nylg = (ny-3)*nproc_y+3

    arpsgrid%first_time = .TRUE.

!-----------------------------------------------------------------------

    terndta = terndtain
    ternfmt = ternfmtin

    ebc = 4      ! initialized for smooth3d, avgsu, avgsv etc.
    wbc = 4
    nbc = 4
    sbc = 4

    IF (loc_x /= 1)       wbc = 0
    IF (loc_x /= nproc_x) ebc = 0
    IF (loc_y /= 1)       sbc = 0
    IF (loc_y /= nproc_y) nbc = 0

    DO i = FINDX_NUM,1,-1
      n = dmp_out_joined/10**(i-1)
      joindmp(i) = MOD(n,2)
    END DO

    runname = runnamein
    dirname = dirnamein
    CALL gtlfnkey(runname, lfnkey)
    CALL strlnth( dirname, ldirnam)

    IF (mp_opt > 0) THEN
      dumpstride = max_fopen
      readstride = nprocs
      IF (ANY(joindmp > 0)) dumpstride = nprocs   ! join and dump
    END IF

!-----------------------------------------------------------------------
!
!  Time conversions.
!  Formats:  timestr='1998-05-25_18:00:00
!
!-----------------------------------------------------------------------
!
    READ(timestr,'(I4.4,a1,I2.2,a1,I2.2,a1,I2.2,a1,I2.2,a1,I2.2)')      &
        year,ach,month,ach,day,ach,hour,ach,minute,ach,second

    IF (lvldbg > 0)  WRITE(6,'(1x,a,I4.4,a,5(I2.2,a))')                 &
        'In ARPS history file, the initial time will be: ',             &
        year,'-',month,'-',day,'_',hour,':',minute,':',second,'.'

    thisdmp = abstimei
    tstart  = 0.0
    tstop   = abstimee-initsec
    latitud = ctrlat

    RETURN
  END SUBROUTINE arpsgrid_meta

  SUBROUTINE arpsgrid_alloc_init(wrfexttrnopt,nx,ny,nz,istatus)

!#######################################################################
    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: wrfexttrnopt
    INTEGER, INTENT(IN)  :: nx, ny, nz
    INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

    INTEGER :: i,j,k
    INTEGER :: nzsoil, nstyps

    REAL    :: latnot(2)

    REAL, ALLOCATABLE :: xs(:)
    REAL, ALLOCATABLE :: ys(:)

    REAL, ALLOCATABLE :: tem1z(:), tem2z(:)
    REAL, ALLOCATABLE :: tem1(:,:,:)

!-----------------------------------------------------------------------

    INCLUDE 'globcst.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    nzsoil  = arpsgrid%nzsoil
    nstyps  = arpsgrid%nstyps

!-----------------------------------------------------------------------
!
! Allocate arrays in the structure "arpsgrid"
!
!-----------------------------------------------------------------------

    ALLOCATE( arpsgrid%x(nx),  STAT = istatus )
    ALLOCATE( arpsgrid%y(ny),  STAT = istatus )
    ALLOCATE( arpsgrid%xs(nx), STAT = istatus )
    ALLOCATE( arpsgrid%ys(ny), STAT = istatus )
    ALLOCATE( arpsgrid%z(nz),  STAT = istatus )
    ALLOCATE( arpsgrid%zp (nx,ny,nz),        STAT = istatus)
    ALLOCATE( arpsgrid%zps(nx,ny,nz),        STAT = istatus)
    ALLOCATE( arpsgrid%zpsoil(nx,ny,nzsoil), STAT = istatus )
    ALLOCATE( arpsgrid%zpu(nx,ny,nz),        STAT = istatus)
    ALLOCATE( arpsgrid%zpv(nx,ny,nz),        STAT = istatus)

    ALLOCATE( arpsgrid%xlat(nx,ny), STAT = istatus )
    ALLOCATE( arpsgrid%ylat(nx,ny), STAT = istatus )
    ALLOCATE( arpsgrid%slat(nx,ny), STAT = istatus )
    ALLOCATE( arpsgrid%xlon(nx,ny), STAT = istatus )
    ALLOCATE( arpsgrid%ylon(nx,ny), STAT = istatus )
    ALLOCATE( arpsgrid%slon(nx,ny), STAT = istatus )

    ALLOCATE( arpsgrid%pprt (nx,ny,nz), STAT = istatus )
    ALLOCATE( arpsgrid%ptprt(nx,ny,nz), STAT = istatus )
    ALLOCATE( arpsgrid%u (nx,ny,nz), STAT = istatus )
    ALLOCATE( arpsgrid%v (nx,ny,nz), STAT = istatus )
    ALLOCATE( arpsgrid%w (nx,ny,nz), STAT = istatus )
    ALLOCATE( arpsgrid%qv(nx,ny,nz), STAT = istatus )
    ALLOCATE( arpsgrid%qscalar(nx,ny,nz,nscalar), STAT = istatus )
    arpsgrid%qscalar(:,:,:,:) = 0.0

    ALLOCATE( arpsgrid%tke(nx,ny,nz), STAT = istatus )
    arpsgrid%tke(:,:,:) = 0.0

    ALLOCATE( arpsgrid%soiltyp (nx,ny,nstyps), STAT = istatus )
    ALLOCATE( arpsgrid%stypfrct(nx,ny,nstyps), STAT = istatus )
    ALLOCATE( arpsgrid%vegtyp(nx,ny),   STAT = istatus )
    ALLOCATE( arpsgrid%veg   (nx,ny),   STAT = istatus )
    ALLOCATE( arpsgrid%lai   (nx,ny),   STAT = istatus )
    ALLOCATE( arpsgrid%roufns(nx,ny),   STAT = istatus )
    ALLOCATE( arpsgrid%tsoil(nx,ny,nzsoil,0:nstyps), STAT = istatus )
    ALLOCATE( arpsgrid%qsoil(nx,ny,nzsoil,0:nstyps), STAT = istatus )
    ALLOCATE( arpsgrid%wetcanp     (nx,ny,0:nstyps), STAT = istatus )
    ALLOCATE( arpsgrid%snowdpth(nx,ny), STAT = istatus )

    ALLOCATE( arpsgrid%pbar(nx,ny,nz),  STAT = istatus )
    ALLOCATE( arpsgrid%ptbar(nx,ny,nz), STAT = istatus )
    ALLOCATE( arpsgrid%qvbar(nx,ny,nz), STAT = istatus )

    ALLOCATE( arpsgrid%ubar(nx,ny,nz), STAT = istatus )
    ALLOCATE( arpsgrid%vbar(nx,ny,nz), STAT = istatus )
    ALLOCATE( arpsgrid%wbar(nx,ny,nz), STAT = istatus )

    ALLOCATE( arpsgrid%rhobar(nx,ny,nz), STAT = istatus )
    ALLOCATE( arpsgrid%rhostr(nx,ny,nz), STAT = istatus )
    ALLOCATE( arpsgrid%wcont(nx,ny,nz),  STAT = istatus )

    ALLOCATE( arpsgrid%hterain(nx,ny),   STAT = istatus)
    ALLOCATE( arpsgrid%mapfct(nx,ny,8),  STAT = istatus)
    ALLOCATE( arpsgrid%j1(nx,ny,nz),            STAT = istatus)
    ALLOCATE( arpsgrid%j2(nx,ny,nz),            STAT = istatus)
    ALLOCATE( arpsgrid%j3(nx,ny,nz),            STAT = istatus)
    ALLOCATE( arpsgrid%j3soil(nx,ny,nzsoil),    STAT = istatus)
    ALLOCATE( arpsgrid%j3soilinv(nx,ny,nzsoil), STAT = istatus)
    ALLOCATE( arpsgrid%aj3z(nx,ny,nz),          STAT = istatus)

    ALLOCATE( arpsgrid%raing(nx,ny),    STAT = istatus)
    ALLOCATE( arpsgrid%rainc(nx,ny),    STAT = istatus)
    ALLOCATE( arpsgrid%t2m(nx,ny),      STAT = istatus)
    arpsgrid%t2m(:,:) = 0.0
    ALLOCATE( arpsgrid%th2m(nx,ny),     STAT = istatus)
    arpsgrid%th2m(:,:) = 0.0
    ALLOCATE( arpsgrid%qv2m(nx,ny),     STAT = istatus)
    arpsgrid%qv2m = 0.0
    ALLOCATE( arpsgrid%u10m(nx,ny),     STAT = istatus)
    ALLOCATE( arpsgrid%v10m(nx,ny),     STAT = istatus)

    ALLOCATE( arpsgrid%raddn(nx,ny),     STAT = istatus)

    ALLOCATE( arpsgrid%wspd10max(nx,ny),  STAT = istatus)
    ALLOCATE( arpsgrid%w_up_max(nx,ny),   STAT = istatus)
    ALLOCATE( arpsgrid%w_dn_max(nx,ny),   STAT = istatus)
    ALLOCATE( arpsgrid%refd_max(nx,ny),   STAT = istatus)
    ALLOCATE( arpsgrid%up_heli_max(nx,ny),STAT = istatus)
    ALLOCATE( arpsgrid%grpl_max(nx,ny),   STAT = istatus)
!    arpsgrid%wspd10max(:,:)   = ARPS_MISSING
!    arpsgrid%w_up_max(:,:)    = ARPS_MISSING
!    arpsgrid%w_dn_max(:,:)    = ARPS_MISSING
!    arpsgrid%refd_max(:,:)    = ARPS_MISSING
!    arpsgrid%up_heli_max(:,:) = ARPS_MISSING
    arpsgrid%grpl_max(:,:)    = ARPS_MISSING

    arps_allocated = .TRUE.

    arpsgrid%soiltyp (:,:,:) = 0
    arpsgrid%stypfrct(:,:,:) = 0.0

!-----------------------------------------------------------------------
!
!  Set up ARPS map projection
!
!-----------------------------------------------------------------------
!
    latnot(1)=arpsgrid%trulat1
    latnot(2)=arpsgrid%trulat2
    CALL setmapr(arpsgrid%mapproj,1.0,latnot,arpsgrid%trulon)

!-----------------------------------------------------------------------
!
!  Set up ARPS grid
!
!-----------------------------------------------------------------------
!
    !
    ! Working arrays
    !
    ALLOCATE(tem1z(nz), STAT = istatus)
    ALLOCATE(tem2z(nz), STAT = istatus)

    ALLOCATE(tem1(nx,ny,nz), STAT = istatus)

    IF (wrfexttrnopt == 1) THEN ! don't really do grid here, just set map proj
      ternopt = -1
      arpsgrid%hterain(:,:) = 0.0
    ELSE IF(wrfexttrnopt == 2 .OR. wrfexttrnopt == 3) THEN
      ternopt = 2
    ELSE
      ternopt = 0
    END IF

    CALL inigrd(nx,ny,nz,nzsoil,arpsgrid%x,arpsgrid%y,arpsgrid%z,       &
                arpsgrid%zp,arpsgrid%zpsoil,arpsgrid%hterain,           &
                arpsgrid%mapfct,arpsgrid%j1,arpsgrid%j2,arpsgrid%j3,    &
                arpsgrid%j3soil,arpsgrid%j3soilinv,                     &
                tem1z,tem2z,tem1)

    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          arpsgrid%aj3z(i,j,k)=0.5*(arpsgrid%j3(i,j,k)+arpsgrid%j3(i,j,k-1))
        END DO
      END DO
    END DO

    DEALLOCATE(tem1z,tem2z)
    DEALLOCATE(tem1)

!-----------------------------------------------------------------------
!
!   Initialize arpsgrid variables and latitude & longitude
!
!-----------------------------------------------------------------------
!
!    ALLOCATE(xs(nx), STAT = istatus)
!    ALLOCATE(ys(ny), STAT = istatus)

    DO i = 1,nx-1
      arpsgrid%xs(i) = 0.5*( arpsgrid%x(i) + arpsgrid%x(i+1) )
    END DO
    arpsgrid%xs(nx) = arpsgrid%xs(nx-1) + arpsgrid%dx

    DO j = 1,ny-1
      arpsgrid%ys(j) = 0.5*( arpsgrid%y(j) + arpsgrid%y(j+1) )
    END DO
    arpsgrid%ys(ny) = arpsgrid%ys(ny-1) + arpsgrid%dy

    CALL xytoll( nx,ny,arpsgrid%x, arpsgrid%ys, arpsgrid%xlat, arpsgrid%xlon )
    CALL xytoll( nx,ny,arpsgrid%xs,arpsgrid%y,  arpsgrid%ylat, arpsgrid%ylon )
    CALL xytoll( nx,ny,arpsgrid%xs,arpsgrid%ys, arpsgrid%slat, arpsgrid%slon )

!    DEALLOCATE( xs, ys )

    RETURN
  END SUBROUTINE arpsgrid_alloc_init

  SUBROUTINE arpsgrid_dealloc(istatus)

!#######################################################################

    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: istatus

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    DEALLOCATE( arpsgrid%xlat, arpsgrid%ylat, arpsgrid%slat )
    DEALLOCATE( arpsgrid%xlon, arpsgrid%ylon, arpsgrid%slon )

    DEALLOCATE( arpsgrid%x )
    DEALLOCATE( arpsgrid%y )
    DEALLOCATE( arpsgrid%xs, arpsgrid%ys )
    DEALLOCATE( arpsgrid%z )
    DEALLOCATE( arpsgrid%zp  )
    DEALLOCATE( arpsgrid%zps )
    DEALLOCATE( arpsgrid%zpsoil )
    DEALLOCATE( arpsgrid%zpu )
    DEALLOCATE( arpsgrid%zpv )

    DEALLOCATE( arpsgrid%pprt )
    DEALLOCATE( arpsgrid%ptprt )
    DEALLOCATE( arpsgrid%u )
    DEALLOCATE( arpsgrid%v )
    DEALLOCATE( arpsgrid%w )
    DEALLOCATE( arpsgrid%qv )
    DEALLOCATE( arpsgrid%qscalar )

    DEALLOCATE( arpsgrid%tke )

    DEALLOCATE( arpsgrid%soiltyp )
    DEALLOCATE( arpsgrid%stypfrct )
    DEALLOCATE( arpsgrid%vegtyp )
    DEALLOCATE( arpsgrid%veg )
    DEALLOCATE( arpsgrid%lai )
    DEALLOCATE( arpsgrid%roufns )
    DEALLOCATE( arpsgrid%tsoil )
    DEALLOCATE( arpsgrid%qsoil )
    DEALLOCATE( arpsgrid%wetcanp )
    DEALLOCATE( arpsgrid%snowdpth )

    DEALLOCATE( arpsgrid%pbar )
    DEALLOCATE( arpsgrid%ptbar )
    DEALLOCATE( arpsgrid%qvbar )
    DEALLOCATE( arpsgrid%ubar )
    DEALLOCATE( arpsgrid%vbar )
    DEALLOCATE( arpsgrid%wbar )
    DEALLOCATE( arpsgrid%rhobar )
    DEALLOCATE( arpsgrid%rhostr )
    DEALLOCATE( arpsgrid%wcont )

    DEALLOCATE( arpsgrid%hterain, STAT = istatus)
    DEALLOCATE( arpsgrid%mapfct,  STAT = istatus)
    DEALLOCATE( arpsgrid%j1,   STAT = istatus)
    DEALLOCATE( arpsgrid%j2,   STAT = istatus)
    DEALLOCATE( arpsgrid%j3,   STAT = istatus)
    DEALLOCATE( arpsgrid%j3soil,    STAT = istatus)
    DEALLOCATE( arpsgrid%j3soilinv, STAT = istatus)
    DEALLOCATE( arpsgrid%aj3z, STAT = istatus)

    DEALLOCATE( arpsgrid%raing, arpsgrid%rainc )
    DEALLOCATE( arpsgrid%t2m, arpsgrid%th2m, arpsgrid%qv2m )
    DEALLOCATE( arpsgrid%u10m, arpsgrid%v10m, arpsgrid%raddn )
    DEALLOCATE( arpsgrid%wspd10max, arpsgrid%w_up_max, arpsgrid%w_dn_max )
    DEALLOCATE( arpsgrid%refd_max, arpsgrid%up_heli_max, arpsgrid%grpl_max)

    arps_allocated = .FALSE.

    RETURN
  END SUBROUTINE arpsgrid_dealloc

  SUBROUTINE arpsgrid_terrain_heights( wrfexttrnopt, extntmrg, intrpopt,&
                 ims,ime,jms,jme, trn_ext, dbglvl, istatus )

!#######################################################################

    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: wrfexttrnopt
    INTEGER, INTENT(IN)  :: extntmrg
    INTEGER, INTENT(IN)  :: intrpopt
    INTEGER, INTENT(IN)  :: ims, ime, jms, jme
    REAL,    INTENT(IN)  :: trn_ext(ims:ime,jms:jme)
!    REAL,    INTENT(IN)  :: iscl(arpsgrid%nx,arpsgrid%ny)
!    REAL,    INTENT(IN)  :: jscl(arpsgrid%nx,arpsgrid%ny)
    INTEGER, INTENT(IN)  :: dbglvl
    INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

    INCLUDE 'globcst.inc'
    INCLUDE 'mp.inc'

    DOUBLE PRECISION   :: ntmergeinv, mfac

    INTEGER :: onvf, idist
    INTEGER :: nx, ny, nz, nzsoil
    INTEGER :: i, j, k

    REAL, ALLOCATABLE :: htrn_tmp(:,:)
    REAL, ALLOCATABLE :: tem1(:,:,:)
    REAL, ALLOCATABLE :: tem1z(:), tem2z(:)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    nx = arpsgrid%nx
    ny = arpsgrid%ny
    nz = arpsgrid%nz
    nzsoil = arpsgrid%nzsoil
!
!-----------------------------------------------------------------------
!
!  Interpolate external terrain to arps grid (if wrfexttrnopt=1 or 3).
!
!-----------------------------------------------------------------------
!
    ALLOCATE(tem1(nx,ny,nz),  STAT = istatus)

    IF ( wrfexttrnopt == 1 .OR. wrfexttrnopt == 3) THEN
                         ! set arps terrain to match external terrain

      ALLOCATE(htrn_tmp(nx,ny), STAT = istatus)
      ALLOCATE(tem1z(nz),       STAT = istatus)
      ALLOCATE(tem2z(nz),       STAT = istatus)

      CALL mkarps2d(nx,ny,intrpopt,1,ims,ime,jms,jme,                   &
                    trn_ext,htrn_tmp,dbglvl,istatus )

      IF (wrfexttrnopt == 1) THEN
        arpsgrid%hterain(:,:) = htrn_tmp(:,:)
      ELSE
        ntmergeinv = 1.d0/extntmrg
        DO j = 1,ny
          DO i = 1,nx
            idist = max(0,min(extntmrg,i-2,nx-2-i,j-2,ny-2-j))
            mfac = idist*ntmergeinv
            arpsgrid%hterain(i,j) = (1.d0-mfac)*htrn_tmp(i,j)           &
                                  + mfac*arpsgrid%hterain(i,j)
          END DO
        END DO
      END IF

      IF( mp_opt > 0) THEN
        CALL mpsendrecv1dew(arpsgrid%hterain,nx,ny,4,4,0,tem1)
        CALL mpsendrecv1dns(arpsgrid%hterain,nx,ny,4,4,0,tem1)
      END IF
      ternopt = -1
      CALL inigrd(nx,ny,nz,nzsoil,arpsgrid%x,arpsgrid%y,arpsgrid%z,     &
                  arpsgrid%zp,arpsgrid%zpsoil,arpsgrid%hterain,         &
                  arpsgrid%mapfct,arpsgrid%j1,arpsgrid%j2,arpsgrid%j3,  &
                  arpsgrid%j3soil,arpsgrid%j3soilinv,                   &
                  tem1z,tem2z,tem1)

      DO k=2,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            arpsgrid%aj3z(i,j,k)=0.5*(arpsgrid%j3(i,j,k)+arpsgrid%j3(i,j,k-1))
          END DO
        END DO
      END DO

      DEALLOCATE(htrn_tmp)
      DEALLOCATE(tem1z, tem2z)

    END IF
!
!-----------------------------------------------------------------------
!
!  Set up z grid at scalar vertical levels.
!
!-----------------------------------------------------------------------
!
    onvf = 0
    CALL avgz(arpsgrid%zp, onvf, nx, ny, nz, 1,nx, 1,ny, 1,nz-1,        &
              arpsgrid%zps)

    DO j = 1, ny
      DO i = 1, nx
        arpsgrid%zps(i,j,arpsgrid%nz)=2.*arpsgrid%zps(i,j,nz-1)         &
                                       - arpsgrid%zps(i,j,nz-2)
      END DO
    END DO

    CALL edgfill(arpsgrid%zps, 1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)
    CALL edgfill(arpsgrid%zp,  1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)
!
!-----------------------------------------------------------------------
!
!  Set up z grid at U stagger vertical levels.
!
!-----------------------------------------------------------------------
!
    CALL avgsu(arpsgrid%zps,nx,ny,nz,1,ny,1,nz,arpsgrid%zpu,tem1)

    IF (loc_x == 1) THEN  ! West boundary of the whole domain
                          ! avgsu calls bcsu and set it differently
                          ! as what we need.
      DO k=1,nz
        DO j=1,ny
          arpsgrid%zpu(1,j,k) = arpsgrid%zps(1,j,k)
        END DO
      END DO
    END IF

    IF (loc_x == nproc_x) THEN  ! east boundary of the whold domain
                                ! avgsu calls bcsu and set it differently
                                ! as what we need.
      DO k=1,nz
        DO j=1,ny
          arpsgrid%zpu(nx,j,k)= arpsgrid%zps(nx-1,j,k)
        END DO
      END DO
    END IF
!
!-----------------------------------------------------------------------
!
!  Set up z grid at V stagger vertical levels.
!
!-----------------------------------------------------------------------
!
    CALL avgsv(arpsgrid%zps,nx,ny,nz,1,nx,1,nz,arpsgrid%zpv,tem1)

    IF (loc_y == 1) THEN     ! south boundary
      DO k=1,nz
        DO i=1,nx
          arpsgrid%zpv(i,1,k)=arpsgrid%zps(i,1,k)
        END DO
      END DO
    END IF

    IF (loc_y == nproc_y) THEN     ! north boundary
      DO k=1,nz
        DO i=1,nx
          arpsgrid%zpv(i,ny,k) = arpsgrid%zps(i,ny-1,k)
        END DO
      END DO
    END IF

!------------------ Before Return --------------------------------------

    DEALLOCATE(tem1)

    RETURN
  END SUBROUTINE arpsgrid_terrain_heights

  SUBROUTINE arpsgrid_qv2rhstar( iamroot, direction, dbglvl, istatus )

!#######################################################################

    IMPLICIT NONE
    LOGICAL, INTENT(IN)  :: iamroot
    INTEGER, INTENT(IN)  :: direction
    INTEGER, INTENT(IN)  :: dbglvl

    INTEGER, INTENT(OUT) :: istatus

!----------------------------------------------------------------------

    REAL, PARAMETER :: rhmax = 1.0

    INTEGER :: i,j,k
    REAL    :: qvmin, qvmax
    REAL    :: pres, temp
    REAL    :: rh, qvsat

    REAL    :: f_qvsat

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    qvmax=-999.
    qvmin=999.

    IF (direction > 0) THEN  ! Convert from qv to RHstar
      ! Not needed for ARPS grid at present
    ELSE                     ! Convert from RHstar to QV

      DO k=1,arpsgrid%nz
        DO j=1,arpsgrid%ny
          DO i=1,arpsgrid%nx
            pres =  arpsgrid%pbar(i,j,k) +arpsgrid%pprt(i,j,k)
            temp = (arpsgrid%ptbar(i,j,k)+arpsgrid%ptprt(i,j,k))*((pres/p0)**rddcp)
            qvsat=f_qvsat( pres, temp )
            IF(qvsat > 1.) THEN
              PRINT *, ' qvsat trouble at: ',i,j,k
              PRINT *, ' temp,press,qvsat: ',temp,pres,qvsat
            END IF
            rh=AMAX1(0.01,rhmax-(arpsgrid%qv(i,j,k)*arpsgrid%qv(i,j,k)))
            rh=AMIN1(rh,rhmax)
            arpsgrid%qv(i,j,k)=rh*qvsat
            IF(arpsgrid%qv(i,j,k) > 0.5) THEN
              PRINT *, ' qvsat trouble at: ',i,j,k
              PRINT *, ' temp,press,qvsat: ',temp,pres,arpsgrid%qv(i,j,k)
              PRINT *, ' rh= ',rh
            END IF
            qvmax=AMAX1(qvmax,arpsgrid%qv(i,j,k))
            qvmin=AMIN1(qvmin,arpsgrid%qv(i,j,k))
            pres = arpsgrid%pbar(i,j,k)
            temp = arpsgrid%ptbar(i,j,k)*((pres/p0)**rddcp)
            qvsat=f_qvsat( pres, temp )
            rh=AMAX1(0.01,rhmax-(arpsgrid%qvbar(i,j,k)*arpsgrid%qvbar(i,j,k)))
            rh=AMIN1(rh,rhmax)
            arpsgrid%qvbar(i,j,k)=rh*qvsat
          END DO
        END DO
      END DO

      CALL mpmax0(qvmax,qvmin)
      IF ( IAMROOT ) WRITE(6,'(3x,a,2(a,f12.8))')                       &
          'ARPS     QV diagnostic information: ',                       &
          'qv     min = ',qvmin, ', max = ',qvmax

    END IF

    RETURN
  END SUBROUTINE arpsgrid_qv2rhstar

  SUBROUTINE arpsgrid_adjustdata(iamroot, earthwind, nscalarq,          &
                            csopt, csfactr, csound,                     &
                            hydradj,dz, wndadj, obropt, obrzero,        &
                            ext_lbc, ext_vbc, tbc, bbc, dbglvl, istatus)
!#######################################################################

    IMPLICIT NONE

    LOGICAL, INTENT(IN)  :: iamroot
    LOGICAL, INTENT(IN)  :: earthwind
    INTEGER, INTENT(IN)  :: nscalarq
    INTEGER, INTENT(IN)  :: csopt, hydradj
    REAL,    INTENT(IN)  :: csfactr, csound, dz
    INTEGER, INTENT(IN)  :: wndadj, obropt
    REAL,    INTENT(IN)  :: obrzero
    INTEGER, INTENT(IN)  :: ext_lbc, ext_vbc, tbc, bbc
    INTEGER, INTENT(IN)  :: dbglvl
    INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

    REAL, POINTER     :: uatv(:,:,:), vatu(:,:,:)
    REAL, ALLOCATABLE :: tem1(:,:,:), tem2(:,:,:)
    REAL, POINTER     :: tem3(:,:,:), tem4(:,:,:)
    REAL, ALLOCATABLE :: tem5(:,:,:), tem6(:,:,:)
    REAL, ALLOCATABLE :: tem7(:,:,:), tem8(:,:,:)
    REAL, POINTER     :: csndsq(:,:,:)

    INTEGER :: nx, ny, nz, nq
    INTEGER :: i,j,k

    REAL    :: csconst, pconst
    REAL    :: tvbar, temp
    REAL    :: qvprt, qtot

!-----------------------------------------------------------------------

    INCLUDE 'mp.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    nx = arpsgrid%nx
    ny = arpsgrid%ny
    nz = arpsgrid%nz

    ALLOCATE(tem1(nx,ny,nz), STAT = istatus)
    ALLOCATE(tem2(nx,ny,nz), STAT = istatus)
    ALLOCATE(tem3(nx,ny,nz), STAT = istatus)
    ALLOCATE(tem4(nx,ny,nz), STAT = istatus)

!-----------------------------------------------------------------------
!
!  External data comes in oriented to true north.
!  Change that to u,v in the new coordinate system.
!
!-----------------------------------------------------------------------
!
    IF (earthwind) THEN
      IF (dbglvl > 3) WRITE(6,'(1x,a)') '--- Rotate wind to ARPS map projection ...'

      uatv => tem3
      vatu => tem4

      DO k = 1, nz
        DO j = 2, ny
          DO i = 1, nx-1
            uatv(i,j,k) = 0.25*( arpsgrid%u(i,  j,k) + arpsgrid%u(i+1,j,  k) + &
                                 arpsgrid%u(i,j-1,k) + arpsgrid%u(i+1,j-1,k) )
          END DO
          uatv(nx,j,k) = 0.5* ( arpsgrid%u(nx,j,k) + arpsgrid%u(nx,j-1,k) )
        END DO

        DO i = 1, nx-1
          uatv(i,1,k)  = 0.5* ( arpsgrid%u(i,1,k) + arpsgrid%u(i+1,1,k) )
        END DO
        uatv(nx,1,k) = arpsgrid%u(nx,1,k)
      END DO

      DO k = 1, nz
        DO j = 1, ny-1
          DO i = 2, nx
            vatu(i,j,k) = 0.25*( arpsgrid%v(i-1,j+1,k) + arpsgrid%v(i,j+1,k) + &
                                 arpsgrid%v(i-1,j,  k) + arpsgrid%v(i,j,  k) )
          END DO
          vatu(1,j,k) = 0.5 * ( arpsgrid%v(1,j,k) + arpsgrid%v(1,j+1,k) )
        END DO

        DO i = 2, nx
          vatu(i,ny,k)  = 0.5* ( arpsgrid%v(i-1,ny,k) + arpsgrid%v(i,ny,k) )
        END DO
        vatu(1,ny,k) = arpsgrid%v(1,ny,k)
      END DO

      IF (mp_opt > 0) THEN
        CALL mpsendrecv2dew(uatv,nx,ny,nz,0,0,2,tem1)
        CALL mpsendrecv2dns(vatu,nx,ny,nz,0,0,1,tem1)
      END IF

      DO k=1,arpsgrid%nz
        CALL uvetomp(nx,ny,arpsgrid%u(:,:,k),vatu(:,:,k),arpsgrid%xlon,   &
                     tem1(:,:,k),tem2(:,:,k))
      END DO
      arpsgrid%u(:,:,:) = tem1(:,:,:)

      DO k=1,arpsgrid%nz
        CALL uvetomp(nx,ny,uatv(:,:,k),arpsgrid%v(:,:,k),arpsgrid%ylon,   &
                     tem1(:,:,k),tem2(:,:,k))
      END DO
      arpsgrid%v(:,:,:) = tem2(:,:,:)
    END IF

!
!-----------------------------------------------------------------------
!
!  Calculate and store the sound wave speed squared in csndsq.
!
!-----------------------------------------------------------------------
!
    csndsq => tem4

    IF (dbglvl > 3) WRITE(6,'(1x,a,I2)') '--- Sound wave adjustment with csopt = ',csopt

    IF(csopt == 1) THEN       ! Original definition of sound speed.
      DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            csndsq(i,j,k)= cpdcv*arpsgrid%pbar(i,j,k)/arpsgrid%rhobar(i,j,k)
          END DO
        END DO
      END DO
    ELSE IF(csopt == 2) THEN   ! Original sound speed multiplied
                               ! by a factor
      csconst = csfactr**2*cpdcv
      DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            csndsq(i,j,k)= csconst * arpsgrid%pbar(i,j,k)/arpsgrid%rhobar(i,j,k)
          END DO
        END DO
      END DO
    ELSE                      ! Specified constant sound speed.
      DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            csndsq(i,j,k)= csound * csound
          END DO
        END DO
      END DO
    END IF

    IF (dbglvl > 3) WRITE(6,'(1x,a,I2)') '--- Hydrostatic adjustment with hydradj = ',hydradj

    IF(hydradj == 1 .OR. hydradj == 2) THEN
      pconst=0.5*g*dz
!
!-----------------------------------------------------------------------
!
!  Create thermal buoyancy at each scalar point,
!  which is stored in tem2
!
!-----------------------------------------------------------------------
!
      DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            qvprt=arpsgrid%qv(i,j,k)-arpsgrid%qvbar(i,j,k)
!            qtot =arpsgrid%qc(i,j,k)+arpsgrid%qr(i,j,k)                 &
!                 +arpsgrid%qi(i,j,k)+arpsgrid%qs(i,j,k)                 &
!                 +arpsgrid%qh(i,j,k)
            qtot=0.0
            DO nq=1,nscalarq
             qtot=qtot+arpsgrid%qscalar(i,j,k,nq)
            END DO
            tem2(i,j,k) = arpsgrid%j3(i,j,k)*arpsgrid%rhobar(i,j,k)*g*  &
                     ( (arpsgrid%ptprt(i,j,k)/arpsgrid%ptbar(i,j,k))    &
                      +(qvprt/(rddrv+arpsgrid%qvbar(i,j,k)))            &
                      -((qvprt+qtot)/(1.+arpsgrid%qvbar(i,j,k))) )
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!
!  Average thermal buoyancy to w points
!  which is stored in tem3
!
!-----------------------------------------------------------------------
!
      CALL avgsw(tem2,nx,ny,nz,1,nx,1,ny, tem3)

      IF(hydradj == 1) THEN

        DO i=1,nx
          DO j=1,ny
            tem1(i,j,1)=arpsgrid%pprt(i,j,1)
            DO k=2,nz-2
              tem1(i,j,k)=                                              &
              (  (1.-(pconst*arpsgrid%j3(i,j,k-1)/csndsq(i,j,k-1)))     &
                * tem1(i,j,k-1) + dz*tem3(i,j,k) ) /                    &
              (1. + (pconst*arpsgrid%j3(i,j,k)/csndsq(i,j,k)))
              arpsgrid%pprt(i,j,k)=tem1(i,j,k)
            END DO
            arpsgrid%pprt(i,j,nz-1)=arpsgrid%pprt(i,j,nz-2)
          END DO
        END DO

      ELSE IF(hydradj == 2) THEN

        DO i=1,nx
          DO j=1,ny
            tem1(i,j,nz-1)=arpsgrid%pprt(i,j,nz-1)
            DO k=nz-2,2,-1
              tem1(i,j,k)=                                            &
                    ( (1.+ (pconst*arpsgrid%j3(i,j,k+1)/csndsq(i,j,k+1)) )*    &
                        tem1(i,j,k+1) -                               &
                        dz*tem3(i,j,k+1) ) /                          &
                        (1.- (pconst*arpsgrid%j3(i,j,k  )/csndsq(i,j,k  )) )
              arpsgrid%pprt(i,j,k)=tem1(i,j,k)
            END DO
            arpsgrid%pprt(i,j,1)=arpsgrid%pprt(i,j,2)
          END DO
        END DO
      END IF   ! hydradj = 1 or 2 ends here

    ELSE IF (hydradj == 3) THEN
!
!-----------------------------------------------------------------------
!
!  Calculate total pressure, store in tem1.
!  Calculate virtual temperature, store in tem2.
!
!-----------------------------------------------------------------------
!
      DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            tem1(i,j,k) = arpsgrid%pbar(i,j,k)+arpsgrid%pprt(i,j,k)
            temp = (arpsgrid%ptbar(i,j,k)+arpsgrid%ptprt(i,j,k))*       &
                   ((tem1(i,j,k)/p0)**rddcp)
            tem2(i,j,k) = temp*(1.0+rvdrd*arpsgrid%qv(i,j,k))/          &
                          (1.0+arpsgrid%qv(i,j,k))
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!
!  Integrate hydrostatic equation, begining at k=2
!
!-----------------------------------------------------------------------
!
      pconst=-g/rd
      DO k=3,nz-1
        DO j=1,ny
          DO i=1,nx
            tvbar=0.5*(tem2(i,j,k)+tem2(i,j,k-1))
            tem1(i,j,k)=tem1(i,j,k-1)*                                &
                       EXP(pconst*(arpsgrid%zps(i,j,k)-arpsgrid%zps(i,j,k-1))/tvbar)
            arpsgrid%pprt(i,j,k)=tem1(i,j,k)-arpsgrid%pbar(i,j,k)
          END DO
        END DO
      END DO
      DO j=1,ny
        DO i=1,nx
          tvbar=0.5*(tem2(i,j,2)+tem2(i,j,1))
          tem1(i,j,1)=tem1(i,j,2)*                                    &
                     EXP(pconst*(arpsgrid%zps(i,j,1)-arpsgrid%zps(i,j,2))/tvbar)
          arpsgrid%pprt(i,j,1)=tem1(i,j,1)-arpsgrid%pbar(i,j,1)
        END DO
      END DO
    END IF
!
!-----------------------------------------------------------------------
!
!  Find vertical velocity and make any u-v adjustments
!  according to wndadj option.
!
!-----------------------------------------------------------------------
!

    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          arpsgrid%rhostr(i,j,k)=ABS(arpsgrid%j3(i,j,k))*arpsgrid%rhobar(i,j,k)
        END DO
      END DO
    END DO
!    IF (mp_opt > 0) THEN
!      CALL acct_interrupt(mp_acct)
!      CALL mpsendrecv2dew(rhostr,nx,ny,nz,4,4,0,tem1)
!      CALL mpsendrecv2dns(rhostr,nx,ny,nz,4,4,0,tem1)
!    END IF

    IF (dbglvl > 3) WRITE(6,'(1x,a,I2)') '--- Wind adjustment with wndadj = ',wndadj

    IF (wndadj > 0) THEN

      ALLOCATE(tem5(nx,ny,nz), STAT = istatus)
      ALLOCATE(tem6(nx,ny,nz), STAT = istatus)
      ALLOCATE(tem7(nx,ny,nz), STAT = istatus)
      ALLOCATE(tem8(nx,ny,nz), STAT = istatus)

      CALL adjuvw( nx,ny,nz, arpsgrid%u,arpsgrid%v,arpsgrid%w,arpsgrid%wcont, &
                   arpsgrid%ptprt,arpsgrid%ptbar,                             &
                   arpsgrid%zp,arpsgrid%j1,arpsgrid%j2,arpsgrid%j3,           &
                   arpsgrid%aj3z,arpsgrid%mapfct,arpsgrid%rhostr,tem1,        &
                   wndadj,obropt,obrzero,0,                                   &
                   tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8)

      DEALLOCATE(tem5, tem6, tem7, tem8 )

    END IF
!
!-----------------------------------------------------------------------
!
!  Enforce vertical w boundary conditions
!
!-----------------------------------------------------------------------
!
    IF (dbglvl > 3) WRITE(6,'(1x,a)') '--- Enforce vertical W boundary conditions ...'

    IF (ext_lbc == 1 .or. ext_vbc == 1) THEN
      CALL rhouvw(nx,ny,nz,arpsgrid%rhostr,tem1,tem2,tem3)
    END IF


    IF (ext_vbc == 1) THEN
      IF (mp_opt > 0) THEN
        CALL mpsendrecv2dew(arpsgrid%u,nx,ny,nz,4,4,1,tem4)
        CALL mpsendrecv2dns(arpsgrid%v,nx,ny,nz,4,4,2,tem4)
        CALL mpsendrecv2dew(arpsgrid%j1,nx,ny,nz,4,4,0,tem4)
        CALL mpsendrecv2dns(arpsgrid%j2,nx,ny,nz,4,4,0,tem4)
      END IF
      CALL vbcw(nx,ny,nz,arpsgrid%w,arpsgrid%wcont,tbc,bbc,             &
                arpsgrid%u,arpsgrid%v,arpsgrid%rhostr,                  &
                tem1,tem2,tem3,arpsgrid%j1,arpsgrid%j2,arpsgrid%j3)
    END IF
!
!
!-----------------------------------------------------------------------
!
!  Assign zero-gradient horizontal boundary conditions
!  to the horizontal & vertical winds.
!
!-----------------------------------------------------------------------
!
    IF (ext_lbc == 1) THEN
      CALL lbcw(nx,ny,nz,1.0, arpsgrid%w,arpsgrid%wcont,                &
                tem1,tem2,tem3,tem4,3,3,3,3,3,3,3,3)
      CALL bcu (nx,ny,nz,1.0, arpsgrid%u,                               &
                tem1,tem2,tem3,tem4,3,3,3,3,1,1,3,3,3,3)
      CALL bcv (nx,ny,nz,1.0, arpsgrid%v,                               &
                tem1,tem2,tem3,tem4, 3,3,3,3,1,1,3,3,3,3)
    END IF

!-----------------------------------------------------------------------
!
! Make sure precipiation are in good range
!
!-----------------------------------------------------------------------

    WHERE(arpsgrid%raing < 0.0) arpsgrid%raing = 0.0
    WHERE(arpsgrid%rainc < 0.0) arpsgrid%rainc = 0.0

!----------------------------------- Before Return ---------------------

    DEALLOCATE(tem1, tem2)
    DEALLOCATE(tem3, tem4)

    RETURN
  END SUBROUTINE arpsgrid_adjustdata

  SUBROUTINE arpsgrid_ternout( dbglvl, istatus )
!#######################################################################

    IMPLICIT NONE

    INTEGER, INTENT(IN)  :: dbglvl
    INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
    INCLUDE 'globcst.inc'
    INCLUDE 'grid.inc'
    INCLUDE 'mp.inc'

    CHARACTER(LEN=MAXFILELEN) :: ternfn, tmpstr
    INTEGER                   :: lternfn

    INTEGER :: i

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    ternfn = runname(1:lfnkey)//".trndata"
    lternfn = lfnkey + 8

    IF( dirname /= ' ' ) THEN
      tmpstr = ternfn
      ternfn = dirname(1:ldirnam)//tmpstr
      lternfn  = lternfn + ldirnam
    END IF

    CALL fnversn(ternfn, lternfn )

    IF(myproc == 0) PRINT *, 'Write terrain data to ',ternfn(1:lternfn)

    IF(mp_opt >0 .AND. joindmp(FINDX_T) > 0) THEN
      CALL writjointrn(arpsgrid%nx,arpsgrid%ny,ternfn(1:lternfn),dx,dy, &
                   mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon, &
                   arpsgrid%hterain)
    ELSE
      IF (mp_opt > 0) THEN
        tmpstr = ternfn
        CALL gtsplitfn(tmpstr,1,1,loc_x,loc_y,1,1,                      &
                       0,0,0,lvldbg,ternfn,istatus)
        lternfn = LEN_TRIM(ternfn)
      END IF

      ! blocking inserted for ordering i/o for message passing
      DO i=0,nprocs-1, dumpstride
        IF ( myproc >= i .AND. myproc <= i+dumpstride-1 ) THEN
          CALL writtrn(arpsgrid%nx,arpsgrid%ny,ternfn(1:lternfn),dx,dy, &
                  mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,  &
                  arpsgrid%hterain)
          IF(terndmp == 1)  &
            CALL trncntl(arpsgrid%nx,arpsgrid%ny,ternfn(1:lternfn),     &
                  arpsgrid%x,arpsgrid%y)
        END IF
        IF (mp_opt > 0) CALL mpbarrier
      END DO

    END IF

    RETURN
  END SUBROUTINE arpsgrid_ternout

  SUBROUTINE arpsgrid_soilout( curtime, dbglvl, istatus )
  !#######################################################################

    IMPLICIT NONE
    REAL,    INTENT(IN)  :: curtime
    INTEGER, INTENT(IN)  :: dbglvl
    INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
    INCLUDE 'globcst.inc'
    INCLUDE 'grid.inc'
    INCLUDE 'mp.inc'

    CHARACTER(LEN=MAXFILELEN) :: soiloutfl, tmpstr
    INTEGER                   :: lfn, tmstrln

    CHARACTER (LEN=80) :: timsnd

    INTEGER ::  tsoilout,qsoilout,zpsoilout, wetcout, snowdout

    INTEGER :: i

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    zpsoilout = 1
    tsoilout  = 1
    qsoilout  = 1
    wetcout   = 1
    snowdout  = 1

    curtim    = curtime
    CALL cvttsnd( curtim, timsnd, tmstrln )

    soiloutfl = runname(1:lfnkey)//".soilvar."//timsnd(1:tmstrln)
    lfn = lfnkey + 9 + tmstrln

    IF( dirname /= ' ' ) THEN
      tmpstr = soiloutfl
      soiloutfl = dirname(1:ldirnam)//tmpstr
      lfn  = lfn + ldirnam
    END IF

    !CALL fnversn(soiloutfl, lfn)

    !IF(myproc == 0) WRITE(6,'(1x,2a)') 'Write soil initial data to ',soiloutfl(1:lfn)

    IF(mp_opt >0 .AND. joindmp(FINDX_S) > 0) THEN
      CALL wrtjoinsoil(arpsgrid%nx,arpsgrid%ny,arpsgrid%nzsoil,arpsgrid%nstyps, &
             soiloutfl,dx,dy,arpsgrid%zpsoil,                           &
             mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,       &
             zpsoilout, tsoilout,qsoilout, wetcout,snowdout,            &
             arpsgrid%tsoil,arpsgrid%qsoil,arpsgrid%wetcanp,            &
             arpsgrid%snowdpth, arpsgrid%soiltyp)
    ELSE
      !IF (mp_opt > 0) THEN
      !  tmpstr = soiloutfl
      !  CALL gtsplitfn(tmpstr,1,1,loc_x,loc_y,1,1,                      &
      !                 0,0,0,lvldbg,soiloutfl,istatus)
      !END IF

      ! blocking inserted for ordering i/o for message passing
      DO i=0, nprocs-1, dumpstride
        IF (myproc >= i .AND. myproc <= i+dumpstride-1) THEN

          CALL wrtsoil(arpsgrid%nx,arpsgrid%ny,arpsgrid%nzsoil,arpsgrid%nstyps, &
             soiloutfl,dx,dy,arpsgrid%zpsoil,                           &
             mapproj,trulat1,trulat2,trulon,sclfct,ctrlat,ctrlon,       &
             zpsoilout, tsoilout,qsoilout, wetcout,snowdout,            &
             arpsgrid%tsoil,arpsgrid%qsoil,arpsgrid%wetcanp,            &
             arpsgrid%snowdpth, arpsgrid%soiltyp)

          IF (soildmp == 1)            &
            CALL soilcntl(arpsgrid%nx,arpsgrid%ny,arpsgrid%nzsoil,      &
                  arpsgrid%zpsoil,soiloutfl,                            &
                  zpsoilout,tsoilout,qsoilout, wetcout,snowdout,        &
                  arpsgrid%x,arpsgrid%y)

        END IF
        IF (mp_opt > 0) CALL mpbarrier
      END DO

    END IF

    RETURN
  END SUBROUTINE arpsgrid_soilout

  SUBROUTINE arpsgrid_histout( IAMROOT, curtime, dbglvl, istatus )
!#######################################################################

    IMPLICIT NONE

    LOGICAL, INTENT(IN)  :: IAMROOT
    REAL,    INTENT(IN)  :: curtime
    INTEGER, INTENT(IN)  :: dbglvl
    INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
    INCLUDE 'globcst.inc'
    INCLUDE 'grid.inc'
    INCLUDE 'mp.inc'

    CHARACTER (LEN=MAXFILELEN) :: basdmpfn
    INTEGER                    :: lbasdmpf

    INTEGER :: grdbas
    INTEGER :: mstout0, rainout0, prcout0, trbout0

    INTEGER :: iniotfu = 21

    INTEGER :: i

    REAL    :: amax, amin
    INTEGER :: nx, ny, nz

    REAL, ALLOCATABLE :: tem1(:,:,:),tem2(:,:,:),tem3(:,:,:),tem4(:,:,:)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    curtim = curtime

    IF (dbglvl > 0 .AND. IAMROOT ) WRITE(6,'(1x,a,f12.2,a)')   &
      'Writing ARPS history file at forecast time: ',curtim,' seconds.'

!
!-----------------------------------------------------------------------
!
!  Print out the max/min of output variables.
!
!-----------------------------------------------------------------------
!
    nx = arpsgrid%nx
    ny = arpsgrid%ny
    nz = arpsgrid%nz

    IF (dbglvl > 1 .AND. mp_opt == 0) THEN

      IF(IAMROOT) WRITE(6,'(/1x,a/)')                                   &
          'Min and max of External data interpolated to the ARPS grid:'

      CALL a3dmax0(arpsgrid%x,1,nx,1,nx,1,1,1,1,1,1,1,1, amax,amin)
      IF(IAMROOT) WRITE(6,'(/1x,2(a,e13.6))')                           &
              'xmin      = ', amin,', xmax      =',amax

      CALL a3dmax0(arpsgrid%y,1,ny,1,ny,1,1,1,1,1,1,1,1, amax,amin)
      IF(IAMROOT) WRITE(6,'(1x,2(a,e13.6))')                            &
              'ymin      = ', amin,', ymax      =',amax

      CALL a3dmax0(arpsgrid%z,1,nz,1,nz,1,1,1,1,1,1,1,1, amax,amin)
      IF(IAMROOT) WRITE(6,'(1x,2(a,e13.6))')                            &
              'zmin      = ', amin,', zmax      =',amax

      CALL a3dmax0(arpsgrid%zp,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz,       &
                amax,amin)
      IF(IAMROOT) WRITE(6,'(1x,2(a,e13.6))')                            &
              'zpmin     = ', amin,', zpmax     =',amax

      CALL a3dmax0(arpsgrid%rhobar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1, &
                amax,amin)
      IF(IAMROOT) WRITE(6,'(1x,2(a,e13.6))')                            &
              'rhobarmin = ', amin,', rhobarmax =',amax

      CALL a3dmax0(arpsgrid%u,1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1,        &
                amax,amin)
      IF(IAMROOT) WRITE(6,'(1x,2(a,e13.6))')                            &
              'umin      = ', amin,', umax      =',amax

      CALL a3dmax0(arpsgrid%v,1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1,        &
                amax,amin)
      IF(IAMROOT) WRITE(6,'(1x,2(a,e13.6))')                            &
              'vmin      = ', amin,', vmax      =',amax

      CALL a3dmax0(arpsgrid%w,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz,        &
                amax,amin)
      IF(IAMROOT) WRITE(6,'(1x,2(a,e13.6))')                            &
              'wmin      = ', amin,', wmax      =',amax

      CALL a3dmax0(arpsgrid%ptbar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,  &
                amax,amin)
      IF(IAMROOT) WRITE(6,'(1x,2(a,e13.6))')                            &
              'ptbarmin  = ', amin,', ptbarmax  =',amax

      CALL a3dmax0(arpsgrid%ptprt,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,  &
                amax,amin)
      IF(IAMROOT) WRITE(6,'(1x,2(a,e13.6))')                            &
              'ptprtmin  = ', amin,', ptprtmax  =',amax

      CALL a3dmax0(arpsgrid%pbar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,   &
                amax,amin)
      IF(IAMROOT) WRITE(6,'(1x,2(a,e13.6))')                            &
              'pbarmin   = ', amin,', pbarmax   =',amax

      CALL a3dmax0(arpsgrid%pprt,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,   &
                amax,amin)
      IF(IAMROOT) WRITE(6,'(1x,2(a,e13.6))')                            &
              'pprtmin   = ', amin,', pprtmax   =',amax

      CALL a3dmax0(arpsgrid%qv,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,     &
                amax,amin)
      IF(IAMROOT) WRITE(6,'(1x,2(a,e13.6))')                            &
              'qvmin     = ', amin,', qvmax     =',amax

      CALL a3dmax0(arpsgrid%raing  ,1,nx,1,nx-1,1,ny,1,ny-1,            &
                   1,1,1,1,amax,amin)
      IF(IAMROOT) WRITE(6,'(1x,2(a,e13.6))')                            &
                   'raing_min = ', amin,', raing_max =',amax
      CALL a3dmax0(arpsgrid%rainc  ,1,nx,1,nx-1,1,ny,1,ny-1,            &
                   1,1,1,1,amax,amin)
      IF(IAMROOT) WRITE(6,'(1x,2(a,e13.6))')                            &
                   'rainc_min = ', amin,', rainc_max =',amax
      CALL a3dmax0(arpsgrid%t2m  ,1,nx,1,nx-1,1,ny,1,ny-1,              &
                   1,1,1,1,amax,amin)
      IF(IAMROOT) WRITE(6,'(1x,2(a,e13.6))')                            &
                   't2m_min   = ', amin,', t2m_max   =',amax
      CALL a3dmax0(arpsgrid%th2m  ,1,nx,1,nx-1,1,ny,1,ny-1,             &
                   1,1,1,1,amax,amin)
      IF(IAMROOT) WRITE(6,'(1x,2(a,e13.6))')                            &
                   'th2m_min  = ', amin,', th2m_max  =',amax
      CALL a3dmax0(arpsgrid%qv2m  ,1,nx,1,nx-1,1,ny,1,ny-1,             &
                   1,1,1,1,amax,amin)
      IF(IAMROOT) WRITE(6,'(1x,2(a,e13.6))')                            &
                   'qv2m_min  = ', amin,', qv2m_max  =',amax
      CALL a3dmax0(arpsgrid%u10m  ,1,nx,1,nx-1,1,ny,1,ny-1,             &
                   1,1,1,1,amax,amin)
      IF(IAMROOT) WRITE(6,'(1x,2(a,e13.6))')                            &
                   'u10m_min  = ', amin,', u10m_max  =',amax
      CALL a3dmax0(arpsgrid%v10m  ,1,nx,1,nx-1,1,ny,1,ny-1,             &
                   1,1,1,1,amax,amin)
      IF(IAMROOT) WRITE(6,'(1x,2(a,e13.6))')                            &
                     'v10m_min  = ', amin,', v10m_max  =',amax
      CALL a3dmax0(arpsgrid%wspd10max  ,1,nx,1,nx-1,1,ny,1,ny-1,        &
                   1,1,1,1,amax,amin)
      IF(IAMROOT) WRITE(6,'(1x,2(a,e13.6))')                            &
                     'wspd10max_min  = ', amin,', wspd10max_max  =',amax
      CALL a3dmax0(arpsgrid%up_heli_max  ,1,nx,1,nx-1,1,ny,1,ny-1,      &
                   1,1,1,1,amax,amin)
      IF(IAMROOT) WRITE(6,'(1x,2(a,e13.6))')                            &
                     'up_heli_max_min  = ', amin,', up_heli_max_max  =',amax
      CALL a3dmax0(arpsgrid%refd_max  ,1,nx,1,nx-1,1,ny,1,ny-1,         &
                   1,1,1,1,amax,amin)
      IF(IAMROOT) WRITE(6,'(1x,2(a,e13.6))')                            &
                     'refd_max_min  = ', amin,', refd_max_max  =',amax
      CALL a3dmax0(arpsgrid%grpl_max  ,1,nx,1,nx-1,1,ny,1,ny-1,         &
                   1,1,1,1,amax,amin)
      IF(IAMROOT) WRITE(6,'(1x,2(a,e13.6))')                            &
                     'grpl_max_min  = ', amin,', grpl_max_max  =',amax

      IF(sfcout == 1) THEN

        CALL a3dmax0(arpsgrid%tsoil,1,nx,1,nx-1,1,ny,1,ny-1,            &
                  1,arpsgrid%nzsoil,1,arpsgrid%nzsoil,amax,amin)
        IF(IAMROOT) WRITE(6,'(1x,2(a,e13.6))')                          &
                'tsoilmin  = ', amin,', tsoilmax  =',amax

        CALL a3dmax0(arpsgrid%qsoil,1,nx,1,nx-1,1,ny,1,ny-1,            &
                  1,arpsgrid%nzsoil,1,arpsgrid%nzsoil,amax,amin)
        IF(IAMROOT) WRITE(6,'(1x,2(a,e13.6))')                          &
                'qsoilmin  = ', amin,', qsoilmax  =',amax

        CALL a3dmax0(arpsgrid%wetcanp,1,nx,1,nx-1,1,ny,1,ny-1,          &
                  1,1,1,1,amax,amin)
        IF(IAMROOT) WRITE(6,'(1x,2(a,e13.6))')                          &
                  'wetcanpmin= ', amin,', wetcanpmax=',amax

      END IF

    END IF

!-----------------------------------------------------------------------
!
!  Data dump of the model grid and base state arrays:
!
!  First find a unique name basdmpfn(1:lbasdmpf) for the grid and
!  base state array dump file
!
!-----------------------------------------------------------------------
!
    CALL gtdmpfn(runname(1:lfnkey),dirname,ldirnam,curtim,hdmpfmt,      &
                 mgrid,nestgrd, hdmpfn, ldmpf)

!    ldirnam=LEN(dirname)
!    CALL strlnth( dirname, ldirnam)

    IF ( hdmpfmt == 5 ) THEN
      WRITE (6,'(a/a)')                                                 &
          'Program wrf2arps does not support Savi3D format.',           &
          'Reset to binary format, hdmpfmt=1.'
      hdmpfmt = 1
    END IF

    ALLOCATE(tem1(nx,ny,nz), STAT = istatus)
    tem1(:,:,:) = 0.
    ALLOCATE(tem2(nx,ny,nz), STAT = istatus)
    ALLOCATE(tem3(nx,ny,nz), STAT = istatus)
    ALLOCATE(tem4(nx,ny,nz), STAT = istatus)

    IF ( hdmpfmt == 9 ) GO TO 450

    IF(arpsgrid%first_time .AND. ABS(curtim) < 1.0E-5) THEN
      ! First run this program AND it is the initial forecast time

      CALL gtbasfn(runname(1:lfnkey),dirname,ldirnam,hdmpfmt,           &
                   mgrid,nestgrd,basdmpfn,lbasdmpf)

      IF(myproc == 0)        &
        PRINT*,'Writing base state history dump ',basdmpfn(1:lbasdmpf)

      mstout0 = mstout
      rainout0= rainout
      prcout0 = prcout
      trbout0 = trbout

      grdbas =1
      mstout =1
      rainout=0
      prcout =0
      trbout =0

      ! blocking inserted for ordering i/o for message passing
      DO i = 0, nprocs-1 ,dumpstride
        IF(myproc >= i .AND. myproc <= i+dumpstride-1)THEN

          CALL dtadump(arpsgrid%nx,arpsgrid%ny,arpsgrid%nz,             &
                   arpsgrid%nzsoil,arpsgrid%nstyps,                     &
                   hdmpfmt,iniotfu,basdmpfn(1:lbasdmpf),grdbas,filcmprs,&
                   arpsgrid%u,arpsgrid%v,arpsgrid%w,                    &
                   arpsgrid%ptprt,arpsgrid%pprt,arpsgrid%qv,            &
                   arpsgrid%qscalar,arpsgrid%tke,tem1,tem1,             &
                   arpsgrid%ubar,arpsgrid%vbar,arpsgrid%wbar,           &
                   arpsgrid%ptbar,arpsgrid%pbar,arpsgrid%rhobar,arpsgrid%qvbar, &
                   arpsgrid%x,arpsgrid%y,arpsgrid%z,arpsgrid%zp,arpsgrid%zpsoil,&
                   arpsgrid%soiltyp,arpsgrid%stypfrct,arpsgrid%vegtyp,  &
                   arpsgrid%lai,arpsgrid%roufns,arpsgrid%veg,           &
                   arpsgrid%tsoil,arpsgrid%qsoil,arpsgrid%wetcanp,      &
                   arpsgrid%snowdpth,arpsgrid%raing,arpsgrid%rainc,tem1,&
                   tem1,tem1,tem1,                                      &
                   tem1,tem1,                                           &
                   tem1,tem1,tem1,tem1,                                 &
                   tem2,tem3,tem4)
        END IF
        IF (mp_opt > 0) CALL mpbarrier
      END DO

      mstout = mstout0
      rainout= rainout0
      prcout = prcout0
      trbout = trbout0

    END IF  ! first_time

    450 CONTINUE

    IF(myproc == 0) THEN
      PRINT*,'curtim = ', curtim
      WRITE(6,'(1x,a,a)') 'History data dump in file ',hdmpfn(1:ldmpf)
    END IF

    grdbas=0

    ! blocking inserted for ordering i/o for message passing
    DO i = 0, nprocs-1, dumpstride
      IF ( myproc >= i .AND. myproc <= i+dumpstride-1 ) THEN
        CALL dtadump(arpsgrid%nx,arpsgrid%ny,arpsgrid%nz,               &
                 arpsgrid%nzsoil,arpsgrid%nstyps,                       &
                 hdmpfmt,iniotfu,hdmpfn(1:ldmpf),grdbas,filcmprs,       &
                 arpsgrid%u,arpsgrid%v,arpsgrid%w,                      &
                 arpsgrid%ptprt,arpsgrid%pprt,arpsgrid%qv,              &
                 arpsgrid%qscalar,arpsgrid%tke,tem1,tem1,               &
                 arpsgrid%ubar,arpsgrid%vbar,tem1,                      &
                 arpsgrid%ptbar,arpsgrid%pbar,arpsgrid%rhobar,arpsgrid%qvbar, &
                 arpsgrid%x,arpsgrid%y,arpsgrid%z,arpsgrid%zp,arpsgrid%zpsoil,&
                 arpsgrid%soiltyp,arpsgrid%stypfrct,arpsgrid%vegtyp,    &
                 arpsgrid%lai,arpsgrid%roufns,arpsgrid%veg,             &
                 arpsgrid%tsoil,arpsgrid%qsoil,arpsgrid%wetcanp,        &
                 arpsgrid%snowdpth,arpsgrid%raing,arpsgrid%rainc,tem1,  &
                 tem1,tem1,tem1,                                        &
                 tem1,tem1,                                             &
                 tem1,tem1,tem1,tem1,                                   &
                 tem2,tem3,tem4)
      END IF
      IF (mp_opt > 0) CALL mpbarrier
    END DO

    DEALLOCATE(tem1)
    DEALLOCATE(tem2, tem3, tem4)

    RETURN
  END SUBROUTINE arpsgrid_histout

  SUBROUTINE arpsgrid_diagout( latdiag, londiag, dbglvl, istatus )
!#######################################################################

    IMPLICIT NONE
    REAL,    INTENT(IN)  :: latdiag, londiag
    INTEGER, INTENT(IN)  :: dbglvl
    INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

    INTEGER :: idiag, jdiag
    REAL    :: xdiag, ydiag
    REAL    :: dd   , latd, lond
    REAL    :: dmin

    REAL    :: ppasc, pmb, theta, tc, smix
    REAL    :: e, bige, alge, tdc, dir, spd

    INTEGER :: i,j,k

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0
!
!-----------------------------------------------------------------------
!
!  Test code, for diagnostic testing.
!  Find x,y of diagnostic sounding location in ARPS grid.
!
!-----------------------------------------------------------------------
!
    CALL lltoxy(1,1,latdiag,londiag,xdiag,ydiag)

    IF (xdiag > arpsgrid%xs(2) .AND. xdiag < arpsgrid%xs(arpsgrid%nx-2) .AND.  &
        ydiag > arpsgrid%ys(2) .AND. ydiag < arpsgrid%ys(arpsgrid%ny-2) ) THEN

      dmin=((xdiag-arpsgrid%xs(1))*(xdiag-arpsgrid%xs(1))+              &
            (ydiag-arpsgrid%ys(1))*(ydiag-arpsgrid%ys(1)))
      idiag=1
      jdiag=1

      DO j=2,arpsgrid%ny-2
        DO i=2,arpsgrid%nx-2
          dd=((xdiag-arpsgrid%xs(i))*(xdiag-arpsgrid%xs(i))+            &
              (ydiag-arpsgrid%ys(j))*(ydiag-arpsgrid%ys(j)))
          IF(dd < dmin) THEN
            dmin=dd
            idiag=i
            jdiag=j
          END IF
        END DO
      END DO
      CALL xytoll(1,1,arpsgrid%xs(idiag),arpsgrid%ys(jdiag),            &
                  latd,lond)
      WRITE(6,'(a,f10.4,f10.4,/a,i5,i5,a,f10.4,f10.4)')                 &
            ' Nearest ARPS pt to diagnostic lat,lon: ',                 &
            latdiag,londiag,                                            &
            ' Diagnostic i,j: ',                                        &
            idiag,jdiag,' lat,lon= ',latd,lond
      WRITE(6,'(///a,/2x,a)')                                           &
          ' ARPS extracted sounding at idiag,jdiag',                    &
          'k   pres   hgt   temp   theta   dewp     u     v     dir    spd'
!
!-----------------------------------------------------------------------
!
!  Convert units of ARPS data and write as a sounding.
!
!-----------------------------------------------------------------------
!
      DO k=arpsgrid%nz-2,1,-1
        ppasc=arpsgrid%pbar(idiag,jdiag,k)+arpsgrid%pprt(idiag,jdiag,k)
        pmb=.01*(arpsgrid%pbar(idiag,jdiag,k)+arpsgrid%pprt(idiag,jdiag,k))
        theta=arpsgrid%ptbar(idiag,jdiag,k)+arpsgrid%ptprt(idiag,jdiag,k)
        tc=(theta*((ppasc/p0)**rddcp))-273.15
        IF( arpsgrid%qv(idiag,jdiag,k) > 0.) THEN
          smix=arpsgrid%qv(idiag,jdiag,k)/(1.-arpsgrid%qv(idiag,jdiag,k))
          e=(pmb*smix)/(0.62197 + smix)
          bige=e/( 1.001 + ( (pmb - 100.) / 900.) * 0.0034)
          alge = ALOG(bige/6.112)
          tdc = (alge * 243.5) / (17.67 - alge)
        ELSE
          tdc = tc-30.
        END IF

        CALL uvrotdd(1,1,londiag,arpsgrid%u(idiag,jdiag,k),             &
                     arpsgrid%v(idiag,jdiag,k),dir,spd)

        WRITE(6,'(i4,f6.0,f9.0,f7.1,f7.1,f7.1,f7.1,f7.1,f7.1,f7.1)')    &
                k,pmb,                                                  &
                arpsgrid%zps(idiag,jdiag,k),                            &
                tc,theta,tdc,                                           &
                arpsgrid%u(idiag,jdiag,k),                              &
                arpsgrid%v(idiag,jdiag,k),                              &
                dir,spd
      END DO
    ELSE
      WRITE(6,'(1x,a)') 'WARNING: Diagnostic lat,lon is outside of the ARPS domain.'
    END IF

    RETURN
  END SUBROUTINE arpsgrid_diagout

  SUBROUTINE arpsgrid_2dout( IAMROOT, curtime, mp_physics,               &
             outheader, gemoutheader,icape,iaccu,iascii,i2dfmt,igempak,  &
             ilite,iltgci,icrtm,isatid,chbgn,chend,icitm,user_emis,      &
             ibeg_offset,iend_offset,jbeg_offset,jend_offset,            &
             dbglvl, istatus )

!#######################################################################

    IMPLICIT NONE
    LOGICAL, INTENT(IN)  :: IAMROOT
    INTEGER, INTENT(IN)  :: mp_physics
    CHARACTER(LEN=80), INTENT(IN) :: outheader, gemoutheader
    INTEGER, INTENT(IN)  :: icape, iaccu, iascii, i2dfmt, igempak
    INTEGER, INTENT(IN)  :: ilite,iltgci,icrtm,isatid,chbgn,chend,icitm,user_emis
    INTEGER, INTENT(IN)  :: ibeg_offset,iend_offset,jbeg_offset,jend_offset
    REAL,    INTENT(IN)  :: curtime

    INTEGER, INTENT(IN)  :: dbglvl
    INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
    INCLUDE 'globcst.inc'

    INTEGER :: lenbin, lengem

    REAL, ALLOCATABLE :: tem1(:,:,:),tem2(:,:,:),tem3(:,:,:)
    REAL, ALLOCATABLE :: tem4(:,:,:),tem5(:,:,:),tem6(:,:,:)

    REAL, ALLOCATABLE :: temh2d(:,:)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    CALL gtdmpfn(runname(1:lfnkey),dirname,ldirnam,curtim,hdmpfmt,      &
                 mgrid,nestgrd, hdmpfn, ldmpf)

    ALLOCATE(tem1(arpsgrid%nx,arpsgrid%ny,arpsgrid%nz), STAT = istatus)
    ALLOCATE(tem2(arpsgrid%nx,arpsgrid%ny,arpsgrid%nz), STAT = istatus)
    ALLOCATE(tem3(arpsgrid%nx,arpsgrid%ny,arpsgrid%nz), STAT = istatus)
    ALLOCATE(tem4(arpsgrid%nx,arpsgrid%ny,arpsgrid%nz), STAT = istatus)
    ALLOCATE(tem5(arpsgrid%nx,arpsgrid%ny,arpsgrid%nz), STAT = istatus)
    ALLOCATE(tem6(arpsgrid%nx,arpsgrid%ny,arpsgrid%nz), STAT = istatus)

    ALLOCATE(temh2d (arpsgrid%nx,arpsgrid%ny), STAT = istatus)
    temh2d(:,:) = ARPS_MISSING

    arpsgrid%u  = arpsgrid%u  - arpsgrid%ubar
    arpsgrid%v  = arpsgrid%v  - arpsgrid%vbar
    arpsgrid%qv = arpsgrid%qv - arpsgrid%qvbar

    lenbin = len_trim(outheader)
    lengem = len_trim(gemoutheader)
    IF(hdmpfn(ldmpf-2:ldmpf-2) == '.') ldmpf = ldmpf - 3

    CALL postcore(arpsgrid%nx,arpsgrid%ny,arpsgrid%nz,                  &
            arpsgrid%nzsoil,arpsgrid%nstyps,curtime,mp_physics,         &
            hdmpfn(1:ldmpf),ldmpf,outheader(1:lenbin),                  &
            lenbin,gemoutheader(1:lengem),lengem,                       &
            icape,iaccu,iascii,i2dfmt,igempak,                          &
            ilite,iltgci,icrtm,isatid,chbgn,chend,icitm,user_emis,      &
            ibeg_offset,iend_offset,jbeg_offset,jend_offset,            &
            arpsgrid%x,arpsgrid%y,arpsgrid%z,arpsgrid%zp,               &
            arpsgrid%u, arpsgrid%v, arpsgrid%w, arpsgrid%ptprt,         &
            arpsgrid%pprt, arpsgrid%qv, arpsgrid%qscalar,               &
            arpsgrid%ubar, arpsgrid%vbar, arpsgrid%wbar, arpsgrid%ptbar,&
            arpsgrid%pbar, arpsgrid%rhobar, arpsgrid%qvbar,             &
            arpsgrid%soiltyp,arpsgrid%stypfrct,arpsgrid%vegtyp,         &
            arpsgrid%lai,arpsgrid%roufns,arpsgrid%veg,                  &
            arpsgrid%tsoil,arpsgrid%qsoil,arpsgrid%wetcanp,arpsgrid%snowdpth, &
            arpsgrid%raing,arpsgrid%rainc,temh2d,temh2d,temh2d,         &
            arpsgrid%t2m,arpsgrid%th2m,                                 &
            arpsgrid%qv2m,arpsgrid%u10m,arpsgrid%v10m,arpsgrid%raddn,   &
            arpsgrid%wspd10max,arpsgrid%w_up_max,arpsgrid%w_dn_max,     &
            arpsgrid%refd_max,arpsgrid%up_heli_max,arpsgrid%grpl_max,   &
            temh2d,temh2d,temh2d,                                       &
            temh2d,temh2d,temh2d,temh2d,temh2d,                         &
            temh2d,temh2d,temh2d,temh2d,                                &
            tem1,tem2,tem3,tem4,tem5,tem6)

    DEALLOCATE( tem1, tem2, tem3, tem4, tem5, tem6 )
    DEALLOCATE( temh2d )

    RETURN
  END SUBROUTINE arpsgrid_2dout

  SUBROUTINE arpsgrid_read_soiltyp_vegtyp( filename, filefmt, dbglvl, istatus )

  !#####################################################################

    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN)  :: filename
    INTEGER, INTENT(IN)  :: filefmt
    INTEGER, INTENT(IN)  :: dbglvl
    INTEGER, INTENT(OUT) :: istatus

  !---------------------------------------------------------------------

    INCLUDE 'mp.inc'

  !---------------------------------------------------------------------

    CHARACTER(LEN=256) :: tmpfile

    INTEGER :: sd_id

    REAL, ALLOCATABLE :: int3d(:,:,:), int2d(:,:)

    LOGICAL :: iread
    INTEGER :: nxdta, nydta

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    istatus = 0

    iread  = .FALSE.
    IF (filefmt == 3) THEN
      IF (myproc == 0) iread = .TRUE.
      nxdta = arpsgrid%nxlg
      nydta = arpsgrid%nylg
      tmpfile = filename
    ELSE IF (filefmt == 103) THEN
      iread = .TRUE.
      nxdta = arpsgrid%nx
      nydta = arpsgrid%ny
      CALL gtsplitfn(filename,1,1,loc_x,loc_y,1,1, &
                         0,0,1,0,tmpfile,istatus)
    ELSE
      WRITE(*,'(1x,a,i0)') 'ERROR: unsupported file format = ',filefmt
      istatus = -1
      RETURN
    END IF

    ALLOCATE(int3d(nxdta,nydta,arpsgrid%nstyps), STAT = istatus)
    CALL check_alloc_status(istatus, "arpsgrid_read_soiltyp_vegtyp:int3d")

    ALLOCATE(int2d(nxdta,nydta),                 STAT = istatus)
    CALL check_alloc_status(istatus, "arpsgrid_read_soiltyp_vegtyp:int2d")

    IF (iread) THEN
      CALL hdfopen(TRIM(tmpfile),1,sd_id)
      IF (sd_id < 0) THEN
        WRITE (6,*) "arpsgrid_read_soiltyp_vegtyp: ERROR opening ",     &
                    trim(tmpfile)," for reading."
        GO TO 110
      END IF

      CALL hdfrd3di(sd_id,"soiltyp",nxdta,nydta,arpsgrid%nstyps,int3d,istatus)
      IF (istatus /= 0) GO TO 110

      CALL hdfrd2di(sd_id,"vegtyp",nxdta,nydta,int2d,istatus)
      IF (istatus /= 0) GO TO 110

      CALL hdfclose(sd_id,istatus)
      IF (istatus == 0) THEN
        WRITE(*,'(1x,2a)')                                              &
          "arpsgrid_read_soiltyp_vegtyp: Successfully read ", trim(tmpfile)
      ELSE
        WRITE(*,'(1x,a,i0,2a)')                                         &
          "arpsgrid_read_soiltyp_vegtyp: ERROR (status=", istatus, ") closing ", trim(tmpfile)
      END IF
    END IF

    110 CONTINUE

    CALL mpupdatei(istatus,1)
    IF (istatus /= 0) RETURN

    IF (filefmt > 100) THEN
      arpsgrid%soiltyp(:,:,:) = int3d(:,:,:)
      arpsgrid%vegtyp(:,:)    = int2d(:,:)
    ELSE
      CALL mpisplit3di(int3d,arpsgrid%nx,arpsgrid%ny,arpsgrid%nstyps,arpsgrid%soiltyp)
      CALL mpisplit3di(int2d,arpsgrid%nx,arpsgrid%ny,1,              arpsgrid%vegtyp)
    END IF

!    print *, "soiltype - ", MAXVAL(arpsgrid%soiltyp), MAXVAL(arpsgrid%vegtyp)

    DEALLOCATE(int3d,int2d)

    RETURN
  END SUBROUTINE arpsgrid_read_soiltyp_vegtyp

END MODULE module_output_arpsgrid
