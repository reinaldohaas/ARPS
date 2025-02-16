SUBROUTINE getcoamps(nx_ext, ny_ext, nz_ext, nzsoil_ext,                &
           dir_extd,extdinit,extdfcst,                                  &
           iproj_ext,scale_ext,                                         &
           trlon_ext,latnot_ext,x0_ext,y0_ext,                          &
           lat_ext,lon_ext,latu_ext,lonu_ext,latv_ext,lonv_ext,         &
           p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                      &
           qc_ext,qr_ext,qi_ext,qs_ext,qh_ext,                          &
           tsoil_ext, qsoil_ext,wetcanp_ext,                            &
           snowdpth_ext,trn_ext,psfc_ext,                               &
           t_2m_ext,qv_2m_ext,u_10m_ext,v_10m_ext,rain_ext,             &
           istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    This programs reads coamps 4.2 and later sigma level output
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    dir_extd      directory where the coamps files are located
!    extdinit      Initialized time in mm-dd-yyyy:hh:mm:ss format
!    extdfcst      Forecast hour in HHH:MM:SS format
!
!  OUTPUT:
!
!    iproj_ext     Map projection number of external data
!    scale_ext     Scale factor of external data
!    trlon_ext     True longitude of external data (degrees E)
!    latnot_ext(2) True latitude(s) of external data (degrees N)
!    x0_ext        x coordinate of origin of external data
!    y0_ext        y coordinate of origin of external data
!    lat_ext       latitude of external data points (degrees N)
!    lon_ext       longitude of external data points (degrees E)
!
!    p_ext         pressure (Pascal)
!    hgt_ext       height (m)
!    t_ext         temperature (K)
!    tsfc_ext      ground/sea-surface temperature (K)
!    tsoil_ext     Deep soil temperature (K) (in deep 1 m layer)
!    wetsfc_ext    Surface soil moisture in the top 1 cm layer
!                                              (fraction,0--1.0))
!    wetdp_ext     Deep soil moisture in the deep 1 m layer
!    wetcanp_ext   Canopy water amount
!    qv_ext        specific humidity(mixing ratio used) (kg/kg)
!    qc_ext        Cloud water mixing ratio (kg/kg)
!    qr_ext        Rain  water mixing ratio (kg/kg)
!    qi_ext        Ice         mixing ratio (kg/kg)
!    qs_ext        Snow        mixing ratio (kg/kg)
!    qh_ext        Hail        mixing ratio (kg/kg)
!    u_ext         u wind component (m/s)
!    v_ext         v wind component (m/s)
!    istatus       status indicator
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER,            INTENT(IN) :: nx_ext, ny_ext, nz_ext, nzsoil_ext

  CHARACTER (LEN=*) , INTENT(IN) :: dir_extd
  CHARACTER (LEN=19), INTENT(IN) :: extdinit
  CHARACTER (LEN=9) , INTENT(IN) :: extdfcst

!
!  Output external grid variables
!
  INTEGER, INTENT(OUT) :: iproj_ext
  REAL   , INTENT(OUT) :: scale_ext,trlon_ext
  REAL   , INTENT(OUT) :: latnot_ext(2)
  REAL   , INTENT(OUT) :: x0_ext,y0_ext

!
!  Output external variable arrays
!
  REAL,    INTENT(OUT) :: lat_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: lon_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: latu_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: lonu_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: latv_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: lonv_ext(nx_ext,ny_ext)

  REAL,    INTENT(OUT) :: p_ext(nx_ext,ny_ext,nz_ext)     ! (Pa)
  REAL,    INTENT(OUT) :: hgt_ext(nx_ext,ny_ext,nz_ext)   ! (m)
  REAL,    INTENT(OUT) :: t_ext(nx_ext,ny_ext,nz_ext)     ! (K)
  REAL,    INTENT(OUT) :: qv_ext(nx_ext,ny_ext,nz_ext)    ! (kg/kg)
  REAL,    INTENT(OUT) :: u_ext(nx_ext,ny_ext,nz_ext)     ! (m/s)
  REAL,    INTENT(OUT) :: v_ext(nx_ext,ny_ext,nz_ext)     ! (m/s)
  REAL,    INTENT(OUT) :: qc_ext(nx_ext,ny_ext,nz_ext)    ! Cloud H2O mixing ratio (kg/kg)
  REAL,    INTENT(OUT) :: qr_ext(nx_ext,ny_ext,nz_ext)    ! Rain  H2O mixing ratio (kg/kg)
  REAL,    INTENT(OUT) :: qi_ext(nx_ext,ny_ext,nz_ext)    ! Ice   H2O mixing ratio (kg/kg)
  REAL,    INTENT(OUT) :: qs_ext(nx_ext,ny_ext,nz_ext)    ! Snow  H2O mixing ratio (kg/kg)
  REAL,    INTENT(OUT) :: qh_ext(nx_ext,ny_ext,nz_ext)    ! Hail  H2O mixing ratio (kg/kg)

  REAL,    INTENT(OUT) :: tsoil_ext (nx_ext,ny_ext,nzsoil_ext)      ! soil temperature (K)
  REAL,    INTENT(OUT) :: qsoil_ext (nx_ext,ny_ext,nzsoil_ext)      ! soil moisture (0-1)
  REAL,    INTENT(OUT) :: wetcanp_ext(nx_ext,ny_ext)      ! Canopy water amount (0-1)

  REAL,    INTENT(OUT) :: snowdpth_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: trn_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: psfc_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: t_2m_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: qv_2m_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: u_10m_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: v_10m_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: rain_ext(nx_ext,ny_ext)

  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'mp.inc'

  REAL, PARAMETER :: gamma  = 0.0068         ! 6.8 degrees per km
  REAL, PARAMETER :: pconst = (g/rd)

  CHARACTER(LEN=6) :: fldnam
  CHARACTER(LEN=3) :: cltyp
  CHARACTER(LEN=1) :: cfluid

  INTEGER          :: inest, m, n
  REAL             :: rlev1, rlev2

  CHARACTER(LEN=10) :: cdtg
  INTEGER           :: ihour, iminute, isecond
  INTEGER           :: iyear, imon, iday, ihr, imin, isec

  INTEGER :: lend, np

  INTEGER :: kka
  INTEGER :: i,j,k,kk
  INTEGER :: nxy_ext, nxyz_ext, nxys_ext

!
!  Original grid variables
!
  INTEGER :: iproj
  REAL :: rscale,trlon,x0,y0
  REAL :: latnot(2)

  REAL :: sigmwa      ! Height of the atmosphere
  REAL :: swlat_ext, swlon_ext

  REAL :: dx_ext, dy_ext
  REAL :: p_2m_ext, deltaz
  LOGICAL :: psfc_good
  !
  ! Working arrays
  !
  REAL, ALLOCATABLE :: x_ext(:),y_ext(:)
  REAL, ALLOCATABLE :: xu_ext(:),yv_ext(:)

  REAL, ALLOCATABLE :: rdata(:)
  REAL, ALLOCATABLE :: sigmma(:)
  REAL, ALLOCATABLE :: dum2d(:,:), dum3d(:,:,:)
  REAL, ALLOCATABLE :: utmp(:,:), vtmp(:,:)

  REAL :: f_qvsatl
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (myproc == 0) THEN
    READ (extdinit,'(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)')                &
                    iyear,  imon, iday,  ihr, imin, isec

    WRITE(cdtg,'(I4.4,3I2.2)') iyear, imon,iday,ihr

    inest  = 1
    cfluid = 'a'
!
!***********************************************************************
! Read flat header file
!***********************************************************************
!
    lend = 2000
    ALLOCATE(rdata(lend), STAT = istatus)
    rdata = 0.0

    fldnam = 'datahd'
    cltyp  = 'sfc'
    rlev1  = 0.0
    rlev2  = 0.0

    m = lend
    n = 1
    ihour   = 0
    iminute  = 0
    isecond = 0

    CALL readff(rdata,lend,fldnam,inest,ihour,iminute,isec,cdtg,        &
                cfluid,cltyp,rlev1,rlev2,dir_extd,m,n,.FALSE.,istatus)

    IF (istatus /= 0) RETURN

!-----------------------------------------------------------------------
!
! COAMPS: map projection. 1: mercator,2: lambert,
!                         3: polar, 4:cartesian, 5: spherical
!
!-----------------------------------------------------------------------

    iproj = INT(rdata(3))

    SELECT CASE (iproj)
    CASE (1)
      iproj_ext = 3
    CASE (2)
      iproj_ext = 2
    CASE (3)
      iproj_ext = 1
    CASE (4)
      iproj_ext = 0
    CASE (5)
      iproj_ext = 5
    CASE DEFAULT
      PRINT *, 'can not handle this projection', iproj
      istatus = -1
      RETURN
    END SELECT

    latnot_ext(1) = rdata(4)
    latnot_ext(2) = rdata(5)
    trlon_ext = rdata(6)
    IF (trlon_ext > 180.0) trlon_ext = trlon_ext - 360.0
    scale_ext = 1.0
    dx_ext = rdata(9)
    dy_ext = rdata(10)

    np = 30 + (inest-1)*30
    swlat_ext = rdata(np+9)
    swlon_ext = rdata(np+10)
    IF (swlon_ext > 180.0) swlon_ext = swlon_ext - 360.0

!-----------------------------------------------------------------------
!
! Sigmma (m)
!
!-----------------------------------------------------------------------
    ALLOCATE(sigmma(nz_ext), STAT = istatus)

    kka = nz_ext-1
    DO k = 1,kka    ! note that nz_ext is kka+1
      kk = nz_ext - k
      sigmma(kk) = rdata(800+k)
    END DO

    sigmwa = sigmma(kka) + rdata(501)/2.0
    sigmma(nz_ext) = sigmma(kka) + rdata(501)

    DEALLOCATE(rdata)

!-----------------------------------------------------------------------
!
!  Get the lat,lon of the COAMPS grid points
!
!-----------------------------------------------------------------------
!
    ALLOCATE(x_ext(nx_ext),  STAT = istatus)
    ALLOCATE(y_ext(ny_ext),  STAT = istatus)
    ALLOCATE(xu_ext(nx_ext), STAT = istatus)
    ALLOCATE(yv_ext(ny_ext), STAT = istatus)

    CALL getmapr(iproj,rscale,latnot,trlon,x0,y0)
    CALL setmapr(iproj_ext,scale_ext,latnot_ext,trlon_ext)
    CALL lltoxy(1,1,swlat_ext,swlon_ext,x0_ext,y0_ext)

    DO i=1,nx_ext
      x_ext(i)=x0_ext+(i-1)*dx_ext
    END DO
    DO j=1,ny_ext
      y_ext(j)=y0_ext+(j-1)*dy_ext
    END DO

    CALL xytoll(nx_ext,ny_ext,x_ext, y_ext, lat_ext,lon_ext)

    xu_ext(:) = x_ext(:) + 0.5*dx_ext
    yv_ext(:) = y_ext(:) + 0.5*dy_ext

    CALL xytoll(nx_ext,ny_ext,xu_ext,y_ext, latu_ext,lonu_ext)
    CALL xytoll(nx_ext,ny_ext,x_ext, yv_ext,latv_ext,lonv_ext)

    DEALLOCATE(x_ext, y_ext, xu_ext, yv_ext)
!
!-----------------------------------------------------------------------
! Read in 2D static files
!-----------------------------------------------------------------------
!
    lend = nx_ext*ny_ext
    m = nx_ext
    n = ny_ext

    fldnam = 'terrht'
    cltyp  = 'sfc'
    rlev1  = 0.0
    rlev2  = 0.0

    trn_ext = -999.9
    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Reading terrain height ...'
    CALL readff(trn_ext,lend,fldnam,inest,ihour,iminute,isecond,cdtg,     &
                cfluid,cltyp,rlev1,rlev2,dir_extd,m,n,.TRUE.,             &
                istatus)

    DO k = 1,nz_ext
      DO j = 1, ny_ext
        DO i = 1, nx_ext
          hgt_ext(i,j,k) = sigmma(k)*(sigmwa-trn_ext(i,j))/sigmwa+trn_ext(i,j)
        END DO
      END DO
    END DO

    ALLOCATE(dum2d(nx_ext,ny_ext), STAT = istatus)

    tsoil_ext = -999.9

    fldnam = 'soltmp'
    cltyp  = 'sfc'
    rlev1  = 0.0
    rlev2  = 0.0

    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Reading soltmp ...'
    CALL readff(dum2d,lend,fldnam,inest,ihour,iminute,isecond,cdtg,       &
                cfluid,cltyp,rlev1,rlev2,dir_extd,m,n,.TRUE.,             &
                istatus)

    IF (istatus == 0) THEN
      tsoil_ext(:,:,2) = dum2d
    END IF

!***********************************************************************
!
! time-dependent variables starts here
!
!***********************************************************************

    READ(extdfcst,'(i3,1x,i2,1x,i2)') ihour,iminute,isecond

!-----------------------------------------------------------------------
!
! Read 2D flat files
!
!-----------------------------------------------------------------------

    lend = nx_ext*ny_ext
    m = nx_ext
    n = ny_ext

    fldnam = 'trpres'
    cltyp  = 'sfc'
    rlev1  = 0.0
    rlev2  = 0.0

    psfc_ext = -999.9
    psfc_good = .FALSE.
    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Reading trpres ...'
    CALL readff(psfc_ext,lend,fldnam,inest,ihour,iminute,isecond,cdtg,    &
                cfluid,cltyp,rlev1,rlev2,dir_extd,m,n,.TRUE.,             &
                istatus)
    IF (istatus == 0) THEN  ! terrain pressure is good
      psfc_good = .TRUE.
      psfc_ext = psfc_ext * 100.0
    END IF

    fldnam = 'airtmp'
    cltyp  = 'zht'
    rlev1  = 2.0
    rlev2  = 0.0

    t_2m_ext = -999.9
    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Reading airtmp at 2m ...'
    CALL readff(t_2m_ext,lend,fldnam,inest,ihour,iminute,isecond,cdtg,    &
                cfluid,cltyp,rlev1,rlev2,dir_extd,m,n,.TRUE.,             &
                istatus)

    fldnam = 'dwptdp'  ! Dew point temperature
    cltyp  = 'zht'
    rlev1  = 2.0
    rlev2  = 0.0

    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Reading dwptdp at 2m ...'
    CALL readff(dum2d,lend,fldnam,inest,ihour,iminute,isecond,cdtg,       &
                cfluid,cltyp,rlev1,rlev2,dir_extd,nx_ext,ny_ext,.TRUE.,   &
                istatus)

    qv_2m_ext(:,:) = -999.9
    IF (istatus == 0 .AND. psfc_good ) THEN

      DO j = 1, ny_ext
        DO i = 1, nx_ext
          dum2d (i,j) = dum2d(i,j)+273.15   !?
          p_2m_ext = psfc_ext(i,j) * EXP( pconst*(-2.0)/t_2m_ext(i,j) )
          qv_2m_ext(i,j) = f_qvsatl( p_2m_ext, dum2d(i,j) )
        END DO
      END DO

    END IF

    fldnam = 'uuwind'
    cltyp  = 'zht'
    rlev1  = 10.0
    rlev2  = 0.0

    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Reading uuwind at 10m ...'
    CALL readff(u_10m_ext,lend,fldnam,inest,ihour,iminute,isecond,cdtg,   &
                cfluid,cltyp,rlev1,rlev2,dir_extd,m,n,.TRUE.,             &
                istatus)

    fldnam = 'vvwind'
    cltyp  = 'zht'
    rlev1  = 10.0
    rlev2  = 0.0

    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Reading vvwind at 10m ...'
    CALL readff(v_10m_ext,lend,fldnam,inest,ihour,iminute,isecond,cdtg,   &
                cfluid,cltyp,rlev1,rlev2,dir_extd,m,n,.TRUE.,             &
                istatus)


    fldnam = 'snowdp'
    cltyp  = 'sfc'
    rlev1  = 0.0
    rlev2  = 0.0

    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Reading snowdp ...'
    CALL readff(snowdpth_ext,lend,fldnam,inest,ihour,iminute,isecond,cdtg,&
                cfluid,cltyp,rlev1,rlev2,dir_extd,m,n,.TRUE.,             &
                istatus)

    fldnam = 'grdwet'
    cltyp  = 'sfc'
    rlev1  = 0.0
    rlev2  = 0.0

    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Reading grdwet ...'
    CALL readff(wetcanp_ext,lend,fldnam,inest,ihour,iminute,isecond,cdtg, &
                cfluid,cltyp,rlev1,rlev2,dir_extd,m,n,.TRUE.,             &
                istatus)

    fldnam = 'conpac'
    cltyp  = 'sfc'
    rlev1  = 0.0
    rlev2  = 0.0

    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Reading conpac ...'
    CALL readff(dum2d,lend,fldnam,inest,ihour,iminute,isecond,cdtg,       &
                cfluid,cltyp,rlev1,rlev2,dir_extd,m,n,.TRUE.,             &
                istatus)

    fldnam = 'stapac'
    cltyp  = 'sfc'
    rlev1  = 0.0
    rlev2  = 0.0

    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Reading stapac ...'
    CALL readff(rain_ext,lend,fldnam,inest,ihour,iminute,isecond,cdtg,    &
                cfluid,cltyp,rlev1,rlev2,dir_extd,m,n,.TRUE.,             &
                istatus)

    rain_ext(:,:) = rain_ext(:,:) + dum2d(:,:)

    fldnam = 'grdtmp'
    cltyp  = 'sfc'
    rlev1  = 0.0
    rlev2  = 0.0

    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Reading grdtmp ...'
    CALL readff(tsoil_ext(:,:,1),lend,fldnam,inest,ihour,iminute,isecond,cdtg,  &
                cfluid,cltyp,rlev1,rlev2,dir_extd,m,n,.TRUE.,             &
                istatus)

    fldnam = 'wvapor'
    cltyp  = 'sfc'
    rlev1  = 0.0
    rlev2  = 0.0

    qsoil_ext = -999.9
    IF (lvldbg > 0) WRITE(*,'(1x,a)') 'Reading wvapor ...'
    CALL readff(qsoil_ext(:,:,1),lend,fldnam,inest,ihour,iminute,isecond,cdtg, &
                cfluid,cltyp,rlev1,rlev2,dir_extd,m,n,.TRUE.,             &
                istatus)

!-----------------------------------------------------------------------
!
! Read 3D flat files
!
!-----------------------------------------------------------------------

    ALLOCATE(dum3d(nx_ext,ny_ext,kka), STAT = istatus)

    lend = nx_ext*ny_ext*(nz_ext-1)
    m = nx_ext
    n = ny_ext


    fldnam = 'ttlprs'
    cltyp  = 'sig'
    rlev2  = sigmma(1)
    rlev1  = sigmma(kka)

    CALL readff(dum3d,lend,fldnam,inest,ihour,iminute,isecond,cdtg,       &
                cfluid,cltyp,rlev1,rlev2,dir_extd,m,n,.TRUE.,             &
                istatus)

    DO k = 1,kka
      DO j = 1, ny_ext
        DO i = 1, nx_ext
          kk = nz_ext - k
          p_ext(i,j,kk) = dum3d(i,j,k)*100.0
        END DO
      END DO
    END DO

    fldnam = 'pottmp'

    CALL readff(dum3d,lend,fldnam,inest,ihour,iminute,isecond,cdtg,       &
                cfluid,cltyp,rlev1,rlev2,dir_extd,m,n,.TRUE.,             &
                istatus)

    DO k = 1,kka
      DO j = 1, ny_ext
        DO i = 1, nx_ext
          kk = nz_ext - k
          t_ext(i,j,kk) = dum3d(i,j,k)*((p_ext(i,j,kk)/p0)**rddcp)
        END DO
      END DO
    END DO

    ! extra-polation, constant lapse rate
    DO j = 1, ny_ext
      DO i = 1, nx_ext
        t_ext(i,j,nz_ext) = t_ext(i,j,kka)-gamma*(hgt_ext(i,j,nz_ext)-hgt_ext(i,j,kka))*0.001
      END DO
    END DO
    !
    ! Hydrostatic interpolation
    !
    deltaz = hgt_ext(1,1,kka)-hgt_ext(1,1,nz_ext)
    DO j = 1, ny_ext
      DO i = 1, nx_ext
        p_ext(i,j,nz_ext) = p_ext(i,j,kka) * EXP ( pconst*deltaz          &
                                                   / (t_ext(i,j,kka) )    &
                                                   )
      END DO
    END DO

!-----------------------------------------------------------------------
!
! Hydrometological variables
!
!-----------------------------------------------------------------------

    fldnam = 'wvapor'

    CALL readff(dum3d,lend,fldnam,inest,ihour,iminute,isecond,cdtg,       &
                cfluid,cltyp,rlev1,rlev2,dir_extd,m,n,.TRUE.,             &
                istatus)
    DO k = 1,kka
      DO j = 1, ny_ext
        DO i = 1, nx_ext
          kk = nz_ext - k
          qv_ext(i,j,kk) = dum3d(i,j,k)
        END DO
      END DO
    END DO

    fldnam = 'cldmix'

    CALL readff(dum3d,lend,fldnam,inest,ihour,iminute,isecond,cdtg,       &
                cfluid,cltyp,rlev1,rlev2,dir_extd,m,n,.TRUE.,             &
                istatus)
    DO k = 1,kka
      DO j = 1, ny_ext
        DO i = 1, nx_ext
          kk = nz_ext - k
          qc_ext(i,j,kk) = dum3d(i,j,k)
        END DO
      END DO
    END DO

    fldnam = 'ranmix'

    CALL readff(dum3d,lend,fldnam,inest,ihour,iminute,isecond,cdtg,       &
                cfluid,cltyp,rlev1,rlev2,dir_extd,m,n,.TRUE.,             &
                istatus)
    DO k = 1,kka
      DO j = 1, ny_ext
        DO i = 1, nx_ext
          kk = nz_ext - k
          qr_ext(i,j,kk) = dum3d(i,j,k)
        END DO
      END DO
    END DO

    fldnam = 'icemix'

    CALL readff(dum3d,lend,fldnam,inest,ihour,iminute,isecond,cdtg,       &
                cfluid,cltyp,rlev1,rlev2,dir_extd,m,n,.TRUE.,             &
                istatus)
    DO k = 1,kka
      DO j = 1, ny_ext
        DO i = 1, nx_ext
          kk = nz_ext - k
          qi_ext(i,j,kk) = dum3d(i,j,k)
        END DO
      END DO
    END DO

    fldnam = 'snomix'

    CALL readff(dum3d,lend,fldnam,inest,ihour,iminute,isecond,cdtg,       &
                cfluid,cltyp,rlev1,rlev2,dir_extd,m,n,.TRUE.,             &
                istatus)
    DO k = 1,kka
      DO j = 1, ny_ext
        DO i = 1, nx_ext
          kk = nz_ext - k
          qs_ext(i,j,kk) = dum3d(i,j,k)
        END DO
      END DO
    END DO

    fldnam = 'grpmix'

    CALL readff(dum3d,lend,fldnam,inest,ihour,iminute,isecond,cdtg,       &
                cfluid,cltyp,rlev1,rlev2,dir_extd,m,n,.TRUE.,             &
                istatus)
    DO k = 1,kka
      DO j = 1, ny_ext
        DO i = 1, nx_ext
          kk = nz_ext - k
          qh_ext(i,j,kk) = dum3d(i,j,k)
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
! Read u-component wind
!
!-----------------------------------------------------------------------

    fldnam = 'uuwind'
    cltyp  = 'sig'
    rlev2  = sigmma(1)
    rlev1  = sigmma(kka)

    CALL readff(dum3d,lend,fldnam,inest,ihour,iminute,isecond,cdtg,       &
                cfluid,cltyp,rlev1,rlev2,dir_extd,m,n,.TRUE.,             &
                istatus)
    DO k = 1,kka
      DO j = 1, ny_ext
        DO i = 1, nx_ext-1
          kk = nz_ext - k
          u_ext(i,j,kk) = dum3d(i,j,k)
        END DO
      END DO
    END DO
    CALL edgfill(u_ext,1,nx_ext,1,nx_ext-1,1,ny_ext,1,ny_ext,             &
                 1,nz_ext,1,nz_ext-1)

!-----------------------------------------------------------------------
!
! Read v-component wind
!
!-----------------------------------------------------------------------

    fldnam = 'vvwind'
    cltyp  = 'sig'
    rlev2  = sigmma(1)
    rlev1  = sigmma(kka)

    CALL readff(dum3d,lend,fldnam,inest,ihour,iminute,isecond,cdtg,       &
                cfluid,cltyp,rlev1,rlev2,dir_extd,m,n,.TRUE.,             &
                istatus)
    DO k = 1,kka
      DO j = 1, ny_ext-1
        DO i = 1, nx_ext
          kk = nz_ext - k
          v_ext(i,j,kk) = dum3d(i,j,k)
        END DO
      END DO
    END DO
    CALL edgfill(v_ext,1,nx_ext,1,nx_ext,1,ny_ext,1,ny_ext-1,             &
                 1,nz_ext,1,nz_ext-1)
!
!-----------------------------------------------------------------------
!
!  Rotate winds to be relative to true north.
!  The COAMPS data are saved as grid-relative.
!
!-----------------------------------------------------------------------
!
    ALLOCATE(utmp(nx_ext,ny_ext), STAT = istatus)
    ALLOCATE(vtmp(nx_ext,ny_ext), STAT = istatus)

    DO k=1,nz_ext
      CALL uvmptoe(nx_ext,ny_ext,u_ext(1,1,k),v_ext(1,1,k),               &
                   lon_ext,utmp,vtmp)
      u_ext(:,:,k) = utmp(:,:)
      v_ext(:,:,k) = vtmp(:,:)
    END DO

    DEALLOCATE(utmp)
    DEALLOCATE(vtmp)
!
!-----------------------------------------------------------------------
!
!  Reset map projection to previous values
!
!-----------------------------------------------------------------------
!
    DEALLOCATE(sigmma)
    DEALLOCATE(dum2d, dum3d)

    CALL setmapr(iproj,rscale,latnot,trlon)
    CALL setorig(1,x0,y0)
  END IF
!
!-----------------------------------------------------------------------
!
!  Set good status
!
!-----------------------------------------------------------------------
!
  nxy_ext  = nx_ext*ny_ext
  nxyz_ext = nx_ext*ny_ext*nz_ext
  nxys_ext = nx_ext*ny_ext*nzsoil_ext

  CALL mpupdatei(iproj_ext,1)
  CALL mpupdater(latnot_ext,2)
  CALL mpupdater(trlon_ext,1)
  CALL mpupdater(scale_ext,1)
  CALL mpupdater(x0_ext,1)
  CALL mpupdater(y0_ext,1)
  CALL mpupdater(dx_ext,1)
  CALL mpupdater(dy_ext,1)
  CALL mpupdater(lat_ext,nxy_ext)
  CALL mpupdater(lon_ext,nxy_ext)
  CALL mpupdater(latv_ext,nxy_ext)
  CALL mpupdater(lonv_ext,nxy_ext)
  CALL mpupdater(latu_ext,nxy_ext)
  CALL mpupdater(lonu_ext,nxy_ext)
  CALL mpupdater(p_ext,nxyz_ext)
  CALL mpupdater(t_ext,nxyz_ext)
  CALL mpupdater(hgt_ext,nxyz_ext)
  CALL mpupdater(qv_ext,nxyz_ext)
  CALL mpupdater(u_ext,nxyz_ext)
  CALL mpupdater(v_ext,nxyz_ext)
  CALL mpupdater(qc_ext,nxyz_ext)
  CALL mpupdater(qr_ext,nxyz_ext)
  CALL mpupdater(qs_ext,nxyz_ext)
  CALL mpupdater(qh_ext,nxyz_ext)

  CALL mpupdater(tsoil_ext,nxys_ext)
  CALL mpupdater(qsoil_ext,nxys_ext)

  CALL mpupdater(wetcanp_ext,nxy_ext)
  CALL mpupdater(snowdpth_ext,nxy_ext)
  CALL mpupdater(trn_ext,nxy_ext)
  CALL mpupdater(psfc_ext,nxy_ext)
  CALL mpupdater(t_2m_ext,nxy_ext)
  CALL mpupdater(qv_2m_ext,nxy_ext)
  CALL mpupdater(u_10m_ext,nxy_ext)
  CALL mpupdater(v_10m_ext,nxy_ext)
  CALL mpupdater(rain_ext,nxy_ext)

  istatus=1

  RETURN
END SUBROUTINE getcoamps

SUBROUTINE getcoampsdims(dir_extd,extdinit,extdfcst,                    &
                         nx_ext, ny_ext, nz_ext,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    This programs reads coamps 4.2 sigma level output header information
!  file to extract dimension size of the external model domain.
!
!-----------------------------------------------------------------------
! MODIFICATION HISTORY:
!
! ORIGINAL AUTHOR (Y. Wang, 01/03/2011)
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    dir_extd      directory where the coamps files are located
!    extdinit      Initialized time in mm-dd-yyyy:hh:mm:ss format
!    extdfcst      Forecast hour in HHH:MM:SS format
!
!  OUTPUT:
!
!    nx_ext     Map projection number of external data
!    ny_ext     Scale factor of external data
!    nz_ext     True longitude of external data (degrees E)
!    istatus    status indicator
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  CHARACTER (LEN=*) , INTENT(IN) :: dir_extd
  CHARACTER (LEN=19), INTENT(IN) :: extdinit
  CHARACTER (LEN=9) , INTENT(IN) :: extdfcst

  INTEGER, INTENT(OUT) :: nx_ext, ny_ext, nz_ext
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  CHARACTER(LEN=6) :: fldnam
  CHARACTER(LEN=3) :: cltyp
  CHARACTER(LEN=1) :: cfluid
  INTEGER          :: inest, m, n
  REAL             :: rlev1, rlev2

  CHARACTER(LEN=10) :: cdtg
  INTEGER           :: ihour, minute, isecond
  INTEGER           :: iyear, imon, iday, ihr, imin, isec

  INTEGER :: lend, np

!
!  COAMPS input fields
!
  REAL, ALLOCATABLE :: rdata(:)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  READ (extdinit,'(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)')                  &
                  iyear,  imon, iday,  ihr, imin, isec

  lend = 2000
  ALLOCATE(rdata(lend), STAT = istatus)
  rdata = 0.0
!
!***********************************************************************
! Read flat header file
!***********************************************************************
!
  fldnam = 'datahd'
  cltyp  = 'sfc'
  rlev1  = 0.0
  rlev2  = 0.0
  inest  = 1
  cfluid = 'a'
  m = lend
  n = 1
  ihour   = 0
  minute  = 0
  isecond = 0

  WRITE(cdtg,'(I4.4,3I2.2)') iyear, imon,iday,ihr

  CALL readff(rdata,lend,fldnam,inest,ihour,minute,isec,cdtg,           &
              cfluid,cltyp,rlev1,rlev2,dir_extd,m,n,.FALSE.,istatus)

  IF (istatus /= 0) RETURN

!-----------------------------------------------------------------------
!
!  Get the COAMPS grid configurations
!
!-----------------------------------------------------------------------
!
  np = 30+(inest-1)*30

  nx_ext = INT( rdata(np+0) )
  ny_ext = INT( rdata(np+1) )
  nz_ext = INT( rdata(2) )

  nz_ext = nz_ext + 1   ! kka+1

!-----------------------------------------------------------------------
!
!  Finialize this subroutine
!
!-----------------------------------------------------------------------
!
  DEALLOCATE(rdata)

  RETURN
END SUBROUTINE getcoampsdims

SUBROUTINE readff(d,lend,fldnam,inest,ihour,minute,isec,cdtg,           &
                  cfluid,cltyp,rlev1,rlev2,dsetff,m,n,lwritu,           &
                  istatus)
!
!***********************************************************************
!          subroutine readff - read a flat file from noraps/coamps
!***********************************************************************
!
  IMPLICIT NONE

!***********************************************************************
!         dimension statements
!***********************************************************************
!
  INTEGER,            INTENT(IN)  :: lend
  REAL,               INTENT(OUT) :: d(lend)

  CHARACTER(LEN=6),   INTENT(IN) :: fldnam
  INTEGER,            INTENT(IN) :: inest
  INTEGER,            INTENT(IN) :: ihour, minute, isec
  CHARACTER(LEN=10),  INTENT(IN) :: cdtg
  CHARACTER(LEN=1  ), INTENT(IN) :: cfluid
  CHARACTER(LEN=3  ), INTENT(IN) :: cltyp
  REAL,               INTENT(IN) :: rlev1, rlev2
  CHARACTER(LEN=256), INTENT(IN) :: dsetff
  LOGICAL,            INTENT(IN) :: lwritu
  INTEGER,            INTENT(IN) :: m, n

  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------

  INTEGER            :: ilev1, ilev2

  CHARACTER(LEN=256) :: cfile
  CHARACTER(LEN=64 ) :: ctemp
  CHARACTER(LEN=7  ) :: outtyp

  INTEGER :: i

  INTEGER :: llen, lsetux
  INTEGER :: funit

  REAL*4, ALLOCATABLE :: temp(:)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = -1

  ilev1=rlev1
  ilev2=rlev2
!
!***********************************************************************
!          build filename
!***********************************************************************
!
  IF(fldnam == 'datahd') THEN
     outtyp='infofld'
  ELSE
     outtyp='fcstfld'
  END IF
  WRITE (ctemp,830) fldnam,cltyp,ilev1,ilev2,inest,cfluid,              &
                        m,n,cdtg,ihour,minute,isec,outtyp
  830 FORMAT (a6,'_',a3,'_',i6.6,'_',i6.6,'_',i1,a1,                    &
              i4.4,'x',i4.4,'_',a10,'_',i4.4,2i2.2,'_',a7)
  llen=lend

  lsetux=LEN_TRIM(dsetff)
  cfile=dsetff(1:lsetux)//ctemp
!
!***********************************************************************
! Reading new 64 characters direct access file
!***********************************************************************

  CALL getunit( funit )

  IF (lwritu) THEN

    ALLOCATE(temp(lend), STAT = istatus)

    llen=4*llen
    OPEN (unit=funit,file=TRIM(cfile),form='unformatted',access='direct',  &
          recl=llen,status='old',err=900)
    READ (funit,rec=1) temp
    DO i=1,lend
      d(i)=temp(i)
    END DO
    CLOSE (funit)

    istatus = 0
  ELSE
    OPEN (unit=funit,file=TRIM(cfile),form='formatted',access='sequential',&
          status='old',err=905)
    READ  (funit,'(5e13.6)',err=910) (d(i),i=1,lend)
    CLOSE (funit)
    istatus = 0
  END IF

!-----------------------------------------------------------------------
!
! Finialize this subroutine
!
!-----------------------------------------------------------------------

  GOTO 1000

  900 CONTINUE
  PRINT *, 'Cannot open file: ',TRIM(cfile),' in readff.'
  istatus = -1
  GOTO 1000

  905 CONTINUE
  PRINT '(4x,2a)', 'ERROR opening flat file in readff, file: ', TRIM(cfile)
  istatus = -2
  GOTO 1000

  910 CONTINUE
  PRINT *,'read error'
  istatus = -3
  CLOSE (funit)

  1000 CONTINUE

  CALL retunit( funit )

  IF (ALLOCATED(temp)) DEALLOCATE(temp)

  RETURN
END SUBROUTINE readff
