!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE HIS2VERGRID               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE his2vergrid(nx,ny,nz,nzsoil,x,y,z,zp,zpsoil,hterain,      &
               uprt,vprt,wprt,ptprt,pprt,qvprt,                      &
               qscalar,tke,kmh,kmv,                                  &
               ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,               &
               nstyps,soiltyp,stypfrct,vegtyp,lai,roufns,            &
               veg,tsoil,qsoil,wetcanp,snowdpth,                     &
               raing,rainc,prcrate,radfrc,radsw,rnflx,               &
               radswnet,radlwin,                                     &
               usflx,vsflx,ptsflx,qvsflx,                            &
               hinfmt,fgrdbasfn,hisfile,nhisfile,                    &
               ibgn,iend,jbgn,jend,                                  &
               mRefopt,mReflist,mRefdir,mRefrunname,                 &
               reffmt,refopt,refcomp,nreflvl,mxscorelvl,threshold)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Performs read in the ARPS' history dump data and
!        the grided observations, and then calclulate
!        ETS and bias
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Hu
!
!  Original Coding: 01/29/2006
!
!  MODIFICATION HISTORY:
!
!  12/01/2011 Y. Wang
!  Upgraded to support ARPS multi-moment microphysics scheme and WRF
!  formular for reflectivity calculation.
!
!-----------------------------------------------------------------------
!
!  DATA ARRAYS READ IN:
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       z coordinate of grid points in physical space (m)
!    zpsoil   z coordinate of soil model
!
!    uprt     Perturbation x component of velocity (m/s)
!    vprt     Perturbation y component of velocity (m/s)
!    wprt     Perturbation z component of velocity (m/s)
!
!    ptprt    Perturbation potential temperature (K)
!    pprt     Perturbation pressure (Pascal)
!
!    qscalar  Hydrometeors scalars
!
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!
!    ubar     Base state x velocity component (m/s)
!    vbar     Base state y velocity component (m/s)
!    wbar     Base state z velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhobar   Base state density (kg/m**3)
!    qvbar    Base state water vapor mixing ratio (kg/kg)
!
!    soiltyp  Soil type
!    vegtyp   Vegetation type
!    lai      Leaf Area Index
!    roufns   Surface roughness
!    veg      Vegetation fraction
!
!    tsoil    soil temperature (K)
!    qsoil    Soil moisture
!    wetcanp  Canopy water amount
!
!    raing    Grid supersaturation rain
!    rainc    Cumulus convective rain
!    prcrate  Precipitation rates
!
!    radfrc   Radiation forcing (K/s)
!    radsw    Solar radiation reaching the surface
!    rnflx    Net radiation flux absorbed by surface
!    radswnet Net shortwave radiation, SWin - SWout
!    radlwin  Incoming longwave radiation
!
!    usflx    Surface flux of u-momentum (kg/(m*s**2))
!    vsflx    Surface flux of v-momentum (kg/(m*s**2))
!    ptsflx   Surface heat flux (K*kg/(m**2 * s ))
!    qvsflx   Surface moisture flux of (kg/(m**2 * s))
!
!  CALCULATED DATA ARRAYS:
!
!    su       Sounding x component of velocity (m/s)
!    sv       Sounding y component of velocity (m/s)
!    stheta   Sounding potential temperature (K)
!    spres    Sounding pressure (mb)
!    stemp    Sounding temperature (C)
!    sdewp    Sounding dew-point (C)
!    sdrct    Sounding wind direction (degrees)
!    ssped    Sounding wind speed (m/s)
!    shght    Sounding height (m)
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!
!   Temporary arrays are defined and used differently by each
!   subroutine.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'vericst.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Arrays to be read in:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx, ny, nz, nzsoil
  REAL :: x     (nx)     ! The x-coord. of the physical and
                                     ! computational grid. Defined at u-point.
  REAL :: y     (ny)     ! The y-coord. of the physical and
                                     ! computational grid. Defined at v-point.
  REAL :: z     (nz)     ! The z-coord. of the computational grid.
                                     ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz) ! The physical height coordinate defined at
                                     ! w-point of the staggered grid.
  REAL :: zpsoil(nx,ny,nzsoil) ! The physical height coordinate of soil model

  REAL :: hterain (nx,ny)       ! Terrain height

  REAL :: uprt   (nx,ny,nz)    ! Perturbation u-velocity (m/s)
  REAL :: vprt   (nx,ny,nz)    ! Perturbation v-velocity (m/s)
  REAL :: wprt   (nx,ny,nz)    ! Perturbation w-velocity (m/s)
  REAL :: ptprt  (nx,ny,nz)    ! Perturbation potential temperature (K)
  REAL :: pprt   (nx,ny,nz)    ! Perturbation pressure (Pascal)
  REAL :: qvprt  (nx,ny,nz)    ! Perturbation water vapor specific
                                         ! humidity (kg/kg)
  REAL :: qscalar(nx,ny,nz,nscalar)

  REAL :: tke    (nx,ny,nz)    ! Turbulent Kinetic Energy ((m/s)**2)
  REAL :: kmh    (nx,ny,nz)    ! Horizontal turb. mixing coef. for
                                         ! momentum. ( m**2/s )
  REAL :: kmv    (nx,ny,nz)    ! Vertical turb. mixing coef. for
                                         ! momentum. ( m**2/s )

  REAL :: ubar   (nx,ny,nz)    ! Base state u-velocity (m/s)
  REAL :: vbar   (nx,ny,nz)    ! Base state v-velocity (m/s)
  REAL :: wbar   (nx,ny,nz)    ! Base state w-velocity (m/s)
  REAL :: ptbar  (nx,ny,nz)    ! Base state potential temperature (K)
  REAL :: pbar   (nx,ny,nz)    ! Base state pressure (Pascal)
  REAL :: rhobar (nx,ny,nz)    ! Base state air density (kg/m**3)
  REAL :: qvbar  (nx,ny,nz)    ! Base state water vapor specific

  INTEGER :: nstyps

  INTEGER :: soiltyp (nx,ny,nstyps)       ! Soil type
  REAL :: stypfrct(nx,ny,nstyps)          ! Soil type
  INTEGER :: vegtyp  (nx,ny)         ! Vegetation type
  REAL :: lai     (nx,ny)            ! Leaf Area Index
  REAL :: roufns  (nx,ny)            ! Surface roughness
  REAL :: veg     (nx,ny)            ! Vegetation fraction

  REAL :: tsoil  (nx,ny,nzsoil,0:nstyps) ! Soil temperature (K)
  REAL :: qsoil  (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny,0:nstyps)       ! Canopy water amount
  REAL :: snowdpth(nx,ny)        ! Snow depth (m)

  REAL :: raing(nx,ny)         ! Grid supersaturation rain
  REAL :: rainc(nx,ny)         ! Cumulus convective rain
  REAL :: prcrate(nx,ny,4)     ! precipitation rate (kg/(m**2*s))
                               !   prcrate(1,1,1) = total precip. rate
                               !   prcrate(1,1,2) = grid scale precip. rate
                               !   prcrate(1,1,3) = cumulative precip. rate
                               !   prcrate(1,1,4) = microphysics precip. rate

  REAL :: radfrc(nx,ny,nz)     ! Radiation forcing (K/s)
  REAL :: radsw (nx,ny)        ! Solar radiation reaching the surface
  REAL :: rnflx (nx,ny)        ! Net radiation flux absorbed by surface
  REAL :: radswnet(nx,ny)      ! Net shortwave radiation
  REAL :: radlwin(nx,ny)       ! Incoming longwave radiation

  REAL :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s)

  INTEGER :: nhisfile

  INTEGER :: ibgn,iend
  INTEGER :: jbgn,jend
!
!
!-----------------------------------------------------------------------
!
!  Map variables, declared here for use in subroutines...
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: xs(:)
  REAL, ALLOCATABLE :: ys(:)
  REAL, ALLOCATABLE :: xmap(:)
  REAL, ALLOCATABLE :: ymap(:)
  REAL, ALLOCATABLE :: latgr(:,:)
  REAL, ALLOCATABLE :: longr(:,:)
!
!-----------------------------------------------------------------------
!
!  Work Arrays
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: tem1(:,:,:)
  REAL, ALLOCATABLE :: tem2(:,:,:)
  REAL, ALLOCATABLE :: tem3(:,:,:)
  REAL, ALLOCATABLE :: tem4(:,:)

  REAL, ALLOCATABLE :: qv(:,:,:)             !zhao kun add qv
  REAL, ALLOCATABLE :: pres(:,:,:)
  REAL, ALLOCATABLE :: temp(:,:,:)

  REAL, ALLOCATABLE :: mosaic3d(:,:,:)
  REAL, ALLOCATABLE :: threatscore(:,:),biasscore(:,:)
!
!-----------------------------------------------------------------------
!
!  options
!
!-----------------------------------------------------------------------
!
  INTEGER :: mRefopt
  CHARACTER (LEN=256) :: mReflist
  CHARACTER (LEN=256) :: mRefdir
  CHARACTER (LEN=256) :: mRefrunname

  INTEGER :: reffmt,refopt
  INTEGER :: refcomp, wrfopt
  INTEGER :: nreflvl

  INTEGER :: mxscorelvl
  REAL :: threshold(20)
!-----------------------------------------------------------------------
!
!  Misc. internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: len1, istatus
  REAL :: time
  INTEGER :: i,j,k
  INTEGER :: ireturn,lengbf,lenfil,nchin,lengbf00Z,lengbf12Z
  INTEGER :: hinfmt
  INTEGER :: nhisfile_max
  PARAMETER (nhisfile_max=200)
  CHARACTER (LEN=256) :: filename
  CHARACTER (LEN=256) :: fgrdbasfn
  CHARACTER (LEN=256) :: grdbasfn
  CHARACTER (LEN=256) :: hisfile(nhisfile_max)
  CHARACTER (LEN=256) :: ftemp
  CHARACTER (LEN=256) :: fmosaic(nhisfile_max)
  REAL :: histime(nhisfile_max)

  CHARACTER (LEN=256) :: refdmpfn
  LOGICAL :: needmRef
  INTEGER :: nf
  CHARACTER (LEN=8) :: the_date
  CHARACTER (LEN=10) :: the_time
  CHARACTER (LEN=5) :: the_zone
  INTEGER :: the_values(8)

  INTEGER :: sd_id,istat
  INTEGER :: nxm,nym,nzm
  INTEGER :: nchanl

  LOGICAL :: fexist

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  nstyp = nstyps

  IF( ibgn == 0 ) ibgn = 2
  IF( iend == 0 ) iend = nx-2
  IF( jbgn == 0 ) jbgn = 2
  IF( jend == 0 ) jend = ny-2

  ibgn=max(2,ibgn)
  iend=min((nx-2),iend)
  jbgn=max(2,jbgn)
  jend=min((ny-2),jend)

!  CALL DATE_AND_TIME(the_date,the_time,the_zone,the_values)

!  DO nf=1,nhisfile
!    lenfil =256
!    CALL slength( hisfile(nf), lenfil)
!    WRITE(6,'(1x,a,a)') 'History file is ',hisfile(nf)(1:lenfil)
!  END DO

  OPEN(12, file=trim(mReflist))
  DO nf=1,nhisfile
    READ(12,*) fmosaic(nf)
    WRITE(*,*) 'mosaic file (',nf,') is:',trim(fmosaic(nf))
  END DO
  CLOSE(12)
!
!  ALLOCATE
!
  ALLOCATE( mosaic3d(nx,ny,nz), STAT=istatus )
  CALL check_alloc_status(istatus, "mosaic3d")
  mosaic3d=0.0
  ALLOCATE( threatscore(mxscorelvl,nhisfile), STAT=istatus )
  CALL check_alloc_status(istatus, "threatscore")
  threatscore=999.0
  ALLOCATE( biasscore(mxscorelvl,nhisfile), STAT=istatus )
  CALL check_alloc_status(istatus, "biasscore")
  biasscore=999.0

  ALLOCATE(xs (nx),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:xs")
  xs = 0.0
  ALLOCATE(ys (ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:ys")
  ys = 0.0
  ALLOCATE(tem1 (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem1")
  tem1 = 0.0

  ALLOCATE(qv (nx,ny,nz),STAT=istatus)             !zhaokun
  CALL check_alloc_status(istatus, "arps:qv")
  qv = 0.0

  ALLOCATE(tem2 (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem2")
  tem2 = 0.0
  ALLOCATE(tem3 (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem3")
  tem3 = 0.0
  ALLOCATE(tem4 (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:tem4")
  tem4 = 0.0
  ALLOCATE(xmap (nx),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:xmap")
  xmap = 0.0
  ALLOCATE(ymap (ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:ymap")
  ymap = 0.0
  ALLOCATE(latgr (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:latgr")
  latgr = 0.0
  ALLOCATE(longr (nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:longr")
  longr = 0.0

  ALLOCATE(pres (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:pres")
  ALLOCATE(temp (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "arps:temp")

!-----------------------------------------------------------------------
!
!  Loop through all history dumps
!
!-----------------------------------------------------------------------
!
  grdbasfn = fgrdbasfn
  lengbf=256
  CALL slength(grdbasfn,lengbf)
  IF ( mp_opt > 0 .AND. readsplit(FINDX_H) <= 0 ) THEN
    CALL gtsplitfn(fgrdbasfn,1,1,loc_x,loc_y,1,1,                      &
                   0,0,1,lvldbg,grdbasfn,ireturn)
    lengbf = LEN_TRIM(grdbasfn)
  END IF

  DO nf = 1, nhisfile

!
!-----------------------------------------------------------------------
!
!    Read all input data arrays
!
!-----------------------------------------------------------------------
!

!zhaokun added to compute accumulated rain at two successive time
!firstly read previous time rain
     IF( mRefopt == 2 .AND. nf>1) then
       filename=hisfile(nf-1)
       lenfil =256
       CALL slength( hisfile(nf-1), lenfil)
       CALL dtaread(nx,ny,nz,nzsoil,nstyps,                             &
                  hinfmt,nchin,grdbasfn(1:lengbf),lengbf,               &
                  filename(1:lenfil),lenfil,time,                       &
                  x,y,z,zp,zpsoil,uprt,vprt,wprt,ptprt,pprt,            &
                  qvprt,qscalar,tke,kmh,kmv,                            &
                  ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,               &
                  soiltyp,stypfrct,vegtyp,lai,roufns,veg,               &
                  tsoil,qsoil,wetcanp,snowdpth,                         &
                  raing,rainc,prcrate,                                  &
                  radfrc,radsw,rnflx,radswnet,radlwin,                  &
                  usflx,vsflx,ptsflx,qvsflx,                            &
                  ireturn, tem1,tem2,tem3)


       DO j=1,ny
         DO i=1,nx
!zh  aokun add to save precipitation
           tem4(i,j)= raing(i,j)+rainc(i,j)
         END DO
       END DO
     ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     filename=hisfile(nf)
     lenfil =256
     CALL slength( hisfile(nf), lenfil)
     IF ( mp_opt > 0 .AND. readsplit(FINDX_H) <= 0 ) THEN
       CALL gtsplitfn(hisfile(nf),1,1,loc_x,loc_y,1,1,                  &
                      0,0,1,lvldbg,filename,ireturn)
       lenfil = LEN_TRIM(filename)
     END IF
     IF ( myproc == 0 )                                                 &
       WRITE(6,'(1x,a,a)') 'History file is ',filename(1:lenfil)

     CALL dtaread(nx,ny,nz,nzsoil,nstyps,                               &
                hinfmt,nchin,grdbasfn(1:lengbf),lengbf,                 &
                filename(1:lenfil),lenfil,time,                         &
                x,y,z,zp,zpsoil,uprt,vprt,wprt,ptprt,pprt,              &
                qvprt,qscalar,tke,kmh,kmv,                              &
                ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,                 &
                soiltyp,stypfrct,vegtyp,lai,roufns,veg,                 &
                tsoil,qsoil,wetcanp,snowdpth,                           &
                raing,rainc,prcrate,                                    &
                radfrc,radsw,rnflx,radswnet,radlwin,                    &
                usflx,vsflx,ptsflx,qvsflx,                              &
                ireturn, tem1,tem2,tem3)

     DO j=1,ny
       DO i=1,nx
         hterain(i,j)=zp(i,j,2)
       END DO
     END DO

    !calculate accumulated rain at two successive time
    IF( mRefopt == 2.and. nf>1) then
      DO j=1,ny
        DO i=1,nx
          !zhaokun add to save precipitation
          tem4(i,j)= raing(i,j)+rainc(i,j)-tem4(i,j)
        END DO
      END DO
    ENDIF
!
!-----------------------------------------------------------------------
!
!    ireturn = 0 for a successful read
!
!-----------------------------------------------------------------------
!
    IF( ireturn == 0 ) THEN   ! successful read

     histime(nf)=time/60.0
     WRITE(6,'(a,f9.1)') ' History file time: ',histime(nf)
!
!-----------------------------------------------------------------------
! get predicted reflectivity
!
!-----------------------------------------------------------------------
      IF( mRefopt == 1 ) then
        !
        !  Calculate temperature (K) at ARPS grid points
        !
        DO k=1,nz
          DO j=1,ny
            DO i=1,nx
              tem2(i,j,k) = ptprt(i,j,k)+ ptbar(i,j,k)
              pres(i,j,k) = pbar(i,j,k )+ pprt (i,j,k)
              qv(i,j,k)   = MAX(0.0,qvprt(i,j,k)+qvbar(i,j,k))
            END DO
          END DO
        END DO
        CALL edgfill(pres,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)

        CALL temper (nx,ny,nz,tem2, pprt ,pbar,temp)
        CALL edgfill(temp,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)
        !
        !  Calculate reflectivity according to precipitation mixing ratio
        !   tem2 = 3d reflectiivty
        !   tem3 = composite reflectivity
        !
        tem2=0
       !tem3=0
       ! IF (refopt == 2) THEN
       !   CALL reflec_ferrier(nx,ny,nz, rhobar, qscalar, tem1, tem2)
       !   !CALL reflec_wrf(nx,ny,nz,qv,qr,qs,qh,1,tem3,tem1,tem2)
       ! ELSE
       !   CALL reflec(nx,ny,nz, rhobar, q, qs, qh, tem2)
       ! ENDIF

        mphyopt = refopt
        ! Milbrandt and Yau 3-moment scheme
        IF (refopt >= 8 .AND. refopt <= 12) THEN
          CALL reflec_MM(nx,ny,nz,rhobar,qscalar,temp,tem2)

        ! Lin,Schultz,Straka Lin,or WSM6 schemes
        ELSE IF (refopt >= 2 .AND. refopt <= 7) THEN
          CALL reflec_ferrier(nx,ny,nz, rhobar, qscalar, temp, tem2)

        ! Warm rain schemes
        ELSE IF (refopt == 1) THEN
          IF(P_QR > 0) THEN
            CALL reflec_wr(nx,ny,nz, rhobar, qscalar(:,:,:,P_QR),tem2)
          ELSE
            tem2 = 0.0
          END IF

        ELSE IF (refopt > 100 .AND. refopt < 200) THEN
          wrfopt = MOD(refopt,100)


          SELECT CASE (wrfopt)

          CASE (1,3)     ! Kessler, WSM 3-class
            print *,' Calling reflec_wrf ...'
            CALL reflec_wrf(nx,ny,nz,qv,                                      &
                 qscalar(:,:,:,P_QR),qscalar(:,:,:,P_QS),qscalar(:,:,:,P_QH), &
                 0,pres,temp,tem2)
          CASE (2,4,6) ! Lin, WSM 5- and 6- classes)
            print *,' Calling reflec_wrf ...'
            CALL reflec_wrf(nx,ny,nz,qv,                                      &
                 qscalar(:,:,:,P_QR),qscalar(:,:,:,P_QS),qscalar(:,:,:,P_QH), &
                 1,pres,temp,tem2)
                                                 ! mixed-phase
          CASE (5)     ! Ferrier microphysics scheme (as in WRF_POST)
            print *,' Calling reflec_ferrier_wrf ...'
            tem1(:,:,:) = 2.0
            CALL reflec_ferrier_wrf(nx,ny,nz,qv,                              &
                 qscalar(:,:,:,P_QC),qscalar(:,:,:,P_QR),qscalar(:,:,:,P_QS), &
                 pres,temp,tem2,tem3,tem1)
          CASE (8) ! Thompson microphysics scheme (old version)
            print *,' Calling CALREF9s ...'
            CALL CALREF9s(nx,ny,nz,                                           &
                 qscalar(:,:,:,P_QR),qscalar(:,:,:,P_QS),qscalar(:,:,:,P_QH), &
                 pres,temp,tem2)
          END SELECT

        ELSE
          print*,'Invalid microphysics option, reflectivity set to zero'
          tem2 = 0.0
        END IF

        tem3 = 0.0
        IF(refcomp == 1) CALL cal_rfc(nx, ny, nz, tem2, tem3)

      END IF
!
!-----------------------------------------------------------------------
!
!      get reflectivity observation
!
!-----------------------------------------------------------------------
!
      print *, ' mRefopt = ',mRefopt
      istatus=0
      IF( mRefopt == 1 ) THEN
        ftemp=trim(mRefdir)//trim(fmosaic(nf))
        WRITE(6,'(a,a)') ' Opening radar file ',trim(ftemp)
        INQUIRE(file=TRIM(ftemp),exist=fexist)

        IF(fexist) THEN

        IF (2==1) THEN
          nchanl=38
          OPEN(UNIT=nchanl,FILE=trim(ftemp),ERR=400,                    &
                  FORM='unformatted',STATUS='old')

          READ(nchanl) nxm,nym,nzm
          IF( (nxm .ne. nx) .or. (nym .ne. ny) .or. (nzm .ne. nz) ) then
            WRITE(*,*) 'Mosaic has a different dimension from the predictio!'
            WRITE(*,*) 'nx,ny,nz=',nx,ny,nz
            WRITE(*,*) 'nxm,nym,nzm=',nxm,nym,nzm
            STOP 123
          END IF
          READ(nchanl) mosaic3d
          CLOSE(nchanl)
        ELSE
          CALL readadcol_Mosaic(nx,ny,nz,reffmt,ftemp,mosaic3d,istatus)
        END IF
        ELSE
          istatus=-1
        END IF
      !zhkun noted for precipitaion ETS
      ELSE if(mRefopt==2) then
        ftemp=trim(mRefdir)//trim(fmosaic(nf))
        WRITE(6,'(a,a)') ' Opening radar file ',trim(ftemp)
        CALL read_Obsrain(nx,ny,nz,ftemp,mosaic3d)
      END IF
!
!-----------------------------------------------------------------------
!
!      calculate ETS
!
!-----------------------------------------------------------------------
      IF( istatus == 0 ) THEN
      IF( mRefopt == 1 ) then

         if(refcomp==1) then
         call ETS_Ref(nx,ny,nz,mosaic3d,tem3,ibgn,iend,jbgn,jend, &
                 refcomp,nreflvl,   &
                 mxscorelvl,threshold,threatscore(:,nf),biasscore(:,nf))
         else

         call ETS_Ref(nx,ny,nz,mosaic3d,tem2,ibgn,iend,jbgn,jend, &
                 refcomp,nreflvl,   &
                 mxscorelvl,threshold,threatscore(:,nf),biasscore(:,nf))
         endif
       else if(mRefopt==2) then
         do i=1,nx
         do j=1,ny
            tem3(i,j,1)=tem4(i,j)
         enddo
         enddo
         call ETS_Ref(nx,ny,nz,mosaic3d,tem3,ibgn,iend,jbgn,jend, &
              refcomp,nreflvl,   &
              mxscorelvl,threshold,threatscore(:,nf),biasscore(:,nf))

      ENDIF
!
! mRefdir,mRefrunname, mReflist
      needmRef=.false.
      IF (needmRef) THEN
        len1 =256
        CALL slength( mRefrunname, len1)
        refdmpfn=mRefrunname(1:len1)//'.ref'//                  &
            filename(lenfil-5:lenfil)
      END IF
      ELSE
        WRITE(6,'(a,a)') 'Error reading mosaic file:',TRIM(ftemp)
      END IF
!
!
!-- ---------------------------------------------------------------------
!
!       Extract model data for MOS.
!
!-- ---------------------------------------------------------------------
!
    ELSE
      WRITE(6,'(a)') ' Error reading data.  HIS2VERGRID ends'
      STOP
    END IF
!
!-----------------------------------------------------------------------
!
!  Go to next history dump
!
!-----------------------------------------------------------------------
!
  END DO  ! nf


!-----------------------------------------------------------------------
!
!  save ETS
!
!-----------------------------------------------------------------------
!
  refdmpfn=trim(mRefrunname)//'_ref.ETS'
  OPEN(12,file=trim(refdmpfn),form='formatted',status='unknown')
  WRITE(12,'(1x,a,2i8)') 'Verif Range i:',ibgn,iend
  WRITE(12,'(1x,a,2i8)') 'Verif Range j:',jbgn,jend
  WRITE(12,'(14x,a,56x,a)') 'Bias Score','Equitable threat score '
  WRITE(12,'(a,8f8.1,8f8.1)') '  Time(min)', &
                              (threshold(k),k=1,mxscorelvl), &
                              (threshold(k),k=1,mxscorelvl)
  DO nf=1,nhisfile
    WRITE(12,'(f9.1,8f8.4,8f8.4)')   histime(nf),&
              (biasscore(k,nf),k=1,mxscorelvl),  &
              (threatscore(k,nf), k=1,mxscorelvl)
  END DO
  CLOSE(12)

  RETURN

  400 continue
  write(*,*) 'Error in open Mosaic file!!!!'
  stop 123

END SUBROUTINE his2vergrid

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE TEMPER                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE temper ( nx,ny,nz,theta, ppert, pbar, t )

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    Using a version of Poisson's formula, calculate temperature.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Joe Bradley
!    12/05/91
!
!  MODIFICATIONS:
!    Modified by Ming Xue so that arrays are only defined at
!             one time level.
!    6/09/92  Added full documentation and phycst include file for
!             rddcp=Rd/Cp  (K. Brewster)
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    theta    Potential temperature (degrees Kelvin)
!    ppert    Perturbation pressure (Pascals)
!    pbar     Base state pressure (Pascals)
!
!  OUTPUT:
!
!    t        Temperature (degrees Kelvin)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx,ny,nz
!
  REAL :: theta(nx,ny,nz)      ! potential temperature (degrees Kelvin)
  REAL :: ppert(nx,ny,nz)      ! perturbation pressure (Pascals)
  REAL :: pbar (nx,ny,nz)      ! base state pressure (Pascals)
!
  REAL :: t    (nx,ny,nz)      ! temperature (degrees Kelvin)
!
!-----------------------------------------------------------------------
!
!  Include file
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
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
!  Calculate the temperature using Poisson's formula.
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1

        t(i,j,k) = theta(i,j,k) *                                       &
             (((ppert(i,j,k) + pbar(i,j,k)) / p0) ** rddcp)

      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE temper

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINES cal_rfc                  ######
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
!  PURPOSE:
!  Calculate rfc value.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Min Zou
!  3/2/98
!
!  MODIFICATION HISTORY:
!
!  Ming Xue (10/16/2001)
!  Now passing in precalculated reflectivity field instead of calculating
!  it inside.
!
!-----------------------------------------------------------------------
!

SUBROUTINE cal_rfc(nx, ny, nz, ref, refc)

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL, INTENT(IN ) :: ref (nx,ny,nz) ! Reflectivity
  REAL, INTENT(OUT) :: refc(nx,ny,nz) ! Composite reflectivity

  INTEGER :: i,j,k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO j=1,ny
    DO i=1,nx
      refc(i,j,1)= ref(i,j,1)
      DO k=2,nz-1
        refc(i,j,1) = MAX(refc(i,j,1),ref(i,j,k))
      END DO
    END DO
  END DO

  DO j=1,ny
    DO i=1,nx
      DO k=2,nz-1
        refc(i,j,k) = refc(i,j,1)
      END DO
    END DO
  END DO


  RETURN
END SUBROUTINE cal_rfc

!
!##################################################################
!##################################################################
!######                                                      ######
!######                   PROGRAM ARPSDAS                    ######
!######             ARPS Data Analysis System                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE ETS_Ref(nx,ny,nz,mosaic3d,fcst3d,           &
                   ibgn,iend,jbgn,jend,                &
                   refcomp,nreflvl,mxscorelvl,lvls,    &
                   threatscore,biasscore)
!
!  PURPOSE: calculate ETS of predicted reflectivity against observations
!
!
!  AUTHOR:
!
!  Ming Hu, CAPS, May, 2006
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

!
  INTEGER :: ibgn,iend,jbgn,jend
  INTEGER :: refcomp,nreflvl
!
! predicted and observed reflectiivty
!
  REAL :: mosaic3d(nx,ny,nz)
  REAL :: fcst3d(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
!  Temporary work arrays
!
!-----------------------------------------------------------------------
!
  REAL :: tem2d(nx,ny)
  REAL :: fcst2d(nx,ny)
!-----------------------------------------------------------------------
!
!  for ETS calculation
!-----------------------------------------------------------------------
!
!
  INTEGER :: mxscorelvl
  REAL :: hits(mxscorelvl), misses(mxscorelvl)
  REAL :: falsealarms(mxscorelvl), corneg(mxscorelvl)
  REAL :: forecastpoints(mxscorelvl), observationpoints(mxscorelvl)
  REAL :: hitsrandom
  REAL :: lvls(20)
  REAL :: threatscore(mxscorelvl),biasscore(mxscorelvl)

!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz
  INTEGER :: nxf,nyf,nzf
!
  integer :: num
  integer :: i,j,k ,ilng,nchanl
  real :: rmax,denom
!
!-----------------------------------------------------------------------
!

  if( mxscorelvl > 20 ) then
    write(*,*) 'Too many threshold for statistic'
    stop 123
  endif

  IF(refcomp == 1 ) then
   tem2d=0.0
   DO j=1,ny
   DO i=1,nx
     rmax=-9999.9
     DO k=2,nz-1
       if( mosaic3d(i,j,k) > rmax) rmax=mosaic3d(i,j,k)
     ENDDO
     tem2d(i,j) = rmax
   ENDDO
   ENDDO

   fcst2d=0.0
   DO j=1,ny
   DO i=1,nx
     rmax=-9999.9
     DO k=2,nz-1
       if( fcst3d(i,j,k) > rmax) rmax=fcst3d(i,j,k)
     ENDDO
     fcst2d(i,j) = rmax
   ENDDO
   ENDDO

  ELSEIF( refcomp == 0 ) then
   DO j=1,ny
   DO i=1,nx
     tem2d(i,j) = mosaic3d(i,j,nreflvl)
     fcst2d(i,j) = fcst3d(i,j,nreflvl)
   ENDDO
   ENDDO
  ELSEIF ( refcomp ==2) then

     DO j=1,ny
     DO i=1,nx
     tem2d(i,j) = mosaic3d(i,j,1)
     fcst2d(i,j) = fcst3d(i,j,1)
     ENDDO
     ENDDO
  ELSE
    write(*,*) ' Unknown choice for refcomp'
    write(*,*) ' refcomp should be 0 or 1'
    stop 123
  ENDIF

!  Equitable threat score and bias score

    num=0
    hits=0
    misses=0
    falsealarms=0
    corneg=0
    forecastpoints=0
    observationpoints=0
    DO i=ibgn,iend
      DO j=jbgn,jend
        num=num+1
        DO k=1,mxscorelvl
          IF(fcst2d(i,j) >= lvls(k) .and. tem2d(i,j) >= lvls(k) ) THEN
            hits(k)=hits(k) + 1
          ELSE IF ( fcst2d(i,j) < lvls(k) .and. tem2d(i,j) >= lvls(k) ) THEN
            misses(k)=misses(k)+1
          ELSE IF ( fcst2d(i,j) >= lvls(k) .and. tem2d(i,j) < lvls(k) ) THEN
            falsealarms(k)=falsealarms(k)+1
          ELSE
            corneg(k)=corneg(k)+1
          END IF

          IF( fcst2d(i,j) >= lvls(k) )  forecastpoints(k) = forecastpoints(k) + 1
          IF( tem2d(i,j) >= lvls(k) )  observationpoints(k) = observationpoints(k) + 1
        ENDDO
      ENDDO
    ENDDO

    IF( num > 0 ) THEN
      write(*,*) ' The total count is:',num
      DO k=1,mxscorelvl
        hitsrandom= forecastpoints(k) * observationpoints(k) / float(num)
        denom= hits(k)+misses(k)+falsealarms(k)-hitsrandom
        IF(denom > 0.0) THEN
          threatscore(k)=(hits(k)-hitsrandom)/denom
        ELSE
          threatscore(k)=0.0
        END IF
        denom=hits(k)+misses(k)
        IF(denom > 0.0) THEN
          biasscore(k)=(hits(k)+falsealarms(k))/denom
        ELSE
          threatscore(k)=0.0
        END IF
      END DO
    ELSE
      write(*,*) ' Error: no data were found in statistic !'
      stop 234
    ENDIF

!  write(*,*) ' Bias Score                     Equitable threat score '
!    write(*,'(f10.2,8f8.4,8f8.4)') (biasscore(k),k=1,mxscorelvl),  &
!                      (threatscore(k), k=1,mxscorelvl)

  return

 400 write(*,*) ' Error in read in Radar file ! '
   stop 8888
 102 write(*,*) ' Error in read in radar data !'


END Subroutine ETS_Ref


!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE READADCOL                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE read_Obsrain(nx,ny,nz,rfname,gridrain)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  03/17/09 Kun Zhao
!  read in hourly precipitation currently support Kun defined format
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    reffmt    file format (1:NSSL,2:CAPS)
!    rfname    radar file name (character*80)
!    gridrain   2d hourly precipitation
!
!  OUTPUT:
!    data are written to file
  INCLUDE 'grid.inc'
!
  INTEGER :: nx,ny,nz
  INTEGER :: reffmt
  CHARACTER (LEN=256) :: rfname
  REAL :: gridrain(nx,ny,nz)
  OPEN(iunit,FILE=TRIM(rfname),STATUS='UNKNOWN',FORM='FORMATTED')
  do j=1,ny
  do i=1,nx
  read(iunit,*)gridrain(i,j,1)
  enddo
  enddo
  CLOSE(iunit)

  return

END SUBROUTINE read_Obsrain


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE READADCOL                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE readadcol_Mosaic(nx,ny,nz,reffmt,rfname,gridref,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  read in gridded radar data to a file as columns with
!  individual lat,lons.
!
!-----------------------------------------------------------------------
!
!  01/28/06 MIng HU
!  Modified to fit the need to write NSSL mosaic reflectiivty
!
!-----------------------------------------------------------------------
!
!  INPUT:
!zhaokun noted
!    reffmt    file format (1:NSSL,2:CAPS)
!    rfname    radar file name (character*80)
!    gridref   save 3d reflectivity file
!
!  OUTPUT:
!    data are written to file
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INCLUDE 'grid.inc'
!
  INTEGER :: nx,ny,nz
  INTEGER :: reffmt
!
  INTEGER :: dmpfmt
  INTEGER :: iradfmt
  INTEGER :: hdf4cmpr
  CHARACTER (LEN=256) :: rfname
  INTEGER :: iyr,imon,iday,ihr,imin,isec
  INTEGER :: isource
!
  REAL :: zs(nz)
  REAL :: zpsc(nx,ny,nz)
  REAL :: gridref(nx,ny,nz)
!
  REAL :: readk(nz)
  REAL :: readhgt(nz)
  REAL :: readref(nz)

  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
!  Radar output descriptors
!
!-----------------------------------------------------------------------
!
!  INTEGER :: mxradvr,nradvr
!  PARAMETER(mxradvr=10,nradvr=6)
!  INTEGER :: iradvr(mxradvr)
!  DATA iradvr /1,2,3,4,5,6,0,0,0,0/
!
!-----------------------------------------------------------------------
!
!  Radar output thresholds
!
!-----------------------------------------------------------------------
!
  REAL :: refmin,refmax,velmin,velmax
  PARAMETER(refmin=-5.0, refmax=100.,                                   &
            velmin=-200.,velmax=200.)
  REAL :: misval
  PARAMETER(misval=-999.0)
!
!-----------------------------------------------------------------------
!
!  Radar output variables
!
!-----------------------------------------------------------------------
!
  REAL :: grdlatc(nx,ny)
  REAL :: grdlonc(nx,ny)
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=4) :: radid
  CHARACTER (LEN=80) :: runname
  CHARACTER (LEN=40) :: varname
  CHARACTER (LEN=20) :: varunits
  CHARACTER (LEN=256) :: dirname

  INTEGER :: iunit,myr,itime
  INTEGER :: i,j,k,klev,kk,kntcol,nn
  INTEGER :: idummy
  INTEGER :: istat,sd_id
  INTEGER :: islash,idot,lrfname
  INTEGER :: ipt
  REAL :: time
  REAL :: gridlat,gridlon,elev,rdummy
  INTEGER(2), allocatable :: itmp(:,:,:) ! Temporary array
  REAL, allocatable :: hmax(:), hmin(:) ! Temporary array
  LOGICAL, PARAMETER :: back = .TRUE.
  LOGICAL, PARAMETER :: noback = .FALSE.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  print *, ' inside readadcol_mosaic'
  print *, ' reffmt =',reffmt
  idummy=-999
  rdummy=-999.
  dmpfmt=1
  iradfmt=1
  hdf4cmpr=1
  time=0.0
  istatus=0
!
  IF (dmpfmt > 1 .AND. hdf4cmpr > 3) THEN
   write(*,*) 'HDF is not available!'
   stop 123
  END IF

  IF (reffmt==3) then
! reffmt is CAPS wrtvar2 format, output of mosaic2arps, composite reflectivity
! for now this is hard-coded for netcdf format, should be a variable
! fill gridref with misval and then composite saved to level 2
     lrfname=LEN_TRIM(rfname)
     islash=INDEX(TRIM(rfname),'/',back)
     dirname=rfname(1:(islash-1))
     idot=INDEX(rfname((islash+1):lrfname),'.')
     runname=rfname((islash+1):(islash+idot-1))
     gridref=misval
     CALL readvar2(nx,ny,1, gridref(1,1,2),'compst',varname,    &
            varunits, time, TRIM(runname), TRIM(dirname), 7, 0, istatus)

  ELSE IF (reffmt==2) then
! reffmt is CAPS format(88d2arps)

  OPEN(iunit,FILE=TRIM(rfname),STATUS='UNKNOWN',FORM='UNFORMATTED')

!kzhao temporatilay add to read arps generated obsref
  read(iunit)nx,ny,nz
  read(iunit)gridref
  CLOSE(iunit)
!    CALL readvar2(nx,ny,nz, gridref,'refmos','ref',    &
!                      'm/s', time, 'refmos', TRIM(rfname),            &
!                      1, 0, 0)
!NSSL format
  else if(reffmt == 1) then

  IF(dmpfmt == 1)THEN
  CALL getunit(iunit)
!
!-----------------------------------------------------------------------
!
!  Open file for output
!
!-----------------------------------------------------------------------
!
  write(*,*) rfname
  OPEN(iunit,FILE=TRIM(rfname),STATUS='UNKNOWN',FORM='UNFORMATTED')
!
!-----------------------------------------------------------------------
!
!  Write radar description variables
!
!-----------------------------------------------------------------------
!
  read(iunit) radid
  read(iunit) idummy,itime,idummy,isource,idummy,                       &
               idummy,idummy,idummy,idummy,idummy
!
  CALL abss2ctim(itime,myr,imon,iday,ihr,imin,isec)
!  write(*,*) 'The time of the data is: ',myr,imon,iday,ihr,imin,isec
!-----------------------------------------------------------------------
!
!  Write grid description variables
!  This should provide enough info to verify that the
!  proper grid has been chosen.  To recreate the grid,
!  icluding elevation information,
!  the reading program should get a grid-base-file
!  named runname.grdbasfil
!
!-----------------------------------------------------------------------
!
  idummy=0
  rdummy=0.
  read(iunit) runname
  read(iunit) iradfmt,strhopt,mapproj,idummy,idummy,                    &
               idummy,idummy,idummy,idummy,idummy
  read(iunit) dx,dy,dz,dzmin,ctrlat,                                    &
               ctrlon,trulat1,trulat2,trulon,sclfct,                    &
               rdummy,rdummy,rdummy,rdummy,rdummy
  read(iunit) idummy,idummy
  ELSE  !HDF4 format
!
   write(*,*) 'HDF is not availble now!'
   stop 123
  ENDIF
!
!-----------------------------------------------------------------------
!
!  For each horizontal grid point form a column of remapped
!  data containing the non-missing grid points
!
!-----------------------------------------------------------------------
!
  IF(dmpfmt==1)THEN

     DO ipt=1,(nx*ny)

       read(iunit,END=51) i,j,rdummy,rdummy,                    &
                   gridlat,gridlon,elev,klev
       read(iunit,END=51) (readk(k),k=1,klev)
       read(iunit,END=51) (readhgt(k),k=1,klev)
       read(iunit,END=51) (readref(k),k=1,klev)

       IF(i <= nx.AND.i >= 1 .AND. j <= ny.AND.j >= 1) THEN
          DO kk=1,klev
             k=nint(readk(kk))
             IF(k <= nz.AND.k >= 1) THEN
                gridref(i,j,k)=readref(kk)
             END IF  ! 1 < k < nz
          END DO  ! kk = 1, klev
       END IF  ! 1 < i < nx  & 1 < j < ny

     END DO  ! ipt = 1, nx*ny
     51  CONTINUE
     ipt=ipt-1
     WRITE(6,'(a,i6,a)') ' End of file reached after reading',          &
                       ipt,' columns'

  ELSE    !HDF4 format

   write(*,*) 'HDF is not availble now!'
   stop 123
  ENDIF


  ENDIF

  RETURN
END SUBROUTINE readadcol_Mosaic
