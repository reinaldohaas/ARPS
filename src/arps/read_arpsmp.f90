!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE READARPSMP                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
SUBROUTINE readarpsmp(ncompressx,ncompressy,nproc_node,ioffset,joffset, &
             hinfmt,grdbasin,lengbf,hisfile,lenhist,nx,ny,nz,nzsoil,nstyps,&
             time,x,y,z,zp,zpsoil,uprt,vprt,wprt,ptprt,pprt,qvprt,      &
             qscalar,tke,kmh,kmv,                                       &
             ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,                    &
             soiltyp,stypfrct,vegtyp,lai,roufns,veg,                    &
             tsoil,qsoil,wetcanp,snowdpth,                              &
             raing,rainc,prcrate, radfrc,radsw,rnflx,radswnet, radlwin, &
             usflx,vsflx,ptsflx,qvsflx,istatus,tem1,tem2,tem3)

!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!    9/18/2008.
!
!  MODIFICATION HISTORY:
!
!  05/08/2012 (Y. Wang)
!  Moved it from arpspltderive.f90 to libarps.a and unified all calls in
!  ARPSPLT, ARPS4WRF, ARPS2COAMPS, ARPSPOST and ARPSENSIC.
!
!-----------------------------------------------------------------------
!
!  DATA ARRAYS READ IN:
!
!    x        x-coordinate of grid points in physical/comp. space (m)
!    y        y-coordinate of grid points in physical/comp. space (m)
!    z        z-coordinate of grid points in computational space (km)
!    zp       z-coordinate of grid points in computational space (m)
!
!    uprt     x-component of perturbation velocity (m/s)
!    vprt     y-component of perturbation velocity (m/s)
!    wprt     vertical component of perturbation velocity in Cartesian
!             coordinates (m/s).
!
!    ptprt    perturbation potential temperature (K)
!    pprt     perturbation pressure (Pascal)
!
!    qvprt    perturbation water vapor mixing ratio (kg/kg)
!    qscalar  Hydrometeor scalars
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!
!    ubar     Base state x-velocity component (m/s)
!    vbar     Base state y-velocity component (m/s)
!    wbar     Base state z-velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhobar   Base state air density (kg/m**3)
!    qvbar    Base state water vapor mixing ratio (kg/kg)
!
!    soiltyp  Soil type
!    stypfrct Soil type fraction
!    vegtyp   Vegetation type
!    lai      Leaf Area Index
!    roufns   Surface roughness
!    veg      Vegetation fraction
!
!    tsoil    soil temperature (K)
!    qsoil    soil moisture
!    wetcanp  Canopy water amount
!    raing    Grid supersaturation rain (mm)
!    rainc    Cumulus convective rain(mm)
!    raint    Total rain (rainc+raing)(mm)
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
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
!    tem4     Temporary work array.
!    tem5     Temporary work array.
!    tem6     Temporary work array.
!    tem7     Temporary work array.
!    tem8     Temporary work array.
!    tem9     Temporary work array.
!
!  (These arrays are defined and used locally (i.e. inside this
!   subroutine), they may also be passed into routines called by
!   this one. Exiting the call to this subroutine, these temporary
!   work arrays may be used for other purposes, and therefore their
!   contents may be overwritten. Please examine the usage of work
!   arrays before you make any change to the code.)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!

  INTEGER, INTENT(IN) :: nx,ny,nz          ! Grid dimensions.

  INTEGER, INTENT(IN) :: nzsoil            ! levels of soil model
  INTEGER, INTENT(IN) :: nstyps            ! Maximum number of soil types.

  INTEGER, INTENT(IN) :: hinfmt
  INTEGER, INTENT(IN) :: lenhist, lengbf

  CHARACTER(LEN=256), INTENT(IN) :: grdbasin
  CHARACTER(LEN=256), INTENT(IN) :: hisfile

  INTEGER, INTENT(OUT) :: istatus
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'indtflg.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'phycst.inc'

  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Arrays to be read in:
!
!-----------------------------------------------------------------------
!
  REAL, INTENT(OUT) :: time

  REAL, DIMENSION(nx), TARGET :: x   ! The x-coord. of the physical and
                                     ! computational grid. Defined at u-point.
  REAL, DIMENSION(ny), TARGET :: y   ! The y-coord. of the physical and
                                     ! computational grid. Defined at v-point.
  REAL, DIMENSION(nz), TARGET :: z   ! The z-coord. of the computational grid.
                                     ! Defined at w-point on the staggered grid.
  REAL, DIMENSION(nx,ny,nz),     TARGET :: zp     ! The physical height coordinate defined at
                                                  ! w-point of the staggered grid.
  REAL, DIMENSION(nx,ny,nzsoil), TARGET :: zpsoil ! The physical height coordinate defined at
                                                  ! w-point of the staggered grid for soil model.

  REAL, DIMENSION(nx,ny,nz), TARGET :: uprt   ! Perturbation u-velocity (m/s)
  REAL, DIMENSION(nx,ny,nz), TARGET :: vprt   ! Perturbation v-velocity (m/s)
  REAL, DIMENSION(nx,ny,nz), TARGET :: wprt   ! Perturbation w-velocity (m/s)
  REAL, DIMENSION(nx,ny,nz), TARGET :: ptprt  ! Perturbation potential temperature
                                            ! from that of base state atmosphere (K)
  REAL, DIMENSION(nx,ny,nz), TARGET :: pprt   ! Perturbation pressure from that
                                            ! of base state atmosphere (Pascal)
  REAL, DIMENSION(nx,ny,nz), TARGET :: qvprt
  REAL, DIMENSION(nx,ny,nz), TARGET :: qv     ! Water vapor specific humidity (kg/kg)
  REAL, DIMENSION(nx,ny,nz,nscalar), TARGET :: qscalar
  REAL, DIMENSION(nx,ny,nz), TARGET :: tke    ! Turbulent Kinetic Energy ((m/s)**2)
  REAL, DIMENSION(nx,ny,nz), TARGET :: kmh    ! Horizontal turb. mixing coef. for
                                            ! momentum. ( m**2/s )
  REAL, DIMENSION(nx,ny,nz), TARGET :: kmv    ! Vertical turb. mixing coef. for
                                            ! momentum. ( m**2/s )
  REAL, DIMENSION(nx,ny,nz), TARGET :: ubar   ! Base state u-velocity (m/s)
  REAL, DIMENSION(nx,ny,nz), TARGET :: vbar   ! Base state v-velocity (m/s)
  REAL, DIMENSION(nx,ny,nz), TARGET :: wbar   ! Base state w-velocity (m/s)
  REAL, DIMENSION(nx,ny,nz), TARGET :: ptbar  ! Base state potential temperature (K)
  REAL, DIMENSION(nx,ny,nz), TARGET :: pbar   ! Base state pressure (Pascal)
  REAL, DIMENSION(nx,ny,nz), TARGET :: rhobar ! Base state density rhobar
  REAL, DIMENSION(nx,ny,nz), TARGET :: qvbar  ! Base state water vapor specific humidity
                                            ! (kg/kg)

  INTEGER, DIMENSION(nx,ny,nstyps), TARGET :: soiltyp ! Soil type
  INTEGER, DIMENSION(nx,ny),        TARGET :: vegtyp  ! Vegetation type
  REAL, DIMENSION(nx,ny,nstyps),    TARGET :: stypfrct  ! Soil type fraction
  REAL, DIMENSION(nx,ny),      TARGET :: lai       ! Leaf Area Index
  REAL, DIMENSION(nx,ny),      TARGET :: roufns   ! Surface roughness
  REAL, DIMENSION(nx,ny),      TARGET :: veg      ! Vegetation fraction

  REAL, DIMENSION(nx,ny,nzsoil,0:nstyps), TARGET :: tsoil   ! soil temperature (K)
  REAL, DIMENSION(nx,ny,nzsoil,0:nstyps), TARGET :: qsoil   ! soil moisture
  REAL, DIMENSION(nx,ny,0:nstyps),        TARGET :: wetcanp   ! Canopy water amount
  REAL, DIMENSION(nx,ny),     TARGET :: snowdpth  ! Snow depth (m)
  REAL, DIMENSION(nx,ny),     TARGET :: raing     ! Grid supersaturation rain
  REAL, DIMENSION(nx,ny),     TARGET :: rainc     ! Cumulus convective rain
  REAL, DIMENSION(nx,ny),     TARGET :: raint     ! Total rain (rainc+raing)
  REAL, DIMENSION(nx,ny,4),  TARGET :: prcrate ! precipitation rate (kg/(m**2*s))
                                             ! prcrate(1,1,1) = total precip. rate
                                             ! prcrate(1,1,2) = grid scale precip. rate
                                             ! prcrate(1,1,3) = cumulus precip. rate
                                             ! prcrate(1,1,4) = microphysics precip. rate

  REAL, DIMENSION(nx,ny,nz), TARGET :: radfrc  ! Radiation forcing (K/s)
  REAL, DIMENSION(nx,ny), TARGET :: radsw     ! Solar radiation reaching the surface
  REAL, DIMENSION(nx,ny), TARGET :: rnflx     ! Net radiation flux absorbed by surface
  REAL, DIMENSION(nx,ny), TARGET :: radswnet  ! Net shortwave radiation
  REAL, DIMENSION(nx,ny), TARGET :: radlwin   ! Incoming longwave radiation


  REAL, DIMENSION(nx,ny), TARGET :: usflx     ! Surface flux of u-momentum (kg/(m*s**2))
  REAL, DIMENSION(nx,ny), TARGET :: vsflx     ! Surface flux of v-momentum (kg/(m*s**2))
  REAL, DIMENSION(nx,ny), TARGET :: ptsflx    ! Surface heat flux (K*kg/(m*s**2))
  REAL, DIMENSION(nx,ny), TARGET :: qvsflx    ! Surface moisture flux (kg/(m**2*s))
!
!-----------------------------------------------------------------------
!
!  Temporary work arrays for general use
!
!-----------------------------------------------------------------------
!
  REAL, INTENT(INOUT) :: tem1(nx,ny,nz)
  REAL, INTENT(INOUT) :: tem2(nx,ny,nz)
  REAL, INTENT(INOUT) :: tem3(nx,ny,nz)
!
!-----------------------------------------------------------------------
!
! Variables for mpi jobs
!
!-----------------------------------------------------------------------
  INTEGER, INTENT(IN) :: ncompressx, ncompressy ! compression in x and y direction:
                                    ! ncompressx=nprocx_in/nproc_x
                                    ! ncompressy=nprocy_in/nproc_y
  INTEGER, INTENT(IN) :: nproc_node
  INTEGER, INTENT(IN) :: ioffset, joffset

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: nxsm, nysm             ! smaller domain

  INTEGER :: n, ii,jj,ia,ja
  INTEGER :: i, j

  REAL, DIMENSION(:),     POINTER :: xsm,ysm
  REAL, DIMENSION(:,:,:), POINTER :: zpsm,zpsoilsm
  REAL, DIMENSION(:,:,:), POINTER :: uprtsm, vprtsm, wprtsm,            &
                                     ptprtsm, pprtsm, qvprtsm
  REAL, DIMENSION(:,:,:,:), POINTER :: qscalarsm
  REAL, DIMENSION(:,:,:), POINTER :: tkesm, kmhsm, kmvsm
  REAL, DIMENSION(:,:,:), POINTER :: ubarsm, vbarsm, wbarsm,            &
                                     ptbarsm,pbarsm, qvbarsm,rhobarsm
  REAL, DIMENSION(:,:,:), POINTER :: prcratesm

  INTEGER, DIMENSION(:,:,:),  POINTER :: soiltypsm
  INTEGER, DIMENSION(:,:),    POINTER :: vegtypsm
  REAL,    DIMENSION(:,:,:),  POINTER :: stypfrctsm
  REAL,    DIMENSION(:,:,:,:),POINTER :: tsoilsm, qsoilsm
  REAL,    DIMENSION(:,:,:),  POINTER :: wetcanpsm
  REAL,    DIMENSION(:,:),    POINTER :: laism, roufnssm, vegsm,        &
                                         snowdpthsm, raingsm, raincsm
  REAL, DIMENSION(:,:,:),  POINTER :: radfrcsm(:,:,:)
  REAL, DIMENSION(:,:),    POINTER :: radswsm, rnflxsm, radswnetsm, radlwinsm
  REAL, DIMENSION(:,:),    POINTER :: usflxsm, vsflxsm, ptsflxsm, qvsflxsm

  CHARACTER(LEN=256) :: grdbasfn, filename
  INTEGER            :: lenbas, lenfil

  INTEGER :: nchin

  INTEGER :: ireturn

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  Set some parameters
!
!-----------------------------------------------------------------------

  nxsm = (nx-3)/ncompressx + 3
  nysm = (ny-3)/ncompressy + 3

  nstyp = nstyps

!-----------------------------------------------------------------------
!
!  Allocate arrays for smaller subdomain if needed
!
!-----------------------------------------------------------------------

  IF(ncompressx > 1 .OR. ncompressy > 1) THEN    ! allocate arrays
                                                 ! otherwise, just link pointers
    ALLOCATE(xsm     (nxsm),            STAT=istatus)
    ALLOCATE(ysm     (nysm),            STAT=istatus)
    ALLOCATE(zpsm    (nxsm,nysm,nz),    STAT=istatus)
    ALLOCATE(zpsoilsm(nxsm,nysm,nzsoil),STAT=istatus)

    ALLOCATE(uprtsm(nxsm,nysm,nz), STAT=istatus)
    ALLOCATE(vprtsm(nxsm,nysm,nz), STAT=istatus)
    ALLOCATE(wprtsm(nxsm,nysm,nz), STAT=istatus)
    ALLOCATE(ptprtsm(nxsm,nysm,nz),STAT=istatus)
    ALLOCATE(pprtsm(nxsm,nysm,nz), STAT=istatus)
    ALLOCATE(qvprtsm(nxsm,nysm,nz),STAT=istatus)

    ALLOCATE(qscalarsm  (nxsm,nysm,nz,nscalar),STAT=istatus)

    ALLOCATE(tkesm (nxsm,nysm,nz),STAT=istatus)
    ALLOCATE(kmhsm (nxsm,nysm,nz),STAT=istatus)
    ALLOCATE(kmvsm (nxsm,nysm,nz),STAT=istatus)

    ALLOCATE(ubarsm  (nxsm,nysm,nz),STAT=istatus)
    ALLOCATE(vbarsm  (nxsm,nysm,nz),STAT=istatus)
    ALLOCATE(wbarsm  (nxsm,nysm,nz),STAT=istatus)
    ALLOCATE(ptbarsm (nxsm,nysm,nz),STAT=istatus)
    ALLOCATE(pbarsm  (nxsm,nysm,nz),STAT=istatus)
    ALLOCATE(rhobarsm(nxsm,nysm,nz),STAT=istatus)
    ALLOCATE(qvbarsm (nxsm,nysm,nz),STAT=istatus)

    ALLOCATE(soiltypsm (nxsm,nysm,nstyps),STAT=istatus)
    ALLOCATE(stypfrctsm(nxsm,nysm,nstyps),STAT=istatus)
    ALLOCATE(vegtypsm  (nxsm,nysm),STAT=istatus)
    ALLOCATE(laism     (nxsm,nysm),STAT=istatus)
    ALLOCATE(roufnssm  (nxsm,nysm),STAT=istatus)
    ALLOCATE(vegsm     (nxsm,nysm),STAT=istatus)

    ALLOCATE(tsoilsm   (nxsm,nysm,nzsoil,0:nstyps),STAT=istatus)
    ALLOCATE(qsoilsm   (nxsm,nysm,nzsoil,0:nstyps),STAT=istatus)
    ALLOCATE(wetcanpsm (nxsm,nysm,0:nstyps),STAT=istatus)
    ALLOCATE(snowdpthsm(nxsm,nysm),STAT=istatus)

    ALLOCATE(raingsm   (nxsm,nysm),STAT=istatus)
    ALLOCATE(raincsm   (nxsm,nysm),STAT=istatus)
    ALLOCATE(prcratesm (nxsm,nysm,4),STAT=istatus)

    ALLOCATE(radfrcsm(nxsm,nysm,nz),STAT=istatus)

    ALLOCATE(radswsm (nxsm,nysm),STAT=istatus)
    ALLOCATE(rnflxsm (nxsm,nysm),STAT=istatus)
    ALLOCATE(radswnetsm (nxsm,nysm),STAT=istatus)
    ALLOCATE(radlwinsm (nxsm,nysm),STAT=istatus)
    ALLOCATE(usflxsm (nxsm,nysm),STAT=istatus)
    ALLOCATE(vsflxsm (nxsm,nysm),STAT=istatus)
    ALLOCATE(ptsflxsm(nxsm,nysm),STAT=istatus)
    ALLOCATE(qvsflxsm(nxsm,nysm),STAT=istatus)
    CALL check_alloc_status(istatus, "arpsplt:qvsflxsm")

    xsm     = 0.0
    ysm     = 0.0
    zpsm    = 0.0
    zpsoilsm= 0.0

    uprtsm  = 0.0
    vprtsm  = 0.0
    wprtsm  = 0.0
    ptprtsm = 0.0
    pprtsm  = 0.0
    qvprtsm = 0.0
    qscalarsm    = 0.0
    tkesm   = 0.0
    kmhsm   = 0.0
    kmvsm   = 0.0

    ubarsm  = 0.0
    vbarsm  = 0.0
    wbarsm  = 0.0
    ptbarsm = 0.0
    pbarsm  = 0.0
    rhobarsm= 0.0
    qvbarsm = 0.0

    soiltypsm = 0.0
    stypfrctsm= 0.0
    vegtypsm  = 0.0
    laism     = 0.0
    roufnssm  = 0.0
    vegsm     = 0.0

    tsoilsm   = 0.0
    qsoilsm   = 0.0
    wetcanpsm = 0.0
    snowdpthsm= 0.0

    raingsm   = 0.0
    raincsm   = 0.0
    prcratesm = 0.0

    radfrcsm  = 0.0

    radswsm    = 0.0
    rnflxsm    = 0.0
    radswnetsm = 0.0
    radlwinsm  = 0.0
    usflxsm    = 0.0
    vsflxsm    = 0.0
    ptsflxsm   = 0.0
    qvsflxsm   = 0.0

  ELSE
    xsm      => x
    ysm      => y
    zpsm     => zp
    zpsoilsm => zpsoil
    uprtsm   => uprt
    vprtsm   => vprt
    wprtsm   => wprt
    ptprtsm  => ptprt
    pprtsm   => pprt
    qvprtsm  => qvprt
    qscalarsm     => qscalar
    tkesm    => tke
    kmhsm    => kmh
    kmvsm    => kmv
    ubarsm   => ubar
    vbarsm   => vbar
    wbarsm   => wbar
    ptbarsm  => ptbar
    pbarsm   => pbar
    rhobarsm => rhobar
    qvbarsm  => qvbar

    soiltypsm  => soiltyp
    stypfrctsm => stypfrct
    vegtypsm   => vegtyp
    laism      => lai
    roufnssm   => roufns
    vegsm      => veg

    tsoilsm    => tsoil
    qsoilsm    => qsoil
    wetcanpsm  => wetcanp
    snowdpthsm => snowdpth

    raingsm    => raing
    raincsm    => rainc
    prcratesm  => prcrate

    radfrcsm   => radfrc
    radswsm    => radsw
    rnflxsm    => rnflx
    radswnetsm => radswnet
    radlwinsm  => radlwin
    usflxsm    => usflx
    vsflxsm    => vsflx
    ptsflxsm   => ptsflx
    qvsflxsm   => qvsflx

  END IF

!-----------------------------------------------------------------------
!
! Reading files
!
!-----------------------------------------------------------------------

  DO jj = 1, ncompressy
    DO ii = 1, ncompressx
      IF ( ncompressx > 1 .OR. ncompressy > 1 .OR.    &
           (mp_opt > 0 .AND. readsplit(FINDX_H) <= 0) ) THEN
        CALL gtsplitfn(grdbasin,ncompressx,ncompressy,loc_x,loc_y,ii,jj, &
                       ioffset-1,joffset-1,1,lvldbg,grdbasfn,istatus)
        lenbas = LEN_TRIM(grdbasfn)

        CALL gtsplitfn(hisfile,ncompressx,ncompressy,loc_x,loc_y,ii,jj,  &
                       ioffset-1,joffset-1,1,lvldbg,filename,istatus)
        lenfil = LEN_TRIM(filename)
      ELSE
        WRITE(grdbasfn,'(a)') grdbasin(1:lengbf)
        WRITE(filename,'(a)') hisfile(1:lenhist)
        lenbas = lengbf
        lenfil = lenhist
      END IF

!        15  CONTINUE         ! also continue to read another time recode
!                             ! from GrADS file
!
!-----------------------------------------------------------------------
!
!  Read all input data arrays
!
!-----------------------------------------------------------------------
!
      CALL setgbrd(0)            ! read grid/base state file

      IF (nproc_node <= 1) THEN  ! the first readstride processes read then
                                 ! the next readstride processes and so on
        DO n = 0, nprocs-1, readstride

          IF(myproc >= n .AND. myproc <= n+readstride-1) THEN

            !IF (myproc == 0 .OR. readsplit(FINDX_H) < 1) WRITE(6,'(1x,a,I4,a,a/)')  &
            IF (myproc == 0) WRITE(6,'(1x,a,I5,a,a/)')                  &
               'process ',myproc,' reading file: ', filename(1:lenfil)

            CALL dtaread(nxsm,nysm,nz,nzsoil, nstyps,                   &
                  hinfmt, nchin,grdbasfn(1:lenbas),lenbas,              &
                  filename(1:lenfil),lenfil,time,                       &
                  xsm,ysm,z,zpsm,zpsoilsm,                              &
                  uprtsm,vprtsm,wprtsm,ptprtsm, pprtsm, qvprtsm,        &
                  qscalarsm,tkesm,kmhsm,kmvsm,                          &
                  ubarsm,vbarsm,wbarsm,ptbarsm,pbarsm,rhobarsm,qvbarsm, &
                  soiltypsm,stypfrctsm,vegtypsm,laism,roufnssm,vegsm,   &
                  tsoilsm,qsoilsm,wetcanpsm,snowdpthsm,                 &
                  raingsm,raincsm,prcratesm,                            &
                  radfrcsm,radswsm,rnflxsm,radswnetsm,radlwinsm,        &
                  usflxsm,vsflxsm,ptsflxsm,qvsflxsm,                    &
                  ireturn, tem1,tem2, tem3)

          END IF
          IF (mp_opt > 0) CALL mpbarrier
        END DO

      ELSE      ! Only one process read at one node

        DO n = 0, nproc_node-1

          IF (MOD(myproc,nproc_node) == n) THEN

            !IF (myproc == 0 .OR. readsplit(FINDX_H) < 1) WRITE(6,'(1x,a,I5,a,a/)')  &
            IF (myproc == 0) WRITE(6,'(1x,a,I5,a,a/)')                  &
              'process ',myproc,' reading file: ', filename(1:lenfil)

            CALL dtaread(nxsm,nysm,nz,nzsoil, nstyps,                   &
                  hinfmt, nchin,grdbasfn(1:lenbas),lenbas,              &
                  filename(1:lenfil),lenfil,time,                       &
                  xsm,ysm,z,zpsm,zpsoilsm,                              &
                  uprtsm,vprtsm,wprtsm,ptprtsm, pprtsm, qvprtsm,        &
                  qscalarsm,tkesm,kmhsm,kmvsm,                          &
                  ubarsm,vbarsm,wbarsm,ptbarsm,pbarsm,rhobarsm,qvbarsm, &
                  soiltypsm,stypfrctsm,vegtypsm,laism,roufnssm,vegsm,   &
                  tsoilsm,qsoilsm,wetcanpsm,snowdpthsm,                 &
                  raingsm,raincsm,prcratesm,                            &
                  radfrcsm,radswsm,rnflxsm,radswnetsm,radlwinsm,        &
                  usflxsm,vsflxsm,ptsflxsm,qvsflxsm,                    &
                  ireturn, tem1,tem2, tem3)

          END IF
          IF (mp_opt > 0) CALL mpbarrier
        END DO

      END IF

      CALL mpsumi(ireturn,1)
!
!-----------------------------------------------------------------------
!
!  ireturn = 0 for a successful read
!  For hinfmt=9, i.e. the GraDs format data, ireturn is used as a
!  flag indicating if there is any data at more time level to be read.
!
!-----------------------------------------------------------------------
!
      IF( ireturn /= 0 .AND. hinfmt /= 9 ) THEN
        WRITE(6,'(1x,a,i3)')                                            &
                'Bad return status from data file reading, ireturn = ', &
                ireturn
        CALL arpsstop('Error inside readarpsmp.',1)
      END IF

      IF(ncompressx > 1 .OR. ncompressy > 1) THEN    ! need join

        DO j = 1, nysm
          ja = (jj-1)*(nysm-3)+j
          DO i = 1, nxsm
            ia = (ii-1)*(nxsm-3)+i

            x(ia) = xsm(i)

            vegtyp(ia,ja) = vegtypsm(i,j)
            lai(ia,ja)    = laism(i,j)
            roufns(ia,ja) = roufnssm(i,j)
            veg(ia,ja)    = vegsm(i,j)
            snowdpth(ia,ja) = snowdpthsm(i,j)
            raing(ia,ja) = raingsm(i,j)
            rainc(ia,ja) = raincsm(i,j)
            radsw(ia,ja)    = radswsm(i,j)
            radswnet(ia,ja) = radswnetsm(i,j)
            radlwin(ia,ja) = radlwinsm(i,j)
            rnflx(ia,ja)    = rnflxsm(i,j)
            usflx(ia,ja)  = usflxsm(i,j)
            vsflx(ia,ja)  = vsflxsm(i,j)
            ptsflx(ia,ja) = ptsflxsm(i,j)
            qvsflx(ia,ja) = qvsflxsm(i,j)

            prcrate(ia,ja,:) = prcratesm(i,j,:)

            zp(ia,ja,:)      = zpsm(i,j,:)
            zpsoil(ia,ja,:)  = zpsoilsm(i,j,:)
            uprt(ia,ja,:)  = uprtsm(i,j,:)
            vprt(ia,ja,:)  = vprtsm(i,j,:)
            wprt(ia,ja,:)  = wprtsm(i,j,:)
            ptprt(ia,ja,:) = ptprtsm(i,j,:)
            pprt(ia,ja,:)  = pprtsm(i,j,:)
            qvprt(ia,ja,:) = qvprtsm(i,j,:)
            qscalar(ia,ja,:,:)  = qscalarsm(i,j,:,:)
            tke(ia,ja,:) = tkesm(i,j,:)
            kmh(ia,ja,:) = kmhsm(i,j,:)
            kmv(ia,ja,:) = kmvsm(i,j,:)
            ubar(ia,ja,:)  = ubarsm(i,j,:)
            vbar(ia,ja,:)  = vbarsm(i,j,:)
            wbar(ia,ja,:)  = wbarsm(i,j,:)
            ptbar(ia,ja,:) = ptbarsm(i,j,:)
            pbar(ia,ja,:)  = pbarsm(i,j,:)
            rhobar(ia,ja,:) = rhobarsm(i,j,:)
            qvbar(ia,ja,:) = qvbarsm(i,j,:)
            radfrc(ia,ja,:) = radfrcsm(i,j,:)

            soiltyp(ia,ja,:)  = soiltypsm(i,j,:)
            stypfrct(ia,ja,:) = stypfrctsm(i,j,:)

            tsoil(ia,ja,:,:) = tsoilsm(i,j,:,:)
            qsoil(ia,ja,:,:) = qsoilsm(i,j,:,:)
            wetcanp(ia,ja,:) = wetcanpsm(i,j,:)
          END DO    !i
          y(ja) = ysm(j)
        END DO      !j

      END IF   ! need join

    END DO   !ii
  END DO   ! jj

  IF(ncompressx > 1 .OR. ncompressy > 1) THEN    ! DEALLOCATE arrays
     DEALLOCATE(xsm, ysm, zpsm, zpsoilsm)

     DEALLOCATE( uprtsm,vprtsm,wprtsm,ptprtsm,pprtsm,qvprtsm )
     DEALLOCATE( qscalarsm )
     DEALLOCATE( tkesm, kmhsm, kmvsm )
     DEALLOCATE( ubarsm, vbarsm, wbarsm, ptbarsm, pbarsm )
     DEALLOCATE( rhobarsm, qvbarsm )

     DEALLOCATE( soiltypsm, stypfrctsm, vegtypsm )
     DEALLOCATE( laism, roufnssm, vegsm )
     DEALLOCATE( tsoilsm, qsoilsm, wetcanpsm, snowdpthsm )

     DEALLOCATE( raingsm, raincsm, prcratesm )
     DEALLOCATE( radfrcsm, radswsm, rnflxsm, radswnetsm, radlwinsm )
     DEALLOCATE( usflxsm, vsflxsm, ptsflxsm, qvsflxsm )
  END IF

!
!-----------------------------------------------------------------------
!
!  Set ice to be icein and then ice is used to control if snow and
!  graupel/hail affect reflectivity calculation in REFLEC.
!
!-----------------------------------------------------------------------
!

  ice = icein

  CALL gtlfnkey( runname, lfnkey )

  mgrid = 1    ! Grid number 1
  nestgrd = 0  ! Not nested grid data

  curtim = time

!-----------------------------------------------------------------------
!
! Return from this subroutine
!
!-----------------------------------------------------------------------

  istatus = ireturn

  RETURN
END SUBROUTINE readarpsmp

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE READARPSVAR                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
SUBROUTINE readarpsvar(ncompressx,ncompressy,nproc_node,                &
             hinfmt,hisfile,lenhist,nx,ny,nz,time,                      &
             varname,varout,tem1,dbglvl, istatus )

!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!    4/29/2010.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Given an ARPS file name "hisfile", read variable "varname" with dimensions
!  nx,ny & nz. The actual file name to be read is construced based on
!  hisfile, but replace with forecast time "time" if time >= 0.0.
!
!  This is a high level driver that calls dtareadvar for actual reading,
!  but doing either patch spliting or joining in this subroutine.
!
!  DATA ARRAYS READ IN:
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!

  INTEGER, INTENT(IN) :: nx,ny,nz          ! Grid dimensions.

  INTEGER, INTENT(IN) :: hinfmt
  INTEGER, INTENT(IN) :: lenhist

  INTEGER, INTENT(IN) :: ncompressx, ncompressy ! compression in x and y direction:
                                    ! ncompressx=nprocx_in/nproc_x
                                    ! ncompressy=nprocy_in/nproc_y
  INTEGER, INTENT(IN) :: nproc_node

  REAL, INTENT(IN) :: time

  CHARACTER(LEN=256), INTENT(IN) :: hisfile
  CHARACTER(LEN=*),   INTENT(IN) :: varname

  INTEGER, INTENT(IN)  :: dbglvl

  INTEGER, INTENT(OUT) :: istatus

!
!-----------------------------------------------------------------------
!
!  Arrays to be read in:
!
!-----------------------------------------------------------------------
!
  REAL, DIMENSION(nx,ny,nz), INTENT(OUT), TARGET :: varout
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Temporary work arrays for general use
!
!-----------------------------------------------------------------------
!
  REAL, INTENT(INOUT) :: tem1(nx,ny,nz)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: nxin, nyin             ! smaller domain

  INTEGER :: n, ii,jj,ia,ja
  INTEGER :: i, j

  REAL, DIMENSION(:,:,:),  POINTER :: varptr(:,:,:)

  CHARACTER(LEN=256) :: filebase, filename
  INTEGER            :: lenfil

  CHARACTER (LEN=80 ) :: timsnd
  INTEGER :: tmstrln, hindex, lindex

  LOGICAL :: ptr_allocated

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  Set some parameters
!
!-----------------------------------------------------------------------

  IF (mp_opt == 0 .OR. readsplit(FINDX_H) > 0) THEN
    nxin = (nx-3)*nproc_x + 3
    nyin = (ny-3)*nproc_y + 3
  ELSE
    nxin = (nx-3)/ncompressx + 3
    nyin = (ny-3)/ncompressy + 3
  END IF

!-----------------------------------------------------------------------
!
!  Allocate arrays for smaller subdomain if needed
!
!-----------------------------------------------------------------------

  ptr_allocated = .FALSE.
  IF( nxin /= nx .OR. nyin /= ny ) THEN    ! allocate arrays
    ptr_allocated = .TRUE.
    ALLOCATE(varptr(nxin,nyin,nz), STAT=istatus)
    varptr  = 0.0
  ELSE                                     ! otherwise, just link pointers
    varptr  => varout
  END IF

  IF (hinfmt == 3) THEN
    hindex = INDEX(hisfile,'hdf',.TRUE.)
    hindex = hindex + 2

    lindex = INDEX(hisfile,'_',.TRUE.)-1
    IF (lindex < hindex) lindex = hindex + 7

  ELSE
    WRITE(6,'(1x,a,I2,a)') 'ERROR: unsupported file format = ',hinfmt,' in readarpsvar.;'
    istatus = -1
    RETURN
  END IF

  IF (time >= 0.0) THEN
    CALL cvttsnd( time , timsnd, tmstrln )
  ELSE
    timsnd = hisfile(hindex+1:lindex)
    tmstrln = lindex-hindex
  END IF

  WRITE(filebase,'(2a)') hisfile(1:hindex),timsnd(1:tmstrln)

!-----------------------------------------------------------------------
!
! Reading files
!
!-----------------------------------------------------------------------

  IF (mp_opt == 0 .OR. readsplit(FINDX_H) > 0) THEN  ! no-MPI or MPI read and split

    filename = filebase
    lenfil   = LEN_TRIM(filename)
    IF (myproc == 0) THEN
      IF (dbglvl > 0) WRITE(6,'(1x,a,I4,a,a/)')                         &
          'process ',myproc,' reading file: ', filename(1:lenfil)
      CALL dtareadvar(nxin,nyin,nz, hinfmt, filename(1:lenfil),lenfil,  &
                      varname,varptr, tem1,  dbglvl, istatus )
    END IF
    CALL mpupdatei(istatus,1)
    IF( istatus /= 0 ) THEN
      WRITE(6,'(1x,a,i3)')                                              &
         'Bad return status from data file reading, ireturn = ', istatus
      RETURN
    END IF

    IF (mp_opt > 0 .AND. readsplit(FINDX_H) > 0) THEN
      CALL mpisplit3d(varptr, nx, ny, nz, varout)
    END IF

  ELSE  ! MPI every PE read one or more files

    DO jj = 1, ncompressy
      DO ii = 1, ncompressx

        CALL gtsplitfn(filebase,ncompressx,ncompressy,loc_x,loc_y,ii,jj,&
                       0,0,1,dbglvl,filename,istatus)
        IF (istatus /= 0) RETURN
        lenfil = LEN_TRIM(filename)

        IF (nproc_node <= 1) THEN  ! the first readstride processes read then
                                   ! the next readstride processes and so on
          DO n = 0, nprocs-1, readstride

            IF(myproc >= n .AND. myproc <= n+readstride-1) THEN

              IF (dbglvl > 0) WRITE(6,'(1x,a,I4,a,a/)')                 &
                 'process ',myproc,' reading file: ', filename(1:lenfil)

              CALL dtareadvar(nxin,nyin,nz, hinfmt, filename(1:lenfil),lenfil, &
                              varname,varptr, tem1,  dbglvl, istatus)

            END IF
            IF (mp_opt > 0) CALL mpbarrier
          END DO

        ELSE      ! Only one process read at one node

          DO n = 0, nproc_node-1

            IF (MOD(myproc,nproc_node) == n) THEN

              IF (dbglvl > 0) WRITE(6,'(1x,a,I5,a,a/)')                 &
                'process ',myproc,' reading file: ', filename(1:lenfil)

              CALL dtareadvar(nxin,nyin,nz, hinfmt, filename(1:lenfil),lenfil, &
                              varname,varptr, tem1,  dbglvl, istatus)

            END IF
            IF (mp_opt > 0) CALL mpbarrier
          END DO

        END IF

        CALL mpsumi(istatus,1)
!
!-----------------------------------------------------------------------
!
!    ireturn = 0 for a successful read
!    For hinfmt=9, i.e. the GraDs format data, ireturn is used as a
!    flag indicating if there is any data at more time level to be read.
!
!-----------------------------------------------------------------------
!
        IF( istatus /= 0 ) THEN
          WRITE(6,'(1x,a,i3)')                                            &
                  'Bad return status from data file reading, ireturn = ', &
                  istatus
          RETURN
        END IF

        IF(ncompressx > 1 .OR. ncompressy > 1) THEN    ! need join

          DO j = 1, nyin
            ja = (jj-1)*(nyin-3)+j
            DO i = 1, nxin
              ia = (ii-1)*(nxin-3)+i
              varout(ia,ja,:) = varptr(i,j,:)
            END DO    !i
          END DO      !j

        END IF   ! need join

      END DO   !ii
    END DO   ! jj
  END IF

!-----------------------------------------------------------------------
!
! Return from this subroutine
!
!-----------------------------------------------------------------------
  IF( ptr_allocated )   DEALLOCATE( varptr )

  RETURN
END SUBROUTINE readarpsvar

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE DTAREADVAR                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
SUBROUTINE dtareadvar(nx,ny,nz,hinfmt, filename, lenfil,                &
             varname,varout,tem1, dbglvl, istatus )

!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!    4/29/2010.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Given an ARPS file name "filename", read variable "varname" with
!  dimensions nx,ny & nz.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!

  INTEGER, INTENT(IN) :: nx,ny,nz          ! Grid dimensions.

  INTEGER, INTENT(IN) :: hinfmt
  INTEGER, INTENT(IN) :: lenfil

  CHARACTER(LEN=*), INTENT(IN) :: filename
  CHARACTER(LEN=*), INTENT(IN) :: varname

  REAL,    INTENT(OUT) :: varout(nx,ny,nz)

  INTEGER, INTENT(IN)  :: dbglvl
  INTEGER, INTENT(OUT) :: istatus

  REAL,  INTENT(INOUT) :: tem1(nx,ny,nz)

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: nxin, nyin,nzin             ! smaller domain

  INTEGER :: i, j

  INTEGER :: sd_id

  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE :: itmp(:,:,:)
                                                ! Temporary array
  REAL, ALLOCATABLE :: hmax(:), hmin(:)         ! Temporary array

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (hinfmt == 3) THEN

    ALLOCATE (itmp(nx,ny,nz),stat=istatus)
    CALL check_alloc_status(istatus, 'DTAREADVAR:itmp')
    ALLOCATE (hmax(nz),stat=istatus)
    CALL check_alloc_status(istatus, 'DTAREADVAR:hmax')
    ALLOCATE (hmin(nz),stat=istatus)
    CALL check_alloc_status(istatus, 'DTAREADVAR:hmin')

    CALL hdfopen(TRIM(filename),1,sd_id)
    IF (sd_id < 0) THEN
      WRITE (6,'(1x,3a)') 'DTAREADVAR-ERROR: opening ',trim(filename),' failed.'
      istatus = -1
      RETURN
    END IF

    CALL hdfrdi(sd_id,'nx',nxin,istatus)
    CALL hdfrdi(sd_id,'ny',nyin,istatus)
    CALL hdfrdi(sd_id,'nz',nzin,istatus)

    IF ( nxin /= nx .OR. nyin /= ny ) THEN
      WRITE(6,'(1x,a)') 'DTAREADVAR-ERROR: Dimensions in DTAREADVAR inconsistent with data.'
      WRITE(6,'(1x,a,3I15)') ' Read were: ', nxin, nyin, nzin
      WRITE(6,'(1x,a,3I15)') ' Expected:  ', nx, ny, nz
      istatus = -2
      RETURN
    END IF

    IF (nz > 1) THEN
      CALL hdfrd3d(sd_id,TRIM(varname),nx,ny,nz,varout,istatus,itmp,hmax,hmin)
    ELSE
      CALL hdfrd2d(sd_id,TRIM(varname),nx,ny,varout,istatus,itmp)
    END IF
    IF (istatus /= 0) THEN
      WRITE(6,'(1x,3a,I3,a)') 'DTAREADVAR-ERROR: Reading ',             &
                   TRIM(varname),' failed with status = ', istatus, '.'
    END IF

    CALL hdfclose(sd_id,istatus)
    IF (istatus == 0) THEN
      IF(dbglvl > 0)  WRITE(6,'(1x,4a)')                                &
         'DTAREADVAR: Successfully read ', trim(varname),' from ', trim(filename)
    ELSE
      WRITE(6,'(1x,a,I3,2a)') 'DTAREADVAR-ERROR:(status = ', istatus,   &
                           ') closing ', trim(filename)
    END IF

    DEALLOCATE( itmp, hmax, hmin )
  ELSE
    WRITE(6,'(1x,a,I2,a)') 'ERROR: unsupported file format = ',hinfmt,' in dtareadvar.;'
    istatus = -1
    RETURN
  END IF

  RETURN
END SUBROUTINE dtareadvar

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GET_EXT_SOIL               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
SUBROUTINE get_ext_soil(ncompressx,ncompressy,nproc_node,               &
             hinfmt,hisfile,lenhist,nx,ny,nzsoil,nstyps,                &
             zpsoil,tsoil,qsoil,istatus)

!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!    5/18/2011.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  DATA ARRAYS READ IN:
!
!    tsoil    soil temperature (K)
!    qsoil    soil moisture
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!

  INTEGER, INTENT(IN) :: nx,ny,nzsoil          ! Grid dimensions.
  INTEGER, INTENT(IN) :: nstyps            ! Maximum number of soil types.

  INTEGER, INTENT(IN) :: hinfmt
  INTEGER, INTENT(IN) :: lenhist

  CHARACTER(LEN=256), INTENT(IN) :: hisfile

  INTEGER, INTENT(OUT) :: istatus
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'

  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Arrays to be read in:
!
!-----------------------------------------------------------------------
!
  REAL, DIMENSION(nx,ny,nzsoil), TARGET :: zpsoil ! The physical height coordinate defined at
                                                  ! w-point of the staggered grid for soil model.
  REAL, DIMENSION(nx,ny,nzsoil,0:nstyps), TARGET :: tsoil   ! soil temperature (K)
  REAL, DIMENSION(nx,ny,nzsoil,0:nstyps), TARGET :: qsoil   ! soil moisture
!
!-----------------------------------------------------------------------
!
!  Temporary work arrays for general use
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
! Variables for mpi jobs
!
!-----------------------------------------------------------------------
  INTEGER, INTENT(IN) :: ncompressx, ncompressy ! compression in x and y direction:
                                    ! ncompressx=nprocx_in/nproc_x
                                    ! ncompressy=nprocy_in/nproc_y
  INTEGER, INTENT(IN) :: nproc_node

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: nxsm, nysm             ! smaller domain

  INTEGER :: n, ii,jj,ia,ja
  INTEGER :: i, j, k

  REAL, DIMENSION(:,:,:),  POINTER :: zpsoilsm
  REAL, DIMENSION(:,:,:,:),POINTER :: tsoilsm, qsoilsm

  CHARACTER(LEN=256) :: filename
  INTEGER            :: lenfil

  INTEGER :: nchin

  INTEGER :: ireturn

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
!  Set some parameters
!
!-----------------------------------------------------------------------

  nxsm = (nx-3)/ncompressx + 3
  nysm = (ny-3)/ncompressy + 3

  nstyp = nstyps

!-----------------------------------------------------------------------
!
!  Allocate arrays for smaller subdomain if needed
!
!-----------------------------------------------------------------------

  IF(ncompressx > 1 .OR. ncompressy > 1) THEN    ! allocate arrays
                                                 ! otherwise, just link pointers
    ALLOCATE(zpsoilsm(nxsm,nysm,nzsoil),STAT=istatus)

    ALLOCATE(tsoilsm   (nxsm,nysm,nzsoil,0:nstyps),STAT=istatus)
    ALLOCATE(qsoilsm   (nxsm,nysm,nzsoil,0:nstyps),STAT=istatus)
    CALL check_alloc_status(istatus, "arpsplt:qvsflxsm")

    zpsoilsm= 0.0

    tsoilsm   = 0.0
    qsoilsm   = 0.0

  ELSE
    zpsoilsm => zpsoil

    tsoilsm    => tsoil
    qsoilsm    => qsoil
  END IF

!-----------------------------------------------------------------------
!
! Reading files
!
!-----------------------------------------------------------------------

  DO jj = 1, ncompressy
    DO ii = 1, ncompressx
!      IF (mp_opt > 0 .AND. readsplit <= 0) THEN
      IF ( ncompressx > 1 .OR. ncompressy > 1 .OR.           &
           (mp_opt > 0 .AND. readsplit(FINDX_H) <= 0) ) THEN

        CALL gtsplitfn(hisfile,ncompressx,ncompressy,loc_x,loc_y,ii,jj, &
                       0,0,1,lvldbg,filename,istatus)
        lenfil = LEN_TRIM(filename)
      ELSE
        WRITE(filename,'(a)') hisfile(1:lenhist)
        lenfil = lenhist
      END IF

!
!-----------------------------------------------------------------------
!
!  Read all input data arrays
!
!-----------------------------------------------------------------------
!
      IF (nproc_node <= 1) THEN  ! the first readstride processes read then
                                 ! the next readstride processes and so on
        DO n = 0, nprocs-1, readstride

          IF(myproc >= n .AND. myproc <= n+readstride-1) THEN

            IF (myproc == 0 .OR. readsplit(FINDX_H) < 1) WRITE(6,'(1x,a,I4,a,a/)')  &
               'process ',myproc,' reading file: ', filename(1:lenfil)

              CALL dtareadvarsoil(nxsm,nysm,nzsoil,nstyps,hinfmt,filename,lenfil, &
                          zpsoilsm,tsoilsm,qsoilsm, 0, istatus )

          END IF
          IF (mp_opt > 0) CALL mpbarrier
        END DO

      ELSE      ! Only one process read at one node

        DO n = 0, nproc_node-1

          IF (MOD(myproc,nproc_node) == n) THEN

            IF (myproc == 0 .OR. readsplit(FINDX_H) < 1) WRITE(6,'(1x,a,I5,a,a/)')  &
              'process ',myproc,' reading file: ', filename(1:lenfil)

              CALL dtareadvarsoil(nxsm,nysm,nzsoil,nstyps,hinfmt,filename,lenfil, &
                          zpsoilsm,tsoilsm,qsoilsm, 0, istatus )

          END IF
          IF (mp_opt > 0) CALL mpbarrier
        END DO

      END IF

      CALL mpsumi(istatus,1)
!
!-----------------------------------------------------------------------
!
!  ireturn = 0 for a successful read
!  For hinfmt=9, i.e. the GraDs format data, ireturn is used as a
!  flag indicating if there is any data at more time level to be read.
!
!-----------------------------------------------------------------------
!
      IF( istatus /= 0 ) THEN
        WRITE(6,'(1x,a,i3)')                                            &
                'Bad return status from data file reading, ireturn = ', &
                ireturn
        CALL arpsstop('Error inside readarpsmp.',1)
      END IF

      IF(ncompressx > 1 .OR. ncompressy > 1) THEN    ! need join

        DO n = 0,nstyps
          DO k = 1,nzsoil
            DO j = 1, nysm
              ja = (jj-1)*(nysm-3)+j
              DO i = 1, nxsm
                ia = (ii-1)*(nxsm-3)+i

                zpsoil(ia,ja,k)  = zpsoilsm(i,j,k)
                tsoil(ia,ja,k,n) = tsoilsm(i,j,k,n)
                qsoil(ia,ja,k,n) = qsoilsm(i,j,k,n)
              END DO    !i
            END DO      !j
          END DO
        END DO

      END IF   ! need join

    END DO   !ii
  END DO   ! jj

  IF(ncompressx > 1 .OR. ncompressy > 1) THEN    ! DEALLOCATE arrays
     DEALLOCATE(zpsoilsm)
     DEALLOCATE(tsoilsm, qsoilsm)
  END IF

  RETURN
END SUBROUTINE get_ext_soil

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE DTAREADVAR                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
SUBROUTINE dtareadvarsoil(nx,ny,nzsoil,nstyps,hinfmt, filename, lenfil, &
                          zpsoil,tsoil,qsoil, dbglvl, istatus )

!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!    5/17/2011.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Given an ARPS file name "filename", read variable "varname" with
!  dimensions nx,ny, nzsoil, & nstyps.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!

  INTEGER, INTENT(IN) :: nx,ny,nzsoil,nstyps         ! Grid dimensions.

  INTEGER, INTENT(IN) :: hinfmt
  INTEGER, INTENT(IN) :: lenfil

  CHARACTER(LEN=*), INTENT(IN) :: filename

  REAL,    INTENT(OUT) :: zpsoil(nx,ny,nzsoil)
  REAL,    INTENT(OUT) :: tsoil(nx,ny,nzsoil,0:nstyps)
  REAL,    INTENT(OUT) :: qsoil(nx,ny,nzsoil,0:nstyps)

  INTEGER, INTENT(IN)  :: dbglvl
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: nxin, nyin,nzsoilin             ! smaller domain

  INTEGER :: i, j

  INTEGER :: sd_id

  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE :: itmp(:,:,:,:)
                                                ! Temporary array
  REAL, ALLOCATABLE :: hmax(:), hmin(:)         ! Temporary array

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (hinfmt == 3) THEN

    ALLOCATE (itmp(nx,ny,nzsoil,0:nstyps),stat=istatus)
    CALL check_alloc_status(istatus, 'DTAREADVAR:itmp')
    ALLOCATE (hmax(nzsoil),stat=istatus)
    CALL check_alloc_status(istatus, 'DTAREADVAR:hmax')
    ALLOCATE (hmin(nzsoil),stat=istatus)
    CALL check_alloc_status(istatus, 'DTAREADVAR:hmin')

    CALL hdfopen(TRIM(filename),1,sd_id)
    IF (sd_id < 0) THEN
      WRITE (6,'(1x,3a)') 'DTAREADVAR-ERROR: opening ',trim(filename),' failed.'
      istatus = -1
      RETURN
    END IF

    CALL hdfrdi(sd_id,'nx',nxin,istatus)
    CALL hdfrdi(sd_id,'ny',nyin,istatus)
    CALL hdfrdi(sd_id,'nzsoil',nzsoilin,istatus)

    IF ( nxin /= nx .OR. nyin /= ny ) THEN
      WRITE(6,'(1x,a)') 'DTAREADVAR-ERROR: Dimensions in DTAREADVAR inconsistent with data.'
      WRITE(6,'(1x,a,3I15)') ' Read were: ', nxin, nyin, nzsoilin
      WRITE(6,'(1x,a,3I15)') ' Expected:  ', nx, ny, nzsoil
      istatus = -2
      RETURN
    END IF

    CALL hdfrd3d(sd_id,'zpsoil',nx,ny,nzsoil,zpsoil,istatus,itmp,hmax,hmin)

    CALL hdfrd4d(sd_id,"tsoil",nx,ny,nzsoil,nstyps+1,tsoil,istatus,     &
                     itmp,hmax,hmin)

    CALL hdfrd4d(sd_id,"qsoil",nx,ny,nzsoil,nstyps+1,qsoil,istatus,     &
                     itmp,hmax,hmin)

    IF (istatus /= 0) THEN
      WRITE(6,'(1x,3a,I3,a)') 'DTAREADVAR-ERROR: Reading '             &
                  // 'soil variables failed with status = ', istatus, '.'
    END IF

    CALL hdfclose(sd_id,istatus)
    IF (istatus == 0) THEN
      IF(dbglvl > 0)  WRITE(6,'(1x,4a)')                                &
         'DTAREADVAR: Successfully read soil variables from ', trim(filename)
    ELSE
      WRITE(6,'(1x,a,I3,2a)') 'DTAREADVAR-ERROR:(status = ', istatus,   &
                           ') closing ', trim(filename)
    END IF

    DEALLOCATE( itmp, hmax, hmin )
  ELSE
    WRITE(6,'(1x,a,I2,a)') 'ERROR: unsupported file format = ',hinfmt,' in dtareadvar.;'
    istatus = -1
    RETURN
  END IF

  RETURN
END SUBROUTINE dtareadvarsoil

SUBROUTINE get_ext_nzsoil(hinfmt, filename, nzsoil, istatus)

!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!    5/17/2011.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Given an ARPS file name "filename", read variable "varname" with
!  dimensions nx,ny, nzsoil, & nstyps.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!

  INTEGER,          INTENT(IN) :: hinfmt
  CHARACTER(LEN=*), INTENT(IN) :: filename

  INTEGER, INTENT(OUT) :: nzsoil
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: sd_id


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (hinfmt == 3) THEN

    CALL hdfopen(TRIM(filename),1,sd_id)
    IF (sd_id < 0) THEN
      WRITE (6,'(1x,3a)') 'DTAREADVAR-ERROR: opening ',trim(filename),' failed.'
      istatus = -1
      RETURN
    END IF

    CALL hdfrdi(sd_id,'nzsoil',nzsoil,istatus)

    IF (istatus /= 0) THEN
      WRITE(6,'(1x,3a,I3,a)') 'DTAREADVAR-ERROR: Reading ',             &
                   'nzsoil',' failed with status = ', istatus, '.'
    END IF

    CALL hdfclose(sd_id,istatus)

    IF (istatus == 0) THEN
      !IF(dbglvl > 0)  WRITE(6,'(1x,3a)')                                &
      !   'DTAREADVAR: Successfully read nzsoil from ', trim(filename)
    ELSE
      WRITE(6,'(1x,a,I3,2a)') 'DTAREADVAR-ERROR:(status = ', istatus,   &
                           ') closing ', trim(filename)
    END IF

  ELSE
    WRITE(6,'(1x,a,I2,a)') 'ERROR: unsupported file format = ',hinfmt,' in dtareadvar.;'
    istatus = -1
    RETURN
  END IF

  RETURN
END SUBROUTINE get_ext_nzsoil
