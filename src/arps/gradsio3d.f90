!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GRADSREAD                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE gradsread(nx,ny,nz,nzsoil,nstyps,                            &
           nchanl, filnam, time,                                        &
           x,y,z,zp,zpsoil, uprt,vprt,w,ptprt,pprt,                     &
           qvprt,qscalar, tke,kmh,kmv,                                  &
           u,v,pt,p,rhobar,qv,                                          &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,                                &
           raing,rainc,prcrate,                                         &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           usflx,vsflx,ptsflx,qvsflx,                                   &
           ireturn, tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read history data into the file "filnam" in the GrADS format.
!
!       X_center(i) = ( X_edge(i) + X_edge(i+1) )/2,  i = 1, n-1, 1
!
!  The last edge value were kept unchanged so that it can be used to
!  restore the stagger grid point values.
!
!       X_edge(n) = X_center(n)
!       X_edge(i) = 2*X_center(i) - X_edge(i+1),      i = n-1, 1, -1
!
!  Therefore, when you display the data by GrADS, set the x, y, and
!  z dimension not to exceed the Max-1.
!
!  Since total and perturbation variables will be dumped, the base
!  state variables can be obtained by Xbar = X - Xprt.
!
!  The smallest time unit of the current version GrADS is minutes.
!  Therefore, the time intervals for history data dumping, thisdmp,
!  should be multiples of 60 seconds.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  10/25/93.
!
!  MODIFICATIONS:
!
!  10/10/94 (Y. Liu)
!  Added an option (flag istgr) so that the history data can be
!  dumped out at stager points.
!
!  02/06/95 (Y. Liu)
!  Added map projection parameters into the GrADS dumping
!
!  03/27/1995 (Yuhe Liu)
!  Added physical vertical coordinates array, zp(nx,ny,nz) into the
!  display data sets in order for the GrADS to display zp and any
!  other variables (interpolated) in the physical coordinates.
!
!  12/09/1998 (Donghai Wang)
!  Added the snow cover.
!
!  05/14/2002 (J. Brotzge)
!  Changed several variables, call routines to allow for multiple soil schemes.
!
!  6 June 2002 Eric Kemp
!  Soil variable updates.
!
!
! 06/28/2002 Zuwen He
!
! This code has not been tested yet.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of grid points in the soil
!
!  OUTPUT:
!
!    nchanl   FORTRAN I/O channel number for history data output.
!    filnam   The name of history data dump file
!    time     Model starting time
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!    zpsoil   Vertical coordinate of grid points in soil (m)
!
!    uprt     x component of velocity (m/s)
!    vprt     y component of velocity (m/s)
!    ptprt    Perturbation potential temperature (K)
!    pprt     Perturbation pressure (Pascal)
!    qvprt    Perturbation water vapor specific humidity (kg/kg)
!
!    qc       Cloud water mixing ratio at a given time level (kg/kg)
!    qr       Rainwater mixing ratio at a given time level (kg/kg)
!    qi       Cloud ice mixing ratio at a given time level (kg/kg)
!    qs       Snow mixing ratio at a given time level (kg/kg)
!    qh       Hail mixing ratio at a given time level (kg/kg)
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!
!    u        x component of velocity (m/s)
!    v        y component of velocity (m/s)
!    w        Vertical component of Cartesian velocity (m/s)
!    pt       Potential temperature (K)
!    p        Pressure (Pascal)
!    rhobar   Base state density (kg/m**3)
!    qv       Water vapor specific humidity (kg/kg)
!
!    soiltyp  Soil type
!    vegtyp   Vegetation type
!    lai      Leaf Area Index
!    roufns   Surface roughness
!    veg      Vegetation fraction
!
!    tsoil    Soil temperature (K)
!    qsoil    Soil moisture (m**3/m**3)
!    wetcanp  Canopy water amount (m)
!
!    raing    Grid supersaturation rain
!    rainc    Cumulus convective rain
!    prcrate  Precipitation rates
!
!    radfrc   Radiation forcing (K/s)
!    radsw    Solar radiation reaching the surface
!    rnflx    Net radiation flux absorbed by surface
!    radswnet Net shortwave radiation
!    radlwin  Incoming longwave radiation
!
!    usflx    Surface flux of u-momentum (kg/(m*s**2))
!    vsflx    Surface flux of v-momentum (kg/(m*s**2))
!    ptsflx   Surface heat flux (K*kg/(m**2 * s ))
!    qvsflx   Surface moisture flux of (kg/(m**2 * s))
!
!    ireturn  Return status indicator
!
!  WORK ARRAY:
!
!    tem1     Temporary work array.
!
!-----------------------------------------------------------------------
!
!  The following parameters are passed into this subroutine through
!  a common block in globcst.inc, and they determine which
!  variables are output.
!
!  varout =0 or 1. If varout=0, dynamical variables are not dumped.
!  mstout =0 or 1. If mstout=0, water variables are not dumped.
!  rainout=0 or 1. If rainout=0, rain variables are not dumped.
!  prcout =0 or 1. If prcout=0, precipitation rates are not dumped.
!  iceout =0 or 1. If iceout=0, qi, qs and qh are not dumped.
!  tkeout =0 or 1. If tkeout=0, tke is not dumped.
!  trbout =0 or 1. If trbout=0, kmh and kmv are not dumped.
!  sfcout =0 or 1. If sfcout=0, surface variables are not dumped.
!  landout=0 or 1. If landout=0, surface property arrays are not dumped.
!  radout =0 or 1. If radout =0, radiation arrays are not dumped.
!  flxout =0 or 1. If flxout =0, surface flux arrays are not dumped.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'phycst.inc'
  INCLUDE 'indtflg.inc'

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of grid points in the soil

  INTEGER :: nchanl            ! FORTRAN I/O channel number for output

  REAL :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz)     ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz)     ! Total w-velocity (m/s)
  REAL :: pt    (nx,ny,nz)     ! Potential temperature (K)
  REAL :: p     (nx,ny,nz)     ! Pressure (Pascal)
  REAL :: qv    (nx,ny,nz)     ! Water vapor specific humidity (kg/kg)

  REAL :: qscalar(nx,ny,nz,nscalar)

  REAL :: uprt  (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL :: vprt  (nx,ny,nz)     ! Total v-velocity (m/s)
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)
  REAL :: qvprt (nx,ny,nz)     ! Water vapor specific humidity (kg/kg)
  REAL :: tke   (nx,ny,nz)     ! Turbulent Kinetic Energy ((m/s)**2)

  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )

  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
  REAL :: zpsoil(nx,ny,nzsoil) ! The physical height coordinate defined at
                               ! w-point of the soil

  INTEGER :: nstyps
  INTEGER :: soiltyp (nx,ny,nstyps)  ! Soil type
  REAL    :: stypfrct(nx,ny,nstyps)  ! Soil type fraction
  INTEGER :: vegtyp (nx,ny)          ! Vegetation type
  REAL :: lai    (nx,ny)          ! Leaf Area Index
  REAL :: roufns (nx,ny)          ! Surface roughness
  REAL :: veg    (nx,ny)          ! Vegetation fraction

  REAL :: tsoil  (nx,ny,nzsoil,0:nstyps) ! Soil temperature (K)
  REAL :: qsoil  (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny,0:nstyps)    ! Canopy water amount (m)
  REAL :: snowdpth(nx,ny)            ! Snow depth (m)

  REAL :: raing(nx,ny)         ! Grid supersaturation rain
  REAL :: rainc(nx,ny)         ! Cumulus convective rain
  REAL :: prcrate(nx,ny,4)     ! precipitation rate (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precip. rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulus precip. rate
                               ! prcrate(1,1,4) = microphysics precip. rate

  REAL :: radfrc(nx,ny,nz)     ! Radiation forcing (K/s)
  REAL :: radsw (nx,ny)        ! Solar radiation reaching the surface
  REAL :: rnflx (nx,ny)        ! Net radiation flux absorbed by surface
  REAL :: radswnet(nx,ny)      ! Net shortwave radiation
  REAL :: radlwin(nx,ny)       ! Incoming longwave radiation

  REAL :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  REAL :: tem1(nx,ny,nz)       ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Parameters describing this routine
!
!-----------------------------------------------------------------------
!
! 06/28/2002 Zuwen He
!
! fmtver??: to label each data a version.
! intver??: an integer to allow faster comparison than fmtver??,
!           which are strings.
!
! Verion 5.00: significant change in soil variables since version 4.10.
!
  CHARACTER (LEN=40) :: fmtver500,fmtver530
  INTEGER  :: intver,intver500,intver530
  PARAMETER (fmtver500='005.00 GrADS Binary Data',intver500=500)
  PARAMETER (fmtver530='005.30 GrADS Binary Data',intver530=530)
  CHARACTER (LEN=40) :: fmtverin

  CHARACTER (LEN=10) :: tmunit
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL    :: time                ! Model starting time
  INTEGER :: nxin,nyin,nzin      ! # of 3-D grid points in data file
  INTEGER :: nzsoilin            ! # of soil levels
  INTEGER :: nstypin             ! # of soil types
  INTEGER :: ireturn, ierr
  INTEGER :: i,j,k,l, istat,is
  INTEGER :: nq, nqin
  INTEGER :: nqscalarin(nscalar)
  INTEGER :: idummy
  REAL    :: rdummy,time0
  INTEGER :: filen
  CHARACTER (LEN=*       ) :: filnam ! File name of data file
  INTEGER :: varnum              ! # of variables in data file
  CHARACTER (LEN=60) :: vartit(50)
  CHARACTER (LEN=8)  :: varnam(50)
  INTEGER :: hdbyte, lev         ! # of bytes in the header of data file
  INTEGER :: istgr
  INTEGER :: nstyp1

  LOGICAL :: firstcall
  SAVE firstcall,time0
  DATA firstcall/.true./

!
!-----------------------------------------------------------------------
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  IF (firstcall) THEN

    CALL getunit( nchanl )

    filen = LEN(filnam)
    CALL strlnth( filnam, filen )

    CALL asnctl ('OLDLOCAL', 1, ierr)
    CALL asnfile(filnam(1:filen), '-F f77 -N ieee', ierr)

    OPEN (UNIT=nchanl,FILE=trim(filnam(1:filen)),STATUS='old',          &
          FORM='unformatted',ACCESS='sequential',IOSTAT= istat )

    hdbyte = 0

    READ (nchanl,ERR=910,END=920) fmtverin
    hdbyte = hdbyte + 8 + 1*40

    IF( fmtverin == fmtver500 ) THEN
      intver = intver500
    ELSE IF( fmtverin == fmtver530 ) THEN
      intver = intver530
    ELSE
      WRITE(6,'(/1x,a,/1x,2a,/1x,3a)')                                  &
          'Data format incompatible with the data reader.',             &
          'Format of data is ',fmtverin,' Format of reader is ',fmtver530, &
          '. Job stopped.'
      CALL arpsstop('arpsstop called from gradsread at initial read',1)
    END IF

    READ (nchanl,ERR=910,END=920) runname
    hdbyte = hdbyte + 8 + 1*80

    READ (nchanl,ERR=910,END=920) nocmnt
    hdbyte = hdbyte + 8 + 1*4

    IF ( nocmnt > 0 ) THEN
      DO l = 1, nocmnt
        READ (nchanl,ERR=910,END=920) cmnt(l)
      END DO
      hdbyte = hdbyte + nocmnt * ( 8 + 80 )
    END IF

    READ (nchanl,ERR=910,END=920) nxin,nyin,nzin,nzsoilin,nstypin
    hdbyte = hdbyte + 8 + 5*4

    IF ( nxin /= nx .OR. nyin /= ny .OR. nzin /= nz .OR. &
         nzsoilin /= nzsoil) THEN
      WRITE (6, '(1x,a/,8(1X,a,i4/),1X,a)' )                         &
          ' Dimensions in GRADSREAD inconsistent with data.',         &
          ' nx        = ',nx,          &
          ' ny        = ',ny,          &
          ' nz        = ',nz,          &
          ' nzsoil    = ',nzsoil,      &
          ' nxin      = ',nxin,        &
          ' nyin      = ',nyin,        &
          ' nzin      = ',nzin,        &
          ' nzsoilin  = ',nzsoilin,    &
          ' Program aborted in GRADSREAD.'
       CALL arpsstop('Stopped at calling gradsread while reading nx..',1)
    END IF

    READ (nchanl,ERR=910,END=920) time0,tmunit
    hdbyte = hdbyte + 8 + 1*4 + 10

    READ (nchanl,ERR=910,END=920)                                       &
                   varin,  mstin,  icein,  trbin,  sfcin,               &
                  rainin, landin, idummy, idummy,  totin,               &
                  tkein ,nscalarin,mapproj, istgr, month,               &
                     day,   year,   hour, minute, second
    hdbyte = hdbyte + 8 + 20*4

    IF (intver >= intver530) THEN

      READ (nchanl,ERR=910,END=920)                                     &
                 p_qcin,p_qrin,p_qiin,p_qsin,p_qgin,p_qhin,             &
                 p_ncin,p_nrin,p_niin,p_nsin,p_ngin,p_nhin,             &
                 p_zrin,p_ziin,p_zsin,p_zgin,p_zhin,idummy,             &
                 idummy,idummy,idummy,idummy,idummy,idummy,             &
                 idummy,idummy,idummy,idummy,idummy,idummy
      hdbyte = hdbyte + 8 + 30*4

    ELSE
      p_qcin=0; p_qrin=0; p_qiin=0; p_qsin=0; p_qgin=0; p_qhin=0
      p_ncin=0; p_nrin=0; p_niin=0; p_nsin=0; p_nhin=0; p_ngin=0
                p_zrin=0; p_ziin=0; p_zsin=0; p_zhin=0; p_zgin=0

      nscalarin = 0
      IF (mstin == 1) THEN
        p_qcin = 1
        p_qrin = 2
        nscalarin = nscalarin + 2
        IF (icein == 1) THEN
          p_qiin = 3
          p_qsin = 4
          p_qhin = 5
          nscalarin = nscalarin + 3
        END IF
      END IF

    END IF

    nqscalarin(:) = 0
    IF (P_QC > 0) nqscalarin(P_QC) = p_qcin
    IF (P_QR > 0) nqscalarin(P_QR) = p_qrin
    IF (P_QI > 0) nqscalarin(P_QI) = p_qiin
    IF (P_QS > 0) nqscalarin(P_QS) = p_qsin
    IF (P_QG > 0) nqscalarin(P_QG) = p_qgin
    IF (P_QH > 0) nqscalarin(P_QH) = p_qhin

    IF (P_NC > 0) nqscalarin(P_NC) = p_ncin
    IF (P_NR > 0) nqscalarin(P_NR) = p_nrin
    IF (P_NI > 0) nqscalarin(P_NI) = p_niin
    IF (P_NS > 0) nqscalarin(P_NS) = p_nsin
    IF (P_NG > 0) nqscalarin(P_NG) = p_ngin
    IF (P_NH > 0) nqscalarin(P_NH) = p_nhin

    IF (P_ZR > 0) nqscalarin(P_ZR) = p_nrin
    IF (P_ZI > 0) nqscalarin(P_ZI) = p_niin
    IF (P_ZS > 0) nqscalarin(P_ZS) = p_nsin
    IF (P_ZG > 0) nqscalarin(P_ZG) = p_ngin
    IF (P_ZH > 0) nqscalarin(P_ZH) = p_nhin

    IF (intver < intver530) THEN
      READ (nchanl,ERR=910,END=920)                                     &
                     umove,   vmove, xgrdorg, ygrdorg, trulat1,         &
                   trulat2,  trulon,  sclfct,  n0rain,  n0snow,         &
                    n0hail, rhosnow, rhohail,  rdummy,  rdummy,         &
                     tstop, thisdmp, latitud,  ctrlat,  ctrlon
      hdbyte = hdbyte + 8 + 20*4
    ELSE
      READ (nchanl,ERR=910,END=920)                                     &
                     umove,   vmove, xgrdorg, ygrdorg, trulat1,         &
                   trulat2,  trulon,  sclfct,  ntcloud, n0rain,         &
                    n0snow,  n0grpl,  n0hail,  rhoice,  rhosnow,        &
                    rhogrpl, rhohail, alpharain, alphaice, alphasnow,   &
                    alphagrpl, alphahail,rdummy, rdummy,  rdummy,       &
                     tstop, thisdmp, latitud,  ctrlat,  ctrlon
      hdbyte = hdbyte + 8 + 30*4
    END IF

    IF ( totin /= 0 ) THEN
      READ (nchanl,ERR=910,END=920)                                     &
                  hdmpopt,nstyp1,  prcin,  radin,  flxin,               &
                  snowcin,snowin, idummy, idummy, idummy,               &
                  idummy, idummy, idummy, idummy, idummy,               &
                  idummy, idummy, idummy, idummy, idummy
      hdbyte = hdbyte + 8 + 20*4
!
! 06/18/2002 Zuwen
!
! NOTE: nstyp1 should be at the same value as nstypin, according to
! the way the data is dumped.
! I prefer nstypin, because read-in of nstyp1 is missed if totin=0.
!
!     IF ( nstyp1 < 1 ) THEN
!       nstyp1 = 1
!     END IF

      READ (nchanl,ERR=910,END=920)                                     &
                  tstrtdmp,rdummy, rdummy, rdummy, rdummy,              &
                  rdummy, rdummy, rdummy, rdummy, rdummy,               &
                  rdummy, rdummy, rdummy, rdummy, rdummy,               &
                  rdummy, rdummy, rdummy, rdummy, rdummy
      hdbyte = hdbyte + 8 + 20*4

    END IF

    IF ( hdmpopt == 2 ) THEN
      READ (nchanl,ERR=910,END=920) numhdmp
      hdbyte = hdbyte + 8 + 4
      IF ( numhdmp > 0 ) THEN
        READ (nchanl,ERR=910,END=920) (hdmptim(i),i=1,numhdmp)
        hdbyte = hdbyte + 8 + numhdmp*4
      END IF
    END IF

    READ (nchanl,ERR=910,END=920) x
    hdbyte = hdbyte + 8 + nx*4

    READ (nchanl,ERR=910,END=920) y
    hdbyte = hdbyte + 8 + ny*4

    READ (nchanl,ERR=910,END=920) z
    hdbyte = hdbyte + 8 + nz*4

    READ (nchanl,ERR=910,END=920) varnum
    hdbyte = hdbyte + 8 + 1*4

    WRITE (6, '(a//a/)' )                                               &
        '     The following variables will be read from the data file', &
        '     Variables      Description'

    DO l = 1, varnum
      READ (nchanl,ERR=910,END=920) varnam(l),lev,vartit(l)
      WRITE (6, '(5x,a8,7x,a60)' ) varnam(l),vartit(l)
    END DO
    hdbyte = hdbyte + varnum*( 8 + 8 + 4 + 60 )

    IF(landin == 1) THEN

      DO is=1,nstypin
        IF (is <= nstyps) THEN
          READ (nchanl,ERR=910,END=920)                                 &
                 ((soiltyp (i,j,is),i=1,nx),j=1,ny)
          READ (nchanl,ERR=910,END=920)                                 &
                 ((stypfrct(i,j,is),i=1,nx),j=1,ny)
        ELSE
          READ (nchanl,ERR=910,END=920)
          READ (nchanl,ERR=910,END=920)
        ENDIF
      END DO

      CALL fix_stypfrct_nstyp(nx,ny,nstypin,nstyp,stypfrct)

      READ (nchanl,ERR=910,END=920) ((vegtyp (i,j),i=1,nx),j=1,ny)
      READ (nchanl,ERR=910,END=920) ((lai    (i,j),i=1,nx),j=1,ny)
      READ (nchanl,ERR=910,END=920) ((roufns (i,j),i=1,nx),j=1,ny)
      READ (nchanl,ERR=910,END=920) ((veg    (i,j),i=1,nx),j=1,ny)
      hdbyte = hdbyte + 5 * ( 8 + nx*ny*4 )
    END IF

    nhisdmp = 1
    firstcall = .false.

    time = time0

  ELSE IF ( hdmpopt == 1 ) THEN     ! not firstcall

    IF ( tstrtdmp > time0 ) THEN
      time = tstrtdmp
    ELSE
      time = tstrtdmp + nhisdmp*thisdmp
    END IF
    nhisdmp = nhisdmp + 1

  ELSE IF ( hdmpopt == 2 ) THEN     ! not firstcall and not hdmpopt=1

    nhisdmp = nhisdmp + 1
    IF ( nhisdmp <= numhdmp ) THEN
      time = hdmptim(nhisdmp)
    ELSE
      WRITE (6,'(a)') 'No more data in this file'
      GO TO 920
    END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  Read in data from the GrADS data file
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Read the physical vertical coordinates into the data set.
!
!-----------------------------------------------------------------------
!
  DO k=1,nz
    READ (nchanl,ERR=910,END=920) ((tem1(i,j,k),i=1,nx),j=1,ny)
  END DO

  IF ( istgr == 0 ) THEN
    DO i=1,nx
      DO j=1,ny
        zp(i,j,nz) = tem1(i,j,nz)
        DO k=nz-1,1,-1
          zp(i,j,k) = 2.*tem1(i,j,k) - zp(i,j,k+1)
        END DO
      END DO
    END DO
  ELSE
    DO k=1,nz
      DO i=1,nx
        DO j=1,ny
          zp(i,j,k) = tem1(i,j,k)
        END DO
      END DO
    END DO
  END IF

  IF ( varin == 1 ) THEN
!
!-----------------------------------------------------------------------
!
!  If varin = 1, read in u, uprt, v, vprt, w, wprt,
!  pt, ptprt, p, pprt, vort, and div.
!
!-----------------------------------------------------------------------
!
    DO k=1,nz
      READ (nchanl,ERR=910,END=920) ((tem1(i,j,k),i=1,nx),j=1,ny)
    END DO

    IF ( istgr == 0 ) THEN
      DO j=1,ny
        DO k=1,nz
          u(nx,j,k) = tem1(nx,j,k)
          DO i=nx-1,1,-1
            u(i,j,k) = 2.*tem1(i,j,k) - u(i+1,j,k)
          END DO
        END DO
      END DO
    ELSE
      DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            u(i,j,k) = tem1(i,j,k)
          END DO
        END DO
      END DO
    END IF

    DO k=1,nz
      READ (nchanl,ERR=910,END=920) ((tem1(i,j,k),i=1,nx),j=1,ny)
    END DO

    IF ( istgr == 0 ) THEN
      DO j=1,ny
        DO k=1,nz
          uprt(nx,j,k) = tem1(nx,j,k)
          DO i=nx-1,1,-1
            uprt(i,j,k) = 2.*tem1(i,j,k) - uprt(i+1,j,k)
          END DO
        END DO
      END DO
    ELSE
      DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            uprt(i,j,k) = tem1(i,j,k)
          END DO
        END DO
      END DO
    END IF

    DO k=1,nz
      READ (nchanl,ERR=910,END=920) ((tem1(i,j,k),i=1,nx),j=1,ny)
    END DO

    IF ( istgr == 0 ) THEN
      DO i=1,nx
        DO k=1,nz
          v(i,ny,k) = tem1(i,ny,k)
          DO j=ny-1,1,-1
            v(i,j,k) = 2.*tem1(i,j,k) - v(i,j+1,k)
          END DO
        END DO
      END DO
    ELSE
      DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            v(i,j,k) = tem1(i,j,k)
          END DO
        END DO
      END DO
    END IF

    DO k=1,nz
      READ (nchanl,ERR=910,END=920) ((tem1(i,j,k),i=1,nx),j=1,ny)
    END DO

    IF ( istgr == 0 ) THEN
      DO i=1,nx
        DO k=1,nz
          vprt(i,ny,k) = tem1(i,ny,k)
          DO j=ny-1,1,-1
            vprt(i,j,k) = 2.*tem1(i,j,k) - vprt(i,j+1,k)
          END DO
        END DO
      END DO
    ELSE
      DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            vprt(i,j,k) = tem1(i,j,k)
          END DO
        END DO
      END DO
    END IF

    DO k=1,nz
      READ (nchanl,ERR=910,END=920) ((tem1(i,j,k),i=1,nx),j=1,ny)
    END DO

    IF ( istgr == 0 ) THEN
      DO i=1,nx
        DO j=1,ny
          w(i,j,nz) = tem1(i,j,nz)
          DO k=nz-1,1,-1
            w(i,j,k) = 2.*tem1(i,j,k) - w(i,j,k+1)
          END DO
        END DO
      END DO
    ELSE
      DO k=1,nz
        DO i=1,nx
          DO j=1,ny
            w(i,j,k) = tem1(i,j,k)
          END DO
        END DO
      END DO
    END IF

    DO k=1,nz
      READ (nchanl,ERR=910,END=920) ((pt  (i,j,k),i=1,nx),j=1,ny)
    END DO
    DO k=1,nz
      READ (nchanl,ERR=910,END=920) ((ptprt(i,j,k),i=1,nx),j=1,ny)
    END DO
    DO k=1,nz
      READ (nchanl,ERR=910,END=920) ((p   (i,j,k),i=1,nx),j=1,ny)
    END DO
    DO k=1,nz
      READ (nchanl,ERR=910,END=920) ((pprt(i,j,k),i=1,nx),j=1,ny)
    END DO

!-----------------------------------------------------------------------
!
!  Vorticity:
!
!-----------------------------------------------------------------------

    DO k=1,nz
      READ (nchanl,ERR=910,END=920) ((tem1(i,j,k),i=1,nx),j=1,ny)
    END DO

!-----------------------------------------------------------------------
!
!  Divergernce:
!
!-----------------------------------------------------------------------

    DO k=1,nz
      READ (nchanl,ERR=910,END=920) ((tem1(i,j,k),i=1,nx),j=1,ny)
    END DO

  END IF       ! End varin

  IF ( mstin == 1 ) THEN
!
!-----------------------------------------------------------------------
!
!  Read in moist variables qv, qvprt, qc, qr, qi, qs, and qh
!
!-----------------------------------------------------------------------
!
    DO k=1,nz
      READ (nchanl,ERR=910,END=920) ((qv   (i,j,k),i=1,nx),j=1,ny)
    END DO
    DO k=1,nz
      READ (nchanl,ERR=910,END=920) ((qvprt(i,j,k),i=1,nx),j=1,ny)
    END DO

    DO nqin = 1,nscalarin
      DO k=1,nz
        READ (nchanl,ERR=910,END=920) ((tem1(i,j,k),i=1,nx),j=1,ny)
      END DO
      DO nq = 1, nscalar
        IF (nqin == nqscalarin(nq)) THEN
          qscalar(:,:,:,nq) = tem1(:,:,:)
          EXIT
        END IF
      END DO
    END DO

    IF ( rainin == 1 ) THEN

      READ (nchanl,ERR=910,END=920) ((raing(i,j),i=1,nx),j=1,ny)
      READ (nchanl,ERR=910,END=920) ((rainc(i,j),i=1,nx),j=1,ny)

    END IF   ! End rainin

    IF ( prcin == 1 ) THEN

      READ (nchanl,ERR=910,END=920) ((prcrate(i,j,1),i=1,nx),j=1,ny)
      READ (nchanl,ERR=910,END=920) ((prcrate(i,j,2),i=1,nx),j=1,ny)
      READ (nchanl,ERR=910,END=920) ((prcrate(i,j,3),i=1,nx),j=1,ny)
      READ (nchanl,ERR=910,END=920) ((prcrate(i,j,4),i=1,nx),j=1,ny)

    END IF   ! End prcin

  END IF       ! End mstin

  IF ( tkein == 1 ) THEN

    DO k=1,nz
      READ (nchanl,ERR=910,END=920) ((tke(i,j,k),i=1,nx),j=1,ny)
    END DO

  END IF

  IF ( trbin == 1 ) THEN
!
!-----------------------------------------------------------------------
!
!  If trbin = 1, read in the turbulence parameter, km.
!
!-----------------------------------------------------------------------
!
    DO k=1,nz
      READ (nchanl,ERR=910,END=920) ((kmh(i,j,k),i=1,nx),j=1,ny)
    END DO
    DO k=1,nz
      READ (nchanl,ERR=910,END=920) ((kmv(i,j,k),i=1,nx),j=1,ny)
    END DO

  END IF       ! trbin

  IF ( sfcin == 1 ) THEN

    DO k=1,nzsoil
      READ (nchanl,ERR=910,END=920) ((zpsoil(i,j,k),i=1,nx),j=1,ny)
    END DO

    DO is=0,nstypin
      IF (is <= nstyps) THEN
        DO k=1,nzsoil
          READ (nchanl,ERR=910,END=920) ((tsoil(i,j,k,is),i=1,nx),j=1,ny)
        END DO
        DO k=1,nzsoil
          READ (nchanl,ERR=910,END=920) ((qsoil(i,j,k,is),i=1,nx),j=1,ny)
        END DO
        READ (nchanl,ERR=910,END=920) ((wetcanp(i,j,is),i=1,nx),j=1,ny)
      ELSE
        DO k=1,nz
          READ (nchanl,ERR=910,END=920)
        END DO
        DO k=1,nz
          READ (nchanl,ERR=910,END=920)
        END DO
        DO k=1,nz
          READ (nchanl,ERR=910,END=920)
        END DO
      ENDIF

    END DO

    CALL fix_soil_nstyp(nx,ny,nzsoil,nstypin,nstyp,tsoil,qsoil,wetcanp)

    IF (snowin == 1) THEN
      READ (nchanl,ERR=910,END=920) ((snowdpth(i,j),i=1,nx),j=1,ny)
    END IF

  END IF

  IF ( radin == 1 ) THEN

    DO k=1,nz
      READ (nchanl,ERR=910,END=920) ((radfrc(i,j,k),i=1,nx),j=1,ny)
    END DO
    READ (nchanl,ERR=910,END=920) ((radsw(i,j),i=1,nx),j=1,ny)
    READ (nchanl,ERR=910,END=920) ((rnflx(i,j),i=1,nx),j=1,ny)
    READ (nchanl,ERR=910,END=920) ((radswnet(i,j),i=1,nx),j=1,ny)
    READ (nchanl,ERR=910,END=920) ((radlwin(i,j),i=1,nx),j=1,ny)

  END IF   ! End radin

  IF ( flxin == 1 ) THEN

    READ (nchanl,ERR=910,END=920) ((usflx(i,j),i=1,nx),j=1,ny)
    READ (nchanl,ERR=910,END=920) ((vsflx(i,j),i=1,nx),j=1,ny)
    READ (nchanl,ERR=910,END=920) ((ptsflx(i,j),i=1,nx),j=1,ny)
    READ (nchanl,ERR=910,END=920) ((qvsflx(i,j),i=1,nx),j=1,ny)

  END IF   ! End flxin

  ireturn=0
  RETURN
!
!-----------------------------------------------------------------------
!
!  Error during read
!
!-----------------------------------------------------------------------
!

  910   CONTINUE
  WRITE(6,'(/a/)') ' Error reading data in GRADSREAD'
  ireturn=1
  RETURN
!
!-----------------------------------------------------------------------
!
!  End-of-file during read
!
!----------------------------------------------------------------------
!

  920   CONTINUE
  WRITE(6,'(/a/)') ' End of file reached in GRADSREAD'
  ireturn=2
END SUBROUTINE gradsread
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GRADSDUMP                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE gradsdump(nx,ny,nz,nzsoil,nstyps,nchanl,filnam,istgr,        &
           u,v,w,ptprt,pprt,qv,qscalar,tke,kmh,kmv,                     &
           ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,                      &
           x,y,z,zp,zpsoil,                                             &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,                                &
           raing,rainc,prcrate,                                         &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           usflx,vsflx,ptsflx,qvsflx,                                   &
           tem1,tem2)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write history data into the file "filnam" in the GrADS format.
!
!  All output data are located at the grid cell centers.
!
!       X_center(i) = ( X_edge(i) + X_edge(i+1) )/2,  i = 1, n-1, 1
!
!  The last edge value were kept unchanged so that it can be used to
!  restore the stagger grid point values.
!
!       X_edge(n) = X_center(n)
!       X_edge(i) = 2*X_center(i) - X_edge(i+1),      i = n-1, 1, -1
!
!  Therefore, when you display the data by GrADS, set the x, y, and
!  z dimension not to exceed the Max-1.
!
!  Since total and perturbation variables will be dumped, the base
!  state variables can be obtained by Xbar = X - Xprt.
!
!  The smallest time unit of the current version GrADS is minutes.
!  Therefore, the time intervals for history data dumping, thisdmp,
!  should be multiples of 60 seconds.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  7/20/93.
!
!
!  MODIFICATION HISTORY:
!
!  9/1/94 (Y. Lu)
!  Cleaned up documentation.
!
!  10/10/94 (Y. Liu)
!  Added an option (flag istgr) so that the history data can be
!  dumped out at stager points.
!
!  02/06/95 (Y. Liu)
!  Added map projection parameters into the GrADS dumping
!
!  03/27/1995 (Yuhe Liu)
!  Added physical vertical coordinates array, zp(nx,ny,nz) into the
!  display data sets in order for the GrADS to display zp and any
!  other variables (interpolated) in the physical coordinates.
!
!  6 June 2002 Eric Kemp
!  Soil variable updates.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of grid points in the vertical
!
!    nchanl   FORTRAN I/O channel number for history data output.
!
!    u        x component of velocity at a given time level (m/s)
!    v        y component of velocity at a given time level (m/s)
!    w        Vertical component of Cartesian velocity at a given
!             time level (m/s)
!    ptprt    Perturbation potential temperature at a given time
!             level (K)
!    pprt     Perturbation pressure at a given time level (Pascal)
!    qv       Water vapor specific humidity at a given time level
!             (kg/kg)
!    qc       Cloud water mixing ratio at a given time level (kg/kg)
!    qr       Rainwater mixing ratio at a given time level (kg/kg)
!    qi       Cloud ice mixing ratio at a given time level (kg/kg)
!    qs       Snow mixing ratio at a given time level (kg/kg)
!    qh       Hail mixing ratio at a given time level (kg/kg)
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!
!    ubar     Base state zonal velocity component (m/s)
!    vbar     Base state meridional velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhobar   Base state density (kg/m**3)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space(m)
!    zpsoil   Vertical coordinate of grid points in physical space(m)
!
!    soiltyp  Soil type
!    vegtyp   Vegetation type
!    lai      Leaf Area Index
!    roufns   Surface roughness
!    veg      Vegetation fraction
!
!    tsoil    Soil temperature (K)
!    qsoil    Soil moisture (m**3/m**3)
!    wetcanp  Canopy water amount (m)
!
!    raing    Grid supersaturation rain
!    rainc    Cumulus convective rain
!    prcrate  Precipitation rates
!
!    radfrc   Radiation forcing (K/s)
!    radsw    Solar radiation reaching the surface
!    rnflx    Net radiation flux absorbed by surface
!
!    usflx    Surface flux of u-momentum (kg/(m*s**2))
!    vsflx    Surface flux of v-momentum (kg/(m*s**2))
!    ptsflx   Surface heat flux (K*kg/(m**2 * s ))
!    qvsflx   Surface moisture flux of (kg/(m**2 * s))
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAY:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!
!
!-----------------------------------------------------------------------
!
!  The following parameters are passed into this subroutine through
!  a common block in globcst.inc, and they determine which
!  variables are output.
!
!  varout =0 or 1. If varout=0, model perturbation variables are not
!                  dumped.
!  mstout =0 or 1. If mstout =0, water variables are not dumped.
!  rainout=0 or 1. If rainout=0, rain variables are not dumped.
!  prcout =0 or 1. If prcout =0, precipitation rates are not dumped.
!  iceout =0 or 1. If iceout =0, qi, qs and qh are not dumped.
!  tkeout =0 or 1. If tkeout =0, tke is not dumped.
!  trbout =0 or 1. If trbout =0, kmh and kmv are not dumped
!  sfcout =0 or 1. If sfcout =0, surface variables are not dumped.
!  landout=0 or 1. If landout=0, surface property arrays are not dumped.
!  radout =0 or 1. If radout =0, radiation arrays are not dumped.
!  flxout =0 or 1. If flxout =0, surface fluxes are not dumped.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'phycst.inc'

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of grid points in the soil

  INTEGER :: istgr             ! Flag for dumping stager point data

  INTEGER :: nchanl            ! FORTRAN I/O channel number for output

  REAL :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz)     ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz)     ! Total w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz)     ! Water vapor specific humidity (kg/kg)

  REAL :: qscalar    (nx,ny,nz,nscalar)

  REAL :: tke   (nx,ny,nz)     ! Turbulent Kinetic Energy ((m/s)**2)
  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: wbar  (nx,ny,nz)     ! Base state w-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific
                               ! humidity(kg/kg)

  REAL :: x     (nx)           ! x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! Physical height coordinate defined at
                               ! w-point of the staggered grid.
  REAL :: zpsoil(nx,ny,nzsoil) ! Physical height coordinate defined at
                               ! w-point of the soil.

  INTEGER :: nstyps
  INTEGER :: soiltyp(nx,ny,nstyps)    ! Soil type
  REAL :: stypfrct(nx,ny,nstyps)   ! Fraction of soil types
  INTEGER :: vegtyp (nx,ny)           ! Vegetation type
  REAL :: lai    (nx,ny)           ! Leaf Area Index
  REAL :: roufns (nx,ny)           ! Surface roughness
  REAL :: veg    (nx,ny)           ! Vegetation fraction

  REAL :: tsoil  (nx,ny,nzsoil,0:nstyps) ! Soil temperature (K)
  REAL :: qsoil  (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny,0:nstyps)     ! Canopy water amount (m)
  REAL :: snowdpth(nx,ny)             ! Snow depth (m)

  REAL :: raing(nx,ny)                ! Grid supersaturation rain
  REAL :: rainc(nx,ny)                ! Cumulus convective rain
  REAL :: prcrate(nx,ny,4)     ! precipitation rate (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precip. rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulus precip. rate
                               ! prcrate(1,1,4) = microphysics precip. rate

  REAL :: radfrc(nx,ny,nz)     ! Radiation forcing (K/s)
  REAL :: radsw (nx,ny)        ! Solar radiation reaching the surface
  REAL :: rnflx (nx,ny)        ! Net radiation flux absorbed by surface
  REAL :: radswnet(nx,ny)      ! Net shortwave radiation
  REAL :: radlwin(nx,ny)       ! Incoming longwave radiation


  REAL :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  REAL :: tem1  (nx,ny,nz)            ! Temporary work array
  REAL :: tem2  (nx,ny,nz)            ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Parameters describing this routine
!
!-----------------------------------------------------------------------
!
! 06/28/2002 Zuwen He
!
! fmtver??: to label each data a version.
! intver??: an integer to allow faster comparison than fmtver??,
!           which are strings.
!
! Verion 5.00: significant change in soil variables since version 4.10.
!
  CHARACTER (LEN=40) :: fmtver,fmtver500,fmtver530
  INTEGER  :: intver,intver500,intver530
  PARAMETER (fmtver500='005.00 GrADS Binary Data',intver500=500)
  PARAMETER (fmtver530='005.30 GrADS Binary Data',intver530=530)

  CHARACTER (LEN=10) :: tmunit
  PARAMETER (tmunit='seconds   ')
  CHARACTER (LEN=2) :: dtunit
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: varnumax
  PARAMETER ( varnumax = 100 )
  CHARACTER (LEN=8) :: varnam(varnumax)
  CHARACTER (LEN=60) :: vartit(varnumax)
  CHARACTER (LEN=10) :: varparam(varnumax)
  INTEGER :: varlev(varnumax)
  INTEGER :: i,j,k,l,m, istat, ierr,is
  INTEGER :: nq
  INTEGER :: varnum
  INTEGER :: tinc,ntm
  INTEGER :: idummy
  REAL :: rdummy
  REAL :: xbgn,ybgn,zbgn, xinc,yinc,zinc
  REAL :: lat11, lon11
  REAL :: latmin, latmax, lonmin, lonmax, latinc, loninc
  CHARACTER (LEN=*) :: filnam
  CHARACTER (LEN=80) :: gradscntl
  INTEGER :: filen, cntlen
  CHARACTER (LEN=20) :: chrstr
  INTEGER :: hdbyte
  INTEGER :: nchout0

  INTEGER :: year1, month1, day1, jday1, loopdy
  INTEGER :: hour1, minute1, second1

  CHARACTER (LEN=3) :: monnam(12)            ! Name of months
  DATA monnam/'jan', 'feb', 'mar', 'apr', 'may', 'jun',                 &
              'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/

  INTEGER :: mndys(12)                 ! days for each months
  DATA mndys/0,31,59,90,120,151,181,212,243,273,304,334/

  LOGICAL :: firstcall
  DATA firstcall/.true./
!
!-----------------------------------------------------------------------
!
  SAVE firstcall
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  intver = intver530  !  for the time being, in the future, we will
                      !  allow to dump data in the different version
                      !  intver will be assigned from input file

  IF (intver == intver500) THEN
    fmtver=fmtver500
  ELSE IF (intver == intver530) THEN
    fmtver=fmtver530
  ELSE
    WRITE(6,'(/1x,a,i10,a/)')                        &
        'Data format, intver=',intver,', not found. The Job stopped.'
    CALL arpsstop('arpstop called from gradsdump.',1)
  END IF

  IF ( firstcall ) THEN

    DO l = 1, varnumax
      varparam(l) = '99'      ! Current version GrADS does not use.
                              ! It can be any integer
    END DO

    CALL getunit( nchanl )

    filen = LEN(filnam)
    CALL strlnth( filnam, filen )

    CALL asnctl ('NEWLOCAL', 1, ierr)
    CALL asnfile(filnam(1:filen), '-F f77 -N ieee', ierr)

    OPEN (UNIT=nchanl,FILE=trim(filnam(1:filen)),STATUS='new',          &
          FORM='unformatted',ACCESS='sequential',IOSTAT= istat )
!
!-----------------------------------------------------------------------
!
!  Write header info. This should be done only one time.
!
!-----------------------------------------------------------------------
!
    xbgn = (x(1) + x(2))/2.
    ybgn = (y(1) + y(2))/2.
    zbgn = (z(1) + z(2))/2.

    xinc = (x(2) - x(1))
    yinc = (y(2) - y(1))
    zinc = (z(2) - z(1))

    CALL xytoll(nx,ny,x,y,tem1(1,1,1),tem1(1,1,2))

    CALL xytoll(1,1,xbgn,ybgn,lat11,lon11)

    CALL a3dmax0lcl(tem1(1,1,1),1,nx,1,nx,1,ny,1,ny-1,1,1,1,1,          &
                 latmax,latmin)
    CALL a3dmax0lcl(tem1(1,1,2),1,nx,1,nx,1,ny,1,ny-1,1,1,1,1,          &
                 lonmax,lonmin)

    latinc = (latmax-latmin)/(ny-1)
    loninc = (lonmax-lonmin)/(nx-1)

    WRITE (6,'(a,f10.4,a,f10.4,a,f10.4)')                               &
             'latmin:latmax:latinc = ',                                 &
              latmin,':',latmax,':',latinc
    WRITE (6,'(a,f10.4,a,f10.4,a,f10.4)')                               &
             'lonmin:lonmax:loninc = ',                                 &
             lonmin,':',lonmax,':',loninc

    IF ( thisdmp <= 0.0 ) THEN
      ntm = 1
    ELSE
      ntm = nint((tstop-tstart)/thisdmp) + 1
    END IF

    IF (thisdmp < 60.) THEN
      WRITE (6, '(/a/a)')                                               &
          'GrADS reqiures the smallest uint minute for time interval.', &
          'Here we use uint MN to represent the second.'
      tinc = nint(thisdmp)
      dtunit = 'MN'
    ELSE IF (thisdmp < 3600.) THEN
      tinc = nint(thisdmp/60.)
      dtunit = 'MN'
    ELSE IF (thisdmp < 86400.) THEN
      tinc = nint(thisdmp/3600.)
      dtunit = 'HR'
    ELSE
      tinc = nint(thisdmp/86400.)
      dtunit = 'DY'
    END IF

    varnum = 0

    varnum = varnum + 1
    varnam(varnum) = 'zp      '
    vartit(varnum) = 'Physical vertical coordinates for the atmos model (m)'
    varlev(varnum) = nz

    IF ( varout == 1 ) THEN

      varnum = varnum + 1
      varnam(varnum) = 'u       '
      vartit(varnum) = 'X-velocity total wind (m/s)'
      varlev(varnum) = nz

      varnum = varnum + 1
      varnam(varnum) = 'uprt    '
      vartit(varnum) = 'X-velocity perturbation (m/s)'
      varlev(varnum) = nz

      varnum = varnum + 1
      varnam(varnum) = 'v       '
      vartit(varnum) = 'Y-velocity total wind (m/s)'
      varlev(varnum) = nz

      varnum = varnum + 1
      varnam(varnum) = 'vprt    '
      vartit(varnum) = 'Y-velocity perturbation (m/s)'
      varlev(varnum) = nz

      varnum = varnum + 1
      varnam(varnum) = 'w       '
      vartit(varnum) = 'Z-velocity total wind (m/s)'
      varlev(varnum) = nz

      varnum = varnum + 1
      varnam(varnum) = 'pt      '
      vartit(varnum) = 'Potential Temperature (K)'
      varlev(varnum) = nz

      varnum = varnum + 1
      varnam(varnum) = 'ptprt   '
      vartit(varnum) = 'Perturbation Potential Temperature (K)'
      varlev(varnum) = nz

      varnum = varnum + 1
      varnam(varnum) = 'p       '
      vartit(varnum) = 'Pressure (pascal)'
      varlev(varnum) = nz

      varnum = varnum + 1
      varnam(varnum) = 'pprt    '
      vartit(varnum) = 'Perturbation Pressure (pascal)'
      varlev(varnum) = nz

      varnum = varnum + 1
      varnam(varnum) = 'vort    '
      vartit(varnum) = 'Vertical vorticity (1/s)'
      varlev(varnum) = nz

      varnum = varnum + 1
      varnam(varnum) = 'div     '
      vartit(varnum) = 'Horizontal divergence (1/s)'
      varlev(varnum) = nz

    END IF

    IF ( mstout == 1 ) THEN

      varnum = varnum + 1
      varnam(varnum) = 'qv      '
      vartit(varnum) = 'Water Vapor Mixing Ratio (g/kg)'
      varlev(varnum) = nz

      varnum = varnum + 1
      varnam(varnum) = 'qvprt   '
      vartit(varnum) = 'Water Vapor Mixing Ratio '//                    &
                       'Perturbation (g/kg)'
      varlev(varnum) = nz

      DO nq = 1,nscalar
        varnum = varnum + 1
        varnam(varnum) = TRIM(qnames(nq))//'      '
        vartit(varnum) = TRIM(qdescp(nq))
        varlev(varnum) = nz
      END DO

      IF ( rainout == 1 ) THEN

        varnum = varnum + 1
        varnam(varnum) = 'raing   '
        vartit(varnum) = 'Grid Supersaturation Rain '
        varlev(varnum) = 0

        varnum = varnum + 1
        varnam(varnum) = 'rainc '
        vartit(varnum) = 'Cumulus Convection Rain '
        varlev(varnum) = 0

      END IF

      IF ( prcout == 1 ) THEN

        varnum = varnum + 1
        varnam(varnum) = 'prcrt1 '
        vartit(varnum) = 'Total precipitation rate '
        varlev(varnum) = 0

        varnum = varnum + 1
        varnam(varnum) = 'prcrt2 '
        vartit(varnum) = 'Grid scale precipitation rate '
        varlev(varnum) = 0

        varnum = varnum + 1
        varnam(varnum) = 'prcrt3 '
        vartit(varnum) = 'Cumulative precipitation rate '
        varlev(varnum) = 0

        varnum = varnum + 1
        varnam(varnum) = 'prcrt4 '
        vartit(varnum) = 'Microphysics precipitation rate '
        varlev(varnum) = 0

      END IF
    END IF

    IF ( tkeout == 1 ) THEN

      varnum = varnum + 1
      varnam(varnum) = 'tke     '
      vartit(varnum) ='Turbulent kinetic energy (m**2/s)'
      varlev(varnum) = nz

    END IF

    IF ( trbout == 1 ) THEN

      varnum = varnum + 1
      varnam(varnum) = 'kmh     '
      vartit(varnum) ='Horizontal Turb. Mixing Coefficient for '        &
                 //'Momentum (m**2/s)'
      varlev(varnum) = nz

      varnum = varnum + 1
      varnam(varnum) = 'kmv     '
      vartit(varnum) = 'Vertical Turb. Mixing Coefficient for '         &
                 //'Momentum (m**2/s)'
      varlev(varnum) = nz

    END IF

    IF ( sfcout == 1 ) THEN
!
! 06/18/2002  Zuwen He
!
! zpsoil is dumped when sfcout is on.
!
!
      varnum = varnum + 1
      varnam(varnum) = 'zpsoil  '
      vartit(varnum) = 'Physical vertical coordinates for the soil model (m)'
      varlev(varnum) = nzsoil
!
! 06/18/2002  Zuwen He
!
! In the previous code, we distinguished nstyp <= 1 with nstyp > 1,
! probably due to disk saving.
! Here we discard that, and dumps everything.
!
      DO is=0,nstyp
        WRITE (chrstr,'(i2.2)') is

        varnum = varnum + 1
        varnam(varnum) = 'tsoil'//chrstr(1:2)//' '
        vartit(varnum) = 'Soil Temperature (K) of soil type '//chrstr(1:2)
        if (is == 0) vartit(varnum) = 'Soil Temperature (K)'
        varlev(varnum) = nzsoil

        varnum = varnum + 1
        varnam(varnum) = 'qsoil'//chrstr(1:2)//' '
        vartit(varnum) = 'Soil moisture (m**3/m**3) of soil type '//chrstr(1:2)
        if (is == 0) vartit(varnum) = 'Soil moisture (m**3/m**3)'
        varlev(varnum) = nzsoil

        varnum = varnum + 1
        varnam(varnum) = 'wr'//chrstr(1:2)//' '
        vartit(varnum) = 'Canopy water amount (m) of soil type '//chrstr(1:2)
        if (is ==0)  vartit(varnum) = 'Canopy water amount (m)'
        varlev(varnum) = 0
      END DO

      IF ( snowout == 1 ) THEN
        varnum = varnum + 1
        varnam(varnum) = 'snowd   '
        vartit(varnum) = 'Snow depth (m)'
        varlev(varnum) = 0
      END IF

    END IF

    IF ( radout == 1 ) THEN

      varnum = varnum + 1
      varnam(varnum) = 'radfrc  '
      vartit(varnum) = 'Total radiation forcing (K/s)'
      varlev(varnum) = nz

      varnum = varnum + 1
      varnam(varnum) = 'radsw   '
      vartit(varnum) = 'Incoming solar rad. at surface (W/m**2)'
      varlev(varnum) = 0

      varnum = varnum + 1
      varnam(varnum) = 'rnflx   '
      vartit(varnum) = 'Surface radiation heat flux (W/m**2)'
      varlev(varnum) = 0

      varnum = varnum + 1
      varnam(varnum) = 'radswnet'
      vartit(varnum) = 'Net solar rad. (W/m**2)'
      varlev(varnum) = 0

      varnum = varnum + 1
      varnam(varnum) = 'radlwin '
      vartit(varnum) = 'Incoming longwave radiation (W/m**2)'
      varlev(varnum) = 0

    END IF

    IF ( flxout == 1 ) THEN

      varnum = varnum + 1
      varnam(varnum) = 'usflx   '
      vartit(varnum) = 'Surface u-momentum flux (kg/m*s**2)'
      varlev(varnum) = 0

      varnum = varnum + 1
      varnam(varnum) = 'vsflx   '
      vartit(varnum) = 'Surface v-momentum flux (kg/m*s**2)'
      varlev(varnum) = 0

      varnum = varnum + 1
      varnam(varnum) = 'ptsflx  '
      vartit(varnum) = 'Surface heat flux (K*kg/s*m**2)'
      varlev(varnum) = 0

      varnum = varnum + 1
      varnam(varnum) = 'qvsflx  '
      vartit(varnum) = 'Surface moisture flux (kg/s*m**2)'
      varlev(varnum) = 0

    END IF

! 06/18/2002 Zuwen He
!
! Writing fileheader
!
!
    hdbyte = 0

    WRITE(nchanl) fmtver
    hdbyte = hdbyte + 8 + 1*40

    WRITE(nchanl) runname
    hdbyte = hdbyte + 8 + 1*80

    WRITE(nchanl) nocmnt
    hdbyte = hdbyte + 8 + 1*4

    IF ( nocmnt > 0 ) THEN
      DO l = 1, nocmnt
        WRITE(nchanl) cmnt(l)
      END DO
      hdbyte = hdbyte + nocmnt * ( 8 + 80 )
    END IF

    WRITE(nchanl) nx,ny,nz,nzsoil,nstyp
    hdbyte = hdbyte + 8 + 5*4

    WRITE(nchanl) curtim,tmunit
    hdbyte = hdbyte + 8 + 4 + 10

    idummy = 0
    rdummy = 0.0

    WRITE(nchanl) varout, mstout, iceout, trbout, sfcout,               &
                 rainout, landout, idummy, idummy, totout,              &
                  tkeout, nscalar, mapproj, istgr, month,               &
                     day,   year,   hour, minute, second
    hdbyte = hdbyte + 8 + 20*4

    WRITE(nchanl) P_QC,  P_QR,  P_QI,  P_QS,  P_QG,  P_QH,              &
                  P_NC,  P_NR,  P_NI,  P_NS,  P_NG,  P_NH,              &
                  P_ZR,  P_ZI,  P_ZS,  P_ZG,  P_ZH,idummy,              &
                idummy,idummy,idummy,idummy,idummy,idummy,              &
                idummy,idummy,idummy,idummy,idummy,idummy
    hdbyte = hdbyte + 8 + 30*4

    WRITE(nchanl)  umove,   vmove, xgrdorg, ygrdorg,  trulat1,          &
                 trulat2,  trulon,  sclfct,  ntcloud, n0rain,           &
                  n0snow, n0grpl, n0hail, rhoice, rhosnow,              &
                   rhogrpl, rhohail,  alpharain, alphaice,alphasnow,    &
                   alphagrpl, alphahail, rdummy, rdummy,  rdummy,       &
                   tstop, thisdmp, latitud, ctrlat,  ctrlon
    hdbyte = hdbyte + 8 + 30*4

    IF ( totout /= 0 ) THEN
      WRITE(nchanl) hdmpopt, nstyp, prcout, radout, flxout,             &
                          0,snowout,idummy, idummy, idummy, & ! 0 for snowcvr
                     idummy, idummy, idummy, idummy, idummy,            &
                    idummy, idummy, idummy, idummy, idummy
      hdbyte = hdbyte + 8 + 20*4

      WRITE(nchanl) tstrtdmp,rdummy, rdummy, rdummy, rdummy,            &
                    rdummy, rdummy, rdummy, rdummy, rdummy,             &
                    rdummy, rdummy, rdummy, rdummy, rdummy,             &
                    rdummy, rdummy, rdummy, rdummy, rdummy
      hdbyte = hdbyte + 8 + 20*4
    END IF

    IF ( hdmpopt == 2 ) THEN
      WRITE (nchanl) numhdmp
      hdbyte = hdbyte + 8 + 4
      IF ( numhdmp > 0 ) THEN
        WRITE (nchanl) (hdmptim(i),i=1,numhdmp)
        hdbyte = hdbyte + 8 + numhdmp*4
      END IF
    END IF

    WRITE(nchanl) x
    hdbyte = hdbyte + 8 + nx*4

    WRITE(nchanl) y
    hdbyte = hdbyte + 8 + ny*4

    WRITE(nchanl) z
    hdbyte = hdbyte + 8 + nz*4

    WRITE(nchanl) varnum
    hdbyte = hdbyte + 8 + 1*4

    DO l = 1, varnum
      WRITE(nchanl) varnam(l),varlev(l),vartit(l)
    END DO
    hdbyte = hdbyte + varnum*( 8 + 8 + 4 + 60 )

    IF(landout == 1) THEN

      DO is=1,nstyp
        CALL iedgfill(soiltyp(1,1,is),1,nx,1,nx-1, 1,ny,1,ny-1,       &
                      1,1,1,1)
        WRITE (nchanl) ((soiltyp(i,j,is),i=1,nx),j=1,ny)
        hdbyte = hdbyte + ( 8 + nx*ny*4 )

        CALL edgfill(stypfrct(1,1,is),1,nx,1,nx-1, 1,ny,1,ny-1,       &
                     1,1,1,1)
        WRITE (nchanl) ((stypfrct(i,j,is),i=1,nx),j=1,ny)
        hdbyte = hdbyte + ( 8 + nx*ny*4 )
      END DO

      CALL iedgfill(vegtyp ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      WRITE (nchanl) ((vegtyp (i,j),i=1,nx),j=1,ny)
      hdbyte = hdbyte + ( 8 + nx*ny*4 )

      CALL edgfill(lai    ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      WRITE (nchanl) ((lai    (i,j),i=1,nx),j=1,ny)
      hdbyte = hdbyte + ( 8 + nx*ny*4 )

      CALL edgfill(roufns ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      WRITE (nchanl) ((roufns (i,j),i=1,nx),j=1,ny)
      hdbyte = hdbyte + ( 8 + nx*ny*4 )

      CALL edgfill(veg    ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      WRITE (nchanl) ((veg    (i,j),i=1,nx),j=1,ny)
      hdbyte = hdbyte + ( 8 + nx*ny*4 )

    END IF
!
!-----------------------------------------------------------------------
!
!  Open the GrADS data control file: runname.ctl
!
!-----------------------------------------------------------------------
!
    cntlen = lfnkey + 10
    gradscntl(1:cntlen) = runname(1:lfnkey)//'.gradscntl'
    CALL fnversn( gradscntl, cntlen )
    CALL getunit (nchout0)

    WRITE (6,'(a)') 'The GrADS data control file is '                   &
                    //gradscntl(1:cntlen)

    OPEN (nchout0, FILE = gradscntl(1:cntlen), STATUS = 'unknown')

    WRITE (nchout0,'(a,a)')                                             &
        'DSET    ',filnam(1:filen)

    WRITE (nchout0,'(a/a)')                                             &
        'TITLE   ARPS model 5.0 output for '//runname(1:lfnkey),'*'

    WRITE (nchout0,'(a)')                                               &
        'OPTIONS sequential cray_32bit_ieee big_endian'

    WRITE (nchout0,'(a,i10)')                                           &
        'FILEHEADER ',hdbyte

    WRITE (nchout0,'(a/a)')                                             &
        'UNDEF   -9.e+33','*'

    IF ( mapproj == 2 ) THEN

!      WRITE (nchout0,'(a)')                                             &
!          '* For lat-lon-lev display, uncomment the following 4 lines.'

      WRITE (nchout0,'(a,1x,i8,1x,i3,a,2f12.6,2i3,3f12.6,2f12.2)')      &
          'PDEF',nx,ny,' LCC',lat11,lon11,1,1,                          &
              trulat1,trulat2,trulon,xinc,yinc

      WRITE (nchout0,'(a,1x,i8,a,f10.4,1x,f10.4)')                      &
          'XDEF',nx,'  LINEAR  ',lonmin,loninc

      WRITE (nchout0,'(a,1x,i8,a,f10.4,1x,f10.4)')                      &
          'YDEF',ny,'  LINEAR  ',latmin,latinc

    ELSE

!
! 06/18/2002 Zuwen He
!
! ZDEF is max(nz,nzsoil), because we are sharing the soil model
! and the atmosphere model in the same dataset.
!
      WRITE (nchout0,'(a)')                                             &
          '* For i-j-k display, uncomment the following 3 lines.'

      WRITE (nchout0,'(a,1x,i8,a,2i10)')                                &
          '* XDEF',nx,'  LINEAR  ',1,1

      WRITE (nchout0,'(a,1x,i8,a,2i10)')                                &
          '* YDEF',ny,'  LINEAR  ',1,1

      WRITE (nchout0,'(a,1x,i8,a,2i10/a)')                              &
          '* ZDEF',max(nz,nzsoil),'  LINEAR  ',1,1,'*'

      WRITE (nchout0,'(a)')                                             &
          '* For x-y-z display, uncomment the following 3 lines.'

      WRITE (nchout0,'(a,1x,i8,a,2f15.4)')                              &
          'XDEF',nx,'  LINEAR  ',xbgn/1000.,xinc/1000.

      WRITE (nchout0,'(a,1x,i8,a,2f15.4)')                              &
          'YDEF',ny,'  LINEAR  ',ybgn/1000.,yinc/1000.

    END IF
!
! 06/18/2002, Zuwen He
!
! adding comments for ZDEF
!
      WRITE (nchout0,'(18(a/),a)')                                      &
          '*          ',                                                &
          '* WARNING: ',                                                &
          '* The associated GrADS dataset includes both the ',          &
          '* atmosphere and soil variables. ',                          &
          '* ZDEF is the larger one between the level numbers ',        &
          '* for the atmosphere and soil models. ',                     &
          '* The actual number of vertical levels for a variable ',     &
          '* can be obtained from the definition of the variable ',     &
          '* located between VARS and ENDVARS.                    ',    &
          '* Vertically, the soil levels increase downward under ',     &
          '* the ground, while atmosphere levels upward above ',        &
          '* the ground. ',                                             &
          '* Since we make the grid locations of the atmospheric ',     &
          '* model as the primary vertical coordinate, when ',          &
          '* plotting the soil variables, the vertical axis is ',       &
          '* actually the atmospheric axis. ',                          &
          '* Thus plotting soil variables is for debug purpose. ',      &
          '*          '

    WRITE (nchout0,'(a,1x,i8,a)')                                     &
          'ZDEF',max(nz,nzsoil),'  LEVELS '
    IF ( ternopt == 0 ) THEN
      WRITE (nchout0,'(8f10.2)')                                        &
          ((zp(2,2,k)+zp(2,2,k+1))/2.,k=1,nz-1),zp(2,2,nz)
    ELSE
      !
      ! Use AGL height as vertical levels
      !
      WRITE (nchout0,'(8f10.2)')                                        &
          ((zp(2,2,k)+zp(2,2,k+1))/2.-zp(2,2,2),k=1,nz-1),zp(2,2,nz)
!      WRITE (nchout0,'(8f10.2)')                                        &
!          ((z(k)+z(k+1))/2.,k=1,nz-1),z(nz)
!      WRITE (nchout0,'(a/a/a)')                                         &
!      '* WARNING  The vertical levels were set to computational ',      &
!      '*          coordinates, z(k), because zp is not uniform in ',    &
!      '*          horizontal when terrain option was turned on'
    END IF

    IF ( initopt /= 2 ) THEN
      WRITE (chrstr,'(i2.2,a,i2.2,a,i2.2,a3,i4.4)')                     &
            hour,':',minute,'Z',day,monnam(month),year
    ELSE
      second1 = MOD( second + nint(tstart), 60 )
      minute1 = ( second + nint(tstart) ) / 60
      minute1 = MOD( minute + minute1, 60 )
      hour1   = ( minute + ( second + nint(tstart) ) / 60 ) /60
      hour1   = MOD( hour + hour1, 24 )
      day1    = ( hour + ( minute                                       &
              + ( second + nint(tstart) ) / 60 ) /60 ) / 24
      jday1   = jday + day1

      loopdy  = 0
      IF ( MOD( year, 4 ) == 0 ) loopdy = 1
      year1 = year + jday1 / ( 365 + loopdy )
      jday1 = MOD( jday1, 365 + loopdy )

      month1 = 1

      DO m = 2, 11
        IF (jday1>mndys(m) .AND. jday1<=mndys(m+1)+loopdy) month1 = m
      END DO
      day1 = jday1 - mndys(month1)

      WRITE (chrstr,'(i2.2,a,i2.2,a,i2.2,a3,i4.4)')                     &
            hour1,':',minute1,'Z',day1,monnam(month1),year1

    END IF

    IF ( hdmpopt == 1 ) THEN
      WRITE (nchout0,'(a/a,f10.2,a/a)')                                 &
          '* WARNING: The time interval is applied after the time ',    &
          '*          level at ',tstrtdmp,' model seconds.',            &
          '*          The first time level was the model initial data.'
    ELSE IF ( hdmpopt == 2 ) THEN
      WRITE (nchout0,'(a/a/a/a)')                                       &
      '* WARNING& ! The actual time levels in the data file may not ',  &
      '*          necessarily be in a constant increment. It was ',     &
      '*          written for the model times listed in the ',          &
      '*          following:'
      DO i=1,numhdmp
        WRITE (nchout0,'(a,f10.2)')                                     &
            '*          ',hdmptim(i)
      END DO
    END IF

    WRITE (nchout0,'(a,1x,i8,a,a,3x,i2.2,a/a)')                         &
        'TDEF',ntm,'  LINEAR  ',chrstr,tinc,dtunit,'*'

    WRITE (nchout0,'(a,1x,i3)')                                         &
        'VARS',varnum

    DO l = 1, varnum

      WRITE (nchout0,'(a8,1x,i3,1x,a10,2x,a)')                          &
             varnam(l),varlev(l),varparam(l),vartit(l)

    END DO

    WRITE (nchout0,'(a)')                                               &
        'ENDVARS'

    CLOSE (nchout0)
    CALL retunit(nchout0)

    firstcall = .false.

  END IF
!
!-----------------------------------------------------------------------
!
!  Write data into the GrADS data file
!
!-----------------------------------------------------------------------
!
  CALL edgfill(zp  ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)

  IF ( istgr == 0 ) THEN
    DO j=1,ny
      DO i=1,nx
        DO k=1,nz-1
          tem1(i,j,k) = .5 * ( zp(i,j,k) + zp(i,j,k+1) )
        END DO
        tem1(i,j,nz) = zp(i,j,nz)
      END DO
    END DO
  ELSE
    DO j=1,ny
      DO i=1,nx
        DO k=1,nz
          tem1(i,j,k) = zp(i,j,k)
        END DO
      END DO
    END DO
  END IF

  DO k=1, nz
    WRITE(nchanl) ((tem1(i,j,k),i=1,nx),j=1,ny)
  ENDDO

  IF ( varout == 1 ) THEN
!
!-----------------------------------------------------------------------
!
!  If varout = 1, Write out u, uprt, v, vprt, w, wprt,
!  pt, ptprt, p, pprt, vort, and div.
!
!-----------------------------------------------------------------------
!
    CALL edgfill(u,   1,nx,1,nx, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL edgfill(ubar,1,nx,1,nx, 1,ny,1,ny-1, 1,nz,1,nz-1)

    IF ( istgr == 0 ) THEN
      DO k=1,nz
        DO j=1,ny
          DO i=1,nx-1
            tem1(i,j,k) = .5 * ( u(i,j,k) + u(i+1,j,k) )
            tem2(i,j,k) = tem1(i,j,k)                                   &
                        - .5 * ( ubar(i,j,k) + ubar(i+1,j,k) )
          END DO
          tem1(nx,j,k) = u(nx,j,k)
          tem2(nx,j,k) = u(nx,j,k) - ubar(nx,j,k)
        END DO
      END DO
    ELSE
      DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            tem1(i,j,k) = u(i,j,k)
            tem2(i,j,k) = u(i,j,k) - ubar(i,j,k)
          END DO
        END DO
      END DO
    END IF

    DO k=1, nz
      WRITE(nchanl) ((tem1(i,j,k),i=1,nx),j=1,ny)
    END DO
    DO k=1, nz
      WRITE(nchanl) ((tem2(i,j,k),i=1,nx),j=1,ny)
    END DO

    CALL edgfill(v,   1,nx,1,nx-1, 1,ny,1,ny, 1,nz,1,nz-1)
    CALL edgfill(vbar,1,nx,1,nx-1, 1,ny,1,ny, 1,nz,1,nz-1)

    IF ( istgr == 0 ) THEN
      DO k=1,nz
        DO i=1,nx
          DO j=1,ny-1
            tem1(i,j,k) = .5 * ( v(i,j,k) + v(i,j+1,k) )
            tem2(i,j,k) = tem1(i,j,k)                                   &
                        - .5 * ( vbar(i,j,k) + vbar(i,j+1,k) )
          END DO
          tem1(i,ny,k) = v(i,ny,k)
          tem2(i,ny,k) = v(i,ny,k) - vbar(i,ny,k)
        END DO
      END DO
    ELSE
      DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            tem1(i,j,k) = v(i,j,k)
            tem2(i,j,k) = v(i,j,k) - vbar(i,j,k)
          END DO
        END DO
      END DO
    END IF

    DO k=1, nz
      WRITE(nchanl) ((tem1(i,j,k),i=1,nx),j=1,ny)
    ENDDO

    DO k=1, nz
      WRITE(nchanl) ((tem2(i,j,k),i=1,nx),j=1,ny)
    ENDDO

    CALL edgfill(w   ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)

    IF ( istgr == 0 ) THEN
      DO j=1,ny
        DO i=1,nx
          DO k=1,nz-1
            tem1(i,j,k) = .5 * ( w(i,j,k) + w(i,j,k+1) )
          END DO
          tem1(i,j,nz) = w(i,j,nz)
        END DO
      END DO
    ELSE
      DO j=1,ny
        DO i=1,nx
          DO k=1,nz
            tem1(i,j,k) = w(i,j,k)
          END DO
        END DO
      END DO
    END IF

    DO k=1, nz
      WRITE(nchanl) ((tem1(i,j,k),i=1,nx),j=1,ny)
    ENDDO

    CALL edgfill(ptprt,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL edgfill(ptbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          tem1(i,j,k) = ptprt(i,j,k) + ptbar(i,j,k)
        END DO
      END DO
    END DO

    DO k=1, nz
      WRITE(nchanl) ((tem1(i,j,k),i=1,nx),j=1,ny)
    ENDDO
    DO k=1, nz
      WRITE(nchanl) ((ptprt(i,j,k),i=1,nx),j=1,ny)
    ENDDO

    CALL edgfill(pprt,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL edgfill(pbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          tem1(i,j,k) = pprt(i,j,k) + pbar(i,j,k)
        END DO
      END DO
    END DO

    DO k=1, nz
      WRITE(nchanl) ((tem1(i,j,k),i=1,nx),j=1,ny)
    END DO
    DO k=1, nz
      WRITE(nchanl) ((pprt(i,j,k),i=1,nx),j=1,ny)
    END DO

!-----------------------------------------------------------------------
!
!  Vorticity:
!
!-----------------------------------------------------------------------

    DO k=2,nz-2
      DO j=2,ny-2
        DO i=2,nx-2
          tem1(i,j,k)=                                                  &
              (v(i+1,j,k)-v(i-1,j,k)+v(i+1,j+1,k)-v(i-1,j+1,k))/        &
              (4*(x(i+1)-x(i)))-                                        &
              (u(i,j+1,k)-u(i,j-1,k)+u(i+1,j+1,k)-u(i+1,j-1,k))/        &
              (4*(y(j+1)-y(j)))
        END DO
      END DO
    END DO

    DO j=2,ny-2
      DO i=2,nx-2
        tem1(i,j,   1)=tem1(i,j,   2)
        tem1(i,j,nz-1)=tem1(i,j,nz-2)
      END DO
    END DO

    DO k=1,nz-1
      DO j=2,ny-2
        tem1(   1,j,k)=tem1(   2,j,k)
        tem1(nx-1,j,k)=tem1(nx-2,j,k)
      END DO
    END DO

    DO k=1,nz-1
      DO i=1,nx-1
        tem1(i,   1,k)=tem1(i,   2,k)
        tem1(i,ny-1,k)=tem1(i,ny-2,k)
      END DO
    END DO

    CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    DO k=1, nz
      WRITE(nchanl) ((tem1(i,j,k),i=1,nx),j=1,ny)
    END DO

!-----------------------------------------------------------------------
!
!  Divergernce:
!
!-----------------------------------------------------------------------

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem1(i,j,k) = (u(i+1,j,k)-u(i,j,k))/(x(i+1)-x(i))             &
                      + (v(i,j+1,k)-v(i,j,k))/(y(j+1)-y(j))
        END DO
      END DO
    END DO

    CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    DO k=1, nz
      WRITE(nchanl) ((tem1(i,j,k),i=1,nx),j=1,ny)
    END DO

  END IF       ! End varout

  IF ( mstout == 1 ) THEN
!
!-----------------------------------------------------------------------
!
!  Write out moist variables qv, qvprt, qc, qr, qi, qs, and qh
!
!-----------------------------------------------------------------------
!
    CALL edgfill(qv,   1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL edgfill(qvbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          tem1(i,j,k) = qv(i,j,k) - qvbar(i,j,k)
        END DO
      END DO
    END DO

    DO k=1, nz
      WRITE(nchanl) ((qv(i,j,k),i=1,nx),j=1,ny)
    END DO
    DO k=1, nz
      WRITE(nchanl) ((tem1(i,j,k),i=1,nx),j=1,ny)
    END DO

    DO nq = 1,nscalar

      CALL edgfill(qscalar(:,:,:,nq),1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      DO k=1, nz
        WRITE(nchanl) ((qscalar(i,j,k,nq),i=1,nx),j=1,ny)
      END DO

    END DO

    IF ( rainout == 1 ) THEN

      CALL edgfill(raing,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      WRITE(nchanl) ((raing(i,j),i=1,nx),j=1,ny)

      CALL edgfill(rainc,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      WRITE(nchanl) ((rainc(i,j),i=1,nx),j=1,ny)

    END IF     ! End rainout

    IF ( prcout == 1 ) THEN

      CALL edgfill(prcrate,1,nx,1,nx-1, 1,ny,1,ny-1, 1,4,1,4)
      WRITE(nchanl) ((prcrate(i,j,1),i=1,nx),j=1,ny)
      WRITE(nchanl) ((prcrate(i,j,2),i=1,nx),j=1,ny)
      WRITE(nchanl) ((prcrate(i,j,3),i=1,nx),j=1,ny)
      WRITE(nchanl) ((prcrate(i,j,4),i=1,nx),j=1,ny)

    END IF     ! End prcout

  END IF       ! End mstout

  IF ( tkeout == 1 ) THEN

    CALL edgfill(tke,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)

    DO k=1, nz
      WRITE(nchanl) ((tke(i,j,k),i=1,nx),j=1,ny)
    END DO

  END IF

  IF ( trbout == 1 ) THEN
!
!-----------------------------------------------------------------------
!
!  If trbout = 1, write out the turbulence parameter, km.
!
!-----------------------------------------------------------------------
!
    CALL edgfill(kmh,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)

    DO k=1, nz
      WRITE(nchanl) ((kmh(i,j,k),i=1,nx),j=1,ny)
    END DO

    CALL edgfill(kmv,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)

    DO k=1, nz
      WRITE(nchanl) ((kmv(i,j,k),i=1,nx),j=1,ny)
    END DO


  END IF       ! trbout

  IF ( sfcout == 1 ) THEN
!
!-----------------------------------------------------------------------
!
!  Write out surface variables tsoil, qsoil, wetcanp
!
!-----------------------------------------------------------------------
!
    CALL edgfill(zpsoil,1,nx,1,nx-1,1,ny,1,ny-1,        &
                 1,nzsoil,1,nzsoil)
    DO k=1, nzsoil
      WRITE(nchanl) ((zpsoil(i,j,k),i=1,nx),j=1,ny)
    END DO

    DO is=0,nstyp
      CALL edgfill(tsoil  (1,1,1,is),1,nx,1,nx-1,1,ny,1,ny-1,        &
                   1,nzsoil,1,nzsoil)
      DO k=1, nzsoil
        WRITE(nchanl) ((tsoil(i,j,k,is),i=1,nx),j=1,ny)
      END DO

      CALL edgfill(qsoil(1,1,1,is),1,nx,1,nx-1,1,ny,1,ny-1,        &
                   1,nzsoil,1,nzsoil)
      DO k=1, nzsoil
        WRITE(nchanl) ((qsoil(i,j,k,is),i=1,nx),j=1,ny)
      END DO

      CALL edgfill(wetcanp(1,1,is),1,nx,1,nx-1,1,ny,1,ny-1,          &
                   1,1,1,1)
      WRITE(nchanl) ((wetcanp(i,j,is),i=1,nx),j=1,ny)
    END DO

    IF (snowout == 1) THEN
      CALL edgfill(snowdpth(1,1),1,nx,1,nx-1,1,ny,1,ny-1,               &
                    1,1,1,1)
      WRITE(nchanl) ((snowdpth(i,j),i=1,nx),j=1,ny)
    END IF

  END IF       ! End sfcout
!
!-----------------------------------------------------------------------
!
!  If radout = 1, write out the radiation arrays
!
!-----------------------------------------------------------------------
!
  IF ( radout == 1 ) THEN

    CALL edgfill(radfrc,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    DO k=1, nz
      WRITE(nchanl) ((radfrc(i,j,k),i=1,nx),j=1,ny)
    END DO

    CALL edgfill(radsw,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    WRITE(nchanl) ((radsw(i,j),i=1,nx),j=1,ny)

    CALL edgfill(rnflx,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    WRITE(nchanl) ((rnflx(i,j),i=1,nx),j=1,ny)

    CALL edgfill(radswnet,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    WRITE(nchanl) ((radswnet(i,j),i=1,nx),j=1,ny)

    CALL edgfill(radlwin,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    WRITE(nchanl) ((radlwin(i,j),i=1,nx),j=1,ny)


  END IF       ! radout
!
!-----------------------------------------------------------------------
!
!  If flxout = 1, write out the surface fluxes
!
!-----------------------------------------------------------------------
!
  IF ( flxout == 1 ) THEN

    CALL edgfill(usflx,1,nx,1,nx, 1,ny,1,ny-1, 1,1,1,1)
    WRITE(nchanl) ((usflx(i,j),i=1,nx),j=1,ny)

    CALL edgfill(vsflx,1,nx,1,nx-1, 1,ny,1,ny, 1,1,1,1)
    WRITE(nchanl) ((vsflx(i,j),i=1,nx),j=1,ny)

    CALL edgfill(ptsflx,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    WRITE(nchanl) ((ptsflx(i,j),i=1,nx),j=1,ny)

    CALL edgfill(qvsflx,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    WRITE(nchanl) ((qvsflx(i,j),i=1,nx),j=1,ny)

  END IF       ! flxout

  RETURN
END SUBROUTINE gradsdump
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GRADSJOINDUMP              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE gradsjoindump(nx,ny,nz,nzsoil,nstyps,nchanl,filnam,istgr,    &
           u,v,w,ptprt,pprt,qv,qscalar,tke,kmh,kmv,                     &
           ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,                      &
           x,y,z,zp,zpsoil,                                             &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,                                &
           raing,rainc,prcrate,                                         &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           usflx,vsflx,ptsflx,qvsflx,                                   &
           tem1,tem2)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write combined history data into the file "filnam" in the GrADS format
!  for MPI jobs.
!
!  All output data are located at the grid cell centers.
!
!       X_center(i) = ( X_edge(i) + X_edge(i+1) )/2,  i = 1, n-1, 1
!
!  The last edge value were kept unchanged so that it can be used to
!  restore the stagger grid point values.
!
!       X_edge(n) = X_center(n)
!       X_edge(i) = 2*X_center(i) - X_edge(i+1),      i = n-1, 1, -1
!
!  Therefore, when you display the data by GrADS, set the x, y, and
!  z dimension not to exceed the Max-1.
!
!  Since total and perturbation variables will be dumped, the base
!  state variables can be obtained by Xbar = X - Xprt.
!
!  The smallest time unit of the current version GrADS is minutes.
!  Therefore, the time intervals for history data dumping, thisdmp,
!  should be multiples of 60 seconds.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  8/20/2002.
!  Based on subroutine gradsdump
!
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of grid points in the vertical
!
!    nchanl   FORTRAN I/O channel number for history data output.
!
!    u        x component of velocity at a given time level (m/s)
!    v        y component of velocity at a given time level (m/s)
!    w        Vertical component of Cartesian velocity at a given
!             time level (m/s)
!    ptprt    Perturbation potential temperature at a given time
!             level (K)
!    pprt     Perturbation pressure at a given time level (Pascal)
!    qv       Water vapor specific humidity at a given time level
!             (kg/kg)
!    qc       Cloud water mixing ratio at a given time level (kg/kg)
!    qr       Rainwater mixing ratio at a given time level (kg/kg)
!    qi       Cloud ice mixing ratio at a given time level (kg/kg)
!    qs       Snow mixing ratio at a given time level (kg/kg)
!    qh       Hail mixing ratio at a given time level (kg/kg)
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!
!    ubar     Base state zonal velocity component (m/s)
!    vbar     Base state meridional velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhobar   Base state density (kg/m**3)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space(m)
!    zpsoil   Vertical coordinate of grid points in physical space(m)
!
!    soiltyp  Soil type
!    vegtyp   Vegetation type
!    lai      Leaf Area Index
!    roufns   Surface roughness
!    veg      Vegetation fraction
!
!    tsoil    Soil temperature (K)
!    qsoil    Soil moisture (m**3/m**3)
!    wetcanp  Canopy water amount (m)
!
!    raing    Grid supersaturation rain
!    rainc    Cumulus convective rain
!    prcrate  Precipitation rates
!
!    radfrc   Radiation forcing (K/s)
!    radsw    Solar radiation reaching the surface
!    rnflx    Net radiation flux absorbed by surface
!
!    usflx    Surface flux of u-momentum (kg/(m*s**2))
!    vsflx    Surface flux of v-momentum (kg/(m*s**2))
!    ptsflx   Surface heat flux (K*kg/(m**2 * s ))
!    qvsflx   Surface moisture flux of (kg/(m**2 * s))
!
!  OUTPUT:
!
!    None.
!
!  WORK ARRAY:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    globx    Temporary work array.
!    globy    Temporary work array.
!    out3d    Temporary work array.
!    outtsoil Temporary work array.
!    outqsoil Temporary work array.
!
!
!-----------------------------------------------------------------------
!
!  The following parameters are passed into this subroutine through
!  a common block in globcst.inc, and they determine which
!  variables are output.
!
!  varout =0 or 1. If varout=0, model perturbation variables are not
!                  dumped.
!  mstout =0 or 1. If mstout =0, water variables are not dumped.
!  rainout=0 or 1. If rainout=0, rain variables are not dumped.
!  prcout =0 or 1. If prcout =0, precipitation rates are not dumped.
!  iceout =0 or 1. If iceout =0, qi, qs and qh are not dumped.
!  tkeout =0 or 1. If tkeout =0, tke is not dumped.
!  trbout =0 or 1. If trbout =0, kmh and kmv are not dumped
!  sfcout =0 or 1. If sfcout =0, surface variables are not dumped.
!  landout=0 or 1. If landout=0, surface property arrays are not dumped.
!  radout =0 or 1. If radout =0, radiation arrays are not dumped.
!  flxout =0 or 1. If flxout =0, surface fluxes are not dumped.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of grid points in the soil

  INTEGER :: istgr             ! Flag for dumping stager point data

  INTEGER :: nchanl            ! FORTRAN I/O channel number for output

  REAL :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz)     ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz)     ! Total w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz)     ! Water vapor specific humidity (kg/kg)

  REAL :: qscalar    (nx,ny,nz,nscalar)

  REAL :: tke   (nx,ny,nz)     ! Turbulent Kinetic Energy ((m/s)**2)
  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: wbar  (nx,ny,nz)     ! Base state w-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific
                               ! humidity(kg/kg)

  REAL :: x     (nx)           ! x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! Physical height coordinate defined at
                               ! w-point of the staggered grid.
  REAL :: zpsoil(nx,ny,nzsoil) ! Physical height coordinate defined at
                               ! w-point of the soil.

  INTEGER :: nstyps
  INTEGER :: soiltyp(nx,ny,nstyps)    ! Soil type
  REAL :: stypfrct(nx,ny,nstyps)   ! Fraction of soil types
  INTEGER :: vegtyp (nx,ny)           ! Vegetation type
  REAL :: lai    (nx,ny)           ! Leaf Area Index
  REAL :: roufns (nx,ny)           ! Surface roughness
  REAL :: veg    (nx,ny)           ! Vegetation fraction

  REAL :: tsoil  (nx,ny,nzsoil,0:nstyps) ! Soil temperature (K)
  REAL :: qsoil  (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny,0:nstyps)     ! Canopy water amount (m)
  REAL :: snowdpth(nx,ny)             ! Snow depth (m)

  REAL :: raing(nx,ny)                ! Grid supersaturation rain
  REAL :: rainc(nx,ny)                ! Cumulus convective rain
  REAL :: prcrate(nx,ny,4)     ! precipitation rate (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precip. rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulus precip. rate
                               ! prcrate(1,1,4) = microphysics precip. rate

  REAL :: radfrc(nx,ny,nz)     ! Radiation forcing (K/s)
  REAL :: radsw (nx,ny)        ! Solar radiation reaching the surface
  REAL :: rnflx (nx,ny)        ! Net radiation flux absorbed by surface
  REAL :: radswnet(nx,ny)      ! Net shortwave radiation
  REAL :: radlwin(nx,ny)       ! Incoming longwave radiation


  REAL :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  REAL :: tem1  (nx,ny,nz)            ! Temporary work array
  REAL :: tem2  (nx,ny,nz)            ! Temporary work array

  REAL, ALLOCATABLE :: gtem1  (:,:,:) ! Temporary work array
  REAL, ALLOCATABLE :: gtem2  (:,:,:) ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Parameters describing this routine
!
!-----------------------------------------------------------------------
!
! fmtver??: to label each data a version.
! intver??: an integer to allow faster comparison than fmtver??,
!           which are strings.
!
! Verion 5.00: significant change in soil variables since version 4.10.
!
  CHARACTER (LEN=40) :: fmtver,fmtver500,fmtver530
  INTEGER  :: intver,intver500,intver530
  PARAMETER (fmtver500='005.00 GrADS Binary Data',intver500=500)
  PARAMETER (fmtver530='005.30 GrADS Binary Data',intver530=530)

  CHARACTER (LEN=10) :: tmunit
  PARAMETER (tmunit='seconds   ')
  CHARACTER (LEN=2) :: dtunit
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: varnumax
  PARAMETER ( varnumax = 100 )
  CHARACTER (LEN=8) :: varnam(varnumax)
  CHARACTER (LEN=60) :: vartit(varnumax)
  CHARACTER (LEN=10) :: varparam(varnumax)
  INTEGER :: varlev(varnumax)
  INTEGER :: i,j,k,l,m, istat, ierr,is
  INTEGER :: nq
  INTEGER :: varnum
  INTEGER :: tinc,ntm
  INTEGER :: idummy
  REAL :: rdummy
  REAL :: xbgn,ybgn,zbgn, xinc,yinc,zinc
  REAL :: lat11, lon11
  REAL :: latmin, latmax, lonmin, lonmax, latinc, loninc
  CHARACTER (LEN=*) :: filnam
  CHARACTER (LEN=80) :: gradscntl
  INTEGER :: filen, cntlen
  CHARACTER (LEN=20) :: chrstr
  INTEGER :: hdbyte
  INTEGER :: nchout0

  INTEGER :: year1, month1, day1, jday1, loopdy
  INTEGER :: hour1, minute1, second1

  INTEGER :: nxlg, nylg, nzlg, nzsoillg, n3d
  INTEGER, ALLOCATABLE :: out3di(:,:,:)
  REAL, ALLOCATABLE :: globx(:), globy(:), out3d(:,:,:), out3dt(:,:,:)
  REAL, AlLOCATABLE :: outtsoil(:,:,:,:), outqsoil(:,:,:,:)

  CHARACTER (LEN=3) :: monnam(12)            ! Name of months
  DATA monnam/'jan', 'feb', 'mar', 'apr', 'may', 'jun',                 &
              'jul', 'aug', 'sep', 'oct', 'nov', 'dec'/

  INTEGER :: mndys(12)                 ! days for each months
  DATA mndys/0,31,59,90,120,151,181,212,243,273,304,334/

  LOGICAL :: firstcall
  DATA firstcall/.true./
!
!-----------------------------------------------------------------------
!
  SAVE firstcall
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  nxlg = nproc_x*(nx-3)+3
  nylg = nproc_y*(ny-3)+3
  nzlg = nz
  nzsoillg = nzsoil
  n3d = MAX(nzlg, nzsoillg, nstyps+1, 4)   ! 3rd dimenson for out3d

  intver = intver530  !  for the time being, in the future, we will
                      !  allow to dump data in the different version
                      !  intver will be assigned from input file
  IF (intver == intver500) THEN
    fmtver=fmtver500
  ELSE IF (intver == intver530) THEN
    fmtver=fmtver530
  ELSE
    WRITE(6,'(/1x,a,i10,a/)')                        &
        'Data format, intver=',intver,', not found. The Job stopped.'
    CALL arpsstop('arpstop called from gradsjoindump.',1)
  END IF

  IF (firstcall) THEN
    IF (myproc == 0) THEN
      ALLOCATE (globx( nxlg ),stat=istat)
      CALL check_alloc_status(istat, "gradsjoindump:globx")

      ALLOCATE (globy( nylg ),stat=istat)
      CALL check_alloc_status(istat, "gradsjoindump:globy")
    END IF

    CALL mpimerge1dx(x,nx,globx)
    CALL mpimerge1dy(y,ny,globy)
  END IF

  IF (myproc == 0) THEN

    ALLOCATE (out3d( nxlg,nylg, n3d ),stat=istat)
    CALL check_alloc_status(istat, "gradsjoindump:out3d")

    ALLOCATE (out3dt( nxlg,nylg, n3d ),stat=istat)
    CALL check_alloc_status(istat, "gradsjoindump:out3dt")

    ALLOCATE (out3di( nxlg,nylg, nstyps ),stat=istat)
    CALL check_alloc_status(istat, "gradsjoindump:out3di")

    ALLOCATE (gtem1( nxlg,nylg, n3d ),stat=istat)
    CALL check_alloc_status(istat, "gradsjoindump:gtem1")

    ALLOCATE (gtem2( nxlg,nylg, n3d ),stat=istat)
    CALL check_alloc_status(istat, "gradsjoindump:gtem2")

    ALLOCATE (outtsoil( nxlg,nylg, nzsoillg, 0:nstyps ),stat=istat)
    CALL check_alloc_status(istat, "gradsjoindump:outtsoil")

    ALLOCATE (outqsoil( nxlg,nylg, nzsoillg, 0:nstyps ),stat=istat)
    CALL check_alloc_status(istat, "gradsjoindump:outqsoil")

    IF ( firstcall ) THEN

      DO l = 1, varnumax
        varparam(l) = '99'      ! Current version GrADS does not use.
                                ! It can be any integer
      END DO

      CALL getunit( nchanl )

      filen = LEN(filnam)
      CALL strlnth( filnam, filen )

      CALL asnctl ('NEWLOCAL', 1, ierr)
      CALL asnfile(filnam(1:filen), '-F f77 -N ieee', ierr)

      OPEN (UNIT=nchanl,FILE=trim(filnam(1:filen)),STATUS='new',          &
            FORM='unformatted',ACCESS='sequential',IOSTAT= istat )
!
!-----------------------------------------------------------------------
!
!  Write header info. This should be done only one time.
!
!-----------------------------------------------------------------------
!
      xbgn = (x(1) + x(2))/2.
      ybgn = (y(1) + y(2))/2.
      zbgn = (z(1) + z(2))/2.

      xinc = (x(2) - x(1))
      yinc = (y(2) - y(1))
      zinc = (z(2) - z(1))

      CALL xytoll(nxlg,nylg,globx,globy,out3d(1,1,1),out3d(1,1,2))

      CALL xytoll(1,1,xbgn,ybgn,lat11,lon11)

      istat = mp_opt
      mp_opt = 0

      CALL a3dmax0(out3d(1,1,1),1,nxlg,1,nxlg,1,nylg,1,nylg-1,1,1,1,1, &
                   latmax,latmin)
      CALL a3dmax0(out3d(1,1,2),1,nxlg,1,nxlg,1,nylg,1,nylg-1,1,1,1,1, &
                   lonmax,lonmin)
      mp_opt = istat

      latinc = (latmax-latmin)/(nylg-1)
      loninc = (lonmax-lonmin)/(nxlg-1)

      WRITE (6,'(a,f10.4,a,f10.4,a,f10.4)')                            &
               'latmin:latmax:latinc = ',                              &
                latmin,':',latmax,':',latinc
      WRITE (6,'(a,f10.4,a,f10.4,a,f10.4)')                            &
               'lonmin:lonmax:loninc = ',                              &
               lonmin,':',lonmax,':',loninc

      IF ( thisdmp <= 0.0 ) THEN
        ntm = 1
      ELSE
        ntm = nint((tstop-tstart)/thisdmp) + 1
      END IF

      IF (thisdmp < 60.) THEN
        WRITE (6, '(/a/a)')                                               &
            'GrADS reqiures the smallest uint minute for time interval.', &
            'Here we use uint MN to represent the second.'
        tinc = nint(thisdmp)
        dtunit = 'MN'
      ELSE IF (thisdmp < 3600.) THEN
        tinc = nint(thisdmp/60.)
        dtunit = 'MN'
      ELSE IF (thisdmp < 86400.) THEN
        tinc = nint(thisdmp/3600.)
        dtunit = 'HR'
      ELSE
        tinc = nint(thisdmp/86400.)
        dtunit = 'DY'
      END IF

      varnum = 0

      varnum = varnum + 1
      varnam(varnum) = 'zp      '
      vartit(varnum) = 'Physical vertical coordinates for the atmos model (m)'
      varlev(varnum) = nzlg

      IF ( varout == 1 ) THEN

        varnum = varnum + 1
        varnam(varnum) = 'u       '
        vartit(varnum) = 'X-velocity total wind (m/s)'
        varlev(varnum) = nzlg

        varnum = varnum + 1
        varnam(varnum) = 'uprt    '
        vartit(varnum) = 'X-velocity perturbation (m/s)'
        varlev(varnum) = nzlg

        varnum = varnum + 1
        varnam(varnum) = 'v       '
        vartit(varnum) = 'Y-velocity total wind (m/s)'
        varlev(varnum) = nzlg

        varnum = varnum + 1
        varnam(varnum) = 'vprt    '
        vartit(varnum) = 'Y-velocity perturbation (m/s)'
        varlev(varnum) = nzlg

        varnum = varnum + 1
        varnam(varnum) = 'w       '
        vartit(varnum) = 'Z-velocity total wind (m/s)'
        varlev(varnum) = nzlg

        varnum = varnum + 1
        varnam(varnum) = 'pt      '
        vartit(varnum) = 'Potential Temperature (K)'
        varlev(varnum) = nzlg

        varnum = varnum + 1
        varnam(varnum) = 'ptprt   '
        vartit(varnum) = 'Perturbation Potential Temperature (K)'
        varlev(varnum) = nzlg

        varnum = varnum + 1
        varnam(varnum) = 'p       '
        vartit(varnum) = 'Pressure (pascal)'
        varlev(varnum) = nzlg

        varnum = varnum + 1
        varnam(varnum) = 'pprt    '
        vartit(varnum) = 'Perturbation Pressure (pascal)'
        varlev(varnum) = nzlg

        varnum = varnum + 1
        varnam(varnum) = 'vort    '
        vartit(varnum) = 'Vertical vorticity (1/s)'
        varlev(varnum) = nzlg

        varnum = varnum + 1
        varnam(varnum) = 'div     '
        vartit(varnum) = 'Horizontal divergence (1/s)'
        varlev(varnum) = nzlg

      END IF

      IF ( mstout == 1 ) THEN

        varnum = varnum + 1
        varnam(varnum) = 'qv      '
        vartit(varnum) = 'Water Vapor Mixing Ratio (g/kg)'
        varlev(varnum) = nzlg

        varnum = varnum + 1
        varnam(varnum) = 'qvprt   '
        vartit(varnum) = 'Water Vapor Mixing Ratio '//                    &
                         'Perturbation (g/kg)'
        varlev(varnum) = nzlg

        DO nq = 1,nscalar
          varnum = varnum + 1
          varnam(varnum) = TRIM(qnames(nq))//'      '
          vartit(varnum) = TRIM(qdescp(nq))
          varlev(varnum) = nzlg
        END DO

        IF ( rainout == 1 ) THEN

          varnum = varnum + 1
          varnam(varnum) = 'raing   '
          vartit(varnum) = 'Grid Supersaturation Rain '
          varlev(varnum) = 0

          varnum = varnum + 1
          varnam(varnum) = 'rainc '
          vartit(varnum) = 'Cumulus Convection Rain '
          varlev(varnum) = 0


        END IF

        IF ( prcout == 1 ) THEN

          varnum = varnum + 1
          varnam(varnum) = 'prcrt1 '
          vartit(varnum) = 'Total precipitation rate '
          varlev(varnum) = 0

          varnum = varnum + 1
          varnam(varnum) = 'prcrt2 '
          vartit(varnum) = 'Grid scale precipitation rate '
          varlev(varnum) = 0

          varnum = varnum + 1
          varnam(varnum) = 'prcrt3 '
          vartit(varnum) = 'Cumulative precipitation rate '
          varlev(varnum) = 0

          varnum = varnum + 1
          varnam(varnum) = 'prcrt4 '
          vartit(varnum) = 'Microphysics precipitation rate '
          varlev(varnum) = 0

        END IF
      END IF

      IF ( tkeout == 1 ) THEN

        varnum = varnum + 1
        varnam(varnum) = 'tke     '
        vartit(varnum) ='Turbulent kinetic energy (m**2/s)'
        varlev(varnum) = nzlg

      END IF

      IF ( trbout == 1 ) THEN

        varnum = varnum + 1
        varnam(varnum) = 'kmh     '
        vartit(varnum) ='Horizontal Turb. Mixing Coefficient for '        &
                   //'Momentum (m**2/s)'
        varlev(varnum) = nzlg

        varnum = varnum + 1
        varnam(varnum) = 'kmv     '
        vartit(varnum) = 'Vertical Turb. Mixing Coefficient for '         &
                   //'Momentum (m**2/s)'
        varlev(varnum) = nzlg

      END IF

      IF ( sfcout == 1 ) THEN
  !
  !
  ! zpsoil is dumped when sfcout is on.
  !
  !
        varnum = varnum + 1
        varnam(varnum) = 'zpsoil  '
        vartit(varnum) = 'Physical vertical coordinates for the soil model (m)'
        varlev(varnum) = nzsoillg
  !
  ! In the previous code, we distinguished nstyp <= 1 with nstyp > 1,
  ! probably due to disk saving.
  ! Here we discard that, and dumps everything.
  !
        DO is=0,nstyp
          WRITE (chrstr,'(i2.2)') is

          varnum = varnum + 1
          varnam(varnum) = 'tsoil'//chrstr(1:2)//' '
          vartit(varnum) = 'Soil Temperature (K) of soil type '//chrstr(1:2)
          if (is == 0) vartit(varnum) = 'Soil Temperature (K)'
          varlev(varnum) = nzsoillg

          varnum = varnum + 1
          varnam(varnum) = 'qsoil'//chrstr(1:2)//' '
          vartit(varnum) = 'Soil moisture (m**3/m**3) of soil type '//chrstr(1:2)
          if (is == 0) vartit(varnum) = 'Soil moisture (m**3/m**3)'
          varlev(varnum) = nzsoillg

          varnum = varnum + 1
          varnam(varnum) = 'wr'//chrstr(1:2)//' '
          vartit(varnum) = 'Canopy water amount (m) of soil type '//chrstr(1:2)
          if (is ==0)  vartit(varnum) = 'Canopy water amount (m)'
          varlev(varnum) = 0
        END DO

        IF ( snowout == 1 ) THEN
          varnum = varnum + 1
          varnam(varnum) = 'snowd   '
          vartit(varnum) = 'Snow depth (m)'
          varlev(varnum) = 0
        END IF

      END IF

      IF ( radout == 1 ) THEN

        varnum = varnum + 1
        varnam(varnum) = 'radfrc  '
        vartit(varnum) = 'Total radiation forcing (K/s)'
        varlev(varnum) = nzlg

        varnum = varnum + 1
        varnam(varnum) = 'radsw   '
        vartit(varnum) = 'Incoming solar rad. at surface (W/m**2)'
        varlev(varnum) = 0

        varnum = varnum + 1
        varnam(varnum) = 'rnflx   '
        vartit(varnum) = 'Surface radiation heat flux (W/m**2)'
        varlev(varnum) = 0

        varnum = varnum + 1
        varnam(varnum) = 'radswnet'
        vartit(varnum) = 'Net solar rad. (W/m**2)'
        varlev(varnum) = 0

        varnum = varnum + 1
        varnam(varnum) = 'radlwin '
        vartit(varnum) = 'Incoming longwave radiation (W/m**2)'
        varlev(varnum) = 0

      END IF

      IF ( flxout == 1 ) THEN

        varnum = varnum + 1
        varnam(varnum) = 'usflx   '
        vartit(varnum) = 'Surface u-momentum flux (kg/m*s**2)'
        varlev(varnum) = 0

        varnum = varnum + 1
        varnam(varnum) = 'vsflx   '
        vartit(varnum) = 'Surface v-momentum flux (kg/m*s**2)'
        varlev(varnum) = 0

        varnum = varnum + 1
        varnam(varnum) = 'ptsflx  '
        vartit(varnum) = 'Surface heat flux (K*kg/s*m**2)'
        varlev(varnum) = 0

        varnum = varnum + 1
        varnam(varnum) = 'qvsflx  '
        vartit(varnum) = 'Surface moisture flux (kg/s*m**2)'
        varlev(varnum) = 0

      END IF

  !
  ! Writing fileheader
  !
  !
      hdbyte = 0

      WRITE(nchanl) fmtver
      hdbyte = hdbyte + 8 + 1*40

      WRITE(nchanl) runname
      hdbyte = hdbyte + 8 + 1*80

      WRITE(nchanl) nocmnt
      hdbyte = hdbyte + 8 + 1*4

      IF ( nocmnt > 0 ) THEN
        DO l = 1, nocmnt
          WRITE(nchanl) cmnt(l)
        END DO
        hdbyte = hdbyte + nocmnt * ( 8 + 80 )
      END IF

      WRITE(nchanl) nxlg,nylg,nzlg,nzsoillg,nstyp
      hdbyte = hdbyte + 8 + 5*4

      WRITE(nchanl) curtim,tmunit
      hdbyte = hdbyte + 8 + 4 + 10

      idummy = 0
      rdummy = 0.0

      WRITE(nchanl) varout, mstout, iceout, trbout, sfcout,             &
                   rainout, landout, idummy, idummy, totout,            &
                    tkeout, nscalar, mapproj, istgr, month,             &
                       day,   year,   hour, minute, second
      hdbyte = hdbyte + 8 + 20*4

      WRITE(nchanl) P_QC,  P_QR,  P_QI,  P_QS,  P_QG,  P_QH,            &
                    P_NC,  P_NR,  P_NI,  P_NS,  P_NG,  P_NH,            &
                    P_ZR,  P_ZI,  P_ZS,  P_ZG,  P_ZH,idummy,            &
                  idummy,idummy,idummy,idummy,idummy,idummy,            &
                  idummy,idummy,idummy,idummy,idummy,idummy
      hdbyte = hdbyte + 8 + 30*4

      WRITE(nchanl)  umove,   vmove, xgrdorg, ygrdorg,  trulat1,        &
                   trulat2,  trulon,  sclfct,  ntcloud, n0rain,         &
                    n0snow, n0grpl, n0hail, rhoice, rhosnow,            &
                     rhogrpl, rhohail, alpharain, alphaice, alphasnow,  &
                     alphagrpl, alphahail,  rdummy,  rdummy, rdummy,    &
                     tstop, thisdmp, latitud, ctrlat,  ctrlon
      hdbyte = hdbyte + 8 + 30*4

      IF ( totout /= 0 ) THEN
        WRITE(nchanl) hdmpopt, nstyp,  prcout, radout, flxout,          &
                            0,snowout, idummy, idummy, idummy, & ! 0 for snowcvr
                       idummy, idummy, idummy, idummy, idummy,          &
                       idummy, idummy, idummy, idummy, idummy
        hdbyte = hdbyte + 8 + 20*4

        WRITE(nchanl) tstrtdmp,rdummy, rdummy, rdummy, rdummy,          &
                      rdummy, rdummy, rdummy, rdummy, rdummy,           &
                      rdummy, rdummy, rdummy, rdummy, rdummy,           &
                      rdummy, rdummy, rdummy, rdummy, rdummy
        hdbyte = hdbyte + 8 + 20*4
      END IF

      IF ( hdmpopt == 2 ) THEN
        WRITE (nchanl) numhdmp
        hdbyte = hdbyte + 8 + 4
        IF ( numhdmp > 0 ) THEN
          WRITE (nchanl) (hdmptim(i),i=1,numhdmp)
          hdbyte = hdbyte + 8 + numhdmp*4
        END IF
      END IF

      WRITE(nchanl) globx
      hdbyte = hdbyte + 8 + nxlg*4

      WRITE(nchanl) globy
      hdbyte = hdbyte + 8 + nylg*4

      WRITE(nchanl) z
      hdbyte = hdbyte + 8 + nzlg*4

      WRITE(nchanl) varnum
      hdbyte = hdbyte + 8 + 1*4

      DO l = 1, varnum
        WRITE(nchanl) varnam(l),varlev(l),vartit(l)
      END DO
      hdbyte = hdbyte + varnum*( 8 + 8 + 4 + 60 )

      DEALLOCATE(globx, globy)

    END IF   ! firstcall
  END IF     ! myproc == 0

  IF (firstcall) THEN

    IF(landout == 1) THEN
      DO is=1,nstyp
        CALL iedgfill(soiltyp(1,1,is),1,nx,1,nx-1, 1,ny,1,ny-1,       &
                      1,1,1,1)
        CALL edgfill(stypfrct(1,1,is),1,nx,1,nx-1, 1,ny,1,ny-1,       &
                     1,1,1,1)
      END DO

      CALL mpimerge3di(soiltyp, nx, ny, nstyp, out3di)
      CALL mpimerge3d (stypfrct,nx, ny, nstyp, out3d)

      IF (myproc == 0) THEN
        DO is=1,nstyp
          WRITE (nchanl) out3d(:,:,is)
          hdbyte = hdbyte + ( 8 + nxlg*nylg*4 )

          WRITE (nchanl) out3d(:,:,is)
          hdbyte = hdbyte + ( 8 + nxlg*nylg*4 )
        END DO
      END IF

      CALL iedgfill(vegtyp ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL mpimerge3di(vegtyp, nx, ny, 1, out3di)

      CALL edgfill(lai    ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL mpimerge3d (lai,    nx, ny, 1, out3d)

      IF (myproc == 0) THEN
        WRITE (nchanl) out3di (:,:,1)
        hdbyte = hdbyte + ( 8 + nxlg*nylg*4 )

        WRITE (nchanl) out3d(:,:,1)
        hdbyte = hdbyte + ( 8 + nxlg*nylg*4 )
      END IF

      CALL edgfill(roufns ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL mpimerge3d (roufns,nx, ny, 1, out3d)
      IF (myproc == 0) THEN
        WRITE (nchanl) out3d (:,:,1)
        hdbyte = hdbyte + ( 8 + nxlg*nylg*4 )
      END IF

      CALL edgfill(veg    ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL mpimerge3d (veg,nx, ny, 1, out3d)
      IF (myproc == 0) THEN
        WRITE (nchanl) out3d (:,:,1)
        hdbyte = hdbyte + ( 8 + nxlg*nylg*4 )
      END IF

    END IF
  !
  !-----------------------------------------------------------------------
  !
  !  Open the GrADS data control file: runname.ctl
  !
  !-----------------------------------------------------------------------
  !
    CALL edgfill(zp  ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)
    CALL mpimerge3d(zp,nx,ny,nz,out3d)

    IF (myproc == 0) THEN

      cntlen = lfnkey + 10
      gradscntl(1:cntlen) = runname(1:lfnkey)//'.gradscntl'
      CALL fnversn( gradscntl, cntlen )
      CALL getunit (nchout0)

      WRITE (6,'(a)') 'The GrADS data control file is '                   &
                      //gradscntl(1:cntlen)

      OPEN (nchout0, FILE = gradscntl(1:cntlen), STATUS = 'unknown')

      WRITE (nchout0,'(a,a)')                                             &
          'DSET    ',filnam(1:filen)

      WRITE (nchout0,'(a/a)')                                             &
          'TITLE   ARPS model 5.0 joined output for '//runname(1:lfnkey),'*'

      WRITE (nchout0,'(a)')                                               &
          'OPTIONS sequential cray_32bit_ieee big_endian'

      WRITE (nchout0,'(a,i10)')                                           &
          'FILEHEADER ',hdbyte

      WRITE (nchout0,'(a/a)')                                             &
          'UNDEF   -9.e+33','*'

      IF ( mapproj == 2 ) THEN

!        WRITE (nchout0,'(a)')                                             &
!            '* For lat-lon-lev display, uncomment the following 4 lines.'

        WRITE (nchout0,'(a,1x,i8,1x,i3,a,2f12.6,2i3,3f12.6,2f12.2)')      &
            'PDEF',nxlg,nylg,' LCC',lat11,lon11,1,1,                      &
                trulat1,trulat2,trulon,xinc,yinc

        WRITE (nchout0,'(a,1x,i8,a,f10.4,1x,f10.4)')                      &
            'XDEF',nxlg,'  LINEAR  ',lonmin,loninc

        WRITE (nchout0,'(a,1x,i8,a,f10.4,1x,f10.4)')                      &
            'YDEF',nylg,'  LINEAR  ',latmin,latinc

      ELSE

  !
  ! ZDEF is max(nz,nzsoil), because we are sharing the soil model
  ! and the atmosphere model in the same dataset.
  !
        WRITE (nchout0,'(a)')                                             &
            '* For i-j-k display, uncomment the following 3 lines.'

        WRITE (nchout0,'(a,1x,i8,a,2i10)')                                &
            '* XDEF',nxlg,'  LINEAR  ',1,1

        WRITE (nchout0,'(a,1x,i8,a,2i10)')                                &
            '* YDEF',nylg,'  LINEAR  ',1,1

        WRITE (nchout0,'(a,1x,i8,a,2i10/a)')                              &
            '* ZDEF',max(nzlg,nzsoillg),'  LINEAR  ',1,1,'*'

        WRITE (nchout0,'(a)')                                             &
            '* For x-y-z display, uncomment the following 3 lines.'

        WRITE (nchout0,'(a,1x,i8,a,2f15.4)')                              &
            'XDEF',nxlg,'  LINEAR  ',xbgn/1000.,xinc/1000.

        WRITE (nchout0,'(a,1x,i8,a,2f15.4)')                              &
            'YDEF',nylg,'  LINEAR  ',ybgn/1000.,yinc/1000.

      END IF
  !
  ! adding comments for ZDEF
  !
        WRITE (nchout0,'(18(a/),a)')                                      &
            '*          ',                                                &
            '* WARNING: ',                                                &
            '* The associated GrADS dataset includes both the ',          &
            '* atmosphere and soil variables. ',                          &
            '* ZDEF is the larger one between the level numbers ',        &
            '* for the atmosphere and soil models. ',                     &
            '* The actual number of vertical levels for a variable ',     &
            '* can be obtained from the definition of the variable ',     &
            '* located between VARS and ENDVARS.                    ',    &
            '* Vertically, the soil levels increase downward under ',     &
            '* the ground, while atmosphere levels upward above ',        &
            '* the ground. ',                                             &
            '* Since we make the grid locations of the atmospheric ',     &
            '* model as the primary vertical coordinate, when ',          &
            '* plotting the soil variables, the vertical axis is ',       &
            '* actually the atmospheric axis. ',                          &
            '* Thus plotting soil variables is for debug purpose. ',      &
            '*          '

      WRITE (nchout0,'(a,1x,i8,a)')                                     &
      'ZDEF',max(nzlg,nzsoillg),'  LEVELS '
      IF ( ternopt == 0 ) THEN
        WRITE (nchout0,'(8f10.2)')                                      &
            ((zp(2,2,k)+zp(2,2,k+1))/2.,k=1,nzlg-1),zp(2,2,nzlg)
      ELSE
        !
        ! Use AGL height as vertical levels
        !
        WRITE (nchout0,'(8f10.2)')                                      &
          ((zp(2,2,k)+zp(2,2,k+1))/2.-zp(2,2,2),k=1,nz-1),zp(2,2,nz)

!        WRITE (nchout0,'(8f10.2)')                                        &
!            ((z(k)+z(k+1))/2.,k=1,nzlg-1),z(nzlg)
!        WRITE (nchout0,'(a/a/a)')                                         &
!        '* WARNING  The vertical levels were set to computational ',      &
!        '*          coordinates, z(k), because zp is not uniform in ',    &
!        '*          horizontal when terrain option was turned on'
      END IF

      IF ( initopt /= 2 ) THEN
        WRITE (chrstr,'(i2.2,a,i2.2,a,i2.2,a3,i4.4)')                     &
              hour,':',minute,'Z',day,monnam(month),year
      ELSE
        second1 = MOD( second + nint(tstart), 60 )
        minute1 = ( second + nint(tstart) ) / 60
        minute1 = MOD( minute + minute1, 60 )
        hour1   = ( minute + ( second + nint(tstart) ) / 60 ) /60
        hour1   = MOD( hour + hour1, 24 )
        day1    = ( hour + ( minute                                       &
                + ( second + nint(tstart) ) / 60 ) /60 ) / 24
        jday1   = jday + day1

        loopdy  = 0
        IF ( MOD( year, 4 ) == 0 ) loopdy = 1
        year1 = year + jday1 / ( 365 + loopdy )
        jday1 = MOD( jday1, 365 + loopdy )

        month1 = 1

        DO m = 2, 11
          IF (jday1>mndys(m) .AND. jday1<=mndys(m+1)+loopdy) month1 = m
        END DO
        day1 = jday1 - mndys(month1)

        WRITE (chrstr,'(i2.2,a,i2.2,a,i2.2,a3,i4.4)')                     &
              hour1,':',minute1,'Z',day1,monnam(month1),year1

      END IF

      IF ( hdmpopt == 1 ) THEN
        WRITE (nchout0,'(a/a,f10.2,a/a)')                              &
            '* WARNING: The time interval is applied after the time ', &
            '*          level at ',tstrtdmp,' model seconds.',         &
            '*          The first time level was the model initial data.'
      ELSE IF ( hdmpopt == 2 ) THEN
        WRITE (nchout0,'(a/a/a/a)')                                    &
        '* WARNING& ! The actual time levels in the data file may not ',  &
        '*          necessarily be in a constant increment. It was ',  &
        '*          written for the model times listed in the ',       &
        '*          following:'
        DO i=1,numhdmp
          WRITE (nchout0,'(a,f10.2)')                                  &
              '*          ',hdmptim(i)
        END DO
      ELSE
        WRITE(6,*) 'hdmpopt = ',hdmpopt
      END IF

      WRITE (nchout0,'(a,1x,i8,a,a,3x,i2.2,a/a)')                         &
          'TDEF',ntm,'  LINEAR  ',chrstr,tinc,dtunit,'*'

      WRITE (nchout0,'(a,1x,i3)')                                         &
          'VARS',varnum

      DO l = 1, varnum

        WRITE (nchout0,'(a8,1x,i3,1x,a10,2x,a)')                          &
               varnam(l),varlev(l),varparam(l),vartit(l)

      END DO

      WRITE (nchout0,'(a)')                                               &
          'ENDVARS'

      CLOSE (nchout0)
      CALL retunit(nchout0)

    END IF ! myproc == 0

    firstcall = .false.

  END IF ! firstcall
!
!-----------------------------------------------------------------------
!
!  Write data into the GrADS data file
!
!-----------------------------------------------------------------------
!
  CALL mpimerge3d(zp,nx,ny,nz,out3dt)
  IF (myproc == 0) THEN

    IF ( istgr == 0 ) THEN
      DO k=1,nzlg-1
        out3d(:,:,k) = .5 * ( out3dt(:,:,k) + out3dt(:,:,k+1) )
      END DO
      out3d(:,:,nzlg) = out3dt(:,:,nzlg)
    ELSE
      out3d(:,:,:) = out3dt(:,:,:)
    END IF

    DO k=1, nzlg
      WRITE(nchanl) out3d(:,:,k)
    ENDDO
  END IF

  IF ( varout == 1 ) THEN
!
!-----------------------------------------------------------------------
!
!  If varout = 1, Write out u, uprt, v, vprt, w, wprt,
!  pt, ptprt, p, pprt, vort, and div.
!
!-----------------------------------------------------------------------
!
    CALL edgfill(u,   1,nx,1,nx, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL mpimerge3d(u,   nx,ny,nz,out3d)

    CALL edgfill(ubar,1,nx,1,nx, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL mpimerge3d(ubar,nx,ny,nz,out3dt)

    IF (myproc == 0) THEN

      IF ( istgr == 0 ) THEN
        DO k=1,nzlg
          DO j=1,nylg
            DO i=1,nxlg-1
              gtem1(i,j,k) = .5 * ( out3d(i,j,k) + out3d(i+1,j,k) )
              gtem2(i,j,k) = gtem1(i,j,k)                                   &
                          - .5 * ( out3dt(i,j,k) + out3dt(i+1,j,k) )
            END DO
            gtem1(nxlg,j,k) = out3d(nxlg,j,k)
            gtem2(nxlg,j,k) = out3d(nxlg,j,k) - out3dt(nxlg,j,k)
          END DO
        END DO
      ELSE
        DO k=1,nzlg
          DO j=1,nylg
            DO i=1,nxlg
              gtem1(i,j,k) = out3d(i,j,k)
              gtem2(i,j,k) = out3d(i,j,k) - out3dt(i,j,k)
            END DO
          END DO
        END DO
      END IF

      DO k=1, nzlg
        WRITE(nchanl) ((gtem1(i,j,k),i=1,nxlg),j=1,nylg)
      END DO
      DO k=1, nzlg
        WRITE(nchanl) ((gtem2(i,j,k),i=1,nxlg),j=1,nylg)
      END DO
    END IF

    CALL edgfill(v,   1,nx,1,nx-1, 1,ny,1,ny, 1,nz,1,nz-1)
    CALL mpimerge3d(v,   nx,ny,nz,out3d)

    CALL edgfill(vbar,1,nx,1,nx-1, 1,ny,1,ny, 1,nz,1,nz-1)
    CALL mpimerge3d(vbar,nx,ny,nz,out3dt)

    IF (myproc == 0) THEN

      IF ( istgr == 0 ) THEN
        DO k=1,nzlg
          DO i=1,nxlg
            DO j=1,nylg-1
              gtem1(i,j,k) = .5 * ( out3d(i,j,k) + out3d(i,j+1,k) )
              gtem2(i,j,k) = gtem1(i,j,k)                                   &
                          - .5 * ( out3dt(i,j,k) + out3dt(i,j+1,k) )
            END DO
            gtem1(i,nylg,k) = out3d(i,nylg,k)
            gtem2(i,nylg,k) = out3d(i,nylg,k) - out3dt(i,nylg,k)
          END DO
        END DO
      ELSE
        DO k=1,nzlg
          DO j=1,nylg
            DO i=1,nxlg
              gtem1(i,j,k) = out3d(i,j,k)
              gtem2(i,j,k) = out3d(i,j,k) - out3dt(i,j,k)
            END DO
          END DO
        END DO
      END IF

      DO k=1, nzlg
        WRITE(nchanl) ((gtem1(i,j,k),i=1,nxlg),j=1,nylg)
      ENDDO

      DO k=1, nzlg
        WRITE(nchanl) ((gtem2(i,j,k),i=1,nxlg),j=1,nylg)
      ENDDO

    END IF

    CALL edgfill(w   ,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz)
    CALL mpimerge3d(w,nx,ny,nz,out3d)

    IF (myproc == 0) THEN

      IF ( istgr == 0 ) THEN
        DO j=1,nylg
          DO i=1,nxlg
            DO k=1,nzlg-1
              gtem1(i,j,k) = .5 * ( out3d(i,j,k) + out3d(i,j,k+1) )
            END DO
            gtem1(i,j,nzlg) = out3d(i,j,nzlg)
          END DO
        END DO
      ELSE
        DO j=1,nylg
          DO i=1,nxlg
            DO k=1,nzlg
              gtem1(i,j,k) = out3d(i,j,k)
            END DO
          END DO
        END DO
      END IF

      DO k=1, nzlg
        WRITE(nchanl) ((gtem1(i,j,k),i=1,nxlg),j=1,nylg)
      ENDDO
    END IF

    CALL edgfill(ptprt,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL mpimerge3d(ptprt,nx,ny,nz,out3d)

    CALL edgfill(ptbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL mpimerge3d(ptbar,nx,ny,nz,out3dt)
    IF (myproc == 0) THEN
      DO k=1,nzlg
        DO j=1,nylg
          DO i=1,nxlg
            gtem1(i,j,k) = out3d(i,j,k) + out3dt(i,j,k)
          END DO
        END DO
      END DO

      DO k=1, nzlg
        WRITE(nchanl) ((gtem1(i,j,k),i=1,nxlg),j=1,nylg)
      ENDDO
      DO k=1, nzlg
        WRITE(nchanl) ((out3d(i,j,k),i=1,nxlg),j=1,nylg)
      ENDDO
    END IF

    CALL edgfill(pprt,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL mpimerge3d(pprt,nx,ny,nz,out3d)

    CALL edgfill(pbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL mpimerge3d(pbar,nx,ny,nz,out3dt)

    IF (myproc == 0) THEN
      DO k=1,nzlg
        DO j=1,nylg
          DO i=1,nxlg
            gtem1(i,j,k) = out3d(i,j,k) + out3dt(i,j,k)
          END DO
        END DO
      END DO

      DO k=1, nzlg
        WRITE(nchanl) ((gtem1(i,j,k),i=1,nxlg),j=1,nylg)
      END DO
      DO k=1, nzlg
        WRITE(nchanl) ((out3d(i,j,k),i=1,nxlg),j=1,nylg)
      END DO
    END IF

!-----------------------------------------------------------------------
!
!  Vorticity:
!
!-----------------------------------------------------------------------

    DO k=2,nz-2
      DO j=2,ny-2
        DO i=2,nx-2
          tem1(i,j,k)=                                                 &
              (v(i+1,j,k)-v(i-1,j,k)+v(i+1,j+1,k)-v(i-1,j+1,k))/       &
              (4*(x(i+1)-x(i)))-                                       &
              (u(i,j+1,k)-u(i,j-1,k)+u(i+1,j+1,k)-u(i+1,j-1,k))/       &
              (4*(y(j+1)-y(j)))
        END DO
      END DO
    END DO

    DO j=2,ny-2
      DO i=2,nx-2
        tem1(i,j,   1)=tem1(i,j,   2)
        tem1(i,j,nz-1)=tem1(i,j,nz-2)
      END DO
    END DO

    DO k=1,nz-1
      DO j=2,ny-2
          tem1(   1,j,k)= tem1(   2,j,k)
          tem1(nx-1,j,k)= tem1(nx-2,j,k)
      END DO
    END DO

    DO k=1,nz-1
      DO i=1,nx-1
          tem1(i,   1,k)=tem1(i,   2,k)
          tem1(i,ny-1,k)=tem1(i,ny-2,k)
      END DO
    END DO

    CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL mpimerge3d(tem1,nx,ny,nz,out3d)
    IF (myproc == 0) THEN
      DO k=1, nzlg
        WRITE(nchanl) ((out3d(i,j,k),i=1,nxlg),j=1,nylg)
      END DO
    END IF

!-----------------------------------------------------------------------
!
!  Divergernce:
!
!-----------------------------------------------------------------------

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem1(i,j,k) = (u(i+1,j,k)-u(i,j,k))/(x(i+1)-x(i))            &
                      + (v(i,j+1,k)-v(i,j,k))/(y(j+1)-y(j))
        END DO
      END DO
    END DO

    CALL edgfill(tem1,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL mpimerge3d(tem1,nx,ny,nz,out3d)
    IF (myproc == 0) THEN
      DO k=1, nzlg
        WRITE(nchanl) ((out3d(i,j,k),i=1,nxlg),j=1,nylg)
      END DO
    END IF

  END IF       ! End varout

  IF ( mstout == 1 ) THEN
!
!-----------------------------------------------------------------------
!
!  Write out moist variables qv, qvprt, qc, qr, qi, qs, and qh
!
!-----------------------------------------------------------------------
!
    CALL edgfill(qv,   1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL mpimerge3d(qv,nx,ny,nz,out3d)

    CALL edgfill(qvbar,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL mpimerge3d(qvbar,nx,ny,nz,out3dt)
    IF (myproc == 0) THEN
      DO k=1,nzlg
        DO j=1,nylg
          DO i=1,nxlg
            gtem1(i,j,k) = out3d(i,j,k) - out3dt(i,j,k)
          END DO
        END DO
      END DO

      DO k=1, nzlg
        WRITE(nchanl) ((out3d(i,j,k),i=1,nxlg),j=1,nylg)
      END DO
      DO k=1, nzlg
        WRITE(nchanl) ((gtem1(i,j,k),i=1,nxlg),j=1,nylg)
      END DO
    END IF

    DO nq = 1,nscalar
      CALL edgfill(qscalar(:,:,:,nq),1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL mpimerge3d(qscalar(:,:,:,nq),nx,ny,nz,out3d)
      IF (myproc == 0) THEN
        DO k=1, nzlg
          WRITE(nchanl) ((out3d(i,j,k),i=1,nxlg),j=1,nylg)
        END DO
      END IF
    END DO

    IF ( rainout == 1 ) THEN

      CALL edgfill(raing,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL mpimerge3d(raing,nx,ny,1,out3d)

      CALL edgfill(rainc,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
      CALL mpimerge3d(rainc,nx,ny,1,out3dt)

      IF (myproc == 0) THEN
        WRITE(nchanl) ((out3d(i,j,1),i=1,nxlg),j=1,nylg)
        WRITE(nchanl) ((out3dt(i,j,1),i=1,nxlg),j=1,nylg)
      END IF

    END IF     ! End rainout

    IF ( prcout == 1 ) THEN

      CALL edgfill(prcrate,1,nx,1,nx-1, 1,ny,1,ny-1, 1,4,1,4)
      CALL mpimerge3d(prcrate,nx,ny,4,out3d)

      IF (myproc == 0) THEN
        WRITE(nchanl) ((out3d(i,j,1),i=1,nxlg),j=1,nylg)
        WRITE(nchanl) ((out3d(i,j,2),i=1,nxlg),j=1,nylg)
        WRITE(nchanl) ((out3d(i,j,3),i=1,nxlg),j=1,nylg)
        WRITE(nchanl) ((out3d(i,j,4),i=1,nxlg),j=1,nylg)
      END IF

    END IF     ! End prcout

  END IF       ! End mstout

  IF ( tkeout == 1 ) THEN

    CALL edgfill(tke,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL mpimerge3d(tke,nx,ny,nz,out3d)

    IF (myproc == 0) THEN
      DO k=1, nzlg
        WRITE(nchanl) ((out3d(i,j,k),i=1,nxlg),j=1,nylg)
      END DO
    END IF

  END IF

  IF ( trbout == 1 ) THEN
!
!-----------------------------------------------------------------------
!
!  If trbout = 1, write out the turbulence parameter, km.
!
!-----------------------------------------------------------------------
!
    CALL edgfill(kmh,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL mpimerge3d(kmh,nx,ny,nz,out3d)

    IF (myproc == 0) THEN
      DO k=1, nzlg
        WRITE(nchanl) ((out3d(i,j,k),i=1,nxlg),j=1,nylg)
      END DO
    END IF

    CALL edgfill(kmv,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL mpimerge3d(kmv,nx,ny,nz,out3d)

    IF (myproc == 0) THEN
      DO k=1, nzlg
        WRITE(nchanl) ((out3d(i,j,k),i=1,nxlg),j=1,nylg)
      END DO
    END IF


  END IF       ! trbout

  IF ( sfcout == 1 ) THEN
!
!-----------------------------------------------------------------------
!
!  Write out surface variables tsoil, qsoil, wetcanp
!
!-----------------------------------------------------------------------
!
    CALL edgfill(zpsoil,1,nx,1,nx-1,1,ny,1,ny-1,        &
                 1,nzsoil,1,nzsoil)
    CALL mpimerge3d(zpsoil,nx,ny,nzsoil,out3d)

    IF (myproc == 0) THEN
      DO k=1, nzsoillg
        WRITE(nchanl) ((out3d(i,j,k),i=1,nxlg),j=1,nylg)
      END DO
    END IF

    DO is=0,nstyp
      CALL edgfill(tsoil  (1,1,1,is),1,nx,1,nx-1,1,ny,1,ny-1,        &
                   1,nzsoil,1,nzsoil)
      CALL edgfill(qsoil(1,1,1,is),1,nx,1,nx-1,1,ny,1,ny-1,        &
                   1,nzsoil,1,nzsoil)
      CALL edgfill(wetcanp(1,1,is),1,nx,1,nx-1,1,ny,1,ny-1,          &
                   1,1,1,1)
    END DO
    CALL mpimerge4d(tsoil,nx,ny,nzsoil,nstyp+1,outtsoil)
    CALL mpimerge4d(qsoil,nx,ny,nzsoil,nstyp+1,outqsoil)
    CALL mpimerge3d(wetcanp,nx,ny,nstyp+1,out3d)

    IF (myproc == 0) THEN
      DO is=0,nstyp
        DO k=1, nzsoillg
          WRITE(nchanl) ((outtsoil(i,j,k,is),i=1,nxlg),j=1,nylg)
        END DO

        DO k=1, nzsoillg
          WRITE(nchanl) ((outqsoil(i,j,k,is),i=1,nxlg),j=1,nylg)
        END DO

        WRITE(nchanl) ((out3d(i,j,is+1),i=1,nxlg),j=1,nylg)
      END DO
    END IF

    IF (snowout == 1) THEN
      CALL edgfill(snowdpth(1,1),1,nx,1,nx-1,1,ny,1,ny-1,               &
                   1,1,1,1)
      CALL mpimerge3d(snowdpth,nxlg,nylg,1,out3d)

      IF (myproc == 0) THEN
        WRITE(nchanl) ((out3d(i,j,1),i=1,nxlg),j=1,nylg)
      END IF
    END IF

  END IF       ! End sfcout
!
!-----------------------------------------------------------------------
!
!  If radout = 1, write out the radiation arrays
!
!-----------------------------------------------------------------------
!
  IF ( radout == 1 ) THEN

    CALL edgfill(radfrc,1,nx,1,nx-1, 1,ny,1,ny-1, 1,nz,1,nz-1)
    CALL mpimerge3d(radfrc,nx,ny,nz,out3d)
    IF (myproc == 0) THEN
      DO k=1, nzlg
        WRITE(nchanl) ((out3d(i,j,k),i=1,nxlg),j=1,nylg)
      END DO
    END IF

    CALL edgfill(radsw,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    CALL mpimerge3d(radsw,nx,ny,1,out3d)
    IF (myproc == 0) THEN
      WRITE(nchanl) ((out3d(i,j,1),i=1,nxlg),j=1,nylg)
    END IF

    CALL edgfill(rnflx,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    CALL mpimerge3d(rnflx,nx,ny,1,out3d)
    IF (myproc == 0) THEN
      WRITE(nchanl) ((out3d(i,j,1),i=1,nxlg),j=1,nylg)
    END IF

    CALL edgfill(radswnet,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    CALL mpimerge3d(radswnet,nx,ny,1,out3d)
    IF (myproc == 0) THEN
      WRITE(nchanl) ((out3d(i,j,1),i=1,nxlg),j=1,nylg)
    END IF

    CALL edgfill(radlwin,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    CALL mpimerge3d(radlwin,nx,ny,1,out3d)
    IF (myproc == 0) THEN
      WRITE(nchanl) ((out3d(i,j,1),i=1,nxlg),j=1,nylg)
    END IF

  END IF       ! radout
!
!-----------------------------------------------------------------------
!
!  If flxout = 1, write out the surface fluxes
!
!-----------------------------------------------------------------------
!
  IF ( flxout == 1 ) THEN

    CALL edgfill(usflx,1,nx,1,nx, 1,ny,1,ny-1, 1,1,1,1)
    CALL mpimerge3d(usflx,nx,ny,1,out3d)

    CALL edgfill(vsflx,1,nx,1,nx-1, 1,ny,1,ny, 1,1,1,1)
    CALL mpimerge3d(vsflx,nx,ny,1,out3dt)

    CALL edgfill(ptsflx,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    CALL mpimerge3d(ptsflx,nx,ny,1,gtem1)

    CALL edgfill(qvsflx,1,nx,1,nx-1, 1,ny,1,ny-1, 1,1,1,1)
    CALL mpimerge3d(qvsflx,nx,ny,1,gtem2)

    IF (myproc == 0) THEN
      WRITE(nchanl) ((out3d(i,j,1),i=1,nxlg),j=1,nylg)
      WRITE(nchanl) ((out3dt(i,j,1),i=1,nxlg),j=1,nylg)
      WRITE(nchanl) ((gtem1(i,j,1),i=1,nxlg),j=1,nylg)
      WRITE(nchanl) ((gtem2(i,j,1),i=1,nxlg),j=1,nylg)
    END IF

  END IF       ! flxout

  IF (myproc == 0) THEN
    DEALLOCATE (out3d, out3dt, out3di)
    DEALLOCATE (outtsoil, outqsoil)
    DEALLOCATE (gtem1, gtem2)
  END IF

  RETURN
END SUBROUTINE gradsjoindump
