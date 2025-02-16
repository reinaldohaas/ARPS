!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RSTOUT                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE rstout(nx,ny,nz,nzsoil,nstyps,exbcbufsz,                     &
           u,v,w,ptprt,pprt,qv,qscalar,tke,                             &
           udteb, udtwb, vdtnb, vdtsb,                                  &
           pdteb ,pdtwb ,pdtnb ,pdtsb,                                  &
           ubar,vbar,ptbar,pbar,rhostr,qvbar,                           &
           x,y,z,zp,zpsoil,hterain, mapfct,                             &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,qvsfc,                          &
           ptcumsrc,qcumsrc,w0avg,nca,kfraincv,                         &
           cldefi,xland,bmjraincv,                                      &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           raing,rainc,prcrate, exbcbuf, tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Dump out a model restart file at a specified model time. Only permanent
!  arrays in the model (which are saved between time steps) need to be
!  dumped for a model restart. For time dependent variables, two time
!  levels (time tpast and tfuture) are needed for a model restart so
!  fields at both time levels are dumped out.
!
!  NOTE: After you make any changes to this subroutine, you should also
!        change the same code in the subroutine RSTJOINOUT below.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  4/01/1992.
!
!  MODIFICATION HISTORY:
!
!  5/06/92 (M. Xue)
!  Added full documentation.
!
!  5/06/92 (M. Xue)
!  Included grid and terrain data in the restart dump.
!
!  6/2/92 (M. Xue and H. Jin)
!  Further facelift.
!
!  2/10/93 (K. Droegemeier)
!  Cleaned up documentation.
!
!  02/07/1995 (Yuhe Liu)
!  Added a new 2-D permanent array, veg(nx,ny), to the argument list
!
!  05/05/1995 (M. Xue)
!  Added rainc and raing into the restart data dump.
!
!  08/22/1995 (M. Xue)
!  Added ptcumsrc and qvcumsrc into the restart data dump.
!
!  08/30/1995 (Yuhe Liu)
!  Added the external boundary data into the restart dump
!
!  2/2/96 (Donghai Wang & Yuhe Liu)
!  Added a 3-D array, mapfct, for map projection factor.
!
!  08/01/97 (Zonghui Huo)
!  Added Kain-fritsch cumulus parameterization scheme.
!
!  11/06/97 (D. Weber)
!  Added three additional levels to the mapfct array.  The three
!  levels (4,5,6) represent the inverse of the first three in order.
!  The inverse map factors are computed to improve efficiency.
!
!  4/15/1998 (Donghai Wang)
!  Added the source terms to the right hand terms of the qc,qr,qi,qs
!  equations due to the K-F cumulus parameterization.
!
!  4/15/1998 (Donghai Wang)
!  Added the running average vertical velocity (array w0avg)
!  for the K-F cumulus parameterization scheme.
!
!  12/09/1998 (Donghai Wang)
!  Added the snow cover.
!
!  03/13/2002 (Eric Kemp)
!  Added arrays for WRF BMJ cumulus scheme.
!
!  April 2002 (Fanyou Kong)
!  Added cnvctopt=5 option for new WRF K-F (KF_ETA) scheme
!  05/14/2002 (J. Brotzge)
!  Added arrays, modified call statements to permit multiple soil schemes
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of grid points in the soil
!
!    u        x component of velocity at times tpast and tpresent (m/s)
!    v        y component of velocity at times tpast and tpresent (m/s)
!    w        Vertical component of Cartesian velocity at times
!    ptprt    Perturbation potential temperature at times tpast and
!             tpresent (K)
!    pprt     Perturbation pressure at times tpast and tpresent (Pascal)
!
!    qv       Water vapor specific humidity at times tpast and tpresent (kg/kg)
!    qc       Cloud water mixing ratio at times tpast and tpresent (kg/kg)
!    qr       Rainwater mixing ratio at times tpast and tpresent (kg/kg)
!    qi       Cloud ice mixing ratio at times tpast and tpresent (kg/kg)
!    qs       Snow mixing ratio at times tpast and tpresent (kg/kg)
!    qh       Hail mixing ratio at times tpast and tpresent (kg/kg)
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!
!    udteb    Time tendency of u field at east boundary (m/s**2)
!    udtwb    Time tendency of u field at west boundary (m/s**2)
!
!    vdtnb    Time tendency of v field at north boundary (m/s**2)
!    vdtsb    Time tendency of v field at south boundary (m/s**2)
!
!    pdteb    Time tendency of pprt field at east boundary (PASCAL/s)
!    pdtwb    Time tendency of pprt field at west boundary (PASCAL/s)
!    pdtnb    Time tendency of pprt field at north boundary (PASCAL/s)
!    pdtsb    Time tendency of pprt field at south boundary (PASCAL/s)
!
!    ubar     Base state zonal velocity component (m/s)
!    vbar     Base state meridional velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhostr   Base state density (kg/m**3)times j3.
!    qvbar    Base state water vapor specific humidity (kg/kg)
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!    zpsoil   Vertical coordinate of grid points in the soil (m)
!    hterain  Terrain height (m)
!
!    mapfct   Map factors at scalar, u and v points
!
!    soiltyp  Soil type
!    vegtyp   Vegetation type
!    lai      Leaf Area Index
!    roufns   Surface roughness
!    veg      Vegetation fraction
!
!    qvsfc    Effective S.H. at sfc.
!    tsoil    Soil temperature (K)
!    qsoil    Soil moisture (m**3/m**3)
!    wetcanp  Canopy water amount
!    ptcumsrc Source term in pt-equation due to cumulus parameterization
!    qcumsrc Source term in water equations due to cumulus parameterization
!    kfraincv   K-F convective rainfall (cm)
!    nca      K-F counter for CAPE release
!    cldefi   BMJ cloud efficiency
!    xland    BMJ land/sea mask
!    bmjraincv   BMJ convective rainfall (cm)
!
!    radfrc   Radiation forcing (K)
!    radsw    Solar radiation reaching the surface
!    rnflx    Net absorbed radiation by the surface
!    radswnet Net solar radiation, SWin - SWout
!    radlwin  Incoming longwave radiation
!
!    raing    Grid scale rainfall
!    rainc    Convective rainfall
!
!  OUTPUT:
!
!    None
!
!  WORK ARRAY:
!
!    tem1    Temporary work array.
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
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'bndry.inc'
  INCLUDE 'exbc.inc'
  INCLUDE 'timelvls.inc'
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of grid points in the soil

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)
  REAL :: qscalar(nx,ny,nz,nt,nscalar)
  REAL :: tke   (nx,ny,nz,nt)  ! Turbulent Kinetic Energy ((m/s)**2)

  REAL :: udteb (ny,nz)        ! T-tendency of u at e-boundary (m/s**2)
  REAL :: udtwb (ny,nz)        ! T-tendency of u at w-boundary (m/s**2)

  REAL :: vdtnb (nx,nz)        ! T-tendency of v at n-boundary (m/s**2)
  REAL :: vdtsb (nx,nz)        ! T-tendency of v at s-boundary (m/s**2)

  REAL :: pdteb (ny,nz)        ! T-tendency of pprt at e-boundary (PASCAL/s)
  REAL :: pdtwb (ny,nz)        ! T-tendency of pprt at w-boundary (PASCAL/s)
  REAL :: pdtnb (nx,nz)        ! T-tendency of pprt at n-boundary (PASCAL/s)
  REAL :: pdtsb (nx,nz)        ! T-tendency of pprt at s-boundary (PASCAL/s)

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL :: rhostr(nx,ny,nz)     ! Base state air density (kg/m**3) times j3.
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                               ! (kg/kg)

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
  REAL :: zpsoil(nx,ny,nzsoil) ! The physical height coordinate defined at
                               ! w-point of the soil.

  REAL :: hterain(nx,ny)       ! Terrain height (m).

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

  INTEGER :: nstyps                   ! Number of soil types
  INTEGER :: soiltyp (nx,ny,nstyps)   ! Soil type
  REAL    :: stypfrct(nx,ny,nstyps)
  INTEGER :: vegtyp(nx,ny)            ! Vegetation type
  REAL    :: lai    (nx,ny)           ! Leaf Area Index
  REAL    :: roufns (nx,ny)           ! Surface roughness
  REAL    :: veg    (nx,ny)           ! Vegetation fraction

  REAL :: tsoil   (nx,ny,nzsoil,0:nstyps) ! Soil temperature(K)
  REAL :: qsoil   (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m**3/m**3)
  REAL :: wetcanp (nx,ny,       0:nstyps) ! Canopy water amount
  REAL :: snowdpth(nx,ny)                 ! Snow depth (m)
  REAL :: qvsfc   (nx,ny,       0:nstyps) ! Effective specific humidity
                                          ! at the surface (kg/kg)

  REAL :: ptcumsrc(nx,ny,nz)   ! Source term in pt-equation due
                               ! to cumulus parameterization
  REAL :: qcumsrc(nx,ny,nz,5)  ! Source term in water equations due
                               ! to cumulus parameterization:
                               ! qcumsrc(1,1,1,1) for qv equation
                               ! qcumsrc(1,1,1,2) for qc equation
                               ! qcumsrc(1,1,1,3) for qr equation
                               ! qcumsrc(1,1,1,4) for qi equation
                               ! qcumsrc(1,1,1,5) for qs equation
  REAL :: w0avg(nx,ny,nz)      ! a closing running average vertical
                               ! velocity in 10min for K-F scheme
  REAL :: kfraincv(nx,ny)      ! K-F convective rainfall (cm)
  INTEGER :: nca(nx,ny)        ! K-F counter for CAPE release

  REAL,INTENT(IN) :: cldefi(nx,ny)    ! BMJ cloud efficiency
  REAL,INTENT(IN) :: xland(nx,ny)     ! BMJ land mask
                                      ! (1.0 = land, 2.0 = sea)
  REAL,INTENT(IN) :: bmjraincv(nx,ny) ! BMJ convective rainfall (cm)

  REAL :: radfrc(nx,ny,nz)     ! Radiation forcing (K/s)
  REAL :: radsw(nx,ny)         ! Solar radiation reacing the surface
  REAL :: rnflx(nx,ny)         ! Net absorbed radiation by the surface
  REAL :: radswnet (nx,ny)     ! Net solar radiation, SWin - SWout
  REAL :: radlwin  (nx,ny)     ! Incoming longwave radiation

  REAL :: raing(nx,ny)         ! Grid supersaturation rain
  REAL :: rainc(nx,ny)         ! Cumulus convective rain
  REAL :: prcrate(nx,ny,4)     ! precipitation rate (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precipitation rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulus precip. rate
                               ! prcrate(1,1,4) = microphysics precip. rate

  INTEGER :: exbcbufsz            ! EXBC buffer size
  REAL    :: exbcbuf( exbcbufsz ) ! EXBC buffer array

  REAL    :: tem1  (nx,ny,nz)     ! Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: basrstout ! Control parameter for the base state
                       ! array output
  INTEGER :: grdrstout ! Control parameter for the grid array output
  INTEGER :: icerstout ! Control parameter for the ice variable output
  INTEGER :: sfcrstout ! Control parameter for the surface variable
                       ! output
  INTEGER :: prcrstout ! Control parameter for the precip. rate and rain output
  INTEGER :: rcumout   ! Control parameter for ptcumsrc and qcumsrc output
  INTEGER :: exbcout   ! Control parameter for external boundary output
  INTEGER :: mapfout   ! Control parameter for map factor output
  INTEGER :: radrstout ! Control parameter for radiation forcing output
  INTEGER :: kfrsout   ! Control parameter for Kain-Fritsch output
  INTEGER :: bmjsout   ! Control parameter for WRF BMJ output

  INTEGER :: idummy
  INTEGER :: istat
  INTEGER :: lrstof
  REAL    :: rdummy

  INTEGER :: i, var, varsize, nq

  INTEGER :: qcbcrd, qrbcrd, qibcrd, qsbcrd, qgbcrd, qhbcrd,            &
             ncbcrd, nrbcrd, nibcrd, nsbcrd, ngbcrd, nhbcrd,            &
                     zrbcrd, zibcrd, zsbcrd, zgbcrd, zhbcrd,            &
             ccbcrd

  CHARACTER(LEN=256) :: filnamr, outdirname
  INTEGER            :: nchout1, oldirnam
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

  CALL get_output_dirname(1,dirname,curtim,1,outdirname,istat)
  oldirnam = LEN_TRIM(outdirname)
!
!-----------------------------------------------------------------------
!
!  Get a name for the restart data file.
!
!-----------------------------------------------------------------------
!

  CALL gtrstfn(runname(1:lfnkey),outdirname,oldirnam,curtim,            &
               mgrid,nestgrd, rstoutf, lrstof )

  CALL getunit( rstount )

  OPEN(UNIT=rstount,FILE=trim(rstoutf(1:lrstof)),FORM='unformatted',    &
       STATUS='new',IOSTAT=istat)

  IF( istat /= 0) THEN

    WRITE(6,'(/a,i2,/a/)')                                              &
        ' Error occured when opening restart output file '              &
        //rstoutf(1:lrstof)//                                           &
        ' using FORTRAN unit ',rstount,' Program stopped in RSTOUT.'
    CALL arpsstop('arpsstop called from RSTOUT problem opening file',1)

  END IF

  WRITE(6,'('' DUMPING OUT RESTART FILE AT TIME '',F10.2,               &
  &       ''(s) in FILE '',a,'' using fortran channel no '', i2)')      &
          curtim, rstoutf(1:lrstof),rstount
!
!-----------------------------------------------------------------------
!
!  Write out the restart data:
!
!-----------------------------------------------------------------------
!
  WRITE(rstount) curtim
  WRITE(rstount) nx,ny,nz, nzsoil

  basrstout = 1
  grdrstout = 1
  icerstout = ice
  mapfout  = 1

  prcrstout = 0
  IF ( moist /= 0 ) prcrstout = 1

  sfcrstout = 0
  IF( sfcphy /= 0 ) sfcrstout = 1

  rcumout=0
  IF ( cnvctopt /= 0 ) rcumout=1

  exbcout = 0
  IF ( lbcopt == 2 ) exbcout = 1

  radrstout = 0
  IF ( radopt > 0 ) radrstout = 1

  kfrsout=0
  IF ( cnvctopt == 3 .OR. cnvctopt == 5) kfrsout=1

  bmjsout=0
  IF ( cnvctopt == 4 ) bmjsout=1

  idummy = 0
  WRITE(rstount) basrstout,grdrstout,icerstout,sfcrstout,prcrstout,     &
                 rcumout,  exbcout,  mapfout,  radrstout,nstyp,         &
                 kfrsout,  rayklow,  bmjsout,  nscalar,  idummy,        &
                 idummy,   idummy,   idummy,   idummy,   idummy,        &
                 idummy

  WRITE(rstount) P_QC,  P_QR,  P_QI,  P_QS,  P_QG,  P_QH,               &
                 P_NC,  P_NR,  P_NI,  P_NS,  P_NG,  P_NH,               &
                 P_ZR,  P_ZI,  P_ZS,  P_ZG,  P_ZH,  P_CC,               &
                 idummy,idummy,idummy,idummy,idummy,idummy,             &
                 idummy,idummy,idummy,idummy,idummy,idummy

  rdummy = 0.0
  WRITE(rstount) dx,dy,dz,umove,vmove,                                  &
                 xgrdorg,ygrdorg,trulat1,trulat2,trulon,                &
                 sclfct,latitud,ctrlat,ctrlon,ntcloud,n0rain,           &
                 n0snow,n0grpl,n0hail,rhoice,rhosnow,rhogrpl,rhohail,   &
                 alpharain,alphaice,alphasnow,alphagrpl,alphahail,rdummy

  IF( grdrstout == 1) THEN

    WRITE(rstount) x
    WRITE(rstount) y
    WRITE(rstount) z
    WRITE(rstount) zp
    WRITE(rstount) zpsoil

  END IF

  IF( basrstout == 1) THEN

    WRITE(rstount) ubar
    WRITE(rstount) vbar
    WRITE(rstount) ptbar
    WRITE(rstount) pbar
    WRITE(rstount) rhostr
    WRITE(rstount) qvbar

  END IF

  CALL cpyary3d(nx,ny,nz,u    (1,1,1,tpast), tem1)
  WRITE(rstount) tem1

  CALL cpyary3d(nx,ny,nz,v    (1,1,1,tpast), tem1)
  WRITE(rstount) tem1

  CALL cpyary3d(nx,ny,nz,w    (1,1,1,tpast), tem1)
  WRITE(rstount) tem1

  CALL cpyary3d(nx,ny,nz,ptprt(1,1,1,tpast), tem1)
  WRITE(rstount) tem1

  CALL cpyary3d(nx,ny,nz,pprt (1,1,1,tpast), tem1)
  WRITE(rstount) tem1

  CALL cpyary3d(nx,ny,nz,qv   (1,1,1,tpast), tem1)
  WRITE(rstount) tem1

  DO nq = 1,nscalar
    CALL cpyary3d(nx,ny,nz,qscalar(1,1,1,tpast,nq), tem1)
    WRITE(rstount) tem1
  END DO

  CALL cpyary3d(nx,ny,nz,tke  (1,1,1,tpast), tem1)
  WRITE(rstount) tem1

  CALL cpyary3d(nx,ny,nz,u    (1,1,1,tpresent), tem1)
  WRITE(rstount) tem1

  CALL cpyary3d(nx,ny,nz,v    (1,1,1,tpresent), tem1)
  WRITE(rstount) tem1

  CALL cpyary3d(nx,ny,nz,w    (1,1,1,tpresent), tem1)
  WRITE(rstount) tem1

  CALL cpyary3d(nx,ny,nz,ptprt(1,1,1,tpresent), tem1)
  WRITE(rstount) tem1

  CALL cpyary3d(nx,ny,nz,pprt (1,1,1,tpresent), tem1)
  WRITE(rstount) tem1

  CALL cpyary3d(nx,ny,nz,qv   (1,1,1,tpresent), tem1)
  WRITE(rstount) tem1

  DO nq = 1,nscalar
    CALL cpyary3d(nx,ny,nz,qscalar(1,1,1,tpresent,nq), tem1)
    WRITE(rstount) tem1
  END DO

  CALL cpyary3d(nx,ny,nz,tke  (1,1,1,tpresent), tem1)
  WRITE(rstount) tem1

  WRITE(rstount) udteb
  WRITE(rstount) udtwb

  WRITE(rstount) vdtnb
  WRITE(rstount) vdtsb

  WRITE(rstount) pdteb
  WRITE(rstount) pdtwb
  WRITE(rstount) pdtnb
  WRITE(rstount) pdtsb

  IF ( sfcrstout /= 0 ) THEN

    PRINT *,'write out sfc/soil variables:'

    WRITE(rstount) soiltyp
    WRITE(rstount) stypfrct
    WRITE(rstount) vegtyp
    WRITE(rstount) lai
    WRITE(rstount) roufns
    WRITE(rstount) veg

    WRITE(rstount) qvsfc
    WRITE(rstount) tsoil
    WRITE(rstount) qsoil
    WRITE(rstount) wetcanp
    WRITE(rstount) snowdpth

  END IF


  IF ( prcrstout /= 0 ) THEN
    WRITE(rstount) raing
    WRITE(rstount) rainc
    WRITE(rstount) prcrate
  END IF

  IF ( rcumout /=  0 ) THEN

    WRITE(rstount) ptcumsrc
    WRITE(rstount) qcumsrc

  END IF

  IF ( exbcout /= 0 ) THEN

    qcbcrd=0; qrbcrd=0; qibcrd=0; qsbcrd=0; qgbcrd=0; qhbcrd=0
    ncbcrd=0; nrbcrd=0; nibcrd=0; nsbcrd=0; ngbcrd=0; nhbcrd=0
              zrbcrd=0; zibcrd=0; zsbcrd=0; zgbcrd=0; zhbcrd=0
    ccbcrd=0

    IF (P_QC > 0) qcbcrd = qscalarbcrd(P_QC)
    IF (P_QR > 0) qrbcrd = qscalarbcrd(P_QR)
    IF (P_QI > 0) qibcrd = qscalarbcrd(P_QI)
    IF (P_QS > 0) qsbcrd = qscalarbcrd(P_QS)
    IF (P_QG > 0) qgbcrd = qscalarbcrd(P_QG)
    IF (P_QH > 0) qhbcrd = qscalarbcrd(P_QH)

    IF (P_NC > 0) ncbcrd = qscalarbcrd(P_NC)
    IF (P_NR > 0) nrbcrd = qscalarbcrd(P_NR)
    IF (P_NI > 0) nibcrd = qscalarbcrd(P_NI)
    IF (P_NS > 0) nsbcrd = qscalarbcrd(P_NS)
    IF (P_NG > 0) ngbcrd = qscalarbcrd(P_NG)
    IF (P_NH > 0) nhbcrd = qscalarbcrd(P_NH)

    IF (P_ZR > 0) zrbcrd = qscalarbcrd(P_ZR)
    IF (P_ZI > 0) zibcrd = qscalarbcrd(P_ZI)
    IF (P_ZS > 0) zsbcrd = qscalarbcrd(P_ZS)
    IF (P_ZG > 0) zgbcrd = qscalarbcrd(P_ZG)
    IF (P_ZH > 0) zhbcrd = qscalarbcrd(P_ZH)

    IF (P_CC > 0) ccbcrd = qscalarbcrd(P_CC)

    WRITE(rstount) abstfcst0, abstfcst,                                 &
                   ubcrd,vbcrd,wbcrd,ptbcrd,prbcrd,                     &
                   qvbcrd,qcbcrd,qrbcrd,qibcrd,qsbcrd,qhbcrd,qgbcrd,    &
                   ncbcrd,nrbcrd,nibcrd,nsbcrd,ngbcrd,nhbcrd,           &
                   zrbcrd,zibcrd,zsbcrd,zgbcrd,zhbcrd,ccbcrd

    varsize = nx*ny*nz
    DO var = 1, exbcbufsz, varsize
      WRITE(rstount) (exbcbuf(i),i=var,var+varsize-1)
    END DO
  END IF

  IF ( mapfout == 1 ) THEN
    WRITE(rstount) mapfct
  END IF

  IF ( radrstout == 1 ) THEN
    WRITE(rstount) radfrc
    WRITE(rstount) radsw
    WRITE(rstount) rnflx
    WRITE(rstount) radswnet
    WRITE(rstount) radlwin
  END IF

  IF ( kfrsout /= 0 ) THEN

    WRITE(rstount) w0avg
    WRITE(rstount) nca
    WRITE(rstount) kfraincv

  END IF

  IF ( bmjsout /= 0 ) THEN

    WRITE(rstount) cldefi
    WRITE(rstount) xland
    WRITE(rstount) bmjraincv

  END IF

  CLOSE (UNIT=rstount)
  CALL retunit( rstount )
!
!-----------------------------------------------------------------------
!
!  Compress the restart file using system command.
!
!-----------------------------------------------------------------------
!
  IF( filcmprs == 1 ) CALL cmprs( rstoutf(1:lrstof) )
!
!-----------------------------------------------------------------------
!
!  Create ready file, indicating restart dump writing is complete
!
!-----------------------------------------------------------------------
!
  IF( readyfl == 1 ) THEN
    WRITE (filnamr,'(a)') trim(rstoutf(1:lrstof)) // "_ready"
    CALL getunit( nchout1 )
    OPEN (UNIT=nchout1,FILE=trim(filnamr))
    WRITE (nchout1,'(a)') trim(rstoutf(1:lrstof))
    CLOSE (nchout1)
    CALL retunit ( nchout1 )
  END IF

  RETURN
END SUBROUTINE rstout

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RSTJOINOUT                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE rstjoinout(nx,ny,nz,nzsoil,nstyps,exbcbufsz,                 &
           u,v,w,ptprt,pprt,qv,qscalar,tke,                             &
           udteb, udtwb, vdtnb, vdtsb,                                  &
           pdteb ,pdtwb ,pdtnb ,pdtsb,                                  &
           ubar,vbar,ptbar,pbar,rhostr,qvbar,                           &
           x,y,z,zp,zpsoil,hterain, mapfct,                             &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,qvsfc,                          &
           ptcumsrc,qcumsrc,w0avg,nca,kfraincv,                         &
           cldefi,xland,bmjraincv,                                      &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           raing,rainc,prcrate, exbcbuf, tem1)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Dump out a model restart file at a specified model time. Only permanent
!  arrays in the model (which are saved between time steps) need to be
!  dumped for a model restart. For time dependent variables, two time
!  levels (time tpast and tfuture) are needed for a model restart so
!  fields at both time levels are dumped out.
!
!  RSTJOINOUT dumps joined restart file for message passing mode.
!
!  NOTE: This suboutine should be consistent with the normal one, RSTOUT.
!        Any changes here should also be copied to subroutine RSTOUT above.
!
!        The parameter list is the same as that of RSTOUT. This will
!        make it easier to call RSTOUT and RSTJOINOUT at the same place
!        of the calling subroutine. It will also be easy to combine
!        these two subroutines into one later if necessary.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  2/25/2003.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of grid points in the soil
!
!    u        x component of velocity at times tpast and tpresent (m/s)
!    v        y component of velocity at times tpast and tpresent (m/s)
!    w        Vertical component of Cartesian velocity at times
!    ptprt    Perturbation potential temperature at times tpast and
!             tpresent (K)
!    pprt     Perturbation pressure at times tpast and tpresent (Pascal)
!
!    qv       Water vapor specific humidity at times tpast and tpresent (kg/kg)
!    qc       Cloud water mixing ratio at times tpast and tpresent (kg/kg)
!    qr       Rainwater mixing ratio at times tpast and tpresent (kg/kg)
!    qi       Cloud ice mixing ratio at times tpast and tpresent (kg/kg)
!    qs       Snow mixing ratio at times tpast and tpresent (kg/kg)
!    qh       Hail mixing ratio at times tpast and tpresent (kg/kg)
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!
!    udteb    Time tendency of u field at east boundary (m/s**2)
!    udtwb    Time tendency of u field at west boundary (m/s**2)
!
!    vdtnb    Time tendency of v field at north boundary (m/s**2)
!    vdtsb    Time tendency of v field at south boundary (m/s**2)
!
!    pdteb    Time tendency of pprt field at east boundary (PASCAL/s)
!    pdtwb    Time tendency of pprt field at west boundary (PASCAL/s)
!    pdtnb    Time tendency of pprt field at north boundary (PASCAL/s)
!    pdtsb    Time tendency of pprt field at south boundary (PASCAL/s)
!
!    ubar     Base state zonal velocity component (m/s)
!    vbar     Base state meridional velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhostr   Base state density (kg/m**3)times j3.
!    qvbar    Base state water vapor specific humidity (kg/kg)
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!    zpsoil   Vertical coordinate of grid points in the soil (m)
!    hterain  Terrain height (m)
!
!    mapfct   Map factors at scalar, u and v points
!
!    soiltyp  Soil type
!    vegtyp   Vegetation type
!    lai      Leaf Area Index
!    roufns   Surface roughness
!    veg      Vegetation fraction
!
!    qvsfc    Effective S.H. at sfc.
!    tsoil    Soil temperature (K)
!    qsoil    Soil moisture (m**3/m**3)
!    wetcanp  Canopy water amount
!    ptcumsrc Source term in pt-equation due to cumulus parameterization
!    qcumsrc Source term in water equations due to cumulus parameterization
!    kfraincv   K-F convective rainfall (cm)
!    nca      K-F counter for CAPE release
!    cldefi   BMJ cloud efficiency
!    xland    BMJ land/sea mask
!    bmjraincv   BMJ convective rainfall (cm)
!
!    radfrc   Radiation forcing (K)
!    radsw    Solar radiation reaching the surface
!    rnflx    Net absorbed radiation by the surface
!    radswnet Net solar radiation, SWin - SWout
!    radlwin  Incoming longwave radiation
!
!    raing    Grid scale rainfall
!    rainc    Convective rainfall
!
!  OUTPUT:
!
!    None
!
!  WORK ARRAY:
!
!    tem1    Temporary work array.
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
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'bndry.inc'
  INCLUDE 'exbc.inc'
  INCLUDE 'mp.inc'
  INCLUDE 'timelvls.inc'
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of grid points in the soil

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)

  REAL :: qscalar(nx,ny,nz,nt,nscalar)

  REAL :: tke   (nx,ny,nz,nt)  ! Turbulent Kinetic Energy ((m/s)**2)

  REAL :: udteb (ny,nz)        ! T-tendency of u at e-boundary (m/s**2)
  REAL :: udtwb (ny,nz)        ! T-tendency of u at w-boundary (m/s**2)

  REAL :: vdtnb (nx,nz)        ! T-tendency of v at n-boundary (m/s**2)
  REAL :: vdtsb (nx,nz)        ! T-tendency of v at s-boundary (m/s**2)

  REAL :: pdteb (ny,nz)        ! T-tendency of pprt at e-boundary (PASCAL/s)
  REAL :: pdtwb (ny,nz)        ! T-tendency of pprt at w-boundary (PASCAL/s)
  REAL :: pdtnb (nx,nz)        ! T-tendency of pprt at n-boundary (PASCAL/s)
  REAL :: pdtsb (nx,nz)        ! T-tendency of pprt at s-boundary (PASCAL/s)

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL :: rhostr(nx,ny,nz)     ! Base state air density (kg/m**3) times j3.
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                               ! (kg/kg)

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
  REAL :: zpsoil(nx,ny,nzsoil) ! The physical height coordinate defined at
                               ! w-point of the soil.

  REAL :: hterain(nx,ny)       ! Terrain height (m).

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

  INTEGER :: nstyps                   ! Number of soil types
  INTEGER :: soiltyp (nx,ny,nstyps)   ! Soil type
  REAL    :: stypfrct(nx,ny,nstyps)
  INTEGER :: vegtyp(nx,ny)            ! Vegetation type
  REAL    :: lai    (nx,ny)           ! Leaf Area Index
  REAL    :: roufns (nx,ny)           ! Surface roughness
  REAL    :: veg    (nx,ny)           ! Vegetation fraction

  REAL :: tsoil  (nx,ny,nzsoil,0:nstyps) ! Soil temperature(K)
  REAL :: qsoil  (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny,       0:nstyps) ! Canopy water amount
  REAL :: snowdpth(nx,ny)                ! Snow depth (m)
  REAL :: qvsfc  (nx,ny,       0:nstyps) ! Effective specific humidity
                                         ! at the surface (kg/kg)

  REAL :: ptcumsrc(nx,ny,nz)   ! Source term in pt-equation due
                               ! to cumulus parameterization
  REAL :: qcumsrc(nx,ny,nz,5)  ! Source term in water equations due
                               ! to cumulus parameterization:
                               ! qcumsrc(1,1,1,1) for qv equation
                               ! qcumsrc(1,1,1,2) for qc equation
                               ! qcumsrc(1,1,1,3) for qr equation
                               ! qcumsrc(1,1,1,4) for qi equation
                               ! qcumsrc(1,1,1,5) for qs equation
  REAL :: w0avg(nx,ny,nz)      ! a closing running average vertical
                               ! velocity in 10min for K-F scheme
  REAL :: kfraincv(nx,ny)      ! K-F convective rainfall (cm)
  INTEGER :: nca(nx,ny)        ! K-F counter for CAPE release

  REAL,INTENT(IN) :: cldefi(nx,ny)    ! BMJ cloud efficiency
  REAL,INTENT(IN) :: xland(nx,ny)     ! BMJ land mask
                                      ! (1.0 = land, 2.0 = sea)
  REAL,INTENT(IN) :: bmjraincv(nx,ny) ! BMJ convective rainfall (cm)

  REAL :: radfrc(nx,ny,nz)     ! Radiation forcing (K/s)
  REAL :: radsw(nx,ny)         ! Solar radiation reacing the surface
  REAL :: rnflx(nx,ny)         ! Net absorbed radiation by the surface
  REAL :: radswnet (nx,ny)     ! Net solar radiation, SWin - SWout
  REAL :: radlwin  (nx,ny)     ! Incoming longwave radiation

  REAL :: raing(nx,ny)         ! Grid supersaturation rain
  REAL :: rainc(nx,ny)         ! Cumulus convective rain
  REAL :: prcrate(nx,ny,4)     ! precipitation rate (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precipitation rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulus precip. rate
                               ! prcrate(1,1,4) = microphysics precip. rate

  INTEGER :: exbcbufsz            ! EXBC buffer size
  REAL    :: exbcbuf( exbcbufsz ) ! EXBC buffer array

  REAL    :: tem1  (nx,ny,nz)     ! Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: basrstout ! Control parameter for the base state
                       ! array output
  INTEGER :: grdrstout ! Control parameter for the grid array output
  INTEGER :: icerstout ! Control parameter for the ice variable output
  INTEGER :: sfcrstout ! Control parameter for the surface variable
                       ! output
  INTEGER :: prcrstout ! Control parameter for the precip. rate and rain output
  INTEGER :: rcumout   ! Control parameter for ptcumsrc and qcumsrc output
  INTEGER :: exbcout   ! Control parameter for external boundary output
  INTEGER :: mapfout   ! Control parameter for map factor output
  INTEGER :: radrstout ! Control parameter for radiation forcing output
  INTEGER :: kfrsout   ! Control parameter for Kain-Fritsch output
  INTEGER :: bmjsout   ! Control parameter for WRF BMJ output

  INTEGER :: idummy
  INTEGER :: istat
  INTEGER :: lrstof
  REAL    :: rdummy

  INTEGER :: qcbcrd, qrbcrd, qibcrd, qsbcrd, qgbcrd, qhbcrd,            &
             ncbcrd, nrbcrd, nibcrd, nsbcrd, ngbcrd, nhbcrd,            &
                     zrbcrd, zibcrd, zsbcrd, zgbcrd, zhbcrd,            &
             ccbcrd

  CHARACTER(LEN=256) :: filnamr, outdirname
  INTEGER            :: nchout1, oldirnam

  INTEGER :: nxlg, nylg, n3rd
  INTEGER :: nq, var
  REAL, ALLOCATABLE :: out1d(:)
  REAL, ALLOCATABLE :: out2d(:,:)
  REAL, ALLOCATABLE :: out3d(:,:,:)
  REAL, ALLOCATABLE :: out4d(:,:,:,:)
  REAL, ALLOCATABLE :: out4dq(:,:,:,:)  ! for qcumsrc(nx,ny,nz,5)

  REAL, ALLOCATABLE :: out2dew(:,:)
  REAL, ALLOCATABLE :: out2dns(:,:)

  INTEGER, ALLOCATABLE :: out2di(:,:)
  INTEGER, ALLOCATABLE :: out3di(:,:,:)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  nxlg = nproc_x*(nx-3) + 3
  nylg = nproc_y*(ny-3) + 3
  n3rd = MAX(nz,nzsoil,nstyps+1,8)

  ALLOCATE (out1d( MAX(nxlg,nylg,nz) ),stat=istat)
  CALL check_alloc_status(istat, "RSTJOINOUT:out1d")

  ALLOCATE (out2d(nxlg,nylg),stat=istat)
  CALL check_alloc_status(istat, "RSTJOINOUT:out2d")

  ALLOCATE (out2di(nxlg,nylg),stat=istat)
  CALL check_alloc_status(istat, "RSTJOINOUT:out2di")

  ALLOCATE (out3d( nxlg,nylg, n3rd ),stat=istat)
  CALL check_alloc_status(istat, "RSTJOINOUT:out3d")

  ALLOCATE (out3di( nxlg,nylg, nstyps ),stat=istat)
  CALL check_alloc_status(istat, "RSTJOINOUT:out3di")

  ALLOCATE (out4d(nxlg, nylg, nzsoil, nstyps+1),stat=istat)
  CALL check_alloc_status(istat, "RSTJOINOUT:out4d")

  ALLOCATE (out4dq(nxlg, nylg, nz, 5),stat=istat)
  CALL check_alloc_status(istat, "RSTJOINOUT:out4dq")

  ALLOCATE (out2dew(nylg,nz),stat=istat)
  CALL check_alloc_status(istat, "RSTJOINOUT:out2dew")

  ALLOCATE (out2dns(nxlg,nz),stat=istat)
  CALL check_alloc_status(istat, "RSTJOINOUT:out2dns")
!
!-----------------------------------------------------------------------
!
!  Get a name for the restart data file.
!
!-----------------------------------------------------------------------
!
  IF(myproc == 0) THEN

    CALL get_output_dirname(0,dirname,curtim,1,outdirname,istat)
    oldirnam = LEN_TRIM(outdirname)

    CALL gtrstfn(runname(1:lfnkey),outdirname,oldirnam,curtim,          &
                 mgrid,nestgrd, rstoutf, lrstof )

    CALL getunit( rstount )

    OPEN(UNIT=rstount,FILE=trim(rstoutf(1:lrstof)),FORM='unformatted',    &
         STATUS='new',IOSTAT=istat)

    IF( istat /= 0) THEN

      WRITE(6,'(/a,i2,/a/)')                                              &
          ' Error occured when opening restart output file '              &
          //rstoutf(1:lrstof)//                                           &
          ' using FORTRAN unit ',rstount,' Program stopped in RSTJOINOUT.'
      CALL arpsstop('arpsstop called from RSTJOINOUT while opening file.',1)

    END IF

    WRITE(6,'('' DUMPING OUT RESTART FILE AT TIME '',F10.2,               &
    &       ''(s) in FILE '',a,'' using fortran channel no '', i2)')      &
            curtim, rstoutf(1:lrstof),rstount
!
!-----------------------------------------------------------------------
!
!  Write out the restart data:
!
!-----------------------------------------------------------------------
!
    WRITE(rstount) curtim
    WRITE(rstount) nxlg,nylg,nz, nzsoil

  END IF  ! myproc == 0

  basrstout = 1
  grdrstout = 1
  icerstout = ice
  mapfout  = 1

  prcrstout = 0
  IF ( moist /= 0 ) prcrstout = 1

  sfcrstout = 0
  IF( sfcphy /= 0 ) sfcrstout = 1

  rcumout=0
  IF ( cnvctopt /= 0 ) rcumout=1

  exbcout = 0
  IF ( lbcopt == 2 ) exbcout = 1

  radrstout = 0
  IF ( radopt > 0 ) radrstout = 1

  kfrsout=0
  IF ( cnvctopt == 3 .OR. cnvctopt == 5) kfrsout=1

  bmjsout=0
  IF ( cnvctopt == 4 ) bmjsout=1

  IF (myproc == 0) THEN
    idummy = 0

    WRITE(rstount) basrstout,grdrstout,icerstout,sfcrstout,prcrstout,   &
                   rcumout,  exbcout,  mapfout,  radrstout,nstyp,       &
                   kfrsout,  rayklow,  bmjsout,  nscalar,  idummy,      &
                   idummy,   idummy,   idummy,   idummy,   idummy,      &
                   idummy

    WRITE(rstount) P_QC,  P_QR,  P_QI,  P_QS,  P_QG,  P_QH,             &
                   P_NC,  P_NR,  P_NI,  P_NS,  P_NG,  P_NH,             &
                   P_ZR,  P_ZI,  P_ZS,  P_ZG,  P_ZH,  P_CC,             &
                   idummy,idummy,idummy,idummy,idummy,idummy,           &
                   idummy,idummy,idummy,idummy,idummy,idummy

    rdummy = 0.0
    WRITE(rstount) dx,dy,dz,umove,vmove,                                &
                   xgrdorg,ygrdorg,trulat1,trulat2,trulon,              &
                   sclfct,latitud,ctrlat,ctrlon,ntcloud,n0rain,         &
                   n0snow,n0grpl,n0hail,rhoice,rhosnow,rhogrpl,rhohail, &
                   alpharain,alphaice,alphasnow,alphagrpl,alphahail,rdummy

  END IF  ! myproc == 0

  IF( grdrstout == 1) THEN

    CALL mpimerge1dx(x,nx,out1d)
    IF(myproc == 0) WRITE(rstount) out1d(1:nxlg)
    CALL mpimerge1dy(y,ny,out1d)
    IF(myproc == 0) WRITE(rstount) out1d(1:nylg)
    IF(myproc == 0) WRITE(rstount) z
    CALL mpimerge3d(zp,nx,ny,nz,out3d)
    IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nz)
    CALL mpimerge3d(zpsoil,nx,ny,nzsoil,out3d)
    IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nzsoil)

  END IF

  IF( basrstout == 1) THEN

    CALL mpimerge3d(ubar,nx,ny,nz,out3d)
    IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nz)
    CALL mpimerge3d(vbar,nx,ny,nz,out3d)
    IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nz)
    CALL mpimerge3d(ptbar,nx,ny,nz,out3d)
    IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nz)
    CALL mpimerge3d(pbar,nx,ny,nz,out3d)
    IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nz)
    CALL mpimerge3d(rhostr,nx,ny,nz,out3d)
    IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nz)
    CALL mpimerge3d(qvbar,nx,ny,nz,out3d)
    IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nz)

  END IF

  CALL mpimerge3d(u(:,:,:,tpast),nx,ny,nz,out3d)
  IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nz)

  CALL mpimerge3d(v(:,:,:,tpast),nx,ny,nz,out3d)
  IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nz)

  CALL mpimerge3d(w(:,:,:,tpast),nx,ny,nz,out3d)
  IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nz)

  CALL mpimerge3d(ptprt(:,:,:,tpast),nx,ny,nz,out3d)
  IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nz)

  CALL mpimerge3d(pprt(:,:,:,tpast),nx,ny,nz,out3d)
  IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nz)

  CALL mpimerge3d(qv(:,:,:,tpast),nx,ny,nz,out3d)
  IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nz)

  DO nq = 1,nscalar
    CALL mpimerge3d(qscalar(:,:,:,tpast,nq),nx,ny,nz,out3d)
    IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nz)
  END DO

  CALL mpimerge3d(tke(:,:,:,tpast),nx,ny,nz,out3d)
  IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nz)

  CALL mpimerge3d(u(:,:,:,tpresent),nx,ny,nz,out3d)
  IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nz)

  CALL mpimerge3d(v(:,:,:,tpresent),nx,ny,nz,out3d)
  IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nz)

  CALL mpimerge3d(w(:,:,:,tpresent),nx,ny,nz,out3d)
  IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nz)

  CALL mpimerge3d(ptprt(:,:,:,tpresent),nx,ny,nz,out3d)
  IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nz)

  CALL mpimerge3d(pprt(:,:,:,tpresent),nx,ny,nz,out3d)
  IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nz)

  CALL mpimerge3d(qv(:,:,:,tpresent),nx,ny,nz,out3d)
  IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nz)

  DO nq = 1,nscalar
    CALL mpimerge3d(qscalar(:,:,:,tpresent,nq),nx,ny,nz,out3d)
    IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nz)
  END DO

  CALL mpimerge3d(tke(:,:,:,tpresent),nx,ny,nz,out3d)
  IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nz)

  CALL mpimerge2dy(udteb,ny,nz,nproc_x,out2dew,istat)
  IF(myproc == 0) WRITE(rstount) out2dew
  CALL mpimerge2dy(udtwb,ny,nz,      1,out2dew,istat)
  IF(myproc == 0) WRITE(rstount) out2dew

  CALL mpimerge2dx(vdtnb,nx,nz,nproc_y,out2dns,istat)
  IF(myproc == 0) WRITE(rstount) out2dns
  CALL mpimerge2dx(vdtsb,nx,nz,      1,out2dns,istat)
  IF(myproc == 0) WRITE(rstount) out2dns

  CALL mpimerge2dy(pdteb,ny,nz,nproc_x,out2dew,istat)
  IF(myproc == 0) WRITE(rstount) out2dew
  CALL mpimerge2dy(pdtwb,ny,nz,      1,out2dew,istat)
  IF(myproc == 0) WRITE(rstount) out2dew

  CALL mpimerge2dx(pdtnb,nx,nz,nproc_y,out2dns,istat)
  IF(myproc == 0) WRITE(rstount) out2dns
  CALL mpimerge2dx(pdtsb,nx,nz,      1,out2dns,istat)
  IF(myproc == 0) WRITE(rstount) out2dns

  IF ( sfcrstout /= 0 ) THEN

    IF(myproc == 0) PRINT *,'write out sfc/soil variables:'

    CALL mpimerge3di(soiltyp,nx,ny,nstyps,out3di)
    IF(myproc == 0) WRITE(rstount) out3di
    CALL mpimerge3d(stypfrct,nx,ny,nstyps,out3d)
    IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nstyps)
    CALL mpimerge2di(vegtyp,nx,ny,out2di)
    IF(myproc == 0) WRITE(rstount) out2di
    CALL mpimerge2d(lai,nx,ny,out2d)
    IF(myproc == 0) WRITE(rstount) out2d
    CALL mpimerge2d(roufns,nx,ny,out2d)
    IF(myproc == 0) WRITE(rstount) out2d
    CALL mpimerge2d(veg,nx,ny,out2d)
    IF(myproc == 0) WRITE(rstount) out2d

    CALL mpimerge3d(qvsfc,nx,ny,nstyps+1,out3d)
    IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nstyps+1)

    CALL mpimerge4d(tsoil,nx,ny,nzsoil,nstyps+1,out4d)
    IF(myproc == 0) WRITE(rstount) out4d
    CALL mpimerge4d(qsoil,nx,ny,nzsoil,nstyps+1,out4d)
    IF(myproc == 0) WRITE(rstount) out4d

    CALL mpimerge3d(wetcanp,nx,ny,nstyps+1,out3d)
    IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nstyps+1)

    CALL mpimerge2d(snowdpth,nx,ny,out2d)
    IF(myproc == 0) WRITE(rstount) out2d

  END IF


  IF ( prcrstout /= 0 ) THEN
    CALL mpimerge2d(raing,nx,ny,out2d)
    IF(myproc == 0) WRITE(rstount) out2d
    CALL mpimerge2d(rainc,nx,ny,out2d)
    IF(myproc == 0) WRITE(rstount) out2d
    CALL mpimerge3d(prcrate,nx,ny,4,out3d)
    IF(myproc == 0) WRITE(rstount) out3d(:,:,1:4)
  END IF

  IF ( rcumout /=  0 ) THEN

    CALL mpimerge3d(ptcumsrc,nx,ny,nz,out3d)
    IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nz)
    CALL mpimerge4d(qcumsrc,nx,ny,nz,5,out4dq)
    IF(myproc == 0) WRITE(rstount) out4dq

  END IF

  IF ( exbcout /= 0 ) THEN

    qcbcrd=0; qrbcrd=0; qibcrd=0; qsbcrd=0; qgbcrd=0; qhbcrd=0
    ncbcrd=0; nrbcrd=0; nibcrd=0; nsbcrd=0; ngbcrd=0; nhbcrd=0
              zrbcrd=0; zibcrd=0; zsbcrd=0; zgbcrd=0; zhbcrd=0
    ccbcrd=0

    IF (P_QC > 0) qcbcrd = qscalarbcrd(P_QC)
    IF (P_QR > 0) qrbcrd = qscalarbcrd(P_QR)
    IF (P_QI > 0) qibcrd = qscalarbcrd(P_QI)
    IF (P_QS > 0) qsbcrd = qscalarbcrd(P_QS)
    IF (P_QG > 0) qgbcrd = qscalarbcrd(P_QG)
    IF (P_QH > 0) qhbcrd = qscalarbcrd(P_QH)

    IF (P_NC > 0) ncbcrd = qscalarbcrd(P_NC)
    IF (P_NR > 0) nrbcrd = qscalarbcrd(P_NR)
    IF (P_NI > 0) nibcrd = qscalarbcrd(P_NI)
    IF (P_NS > 0) nsbcrd = qscalarbcrd(P_NS)
    IF (P_NG > 0) ngbcrd = qscalarbcrd(P_NG)
    IF (P_NH > 0) nhbcrd = qscalarbcrd(P_NH)

    IF (P_ZR > 0) zrbcrd = qscalarbcrd(P_ZR)
    IF (P_ZI > 0) zibcrd = qscalarbcrd(P_ZI)
    IF (P_ZS > 0) zsbcrd = qscalarbcrd(P_ZS)
    IF (P_ZG > 0) zgbcrd = qscalarbcrd(P_ZG)
    IF (P_ZH > 0) zhbcrd = qscalarbcrd(P_ZH)

    IF (P_CC > 0) ccbcrd = qscalarbcrd(P_CC)

    IF (myproc == 0) WRITE(rstount) abstfcst0, abstfcst,                &
                   ubcrd,vbcrd,wbcrd,ptbcrd,prbcrd,                     &
                   qvbcrd,qcbcrd,qrbcrd,qibcrd,qsbcrd,qhbcrd,qgbcrd,    &
                   ncbcrd,nrbcrd,nibcrd,nsbcrd,ngbcrd,nhbcrd,           &
                   zrbcrd,zibcrd,zsbcrd,zgbcrd,zhbcrd,ccbcrd

    ! Assume each variable in exbcbuf is of size (nx*ny*nz).

    DO var = 1, exbcbufsz, nx*ny*nz
      CALL mpimerge3d(exbcbuf(var),nx,ny,nz,out3d)
      IF(myproc ==0) WRITE(rstount) out3d(:,:,1:nz)
    END DO

  END IF

  IF ( mapfout == 1 ) THEN
    CALL mpimerge3d(mapfct,nx,ny,8,out3d)
    IF(myproc == 0) WRITE(rstount) out3d(:,:,1:8)
  END IF

  IF ( radrstout == 1 ) THEN
    CALL mpimerge3d(radfrc,nx,ny,nz,out3d)
    IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nz)
    CALL mpimerge2d(radsw,nx,ny,out2d)
    IF(myproc == 0) WRITE(rstount) out2d
    CALL mpimerge2d(rnflx,nx,ny,out2d)
    IF(myproc == 0) WRITE(rstount) out2d
    CALL mpimerge2d(radswnet,nx,ny,out2d)
    IF(myproc == 0) WRITE(rstount) out2d
    CALL mpimerge2d(radlwin,nx,ny,out2d)
    IF(myproc == 0) WRITE(rstount) out2d
  END IF

  IF ( kfrsout /= 0 ) THEN

    CALL mpimerge3d(w0avg,nx,ny,nz,out3d)
    IF(myproc == 0) WRITE(rstount) out3d(:,:,1:nz)
    CALL mpimerge2di(nca,nx,ny,out2di)
    IF(myproc == 0) WRITE(rstount) out2di
    CALL mpimerge2d(kfraincv,nx,ny,out2d)
    IF(myproc == 0) WRITE(rstount) out2d

  END IF

  IF ( bmjsout /= 0 ) THEN

    CALL mpimerge2d(cldefi,nx,ny,out2d)
    IF(myproc == 0) WRITE(rstount) out2d
    CALL mpimerge2d(xland,nx,ny,out2d)
    IF(myproc == 0) WRITE(rstount) out2d
    CALL mpimerge2d(bmjraincv,nx,ny,out2d)
    IF(myproc == 0) WRITE(rstount) out2d

  END IF

  IF(myproc == 0) THEN
    CLOSE (UNIT=rstount)
    CALL retunit( rstount )
  END IF

  DEALLOCATE(out1d,out2d,out3d)
  DEALLOCATE(out2di,out3di)
  DEALLOCATE(out4d,out4dq)
  DEALLOCATE(out2dew,out2dns)
!
!-----------------------------------------------------------------------
!
!  Compress the restart file using system command.
!
!-----------------------------------------------------------------------
!
  IF(filcmprs == 1 .AND. myproc == 0) CALL cmprs( rstoutf(1:lrstof) )
!
!-----------------------------------------------------------------------
!
!  Create ready file, indicating restart dump writing is complete
!
!-----------------------------------------------------------------------
!
  IF( readyfl == 1 .AND. myproc == 0) THEN
    WRITE (filnamr,'(a)') trim(rstoutf(1:lrstof)) // "_ready"
    CALL getunit( nchout1 )
    OPEN (UNIT=nchout1,FILE=trim(filnamr))
    WRITE (nchout1,'(a)') trim(rstoutf(1:lrstof))
    CLOSE (nchout1)
    CALL retunit ( nchout1 )
  END IF

  RETURN
END SUBROUTINE rstjoinout
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RSTIN                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE rstin(nx,ny,nz,nzsoil,nts,nstyps,exbcbufsz,                  &
           u,v,w,ptprt,pprt,qv,qscalar,tke,                             &
           udteb, udtwb, vdtnb, vdtsb,                                  &
           pdteb ,pdtwb ,pdtnb ,pdtsb,                                  &
           ubar,vbar,ptbar,pbar,rhostr,qvbar,                           &
           x,y,z,zp,zpsoil,hterain,mapfct,j1,j2,j3,j3soil,              &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,qvsfc,                          &
           ptcumsrc,qcumsrc,w0avg,nca,kfraincv,                         &
           cldefi,xland,bmjraincv,                                      &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           raing,rainc,prcrate, exbcbuf, tem1, tem2)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in data from a restart file to initialize u,v,w,prprt,pprt,
!  qv,qc,qr,qi,qs,qh and tke at time tpast and tpresent, the base state
!  variables ubar,vbar,ptbar,pbar,rhostr,qvbar, and the time tendencies
!    of variables at the lateral boundaries.
!
!  Fields at tfuture are set to the values at tpresent.
!
!  This subroutine also sets the value of tstart.
!
!  NOTE: After you make any changes to this subroutine, you should also
!        change the same code in the subroutine RSTJOINOUT below.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  4/01/1992.
!
!  MODIFICATION HISTORY:
!
!  5/06/92 (M. Xue)
!  Added full documentation.
!
!  10/15/1992 (M. Xue)
!  Reading of grid and base state arrays added.
!
!  2/10/93 (K. Droegemeier)
!  Cleaned up documentation.
!
!  9/7/93 (Ming Xue)
!  Changed cpyary to cpyary3d.
!
!  9/7/93 (A. Shapiro & Ming Xue)
!  Adjustment to tpast values after umove and vmove are changed.
!
!  02/07/1995 (Yuhe Liu)
!  Added a new 2-D permanent array, veg(nx,ny), to the argument list
!
!  05/05/1995 (M. Xue)
!  Added rainc and raing into the restart data dump.
!
!  08/22/1995 (M. Xue)
!  Added ptcumsrc and qvcumsrc into the restart data dump.
!
!  08/30/1995 (Yuhe Liu)
!  Added the external boundary data into the restart dump
!
!  9/10/1995 (M. Xue)
!  When umove or vmove in arps40.input is 999.0, (umove,vmove)
!  in the restart data is used. No adjustment will be
!  made to the wind fields in this case.
!
!  2/2/96 (Donghai Wang & Yuhe Liu)
!  Added a 3-D array, mapfct, for map projection factor.
!
!  08/01/97 (Zonghui Huo)
!  Added Kain-fritsch cumulus parameterization scheme.
!
!  11/06/97 (D. Weber)
!  Added three additional levels to the mapfct array.  The three
!  levels (4,5,6) represent the inverse of the first three in order.
!  The inverse map factors are computed to improve efficiency.
!
!  4/15/1998 (Donghai Wang)
!  Added the source terms to the right hand terms of the qc,qr,qi,qs
!  equations due to the K-F cumulus parameterization.
!
!  4/15/1998 (Donghai Wang)
!  Added the running average vertical velocity (array w0avg)
!  for the K-F cumulus parameterization scheme.
!
!  12/09/1998 (Donghai Wang)
!  Added the snow cover.
!
!  03/13/2002 (Eric Kemp)
!  Added arrays for WRF BMJ cumulus scheme.
!
!  05/14/2002  (J. Brotzge)
!  Added arrays, modified call statements to allow for multiple soil schemes
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of grid points in the soil
!    nts      Number of time levels to be initialized.
!
!  OUTPUT:
!
!    u        x component of velocity at times tpast and tpresent (m/s)
!    v        y component of velocity at times tpast and tpresent (m/s)
!    w        Vertical component of Cartesian velocity at times
!             tpast and tpresent (m/s)
!    ptprt    Perturbation potential temperature at times tpast and
!             tpresent (K)
!    pprt     Perturbation pressure at times tpast and tpresent (Pascal)
!
!    qv       Water vapor specific humidity at times tpast and tpresent (kg/kg)
!    qc       Cloud water mixing ratio at times tpast and tpresent (kg/kg)
!    qr       Rainwater mixing ratio at times tpast and tpresent (kg/kg)
!    qi       Cloud ice mixing ratio at times tpast and tpresent (kg/kg)
!    qs       Snow mixing ratio at times tpast and tpresent (kg/kg)
!    qh       Hail mixing ratio at times tpast and tpresent (kg/kg)
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!
!    udteb    Time tendency of u field at east boundary (m/s**2)
!    udtwb    Time tendency of u field at west boundary (m/s**2)
!
!    vdtnb    Time tendency of v field at north boundary (m/s**2)
!    vdtsb    Time tendency of v field at south boundary (m/s**2)
!
!    pdteb    Time tendency of pprt field at east boundary (PASCAL/s)
!    pdtwb    Time tendency of pprt field at west boundary (PASCAL/s)
!    pdtnb    Time tendency of pprt field at north boundary (PASCAL/s)
!    pdtsb    Time tendency of pprt field at south boundary (PASCAL/s)
!
!    ubar     Base state zonal velocity component (m/s)
!    vbar     Base state meridional velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhostr   Base state density (kg/m**3) times j3.
!    qvbar    Base state water vapor specific humidity (kg/kg)
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!    zpsoil   Vertical coordinate of grid points in the soil (m)
!    hterain  Terrain height (m)
!
!    mapfct   Map factors at scalar, u and v points
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    j3soil   Coordinate transformation Jacobian  d(zpsoil)/dz
!
!    soiltyp  Soil type
!    vegtyp   Vegetation type
!    lai      Leaf Area Index
!    roufns   Surface roughness
!    veg      Vegetation fraction
!
!    qvsfc    Effective S.H. at sfc.
!    tsoil    Soil temperature (K)
!    qsoil    Soil moisture (m**3/m**3)
!    wetcanp  Canopy water amount
!    ptcumsrc Source term in pt-equation due to cumulus parameterization
!    qcumsrc Source term in water equations due to cumulus parameterization
!    kfraincv   K-F convective rainfall (cm)
!    nca      K-F counter for CAPE release
!    cldefi   BMJ cloud efficiency
!    xland    BMJ land/sea mask
!    bmjraincv   BMJ convective rainfall (cm)
!
!    radfrc   Radiation forcing (K)
!    radsw    Solar radiation reaching the surface
!    rnflx    Net absorbed radiation by the surface
!    radswnet Net shortwave radiation
!    radlwin  Incoming longwave radiation
!
!    raing    Grid scale rainfall
!    rainc    Convective rainfall
!
!    tstart   The time when the time integration starts, which is set to
!             the time of the restart data
!
!    tem1     Temporary work array
!    tem2     Temporary work array
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
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'bndry.inc'
  INCLUDE 'exbc.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nts          ! Number of time levels to be initialized.
  INTEGER :: tpast        ! Index of time level for the past time.
  INTEGER :: tpresent     ! Index of time level for the present time.
  INTEGER :: tfuture      ! Index of time level for the future time.

  INTEGER :: nx,ny,nz     ! Number of grid points in 3 directions
  INTEGER :: nzsoil       ! Number of grid points in 3 directions

  REAL :: u     (nx,ny,nz,nts) ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nts) ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nts) ! Total w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz,nts) ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nts) ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz,nts) ! Water vapor specific humidity (kg/kg)
  REAL :: qscalar(nx,ny,nz,nts,nscalar)
  REAL :: tke   (nx,ny,nz,nts) ! Turbulent Kinetic Energy ((m/s)**2)

  REAL :: udteb (ny,nz)        ! T-tendency of u at e-boundary (m/s**2)
  REAL :: udtwb (ny,nz)        ! T-tendency of u at w-boundary (m/s**2)

  REAL :: vdtnb (nx,nz)        ! T-tendency of v at n-boundary (m/s**2)
  REAL :: vdtsb (nx,nz)        ! T-tendency of v at s-boundary (m/s**2)

  REAL :: pdteb (ny,nz)        ! T-tendency of pprt at e-boundary (PASCAL/s)
  REAL :: pdtwb (ny,nz)        ! T-tendency of pprt at w-boundary (PASCAL/s)
  REAL :: pdtnb (nx,nz)        ! T-tendency of pprt at n-boundary (PASCAL/s)
  REAL :: pdtsb (nx,nz)        ! T-tendency of pprt at s-boundary (PASCAL/s)

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL :: rhostr(nx,ny,nz)     ! Base state air density (kg/m**3) time j3.
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                               ! (kg/kg)

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
  REAL :: zpsoil(nx,ny,nzsoil) ! The physical height coordinate defined at
                               ! w-point of the soil.

  REAL :: hterain(nx,ny)       ! Terrain height (m).

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian -d(zp)/dx.
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian -d(zp)/dy.
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian  d(zp)/dz.
  REAL :: j3soil(nx,ny,nzsoil) ! Coordinate transformation Jacobian  d(zpsoil)/dz.

  INTEGER :: nstyps                   ! Number of soil types
  INTEGER :: soiltyp (nx,ny,nstyps)   ! Soil type
  REAL    :: stypfrct(nx,ny,nstyps)
  INTEGER :: vegtyp (nx,ny)           ! Vegetation type
  REAL    :: lai    (nx,ny)           ! Leaf Area Index
  REAL    :: roufns (nx,ny)           ! Surface roughness
  REAL    :: veg    (nx,ny)           ! Vegetation fraction

  REAL :: qvsfc  (nx,ny,0:nstyps) ! Effective S. H. at the surface (kg/kg)
  REAL :: tsoil  (nx,ny,nzsoil,0:nstyps) ! Soil temperature(K)
  REAL :: qsoil  (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny,0:nstyps) ! Canopy water amount
  REAL :: snowdpth(nx,ny)         ! Snow depth (m)
  REAL :: ptcumsrc(nx,ny,nz)   ! Source term in pt-equation due
                               ! to cumulus parameterization
  REAL :: qcumsrc(nx,ny,nz,5)  ! Source term in water equations due
                               ! to cumulus parameterization:
                               ! qcumsrc(1,1,1,1) for qv equation
                               ! qcumsrc(1,1,1,2) for qc equation
                               ! qcumsrc(1,1,1,3) for qr equation
                               ! qcumsrc(1,1,1,4) for qi equation
                               ! qcumsrc(1,1,1,5) for qs equation
  REAL :: w0avg(nx,ny,nz)      ! a closing running average vertical
                               ! velocity in 10min for K-F scheme
  REAL :: kfraincv(nx,ny)      ! K-F convective rainfall (cm)
  INTEGER :: nca(nx,ny)        ! K-F counter for CAPE release

  REAL,INTENT(OUT) :: cldefi(nx,ny)    ! BMJ cloud efficiency
  REAL,INTENT(OUT) :: xland(nx,ny)     ! BMJ land mask
  REAL,INTENT(OUT) :: bmjraincv(nx,ny) ! BMJ convective rainfall (cm)
                                       ! (1.0 = land, 2.0 = sea)

  REAL :: radfrc(nx,ny,nz)     ! Radiation forcing (K/s)
  REAL :: radsw(nx,ny)         ! Solar radiation reacing the surface
  REAL :: rnflx(nx,ny)         ! Net absorbed radiation by the surface
  REAL :: radswnet(nx,ny)      ! Net shortwave radiation
  REAL :: radlwin(nx,ny)       ! Incoming longwave radiation

  REAL :: raing(nx,ny)         ! Grid supersaturation rain
  REAL :: rainc(nx,ny)         ! Cumulus convective rain
  REAL :: prcrate(nx,ny,4)     ! precipitation rate (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precipitation rate
                               ! prcrate(1,1,2) = grid scale precip.  rate
                               ! prcrate(1,1,3) = cumulus precip.  rate
                               ! prcrate(1,1,4) = microphysics precip. rate

  INTEGER :: exbcbufsz         ! EXBC buffer size
  REAL :: exbcbuf( exbcbufsz ) ! EXBC buffer array

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array.

  INTEGER :: grdrstin ! Parameter indicating if the restart data contains
                      ! the grid variables.
  INTEGER :: basrstin ! Parameter indicating if the restart data contains
                      ! the base state variables.
  INTEGER :: icerstin ! Parameter indicating if the restart data contains
                      ! the ice variables.
  INTEGER :: sfcrstin ! Parameter indicating if the restart data contains
                      ! the surface variables.
  INTEGER :: prcrsin  ! Parameter indicating if the restart data contains
                      ! precipitation rate and rainfall
  INTEGER :: rcumin   ! Parameter indicating if the cumulus source terms
                      ! data are present.
  INTEGER :: exbcin   ! Parameter indicating if the external boundary
                      ! data are present.
  INTEGER :: mapfin   ! Parameter indicating if the map factor
                      ! data are present.
  INTEGER :: radrstin ! Parameter indicating if the radiation forcing
                      ! arrays are present.
  INTEGER :: kfrsin   ! Parameter indicating if k-f variable exists
  INTEGER :: bmjsin   ! Parameter indicating if BMJ variable exists

  REAL :: umoveold    ! The domain translation speed of in restart data
  REAL :: vmoveold    ! The domain translation speed of in restart data
  REAL :: uchange     ! Change in domain translation speed
  REAL :: vchange     ! Change in domain translation speed
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: tim
  INTEGER :: i, j, k, n
  INTEGER :: nxin, nyin, nzin,nzsoilin
  INTEGER :: istat, idummy
  REAL    :: datatim,rdummy
  INTEGER :: nsclin
  INTEGER :: p_qcin,p_qrin,p_qiin,p_qsin,p_qgin,p_qhin,                 &
             p_ncin,p_nrin,p_niin,p_nsin,p_ngin,p_nhin,                 &
                    p_zrin,p_ziin,p_zsin,p_zgin,p_zhin,                 &
             p_ccin
  INTEGER :: nqscalarin(nscalar)

  REAL    :: amin, amax
  LOGICAL :: fexist,cmprsed
  INTEGER :: lrstfn

  INTEGER :: var,varsize
  INTEGER :: nq,nqin

  INTEGER :: qcbcrd, qrbcrd, qibcrd, qsbcrd, qhbcrd, qgbcrd,            &
             ncbcrd, nrbcrd, nibcrd, nsbcrd, nhbcrd, ngbcrd,            &
                     zrbcrd, zibcrd, zsbcrd, zhbcrd, zgbcrd,            &
             ccbcrd

  CHARACTER (LEN=256) :: savename

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  IF(nts == 3 ) THEN
    tpast    = 1
    tpresent = 2
    tfuture  = 3
  ELSE IF(nts == 2 ) THEN
    tpast    = 1
    tpresent = 2
    tfuture  = 2
  ELSE
    tpast    = 1
    tpresent = 1
    tfuture  = 1
  END IF

!   Added a wrapper for scheduling the number of open files...
!   The wrapper will conclude prior to the maxmin checking....
!   note call to jacob is moved outside the file open control loop.
!   due to message passing code in the call to jacob.

!  blocking inserted for ordering i/o for message passing
  DO n=0,nprocs-1,max_fopen
    IF (myproc >= n .AND. myproc <= n+max_fopen-1) THEN

      CALL getunit( rstiunt )

      lrstfn = 256
      CALL strlnth( rstinf, lrstfn)

      IF (mp_opt > 0) THEN
        savename(1:256) = rstinf(1:256)
        CALL gtsplitfn(savename,1,1,loc_x,loc_y,1,1,0,0,1,lvldbg,rstinf,istat)
        lrstfn = LEN_TRIM(rstinf)
      END IF

      cmprsed = .false.

      INQUIRE(FILE=rstinf(1:lrstfn), EXIST = fexist )
      IF( fexist ) GO TO 100

      INQUIRE(FILE=rstinf(1:lrstfn)//'.Z', EXIST = fexist )
      IF( fexist ) THEN
        cmprsed = .true.
        CALL uncmprs( rstinf(1:lrstfn)//'.Z' )
        GO TO 100
      END IF

      INQUIRE(FILE=rstinf(1:lrstfn)//'.gz', EXIST = fexist )
      IF( fexist ) THEN
        cmprsed = .true.
        CALL uncmprs( rstinf(1:lrstfn)//'.gz' )
        GO TO 100
      END IF

      CALL wrtcomment('File '//rstinf(1:lrstfn)//                       &
          ' or its compressed version not found.',1)
      CALL arpsstop('arpsstop called from RSTIN compressed file not '// &
           'found',1)

      100   CONTINUE

      OPEN(UNIT=rstiunt,FILE=trim(rstinf(1:lrstfn)),                    &
          FORM='unformatted',STATUS='old',IOSTAT=istat)

      IF (mp_opt > 0) THEN
        rstinf(1:256) = savename(1:256)
        lrstfn = lrstfn - 5
      END IF

      IF( istat /= 0) THEN

        WRITE(6,'(/1x,a,i2/)')                                          &
            'Error occured when opening restart input file '//          &
            rstinf(1:lrstfn)//                                          &
            ' using FORTRAN unit ',rstiunt
      CALL arpsstop('arpsstop called from RSTIN restart file not'//     &
           'found',1)

      END IF

      IF(myproc == 0) WRITE(6,'(/1x,a,/1x,a,i2/)')                      &
          'This is a restart run. Input was read from restart file ',   &
          rstinf(1:lrstfn)//' using fortran unit ',rstiunt
!
!
!-----------------------------------------------------------------------
!
!  Read in the restart data:
!
!-----------------------------------------------------------------------
!
      READ(rstiunt,ERR=999) datatim

      tstart = datatim
      IF(myproc == 0)   &
      WRITE(6,'(a,f8.1)') ' Restart data is at time ',  datatim

      READ(rstiunt,ERR=999) nxin,nyin,nzin, nzsoilin

      IF((nx /= nxin).OR.(ny /= nyin).OR.(nz /= nzin) .OR. (nzsoil /= nzsoilin)) THEN

        WRITE(6,'(a,/a,i5,a,i5,a,i5,/a,i5,a,i5,a,i5)')                    &
            ' Array dimension(s) in the restart data inconsistent with ', &
            ' model definitions, dimensions in input data were nx=',nxin, &
            ', ny=',nyin,', nz=',nzin,' the model definitions were nx=',  &
            nx,' ny= ', ny, ' nz= ',nz
        WRITE(6,'(a)') ' Job stopped in subroutine rstin.'
        CALL arpsstop('arpsstop called from RSTIN dimensions '//          &
             'inconsistent',1)

      END IF

      READ(rstiunt) basrstin,grdrstin,icerstin,sfcrstin, prcrsin,       &
                    rcumin,  exbcin,  mapfin,  radrstin, nstyp,         &
                    kfrsin,  rayklow, bmjsin,  nsclin,   idummy,        &
                    idummy,  idummy,  idummy,  idummy,   idummy,        &
                    idummy   ! a bug since arps5.0.0IHOP_0

      IF (nsclin > 0) THEN   ! new version
      READ(rstiunt) p_qcin,p_qrin,p_qiin,p_qsin,p_qgin,p_qhin,          &
                    p_ncin,p_nrin,p_niin,p_nsin,p_ngin,p_nhin,          &
                    p_zrin,p_ziin,p_zsin,p_zgin,p_zhin,p_ccin,          &
                    idummy,idummy,idummy,idummy,idummy,idummy,          &
                    idummy,idummy,idummy,idummy,idummy,idummy
      ELSE                   ! old version
        p_qgin=0;
        p_ncin=0; p_nrin=0; p_niin=0; p_nsin=0; p_ngin=0; p_nhin=0
                  p_zrin=0; p_ziin=0; p_zsin=0; p_zgin=0; p_zhin=0
        p_ccin=0

        IF (icerstin /= 0) THEN
          nsclin = 5
          p_qcin = 1; p_qrin = 2; p_qiin = 3; p_qsin = 4; p_qhin = 5
        ELSE
          nsclin = 2
          p_qcin = 1; p_qrin = 2; p_qiin = 0; p_qsin = 0; p_qhin = 0
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

      IF (P_ZR > 0) nqscalarin(P_ZR) = p_zrin
      IF (P_ZI > 0) nqscalarin(P_ZI) = p_ziin
      IF (P_ZS > 0) nqscalarin(P_ZS) = p_zsin
      IF (P_ZG > 0) nqscalarin(P_ZG) = p_zgin
      IF (P_ZH > 0) nqscalarin(P_ZH) = p_zhin

      IF (P_CC > 0) nqscalarin(P_CC) = p_ccin

      READ(rstiunt) dx,          dy,     dz,umoveold, vmoveold,         &
                    xgrdorg,ygrdorg,trulat1, trulat2,  trulon,          &
                    sclfct, latitud, ctrlat,  ctrlon,  ntcloud, n0rain, &
                    n0snow,n0grpl,n0hail,rhoice,rhosnow,rhogrpl,rhohail,&
                    alpharain,alphaice,alphasnow,alphagrpl,alphahail,rdummy

      IF( grdrstin == 1) THEN

        READ(rstiunt) x
        READ(rstiunt) y
        READ(rstiunt) z
        READ(rstiunt) zp
        READ(rstiunt) zpsoil

        DO i=1,nx
          DO j=1,ny
            hterain(i,j) = zp(i,j,2)
          END DO
        END DO

! let us first set it to 1.0, it may need to be changed later.

        j3soil = 1.0

      END IF

      IF( basrstin == 1) THEN

        READ(rstiunt) ubar
        READ(rstiunt) vbar
        READ(rstiunt) ptbar
        READ(rstiunt) pbar
        READ(rstiunt) rhostr
        READ(rstiunt) qvbar

        IF(myproc == 0)  WRITE(6,'(/1x,a/,1x,a/)')                      &
            'Base state arrays are read in from restart data set',      &
            'the base state set in INIBASE is superceded.'

      END IF


      tim = tpast

      READ(rstiunt,ERR=999) tem1
      CALL cpyary3d(nx,ny,nz,tem1,u    (1,1,1,tim))

      READ(rstiunt,ERR=999) tem1
      CALL cpyary3d(nx,ny,nz,tem1,v    (1,1,1,tim))

      READ(rstiunt,ERR=999) tem1
      CALL cpyary3d(nx,ny,nz,tem1,w    (1,1,1,tim))

      READ(rstiunt,ERR=999) tem1
      CALL cpyary3d(nx,ny,nz,tem1,ptprt(1,1,1,tim))

      READ(rstiunt,ERR=999) tem1
      CALL cpyary3d(nx,ny,nz,tem1,pprt (1,1,1,tim))

      READ(rstiunt,ERR=999) tem1
      CALL cpyary3d(nx,ny,nz,tem1,qv   (1,1,1,tim))

      DO nqin = 1,nsclin
        READ(rstiunt,ERR=999) tem1
        DO nq = 1,nscalar
          IF (nqin == nqscalarin(nq)) THEN
            CALL cpyary3d(nx,ny,nz,tem1,qscalar(1,1,1,tim,nq))
            EXIT
          END IF
        END DO
      END DO

      READ(rstiunt,ERR=999) tem1
      CALL cpyary3d(nx,ny,nz,tem1,tke   (1,1,1,tim))

      tim = tpresent

      READ(rstiunt,ERR=999) tem1
      CALL cpyary3d(nx,ny,nz,tem1,u    (1,1,1,tim))

      READ(rstiunt,ERR=999) tem1
      CALL cpyary3d(nx,ny,nz,tem1,v    (1,1,1,tim))

      READ(rstiunt,ERR=999) tem1
      CALL cpyary3d(nx,ny,nz,tem1,w    (1,1,1,tim))

      READ(rstiunt,ERR=999) tem1
      CALL cpyary3d(nx,ny,nz,tem1,ptprt(1,1,1,tim))

      READ(rstiunt,ERR=999) tem1
      CALL cpyary3d(nx,ny,nz,tem1,pprt (1,1,1,tim))

      READ(rstiunt,ERR=999) tem1
      CALL cpyary3d(nx,ny,nz,tem1,qv   (1,1,1,tim))

      DO nqin = 1,nsclin
        READ(rstiunt,ERR=999) tem1
        DO nq = 1,nscalar
          IF (nqin == nqscalarin(nq)) THEN
            CALL cpyary3d(nx,ny,nz,tem1,qscalar(1,1,1,tim,nq))
            EXIT
          END IF
        END DO
      END DO

      READ(rstiunt,ERR=999) tem1
      CALL cpyary3d(nx,ny,nz,tem1,tke  (1,1,1,tim))

      READ(rstiunt,ERR=999) udteb
      READ(rstiunt,ERR=999) udtwb

      READ(rstiunt,ERR=999) vdtnb
      READ(rstiunt,ERR=999) vdtsb

      READ(rstiunt,ERR=999) pdteb
      READ(rstiunt,ERR=999) pdtwb
      READ(rstiunt,ERR=999) pdtnb
      READ(rstiunt,ERR=999) pdtsb

!-----------------------------------------------------------------------
!
!TINA initialize cc bubble here if using restart file and ccin = 1
!This initializes tpast,present,future
!
!-----------------------------------------------------------------------
      IF (P_CC > 0) THEN
        IF (ccin == 1 .AND. cpoint < 0) THEN
          CALL ccinit(nx,ny,nz,3,x,y,z,zp,qscalar(:,:,:,:,P_CC),tem1)
          IF (myproc == 0) print *,'Initializing cc as bubble'
        END IF
      END IF
!TINA/michi

!-----------------------------------------------------------------------
!
!  Set the future values of variables to their current values.
!  This is done primarily for safety reasons since the arrays at
!  tfuture will be overwritten by the new values during the
!  time integration.
!
!-----------------------------------------------------------------------
!
      CALL cpyary3d(nx,ny,nz,u (1,1,1,tpresent) , u (1,1,1,tfuture))
      CALL cpyary3d(nx,ny,nz,v (1,1,1,tpresent) , v (1,1,1,tfuture))
      CALL cpyary3d(nx,ny,nz,w (1,1,1,tpresent) , w (1,1,1,tfuture))
      CALL cpyary3d(nx,ny,nz,ptprt(1,1,1,tpresent),                     &
                                               ptprt(1,1,1,tfuture))
      CALL cpyary3d(nx,ny,nz,pprt (1,1,1,tpresent),                     &
                                               pprt (1,1,1,tfuture))
      CALL cpyary3d(nx,ny,nz,qv(1,1,1,tpresent) , qv(1,1,1,tfuture))
      DO nq = 1,nscalar
        CALL cpyary3d(nx,ny,nz,qscalar(1,1,1,tpresent,nq) ,             &
                                          qscalar(1,1,1,tfuture,nq))
      END DO
      CALL cpyary3d(nx,ny,nz,tke(1,1,1,tpresent),tke(1,1,1,tfuture))

      IF ( sfcrstin /= 0 ) THEN

        IF(myproc == 0) PRINT *,'read in sfc/soil variables:'

        READ(rstiunt,ERR=999) soiltyp
        READ(rstiunt,ERR=999) stypfrct

        READ(rstiunt,ERR=999) vegtyp
        READ(rstiunt,ERR=999) lai
        READ(rstiunt,ERR=999) roufns
        READ(rstiunt,ERR=999) veg

        READ(rstiunt,ERR=999) qvsfc
        READ(rstiunt,ERR=999) tsoil
        READ(rstiunt,ERR=999) qsoil
        READ(rstiunt,ERR=999) wetcanp
        READ(rstiunt,ERR=999) snowdpth

      END IF

      IF ( prcrsin /= 0 ) THEN
        READ(rstiunt,ERR=999) raing
        READ(rstiunt,ERR=999) rainc
        READ(rstiunt,ERR=999) prcrate
      END IF

      IF ( rcumin /= 0 ) THEN
        READ(rstiunt,ERR=999) ptcumsrc
        READ(rstiunt,ERR=999) qcumsrc
      END IF

      IF ( exbcin /= 0 ) THEN
        IF ( lbcopt == 2 ) THEN
          READ(rstiunt,ERR=999) abstfcst0, abstfcst,                    &
                   ubcrd,vbcrd,wbcrd,ptbcrd,prbcrd,                     &
                   qvbcrd,qcbcrd,qrbcrd,qibcrd,qsbcrd,qhbcrd,qgbcrd,    &
                   ncbcrd,nrbcrd,nibcrd,nsbcrd,ngbcrd,nhbcrd,           &
                   zrbcrd,zibcrd,zsbcrd,zgbcrd,zhbcrd,ccbcrd

          qscalarbcrd(:) = 0
          IF (P_QC >0) qscalarbcrd(P_QC) = qcbcrd
          IF (P_QR >0) qscalarbcrd(P_QR) = qrbcrd
          IF (P_QI >0) qscalarbcrd(P_QI) = qibcrd
          IF (P_QS >0) qscalarbcrd(P_QS) = qsbcrd
          IF (P_QH >0) qscalarbcrd(P_QH) = qhbcrd
          IF (P_QG >0) qscalarbcrd(P_QG) = qgbcrd

          IF (P_NC >0) qscalarbcrd(P_NC) = ncbcrd
          IF (P_NR >0) qscalarbcrd(P_NR) = nrbcrd
          IF (P_NI >0) qscalarbcrd(P_NI) = nibcrd
          IF (P_NS >0) qscalarbcrd(P_NS) = nsbcrd
          IF (P_NH >0) qscalarbcrd(P_NH) = nhbcrd
          IF (P_NG >0) qscalarbcrd(P_NG) = ngbcrd

          IF (P_ZR >0) qscalarbcrd(P_ZR) = zrbcrd
          IF (P_ZI >0) qscalarbcrd(P_ZI) = zibcrd
          IF (P_ZS >0) qscalarbcrd(P_ZS) = zsbcrd
          IF (P_ZH >0) qscalarbcrd(P_ZH) = zhbcrd
          IF (P_ZG >0) qscalarbcrd(P_ZG) = zgbcrd

          IF (P_CC >0) qscalarbcrd(P_CC) = ccbcrd

!   Discretized for the auto split and auto join of the message passing mode

          varsize = nx*ny*nz
          DO var = 1, exbcbufsz, varsize
            READ(rstiunt,ERR=999) (exbcbuf(i),i=var,var+varsize-1)
          END DO

        ELSE
          WRITE(6,'(a/a/a/a)')                                          &
              'WARNING: The restart file contains EXBC arrays, while',  &
              '         the this run does not have EXBC option.',       &
              '         Therefore, the results from restart run may be', &
              '         alterred. The program will continue.'

          READ(rstiunt,ERR=999)

          READ(rstiunt,ERR=999)
          READ(rstiunt,ERR=999)
          READ(rstiunt,ERR=999)
          READ(rstiunt,ERR=999)
          READ(rstiunt,ERR=999)
          READ(rstiunt,ERR=999)
          READ(rstiunt,ERR=999)
          READ(rstiunt,ERR=999)
          READ(rstiunt,ERR=999)
          READ(rstiunt,ERR=999)
          READ(rstiunt,ERR=999)
        END IF
      END IF

      IF ( mapfin == 1 ) THEN
        READ(rstiunt,ERR=999) mapfct
      END IF

      IF ( radrstin == 1 ) THEN
        READ(rstiunt,ERR=999) radfrc
        READ(rstiunt,ERR=999) radsw
        READ(rstiunt,ERR=999) rnflx
        READ(rstiunt,ERR=999) radswnet
        READ(rstiunt,ERR=999) radlwin

      END IF

      IF ( kfrsin /= 0) THEN
        READ(rstiunt,ERR=999) w0avg
        READ(rstiunt,ERR=999) nca
        READ(rstiunt,ERR=999) kfraincv
      END IF

      IF ( bmjsin /= 0) THEN
        READ(rstiunt,ERR=999) cldefi
        READ(rstiunt,ERR=999) xland
        READ(rstiunt,ERR=999) bmjraincv
      END IF

      CLOSE (UNIT=rstiunt)
      CALL retunit( rstiunt )
!
!-----------------------------------------------------------------------
!
!  Reset the model u and v velocity values using the new
!  domain translation speed.
!
!-----------------------------------------------------------------------
!

      IF( nint(umove) == 999 .OR. nint(vmove) == 999 ) THEN

        umove = umoveold
        vmove = vmoveold

      ELSE IF (umoveold /= umove .OR. vmoveold /= vmove ) THEN

        WRITE(6,'(3(/1x,a)/)')                                          &
            'ATTENTION: UMOVE or VMOVE in the input file were different ', &
            'from those in the restart file. Subroutine ADJUVMV is called',&
            'to adjust the time-dependent variables for option grdtrns!=0.'

        IF ( grdtrns /= 0 ) THEN
          uchange = umove - umoveold
          vchange = vmove - vmoveold

          CALL adjuvmv(nx,ny,nz,                                        &
                       ubar,vbar,u,v,w,ptprt,pprt,qv,qvbar,qscalar,     &
                       uchange, vchange, tem1, tem2)

        END IF
      END IF

    END IF  ! end of FOPEN wrapper for file read/write...
    IF (mp_opt > 0) CALL mpbarrier
  END DO

  CALL jacob(nx,ny,nz,x,y,z,zp,j1,j2,j3,tem1)

  IF(myproc == 0) THEN
    WRITE(6,'(/1x,a/,1x,a/)')                                           &
      'Grid definition arrays are read in from initialization data',    &
      'those set in INIGRD are superceded.'

!
!-----------------------------------------------------------------------
!
!  Print out the domain-wide max/min of output variables.
!
!-----------------------------------------------------------------------

    WRITE(6,'(/1x,a/)')                                               &
      'Min. and max. of the data arrays read in from restart data:'

  END IF  ! myproc == 0

  CALL a3dmax0(x,1,nx,1,nx,1,1,1,1, 1,1,1,1, amax,amin)
  IF(myproc == 0) &
  WRITE(6,'(/1x,2(a,e13.6))') 'xmin    = ', amin,',  xmax    =',amax

  CALL a3dmax0(y,1,ny,1,ny,1,1,1,1, 1,1,1,1, amax,amin)
  IF(myproc == 0) &
  WRITE(6,'(1x,2(a,e13.6))') 'ymin    = ', amin,',  ymax    =',amax

  CALL a3dmax0(z,1,nz,1,nz,1,1,1,1, 1,1,1,1, amax,amin)
  IF(myproc == 0) &
  WRITE(6,'(1x,2(a,e13.6))') 'zmin    = ', amin,',  zmax    =',amax

  CALL a3dmax0(zp,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz, amax,amin)
  IF(myproc == 0) &
  WRITE(6,'(1x,2(a,e13.6))') 'zpmin   = ', amin,',  zpmax   =',amax

  CALL a3dmax0(zpsoil,1,nx,1,nx-1,1,ny,1,ny-1,1,nzsoil,1,nzsoil, amax,amin)
  IF(myproc == 0) &
  WRITE(6,'(1x,2(a,e13.6))') 'zpsoilmin   = ', amin,',  zpsoilmax   =',amax

  CALL a3dmax0(hterain,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1, amax,amin)
  IF(myproc == 0) &
  WRITE(6,'(1x,2(a,e13.6))') 'hmin    = ', amin,',  hmax    =',amax

  CALL a3dmax0(ubar,1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1, amax,amin)
  IF(myproc == 0) &
  WRITE(6,'(1x,2(a,e13.6))') 'ubarmin = ', amin,',  ubarmax =',amax

  CALL a3dmax0(vbar,1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1, amax,amin)
  IF(myproc == 0) &
  WRITE(6,'(1x,2(a,e13.6))') 'vbarmin = ', amin,',  vbarmax =',amax

  CALL a3dmax0(ptbar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
  IF(myproc == 0) &
  WRITE(6,'(1x,2(a,e13.6))') 'ptbarmin= ', amin,',  ptbarmax=',amax

  CALL a3dmax0(pbar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1, amax,amin)
  IF(myproc == 0) &
  WRITE(6,'(1x,2(a,e13.6))') 'pbarmin = ', amin,',  pbarmax =',amax

  CALL a3dmax0(rhostr,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
  IF(myproc == 0) &
  WRITE(6,'(1x,2(a,e13.6))') 'rhostrmin=', amin,', rhostrmax=',amax

  CALL a3dmax0(qvbar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
  IF(myproc == 0) &
  WRITE(6,'(1x,2(a,e13.6))') 'qvbarmin= ', amin,',  qvbarmax=',amax


  DO i = 1,2

    IF( i == 1) THEN
      IF(myproc == 0) WRITE(6,'(/1x,a/)') 'Min/max of fields at tpresent:'
      tim = tpresent
    ELSE
      IF(myproc == 0) WRITE(6,'(/1x,a/)') 'Min/max of fields at tpast:'
      tim = tpast
    END IF

    CALL a3dmax0(u(1,1,1,tim),1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1,        &
                 amax,amin)
    IF(myproc == 0) &
    WRITE(6,'(1x,2(a,e13.6))') 'umin    = ', amin,',  umax    =',       &
          amax

    CALL a3dmax0(v(1,1,1,tim),1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1,        &
                 amax,amin)
    IF(myproc == 0) &
    WRITE(6,'(1x,2(a,e13.6))') 'vmin    = ', amin,',  vmax    =',       &
          amax

    CALL a3dmax0(w(1,1,1,tim),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz,        &
                 amax,amin)
    IF(myproc == 0) &
    WRITE(6,'(1x,2(a,e13.6))') 'wmin    = ', amin,',  wmax    =',       &
          amax

    CALL a3dmax0(ptprt(1,1,1,tim),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,       &
                 nz-1,amax,amin)
    IF(myproc == 0) &
    WRITE(6,'(1x,2(a,e13.6))') 'ptprtmin= ', amin,',  ptprtmax=',       &
          amax

    CALL a3dmax0(pprt(1,1,1,tim),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,        &
                 nz-1,amax,amin)
    IF(myproc == 0) &
    WRITE(6,'(1x,2(a,e13.6))') 'pprtmin = ', amin,',  pprtmax =',       &
          amax

    CALL a3dmax0(qv(1,1,1,tim),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,     &
                 amax,amin)
    IF(myproc == 0) &
    WRITE(6,'(1x,2(a,e13.6))') 'qvmin   = ', amin,',  qvmax   =',       &
          amax

    DO nq = 1,nscalar
      CALL a3dmax0(qscalar(1,1,1,tim,nq),1,nx,1,nx-1,1,ny,1,ny-1,       &
                   1,nz,1,nz-1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,a,e13.6,2a,e13.6)')                  &
                   TRIM(qnames(nq))//'min   = ', amin,', ',             &
                   TRIM(qnames(nq))//'max   = ', amax
    END DO

  END DO

  IF(myproc == 0) &
  WRITE(6,'(/1x,a/)')                                                   &
        'Min/max of fields for other one time level arrays:'

  CALL a3dmax0(tsoil(1,1,1,0),1,nx,1,nx-1,1,ny,1,ny-1,1,nzsoil,1,nzsoil, &
               amax,amin)
  IF(myproc == 0) &
  WRITE(6,'(1x,2(a,e13.6))') 'tsoilmin= ', amin,', tsoilmax =',amax

  CALL a3dmax0(qsoil(1,1,1,0),1,nx,1,nx-1,1,ny,1,ny-1,1,nzsoil,1,nzsoil, &
               amax,amin)
  IF(myproc == 0) &
  WRITE(6,'(1x,2(a,e13.6))') 'qsoilmin = ', amin,', qsoilmax  =',amax

  CALL a3dmax0(wetcanp(1,1,0),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,          &
               amax,amin)
  IF(myproc == 0) &
  WRITE(6,'(1x,2(a,e13.6))') 'wetcmin = ', amin,', wetcmax  =',amax

  CALL a3dmax0(raing,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1, amax,amin)
  IF(myproc == 0) &
  WRITE(6,'(1x,2(a,e13.6))') 'raingmin= ', amin,', raingmax =',amax

  CALL a3dmax0(rainc,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1, amax,amin)
  IF(myproc == 0) &
  WRITE(6,'(1x,2(a,e13.6))') 'raincmin= ', amin,', raincmax =',amax

  CALL a3dmax0(ptcumsrc,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,            &
               amax,amin)
  IF(myproc == 0) &
  WRITE(6,'(1x,2(a,e13.6))') 'ptcummin= ', amin,',  ptcummax=',amax

  CALL a3dmax0(qcumsrc(1,1,1,1),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,    &
               amax,amin)
  IF(myproc == 0) &
  WRITE(6,'(1x,2(a,e13.6))') 'qvcummin= ', amin,',  qvcummax=',amax
  CALL a3dmax0(qcumsrc(1,1,1,2),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,    &
               amax,amin)
  IF(myproc == 0) &
  WRITE(6,'(1x,2(a,e13.6))') 'qccummin= ', amin,',  qccummax=',amax
  CALL a3dmax0(qcumsrc(1,1,1,3),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,    &
               amax,amin)
  IF(myproc == 0) &
  WRITE(6,'(1x,2(a,e13.6))') 'qrcummin= ', amin,',  qrcummax=',amax
  CALL a3dmax0(qcumsrc(1,1,1,4),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,    &
               amax,amin)
  IF(myproc == 0) &
  WRITE(6,'(1x,2(a,e13.6))') 'qicummin= ', amin,',  qicummax=',amax
  CALL a3dmax0(qcumsrc(1,1,1,5),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,    &
               amax,amin)
  IF(myproc == 0) &
  WRITE(6,'(1x,2(a,e13.6))') 'qscummin= ', amin,',  qscummax=',amax
!
!-----------------------------------------------------------------------
!
!  Compress the restart file if it was originally compressed.
!
!-----------------------------------------------------------------------
!
  IF( cmprsed .AND. filcmprs == 1 ) THEN
    CALL cmprs( rstinf(1:lrstfn) )
  END IF

  RETURN

  999   CONTINUE
  WRITE(6,'(a)') ' Error reading restart data '//rstinf
  WRITE(6,'(a,i3,a)') ' Fortran channel ',rstiunt,' was used.'
  WRITE(6,'(a)') ' Job stopped in subroutine rstin!'
  CLOSE (UNIT=rstiunt)

  CALL arpsstop('arpsstop called from RSTIN error reading restart'//   &
                 'file ',1)

END SUBROUTINE rstin

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RSTINSPLIT                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE rstinsplit(nx,ny,nz,nzsoil,nts,nstyps,exbcbufsz,             &
           u,v,w,ptprt,pprt,qv,qscalar,tke,                             &
           udteb, udtwb, vdtnb, vdtsb,                                  &
           pdteb ,pdtwb ,pdtnb ,pdtsb,                                  &
           ubar,vbar,ptbar,pbar,rhostr,qvbar,                           &
           x,y,z,zp,zpsoil,hterain,mapfct,j1,j2,j3,j3soil,              &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,qvsfc,                          &
           ptcumsrc,qcumsrc,w0avg,nca,kfraincv,                         &
           cldefi,xland,bmjraincv,                                      &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           raing,rainc,prcrate, exbcbuf, tem1, tem2)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in data from a restart file to initialize u,v,w,prprt,pprt,
!  qv,qc,qr,qi,qs,qh and tke at time tpast and tpresent, the base state
!  variables ubar,vbar,ptbar,pbar,rhostr,qvbar, and the time tendencies
!    of variables at the lateral boundaries.
!
!  Fields at tfuture are set to the values at tpresent.
!
!  This subroutine also sets the value of tstart.
!
!  RSTINSPLIT read in data and split for message passing mode.
!
!  NOTE: This suboutine should be consistent with the normal one, RSTIN.
!        Any changes here should also be copied to subroutine RSTIN above.
!
!        The parameter list is the same as that of RSTIN. This will
!        make it easier to call RSTIN and RSTINSPLIT at the same place
!        of the calling subroutine. It will also be easy to combine
!        these two subroutines into one later if necessary.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  2/25/2003.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of grid points in the soil
!    nts      Number of time levels to be initialized.
!
!  OUTPUT:
!
!    u        x component of velocity at times tpast and tpresent (m/s)
!    v        y component of velocity at times tpast and tpresent (m/s)
!    w        Vertical component of Cartesian velocity at times
!             tpast and tpresent (m/s)
!    ptprt    Perturbation potential temperature at times tpast and
!             tpresent (K)
!    pprt     Perturbation pressure at times tpast and tpresent (Pascal)
!
!    qv       Water vapor specific humidity at times tpast and tpresent (kg/kg)
!    qc       Cloud water mixing ratio at times tpast and tpresent (kg/kg)
!    qr       Rainwater mixing ratio at times tpast and tpresent (kg/kg)
!    qi       Cloud ice mixing ratio at times tpast and tpresent (kg/kg)
!    qs       Snow mixing ratio at times tpast and tpresent (kg/kg)
!    qh       Hail mixing ratio at times tpast and tpresent (kg/kg)
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!
!    udteb    Time tendency of u field at east boundary (m/s**2)
!    udtwb    Time tendency of u field at west boundary (m/s**2)
!
!    vdtnb    Time tendency of v field at north boundary (m/s**2)
!    vdtsb    Time tendency of v field at south boundary (m/s**2)
!
!    pdteb    Time tendency of pprt field at east boundary (PASCAL/s)
!    pdtwb    Time tendency of pprt field at west boundary (PASCAL/s)
!    pdtnb    Time tendency of pprt field at north boundary (PASCAL/s)
!    pdtsb    Time tendency of pprt field at south boundary (PASCAL/s)
!
!    ubar     Base state zonal velocity component (m/s)
!    vbar     Base state meridional velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhostr   Base state density (kg/m**3) times j3.
!    qvbar    Base state water vapor specific humidity (kg/kg)
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!    zpsoil   Vertical coordinate of grid points in the soil (m)
!    hterain  Terrain height (m)
!
!    mapfct   Map factors at scalar, u and v points
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    j3soil   Coordinate transformation Jacobian  d(zpsoil)/dz
!
!    soiltyp  Soil type
!    vegtyp   Vegetation type
!    lai      Leaf Area Index
!    roufns   Surface roughness
!    veg      Vegetation fraction
!
!    qvsfc    Effective S.H. at sfc.
!    tsoil    Soil temperature (K)
!    qsoil    Soil moisture (m**3/m**3)
!    wetcanp  Canopy water amount
!    ptcumsrc Source term in pt-equation due to cumulus parameterization
!    qcumsrc Source term in water equations due to cumulus parameterization
!    kfraincv   K-F convective rainfall (cm)
!    nca      K-F counter for CAPE release
!    cldefi   BMJ cloud efficiency
!    xland    BMJ land/sea mask
!    bmjraincv   BMJ convective rainfall (cm)
!
!    radfrc   Radiation forcing (K)
!    radsw    Solar radiation reaching the surface
!    rnflx    Net absorbed radiation by the surface
!    radswnet Net shortwave radiation
!    radlwin  Incoming longwave radiation
!
!    raing    Grid scale rainfall
!    rainc    Convective rainfall
!
!    tstart   The time when the time integration starts, which is set to
!             the time of the restart data
!
!    tem1     Temporary work array
!    tem2     Temporary work array
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
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'bndry.inc'
  INCLUDE 'exbc.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'            ! Message passing parameters.

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nts          ! Number of time levels to be initialized.
  INTEGER :: tpast        ! Index of time level for the past time.
  INTEGER :: tpresent     ! Index of time level for the present time.
  INTEGER :: tfuture      ! Index of time level for the future time.

  INTEGER :: nx,ny,nz     ! Number of grid points in 3 directions
  INTEGER :: nzsoil       ! Number of grid points in 3 directions

  REAL :: u     (nx,ny,nz,nts) ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nts) ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nts) ! Total w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz,nts) ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nts) ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz,nts) ! Water vapor specific humidity (kg/kg)
  REAL :: qscalar(nx,ny,nz,nts,nscalar)
  REAL :: tke   (nx,ny,nz,nts) ! Turbulent Kinetic Energy ((m/s)**2)

  REAL :: udteb (ny,nz)        ! T-tendency of u at e-boundary (m/s**2)
  REAL :: udtwb (ny,nz)        ! T-tendency of u at w-boundary (m/s**2)

  REAL :: vdtnb (nx,nz)        ! T-tendency of v at n-boundary (m/s**2)
  REAL :: vdtsb (nx,nz)        ! T-tendency of v at s-boundary (m/s**2)

  REAL :: pdteb (ny,nz)        ! T-tendency of pprt at e-boundary (PASCAL/s)
  REAL :: pdtwb (ny,nz)        ! T-tendency of pprt at w-boundary (PASCAL/s)
  REAL :: pdtnb (nx,nz)        ! T-tendency of pprt at n-boundary (PASCAL/s)
  REAL :: pdtsb (nx,nz)        ! T-tendency of pprt at s-boundary (PASCAL/s)

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL :: rhostr(nx,ny,nz)     ! Base state air density (kg/m**3) time j3.
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                               ! (kg/kg)

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
  REAL :: zpsoil(nx,ny,nzsoil) ! The physical height coordinate defined at
                               ! w-point of the soil.

  REAL :: hterain(nx,ny)       ! Terrain height (m).

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian -d(zp)/dx.
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian -d(zp)/dy.
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian  d(zp)/dz.
  REAL :: j3soil(nx,ny,nzsoil) ! Coordinate transformation Jacobian  d(zpsoil)/dz.

  INTEGER :: nstyps                   ! Number of soil types
  INTEGER :: soiltyp (nx,ny,nstyps)   ! Soil type
  REAL    :: stypfrct(nx,ny,nstyps)
  INTEGER :: vegtyp(nx,ny)            ! Vegetation type
  REAL    :: lai    (nx,ny)           ! Leaf Area Index
  REAL    :: roufns (nx,ny)           ! Surface roughness
  REAL    :: veg    (nx,ny)           ! Vegetation fraction

  REAL :: qvsfc  (nx,ny,       0:nstyps) ! Effective S. H. at the surface (kg/kg)
  REAL :: tsoil  (nx,ny,nzsoil,0:nstyps) ! Soil temperature(K)
  REAL :: qsoil  (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny,       0:nstyps) ! Canopy water amount
  REAL :: snowdpth(nx,ny)         ! Snow depth (m)
  REAL :: ptcumsrc(nx,ny,nz)   ! Source term in pt-equation due
                               ! to cumulus parameterization
  REAL :: qcumsrc(nx,ny,nz,5)  ! Source term in water equations due
                               ! to cumulus parameterization:
                               ! qcumsrc(1,1,1,1) for qv equation
                               ! qcumsrc(1,1,1,2) for qc equation
                               ! qcumsrc(1,1,1,3) for qr equation
                               ! qcumsrc(1,1,1,4) for qi equation
                               ! qcumsrc(1,1,1,5) for qs equation
  REAL :: w0avg(nx,ny,nz)      ! a closing running average vertical
                               ! velocity in 10min for K-F scheme
  REAL :: kfraincv(nx,ny)      ! K-F convective rainfall (cm)
  INTEGER :: nca(nx,ny)        ! K-F counter for CAPE release

  REAL,INTENT(OUT) :: cldefi(nx,ny)    ! BMJ cloud efficiency
  REAL,INTENT(OUT) :: xland(nx,ny)     ! BMJ land mask
  REAL,INTENT(OUT) :: bmjraincv(nx,ny) ! BMJ convective rainfall (cm)
                                       ! (1.0 = land, 2.0 = sea)

  REAL :: radfrc(nx,ny,nz)     ! Radiation forcing (K/s)
  REAL :: radsw(nx,ny)         ! Solar radiation reacing the surface
  REAL :: rnflx(nx,ny)         ! Net absorbed radiation by the surface
  REAL :: radswnet(nx,ny)      ! Net shortwave radiation
  REAL :: radlwin(nx,ny)       ! Incoming longwave radiation

  REAL :: raing(nx,ny)         ! Grid supersaturation rain
  REAL :: rainc(nx,ny)         ! Cumulus convective rain
  REAL :: prcrate(nx,ny,4)     ! precipitation rate (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precipitation rate
                               ! prcrate(1,1,2) = grid scale precip.  rate
                               ! prcrate(1,1,3) = cumulus precip.  rate
                               ! prcrate(1,1,4) = microphysics precip. rate

  INTEGER :: exbcbufsz         ! EXBC buffer size
  REAL :: exbcbuf( exbcbufsz ) ! EXBC buffer array

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array.

  INTEGER :: grdrstin ! Parameter indicating if the restart data contains
                      ! the grid variables.
  INTEGER :: basrstin ! Parameter indicating if the restart data contains
                      ! the base state variables.
  INTEGER :: icerstin ! Parameter indicating if the restart data contains
                      ! the ice variables.
  INTEGER :: sfcrstin ! Parameter indicating if the restart data contains
                      ! the surface variables.
  INTEGER :: prcrsin  ! Parameter indicating if the restart data contains
                      ! precipitation rate and rainfall
  INTEGER :: rcumin   ! Parameter indicating if the cumulus source terms
                      ! data are present.
  INTEGER :: exbcin   ! Parameter indicating if the external boundary
                      ! data are present.
  INTEGER :: mapfin   ! Parameter indicating if the map factor
                      ! data are present.
  INTEGER :: radrstin ! Parameter indicating if the radiation forcing
                      ! arrays are present.
  INTEGER :: kfrsin   ! Parameter indicating if k-f variable exists
  INTEGER :: bmjsin   ! Parameter indicating if BMJ variable exists

  REAL :: umoveold    ! The domain translation speed of in restart data
  REAL :: vmoveold    ! The domain translation speed of in restart data
  REAL :: uchange     ! Change in domain translation speed
  REAL :: vchange     ! Change in domain translation speed
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: tim
  INTEGER :: i, j, k, n
  INTEGER :: nxin, nyin, nzin,nzsoilin
  INTEGER :: istat, idummy
  REAL    :: datatim,rdummy
  INTEGER :: nsclin
  INTEGER :: p_qcin,p_qrin,p_qiin,p_qsin,p_qgin,p_qhin,                 &
             p_ncin,p_nrin,p_niin,p_nsin,p_ngin,p_nhin,                 &
                    p_zrin,p_ziin,p_zsin,p_zgin,p_zhin,                 &
             p_ccin
  INTEGER :: nqscalarin(nscalar)

  REAL    :: amin, amax
  LOGICAL :: fexist,cmprsed
  INTEGER :: lrstfn

  INTEGER :: nxlg, nylg, n3rd, var
  INTEGER :: nq, nqin
  INTEGER :: qcbcrd, qrbcrd, qibcrd, qsbcrd, qgbcrd, qhbcrd,            &
             ncbcrd, nrbcrd, nibcrd, nsbcrd, ngbcrd, nhbcrd,            &
                     zrbcrd, zibcrd, zsbcrd, zgbcrd, zhbcrd,            &
             ccbcrd

  REAL, ALLOCATABLE :: in1d(:)
  REAL, ALLOCATABLE :: in2d(:,:)
  REAL, ALLOCATABLE :: in3d(:,:,:)
  REAL, ALLOCATABLE :: in4d(:,:,:,:), in4dq(:,:,:,:)
  REAL, ALLOCATABLE :: in2dew(:,:),in2dns(:,:)

  INTEGER, ALLOCATABLE :: in2di(:,:)
  INTEGER, ALLOCATABLE :: in3di(:,:,:)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  nxlg = nproc_x*(nx-3) + 3
  nylg = nproc_y*(ny-3) + 3
  n3rd = MAX(nz,nzsoil,nstyps+1,8)

  ALLOCATE (in1d( MAX(nxlg,nylg,nz) ),stat=istat)
  CALL check_alloc_status(istat, "rstinsplit:in1d")

  ALLOCATE (in2d(nxlg,nylg),stat=istat)
  CALL check_alloc_status(istat, "rstinsplit:in2d")

  ALLOCATE (in2di(nxlg,nylg),stat=istat)
  CALL check_alloc_status(istat, "rstinsplit:in2di")

  ALLOCATE (in2dns(nxlg,nz),stat=istat)
  CALL check_alloc_status(istat, "rstinsplit:in2dns")

  ALLOCATE (in2dew(nylg,nz),stat=istat)
  CALL check_alloc_status(istat, "rstinsplit:in2dew")

  ALLOCATE (in3d( nxlg,nylg, n3rd ),stat=istat)
  CALL check_alloc_status(istat, "rstinsplit:in3d")

  ALLOCATE (in3di( nxlg,nylg, nstyps ),stat=istat)
  CALL check_alloc_status(istat, "rstinsplit:in3di")

  ALLOCATE (in4d(nxlg, nylg, nzsoil, nstyps+1),stat=istat)
  CALL check_alloc_status(istat, "rstinsplit:in4d")

  ALLOCATE (in4dq(nxlg, nylg, nz, 5),stat=istat)
  CALL check_alloc_status(istat, "rstinsplit:in4dq")

  IF(nts == 3 ) THEN
    tpast    = 1
    tpresent = 2
    tfuture  = 3
  ELSE IF(nts == 2 ) THEN
    tpast    = 1
    tpresent = 2
    tfuture  = 2
  ELSE
    tpast    = 1
    tpresent = 1
    tfuture  = 1
  END IF

!  only processor 0 needs to open the input file

  IF(myproc == 0) THEN

      CALL getunit( rstiunt )

      lrstfn = 256
      CALL strlnth( rstinf, lrstfn)

      cmprsed = .false.

      INQUIRE(FILE=rstinf(1:lrstfn), EXIST = fexist )
      IF( fexist ) GO TO 100

      INQUIRE(FILE=rstinf(1:lrstfn)//'.Z', EXIST = fexist )
      IF( fexist ) THEN
        cmprsed = .true.
        CALL uncmprs( rstinf(1:lrstfn)//'.Z' )
        GO TO 100
      END IF

      INQUIRE(FILE=rstinf(1:lrstfn)//'.gz', EXIST = fexist )
      IF( fexist ) THEN
        cmprsed = .true.
        CALL uncmprs( rstinf(1:lrstfn)//'.gz' )
        GO TO 100
      END IF

      CALL wrtcomment('File '//rstinf(1:lrstfn)//                       &
          ' or its compressed version not found.',1)
      CALL arpsstop('arpsstop called from RSTINSPLIT compressed file not '// &
           'found',1)

      100   CONTINUE

      OPEN(UNIT=rstiunt,FILE=trim(rstinf(1:lrstfn)),                    &
          FORM='unformatted',STATUS='old',IOSTAT=istat)

      IF( istat /= 0) THEN

        WRITE(6,'(/1x,a,i2/)')                                          &
            'Error occured when opening restart input file '//          &
            rstinf(1:lrstfn)//                                          &
            ' using FORTRAN unit ',rstiunt
        CALL arpsstop('arpsstop called from RSTINSPLIT restart file'//  &
           ' cannot be opened',1)

      END IF

      WRITE(6,'(/1x,a,/1x,a,i2/)')                                      &
          'This is a restart run. Input was read from restart file ',   &
          rstinf(1:lrstfn)//' using fortran unit ',rstiunt

!
!
!-----------------------------------------------------------------------
!
!  Read in the restart data:
!
!-----------------------------------------------------------------------
!
      READ(rstiunt,ERR=999) datatim

      tstart = datatim
      WRITE(6,'(a,f8.1)') ' Restart data is at time ',  datatim

      READ(rstiunt,ERR=999) nxin,nyin,nzin, nzsoilin

      IF((nxlg /= nxin) .OR. (nylg   /= nyin) .OR.             &
         (nz   /= nzin) .OR. (nzsoil /= nzsoilin)) THEN

        WRITE(6,'(a,/a,i5,a,i5,a,i5,/a,a,i5,a,i5,a,i5)')                  &
            ' Array dimension(s) in the restart data inconsistent with ', &
            ' model definitions, dimensions in input data were nx=',nxin, &
            ', ny=',nyin,', nz=',nzin,' the model definitions were ',     &
            'nxlg=',nxlg,' nylg= ', nylg, ' nz= ',nz
        WRITE(6,'(a)') ' Job stopped in subroutine rstinsplit.'
        CALL arpsstop('arpsstop called from RSTINSPLIT dimensions '//     &
             'inconsistent',1)

      END IF

      READ(rstiunt) basrstin,grdrstin,icerstin,sfcrstin,prcrsin,        &
                      rcumin,  exbcin,  mapfin,radrstin,  nstyp,        &
                      kfrsin, rayklow,  bmjsin,  nsclin, idummy,        &
                      idummy,  idummy,  idummy,  idummy, idummy,        &
                      idummy

      print*,'nsclin = ',nsclin

      IF (nsclin > 0) THEN   ! new version
      READ(rstiunt) p_qcin,p_qrin,p_qiin,p_qsin,p_qgin,p_qhin,          &
                    p_ncin,p_nrin,p_niin,p_nsin,p_ngin,p_nhin,          &
                    p_zrin,p_ziin,p_zsin,p_zgin,p_zhin,p_ccin,          &
                    idummy,idummy,idummy,idummy,idummy,idummy,          &
                    idummy,idummy,idummy,idummy,idummy,idummy
      ELSE                   ! old version
        p_qgin=0;
        p_ncin=0; p_nrin=0; p_niin=0; p_nsin=0; p_ngin=0; p_nhin=0
                  p_zrin=0; p_ziin=0; p_zsin=0; p_zgin=0; p_zhin=0
        p_ccin = 0

        IF (icerstin /= 0) THEN
          nsclin = 5
          p_qcin = 1; p_qrin = 2; p_qiin = 3; p_qsin = 4; p_qhin = 5
        ELSE
          nsclin = 2
          p_qcin = 1; p_qrin = 2; p_qiin = 0; p_qsin = 0; p_qhin = 0
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

      IF (P_ZR > 0) nqscalarin(P_ZR) = p_zrin
      IF (P_ZI > 0) nqscalarin(P_ZI) = p_ziin
      IF (P_ZS > 0) nqscalarin(P_ZS) = p_zsin
      IF (P_ZG > 0) nqscalarin(P_ZG) = p_zgin
      IF (P_ZH > 0) nqscalarin(P_ZH) = p_zhin

      IF (P_CC > 0) nqscalarin(P_CC) = p_ccin

      READ(rstiunt) dx,dy,dz,umoveold,vmoveold,                         &
                xgrdorg,ygrdorg,trulat1,trulat2,trulon,                 &
                sclfct,latitud,ctrlat,ctrlon,ntcloud,n0rain,            &
                n0snow,n0grpl,n0hail,rhoice,rhosnow,rhogrpl,rhohail,    &
                alpharain,alphaice,alphasnow,alphagrpl,alphahail,rdummy

  END IF   ! myproc == 0

  CALL mpupdater(tstart,1)

  CALL mpupdatei(basrstin,1)
  CALL mpupdatei(grdrstin,1)
  CALL mpupdatei(icerstin,1)
  CALL mpupdatei(sfcrstin,1)
  CALL mpupdatei(prcrsin,1)
  CALL mpupdatei(rcumin,1)
  CALL mpupdatei(exbcin,1)
  CALL mpupdatei(mapfin,1)
  CALL mpupdatei(radrstin,1)
  CALL mpupdatei(nstyp,1)
  CALL mpupdatei(kfrsin,1)
  CALL mpupdatei(rayklow,1)
  CALL mpupdatei(bmjsin,1)
  CALL mpupdatei(nsclin,1)
  CALL mpupdatei(nqscalarin,nscalar)

  CALL mpupdater(dx,1)
  CALL mpupdater(dy,1)
  CALL mpupdater(dz,1)
  CALL mpupdater(umoveold,1)
  CALL mpupdater(vmoveold,1)
  CALL mpupdater(xgrdorg,1)
  CALL mpupdater(ygrdorg,1)
  CALL mpupdater(trulat1,1)
  CALL mpupdater(trulat2,1)
  CALL mpupdater(trulon,1)
  CALL mpupdater(sclfct,1)
  CALL mpupdater(latitud,1)
  CALL mpupdater(ctrlat,1)
  CALL mpupdater(ctrlon,1)
  CALL mpupdater(ntcloud,1)
  CALL mpupdater(n0rain,1)
  CALL mpupdater(n0snow,1)
  CALL mpupdater(n0grpl,1)
  CALL mpupdater(n0hail,1)
  CALL mpupdater(rhoice,1)
  CALL mpupdater(rhosnow,1)
  CALL mpupdater(rhogrpl,1)
  CALL mpupdater(rhohail,1)
  CALL mpupdater(alpharain,1)
  CALL mpupdater(alphaice,1)
  CALL mpupdater(alphasnow,1)
  CALL mpupdater(alphagrpl,1)
  CALL mpupdater(alphahail,1)

  IF( grdrstin == 1) THEN

    IF(myproc == 0) READ(rstiunt,ERR=999) (in1d(i),i=1,nxlg)
    CALL mpisplit1dx(in1d,nx,x)
    IF(myproc == 0) READ(rstiunt,ERR=999) (in1d(j),j=1,nylg)
    CALL mpisplit1dy(in1d,ny,y)
    IF(myproc == 0) READ(rstiunt,ERR=999) z
    CALL mpupdater(z,nz)
    IF(myproc == 0) READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nz)
    CALL mpisplit3d(in3d,nx,ny,nz,zp)
    IF(myproc == 0) READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nzsoil)
    CALL mpisplit3d(in3d,nx,ny,nzsoil,zpsoil)

    DO i=1,nx
      DO j=1,ny
        hterain(i,j) = zp(i,j,2)
      END DO
    END DO

    j3soil = 1.0
    ! let us first set it to 1.0, it may need to be changed later.

  END IF

  IF( basrstin == 1) THEN

    IF(myproc == 0) READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nz)
    CALL mpisplit3d(in3d,nx,ny,nz,ubar)
    IF(myproc == 0) READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nz)
    CALL mpisplit3d(in3d,nx,ny,nz,vbar)
    IF(myproc == 0) READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nz)
    CALL mpisplit3d(in3d,nx,ny,nz,ptbar)
    IF(myproc == 0) READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nz)
    CALL mpisplit3d(in3d,nx,ny,nz,pbar)
    IF(myproc == 0) READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nz)
    CALL mpisplit3d(in3d,nx,ny,nz,rhostr)
    IF(myproc == 0) READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nz)
    CALL mpisplit3d(in3d,nx,ny,nz,qvbar)

    IF(myproc == 0) WRITE(6,'(/1x,a/,1x,a/)')                           &
            'Base state arrays are read in from restart data set',      &
            'the base state set in INIBASE is superceded.'

  END IF

  tim = tpast

  IF(myproc == 0) READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nz)
  CALL mpisplit3d(in3d,nx,ny,nz,u(:,:,:,tim))

  IF(myproc == 0) READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nz)
  CALL mpisplit3d(in3d,nx,ny,nz,v(:,:,:,tim))

  IF(myproc == 0) READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nz)
  CALL mpisplit3d(in3d,nx,ny,nz,w(:,:,:,tim))

  IF(myproc == 0) READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nz)
  CALL mpisplit3d(in3d,nx,ny,nz,ptprt(:,:,:,tim))

  IF(myproc == 0) READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nz)
  CALL mpisplit3d(in3d,nx,ny,nz,pprt(:,:,:,tim))

  IF(myproc == 0) READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nz)
  CALL mpisplit3d(in3d,nx,ny,nz,qv(:,:,:,tim))

  DO nqin = 1,nsclin
    IF(myproc == 0) READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nz)
    DO nq = 1,nscalar
      IF (nqin == nqscalarin(nq)) THEN
        CALL mpisplit3d(in3d,nx,ny,nz,qscalar(:,:,:,tim,nq))
        EXIT
      END IF
    END DO
  END DO

  IF(myproc == 0) READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nz)
  CALL mpisplit3d(in3d,nx,ny,nz,tke(:,:,:,tim))

  tim = tpresent

  IF(myproc == 0) READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nz)
  CALL mpisplit3d(in3d,nx,ny,nz,u(:,:,:,tim))

  IF(myproc == 0) READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nz)
  CALL mpisplit3d(in3d,nx,ny,nz,v(:,:,:,tim))

  IF(myproc == 0) READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nz)
  CALL mpisplit3d(in3d,nx,ny,nz,w(:,:,:,tim))

  IF(myproc == 0) READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nz)
  CALL mpisplit3d(in3d,nx,ny,nz,ptprt(:,:,:,tim))

  IF(myproc == 0) READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nz)
  CALL mpisplit3d(in3d,nx,ny,nz,pprt(:,:,:,tim))

  IF(myproc == 0) READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nz)
  CALL mpisplit3d(in3d,nx,ny,nz,qv(:,:,:,tim))

  DO nqin = 1,nsclin
    IF(myproc == 0) READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nz)
    DO nq = 1,nscalar
      IF (nqin == nqscalarin(nq)) THEN
        CALL mpisplit3d(in3d,nx,ny,nz,qscalar(:,:,:,tim,nq))
        EXIT
      END IF
    END DO
  END DO

  IF(myproc == 0) READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nz)
  CALL mpisplit3d(in3d,nx,ny,nz,tke(:,:,:,tim))

  IF(myproc == 0) READ(rstiunt,ERR=999) in2dew
  CALL mpisplit2dew(in2dew,ny,nz,udteb)
  IF(myproc == 0) READ(rstiunt,ERR=999) in2dew
  CALL mpisplit2dew(in2dew,ny,nz,udtwb)

  IF(myproc == 0) READ(rstiunt,ERR=999) in2dns
  CALL mpisplit2dns(in2dns,nx,nz,vdtnb)
  IF(myproc == 0) READ(rstiunt,ERR=999) in2dns
  CALL mpisplit2dns(in2dns,nx,nz,vdtsb)

  IF(myproc == 0) READ(rstiunt,ERR=999) in2dew
  CALL mpisplit2dew(in2dew,ny,nz,pdteb)
  IF(myproc == 0) READ(rstiunt,ERR=999) in2dew
  CALL mpisplit2dew(in2dew,ny,nz,pdtwb)
  IF(myproc == 0) READ(rstiunt,ERR=999) in2dns
  CALL mpisplit2dns(in2dns,nx,nz,pdtnb)
  IF(myproc == 0) READ(rstiunt,ERR=999) in2dns
  CALL mpisplit2dns(in2dns,nx,nz,pdtsb)

!-----------------------------------------------------------------------
!
!TINA initialize cc bubble here if using restart file and ccin = 1
!This initializes tpast,present,future
!
!-----------------------------------------------------------------------
  IF (P_CC > 0) THEN
    IF (ccin == 1 .AND. cpoint < 0) THEN
      CALL ccinit(nx,ny,nz,3,x,y,z,zp,qscalar(:,:,:,:,P_CC),tem1)
      IF (myproc == 0) print *,'Initializing cc as bubble'
    END IF
  END IF
!TINA/michi

!-----------------------------------------------------------------------
!
!  Set the future values of variables to their current values.
!  This is done primarily for safety reasons since the arrays at
!  tfuture will be overwritten by the new values during the
!  time integration.
!
!-----------------------------------------------------------------------
!
  CALL cpyary3d(nx,ny,nz,u (1,1,1,tpresent) , u (1,1,1,tfuture))
  CALL cpyary3d(nx,ny,nz,v (1,1,1,tpresent) , v (1,1,1,tfuture))
  CALL cpyary3d(nx,ny,nz,w (1,1,1,tpresent) , w (1,1,1,tfuture))
  CALL cpyary3d(nx,ny,nz,ptprt(1,1,1,tpresent),                     &
          ptprt(1,1,1,tfuture))
  CALL cpyary3d(nx,ny,nz,pprt (1,1,1,tpresent),                     &
          pprt (1,1,1,tfuture))
  CALL cpyary3d(nx,ny,nz,qv(1,1,1,tpresent) , qv(1,1,1,tfuture))
  DO nq = 1,nscalar
    CALL cpyary3d(nx,ny,nz,qscalar(1,1,1,tpresent,nq),                  &
                                              qscalar(1,1,1,tfuture,nq))
  END DO
  CALL cpyary3d(nx,ny,nz,tke(1,1,1,tpresent), tke(1,1,1,tfuture))

  IF ( sfcrstin /= 0 ) THEN

    IF(myproc == 0) PRINT *,'read in sfc/soil variables:'

    IF(myproc == 0) READ(rstiunt,ERR=999) in3di
    CALL mpisplit3di(in3di,nx,ny,nstyps,soiltyp)
    IF(myproc == 0)            &
       READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nstyps)
    CALL mpisplit3d(in3d,nx,ny,nstyps,stypfrct)

    IF(myproc == 0) READ(rstiunt,ERR=999) in2di
    CALL mpisplit2di(in2di,nx,ny,vegtyp)
    IF(myproc == 0) READ(rstiunt,ERR=999) in2d
    CALL mpisplit2d(in2d,nx,ny,lai)
    IF(myproc == 0) READ(rstiunt,ERR=999) in2d
    CALL mpisplit2d(in2d,nx,ny,roufns)
    IF(myproc == 0) READ(rstiunt,ERR=999) in2d
    CALL mpisplit2d(in2d,nx,ny,veg)

    IF(myproc == 0)            &
       READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nstyps+1)
    CALL mpisplit3d(in3d,nx,ny,nstyps+1,qvsfc)
    IF(myproc == 0) READ(rstiunt,ERR=999) in4d
    CALL mpisplit4d(in4d,nx,ny,nzsoil,nstyps+1,tsoil)
    IF(myproc == 0) READ(rstiunt,ERR=999) in4d
    CALL mpisplit4d(in4d,nx,ny,nzsoil,nstyps+1,qsoil)
    IF(myproc == 0)            &
       READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nstyps+1)
    CALL mpisplit3d(in3d,nx,ny,nstyps+1,wetcanp)
    IF(myproc == 0) READ(rstiunt,ERR=999) in2d
    CALL mpisplit2d(in2d,nx,ny,snowdpth)

  END IF

  IF ( prcrsin /= 0 ) THEN
    IF(myproc == 0) READ(rstiunt,ERR=999) in2d
    CALL mpisplit2d(in2d,nx,ny,raing)
    IF(myproc == 0) READ(rstiunt,ERR=999) in2d
    CALL mpisplit2d(in2d,nx,ny,rainc)
    IF(myproc == 0) READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,4)
    CALL mpisplit3d(in3d,nx,ny,4,prcrate)
  END IF

  IF ( rcumin /= 0 ) THEN
    IF(myproc == 0)            &
          READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nz)
    CALL mpisplit3d(in3d,nx,ny,nz,ptcumsrc)
    IF(myproc == 0) READ(rstiunt,ERR=999) in4dq
    CALL mpisplit4d(in4dq,nx,ny,nz,5,qcumsrc)
  END IF

  IF ( exbcin /= 0 ) THEN
    IF ( lbcopt == 2 ) THEN
      IF(myproc == 0) THEN

        READ(rstiunt,ERR=999) abstfcst0, abstfcst,                      &
                 ubcrd,vbcrd,wbcrd,ptbcrd,prbcrd,                       &
                 qvbcrd,qcbcrd,qrbcrd,qibcrd,qsbcrd,qhbcrd,qgbcrd,      &
                 ncbcrd,nrbcrd,nibcrd,nsbcrd,ngbcrd,nhbcrd,             &
                 zrbcrd,zibcrd,zsbcrd,zgbcrd,zhbcrd,ccbcrd

        qscalarbcrd(:) = 0
        IF (P_QC >0) qscalarbcrd(P_QC) = qcbcrd
        IF (P_QR >0) qscalarbcrd(P_QR) = qrbcrd
        IF (P_QI >0) qscalarbcrd(P_QI) = qibcrd
        IF (P_QS >0) qscalarbcrd(P_QS) = qsbcrd
        IF (P_QH >0) qscalarbcrd(P_QH) = qhbcrd
        IF (P_QG >0) qscalarbcrd(P_QG) = qgbcrd

        IF (P_NC >0) qscalarbcrd(P_NC) = ncbcrd
        IF (P_NR >0) qscalarbcrd(P_NR) = nrbcrd
        IF (P_NI >0) qscalarbcrd(P_NI) = nibcrd
        IF (P_NS >0) qscalarbcrd(P_NS) = nsbcrd
        IF (P_NH >0) qscalarbcrd(P_NH) = nhbcrd
        IF (P_NG >0) qscalarbcrd(P_NG) = ngbcrd

        IF (P_ZR >0) qscalarbcrd(P_ZR) = zrbcrd
        IF (P_ZI >0) qscalarbcrd(P_ZI) = zibcrd
        IF (P_ZS >0) qscalarbcrd(P_ZS) = zsbcrd
        IF (P_ZH >0) qscalarbcrd(P_ZH) = zhbcrd
        IF (P_ZG >0) qscalarbcrd(P_ZG) = zgbcrd

        IF (P_CC >0) qscalarbcrd(P_CC) = ccbcrd

      END IF

      CALL mpupdatei(abstfcst0,1)
      CALL mpupdatei(abstfcst,1)
      CALL mpupdatei(ubcrd,1)
      CALL mpupdatei(vbcrd,1)
      CALL mpupdatei(wbcrd,1)
      CALL mpupdatei(ptbcrd,1)
      CALL mpupdatei(prbcrd,1)
      CALL mpupdatei(qvbcrd,1)
      CALL mpupdatei(qscalarbcrd,20)

      DO var=1,exbcbufsz,nx*ny*nz
          IF(myproc ==0)          &
            READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nz)
          CALL mpisplit3d(in3d,nx,ny,nz,exbcbuf(var))
      END DO

    ELSE
      IF(myproc == 0) THEN
          WRITE(6,'(a/a/a/a)')                                          &
              'WARNING: The restart file contains EXBC arrays, while',  &
              '         this run does not have EXBC option.',           &
              '         Therefore, the results from restart run may be', &
              '         alterred. The program will continue.'

          READ(rstiunt,ERR=999)

          READ(rstiunt,ERR=999)
          READ(rstiunt,ERR=999)
          READ(rstiunt,ERR=999)
          READ(rstiunt,ERR=999)
          READ(rstiunt,ERR=999)
          READ(rstiunt,ERR=999)
          READ(rstiunt,ERR=999)
          READ(rstiunt,ERR=999)
          READ(rstiunt,ERR=999)
          READ(rstiunt,ERR=999)
          READ(rstiunt,ERR=999)
      END IF   ! myproc == 0
    END IF
  END IF

  IF ( mapfin == 1 ) THEN
    IF(myproc == 0)            &
          READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,8)
    CALL mpisplit3d(in3d,nx,ny,8,mapfct)
  END IF

  IF ( radrstin == 1 ) THEN
    IF(myproc == 0)            &
          READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nz)
    CALL mpisplit3d(in3d,nx,ny,nz,radfrc)
    IF(myproc == 0) READ(rstiunt,ERR=999) in2d
    CALL mpisplit2d(in2d,nx,ny,radsw)
    IF(myproc == 0) READ(rstiunt,ERR=999) in2d
    CALL mpisplit2d(in2d,nx,ny,rnflx)
    IF(myproc == 0) READ(rstiunt,ERR=999) in2d
    CALL mpisplit2d(in2d,nx,ny,radswnet)
    IF(myproc == 0) READ(rstiunt,ERR=999) in2d
    CALL mpisplit2d(in2d,nx,ny,radlwin)

  END IF

  IF ( kfrsin /= 0) THEN
    IF(myproc == 0)            &
          READ(rstiunt,ERR=999) (((in3d(i,j,k),i=1,nxlg),j=1,nylg),k=1,nz)
    CALL mpisplit3d(in3d,nx,ny,nz,w0avg)
    IF(myproc == 0) READ(rstiunt,ERR=999) in2di
    CALL mpisplit2di(in2di,nx,ny,nca)
    IF(myproc == 0) READ(rstiunt,ERR=999) in2d
    CALL mpisplit2d(in2d,nx,ny,kfraincv)
  END IF

  IF ( bmjsin /= 0) THEN
    IF(myproc == 0) READ(rstiunt, ERR=999) in2d
    CALL mpisplit2d(in2d,nx,ny,cldefi)
    IF(myproc == 0) READ(rstiunt, ERR=999) in2d
    CALL mpisplit2d(in2d,nx,ny,xland)
    IF(myproc == 0) READ(rstiunt, ERR=999) in2d
    CALL mpisplit2d(in2d,nx,ny,bmjraincv)
  END IF

  IF(myproc == 0) THEN
    CLOSE (UNIT=rstiunt)
    CALL retunit( rstiunt )
  END IF
!
!-----------------------------------------------------------------------
!
!  Reset the model u and v velocity values using the new
!  domain translation speed.
!
!-----------------------------------------------------------------------
!

  IF( nint(umove) == 999 .OR. nint(vmove) == 999 ) THEN

    umove = umoveold
    vmove = vmoveold

  ELSE IF (umoveold /= umove .OR. vmoveold /= vmove ) THEN

    WRITE(6,'(3(/1x,a)/)')                                          &
        'ATTENTION: UMOVE or VMOVE in the input file were different ', &
        'from those in the restart file. Subroutine ADJUVMV is called', &
        'to adjust the time-dependent variables for option grdtrns!=0.'

    IF ( grdtrns /= 0 ) THEN
      uchange = umove - umoveold
      vchange = vmove - vmoveold

      CALL adjuvmv(nx,ny,nz,                                        &
                   ubar,vbar,u,v,w,ptprt,pprt,qv,qvbar,qscalar,     &
                   uchange, vchange, tem1, tem2)

    END IF
  END IF

  CALL jacob(nx,ny,nz,x,y,z,zp,j1,j2,j3,tem1)

  IF(myproc == 0) WRITE(6,'(/1x,a/,1x,a/)')                             &
      'Grid definition arrays are read in from initialization data',    &
      'those set in INIGRD are superceded.'

!
!-----------------------------------------------------------------------
!
!  Print out the domain-wide max/min of output variables.
!
!-----------------------------------------------------------------------

  IF(myproc ==0) WRITE(6,'(/1x,a/)')                                   &
      'Min. and max. of the data arrays read in from restart data:'

  CALL a3dmax0(x,1,nx,1,nx,1,1,1,1, 1,1,1,1, amax,amin)
  IF(myproc ==0)   &
    WRITE(6,'(/1x,2(a,e13.6))') 'xmin    = ', amin,',  xmax    =',amax

  CALL a3dmax0(y,1,ny,1,ny,1,1,1,1, 1,1,1,1, amax,amin)
  IF(myproc ==0)   &
    WRITE(6,'(1x,2(a,e13.6))') 'ymin    = ', amin,',  ymax    =',amax

  CALL a3dmax0(z,1,nz,1,nz,1,1,1,1, 1,1,1,1, amax,amin)
  IF(myproc ==0)   &
    WRITE(6,'(1x,2(a,e13.6))') 'zmin    = ', amin,',  zmax    =',amax

  CALL a3dmax0(zp,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz, amax,amin)
  IF(myproc ==0)   &
    WRITE(6,'(1x,2(a,e13.6))') 'zpmin   = ', amin,',  zpmax   =',amax

  CALL a3dmax0(zpsoil,1,nx,1,nx-1,1,ny,1,ny-1,1,nzsoil,1,nzsoil, amax,amin)
  IF(myproc ==0)   &
    WRITE(6,'(1x,2(a,e13.6))') 'zpsoilmin   = ', amin,',  zpsoilmax   =',amax

  CALL a3dmax0(hterain,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1, amax,amin)
  IF(myproc ==0)   &
    WRITE(6,'(1x,2(a,e13.6))') 'hmin    = ', amin,',  hmax    =',amax

  CALL a3dmax0(ubar,1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1, amax,amin)
  IF(myproc ==0)   &
    WRITE(6,'(1x,2(a,e13.6))') 'ubarmin = ', amin,',  ubarmax =',amax

  CALL a3dmax0(vbar,1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1, amax,amin)
  IF(myproc ==0)   &
    WRITE(6,'(1x,2(a,e13.6))') 'vbarmin = ', amin,',  vbarmax =',amax

  CALL a3dmax0(ptbar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
  IF(myproc ==0)   &
    WRITE(6,'(1x,2(a,e13.6))') 'ptbarmin= ', amin,',  ptbarmax=',amax

  CALL a3dmax0(pbar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1, amax,amin)
  IF(myproc ==0)   &
    WRITE(6,'(1x,2(a,e13.6))') 'pbarmin = ', amin,',  pbarmax =',amax

  CALL a3dmax0(rhostr,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
  IF(myproc ==0)   &
    WRITE(6,'(1x,2(a,e13.6))') 'rhostrmin=', amin,', rhostrmax=',amax

  CALL a3dmax0(qvbar,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,amax,amin)
  IF(myproc ==0)   &
    WRITE(6,'(1x,2(a,e13.6))') 'qvbarmin= ', amin,',  qvbarmax=',amax


  DO i = 1,2

    IF( i == 1) THEN
      IF(myproc ==0)   &
      WRITE(6,'(/1x,a/)') 'Min/max of fields at tpresent:'
      tim = tpresent
    ELSE
      IF(myproc ==0)   &
      WRITE(6,'(/1x,a/)') 'Min/max of fields at tpast:'
      tim = tpast
    END IF

    CALL a3dmax0(u(1,1,1,tim),1,nx,1,nx,1,ny,1,ny-1,1,nz,1,nz-1,        &
                 amax,amin)
    IF(myproc ==0)   &
    WRITE(6,'(1x,2(a,e13.6))') 'umin    = ', amin,',  umax    =',       &
          amax

    CALL a3dmax0(v(1,1,1,tim),1,nx,1,nx-1,1,ny,1,ny,1,nz,1,nz-1,        &
                 amax,amin)
    IF(myproc ==0)   &
    WRITE(6,'(1x,2(a,e13.6))') 'vmin    = ', amin,',  vmax    =',       &
          amax

    CALL a3dmax0(w(1,1,1,tim),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz,        &
                 amax,amin)
    IF(myproc ==0)   &
    WRITE(6,'(1x,2(a,e13.6))') 'wmin    = ', amin,',  wmax    =',       &
          amax

    CALL a3dmax0(ptprt(1,1,1,tim),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,       &
                 nz-1,amax,amin)
    IF(myproc ==0)   &
    WRITE(6,'(1x,2(a,e13.6))') 'ptprtmin= ', amin,',  ptprtmax=',       &
          amax

    CALL a3dmax0(pprt(1,1,1,tim),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,        &
                 nz-1,amax,amin)
    IF(myproc ==0)   &
    WRITE(6,'(1x,2(a,e13.6))') 'pprtmin = ', amin,',  pprtmax =',       &
          amax

    CALL a3dmax0(qv(1,1,1,tim),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,     &
                 amax,amin)
    IF(myproc ==0)   &
    WRITE(6,'(1x,2(a,e13.6))') 'qvmin   = ', amin,',  qvmax   =',       &
          amax

    DO nq = 1,nscalar
      CALL a3dmax0(qscalar(1,1,1,tim,nq),1,nx,1,nx-1,1,ny,1,ny-1,       &
                   1,nz,1,nz-1,amax,amin)
      IF(myproc == 0) WRITE(6,'(1x,a,e13.6,2a,e13.6)')                  &
                   TRIM(qnames(nq))//'min   = ', amin,', ',             &
                   TRIM(qnames(nq))//'max   = ', amax
    END DO

  END DO

  IF(myproc ==0) WRITE(6,'(/1x,a/)')                                     &
        'Min/max of fields for other one time level arrays:'

  CALL a3dmax0(tsoil(1,1,1,0),1,nx,1,nx-1,1,ny,1,ny-1,1,nzsoil,1,nzsoil, &
               amax,amin)
  IF(myproc ==0)   &
  WRITE(6,'(1x,2(a,e13.6))') 'tsoilmin= ', amin,', tsoilmax =',amax

  CALL a3dmax0(qsoil(1,1,1,0),1,nx,1,nx-1,1,ny,1,ny-1,1,nzsoil,1,nzsoil, &
               amax,amin)
  IF(myproc ==0)   &
  WRITE(6,'(1x,2(a,e13.6))') 'qsoilmin = ', amin,', qsoilmax  =',amax

  CALL a3dmax0(wetcanp(1,1,0),1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1,          &
               amax,amin)
  IF(myproc ==0)   &
  WRITE(6,'(1x,2(a,e13.6))') 'wetcmin = ', amin,', wetcmax  =',amax

  CALL a3dmax0(raing,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1, amax,amin)
  IF(myproc ==0)   &
  WRITE(6,'(1x,2(a,e13.6))') 'raingmin= ', amin,', raingmax =',amax

  CALL a3dmax0(rainc,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1, amax,amin)
  IF(myproc ==0)   &
  WRITE(6,'(1x,2(a,e13.6))') 'raincmin= ', amin,', raincmax =',amax

  CALL a3dmax0(ptcumsrc,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,            &
               amax,amin)
  IF(myproc ==0)   &
  WRITE(6,'(1x,2(a,e13.6))') 'ptcummin= ', amin,',  ptcummax=',amax

  CALL a3dmax0(qcumsrc(1,1,1,1),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,    &
               amax,amin)
  IF(myproc ==0)   &
  WRITE(6,'(1x,2(a,e13.6))') 'qvcummin= ', amin,',  qvcummax=',amax
  CALL a3dmax0(qcumsrc(1,1,1,2),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,    &
               amax,amin)
  IF(myproc ==0)   &
  WRITE(6,'(1x,2(a,e13.6))') 'qccummin= ', amin,',  qccummax=',amax
  CALL a3dmax0(qcumsrc(1,1,1,3),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,    &
               amax,amin)
  IF(myproc ==0)   &
  WRITE(6,'(1x,2(a,e13.6))') 'qrcummin= ', amin,',  qrcummax=',amax
  CALL a3dmax0(qcumsrc(1,1,1,4),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,    &
               amax,amin)
  IF(myproc ==0)   &
  WRITE(6,'(1x,2(a,e13.6))') 'qicummin= ', amin,',  qicummax=',amax
  CALL a3dmax0(qcumsrc(1,1,1,5),1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1,    &
               amax,amin)
  IF(myproc ==0)   &
  WRITE(6,'(1x,2(a,e13.6))') 'qscummin= ', amin,',  qscummax=',amax
!
!-----------------------------------------------------------------------
!
!  Compress the restart file if it was originally compressed.
!
!-----------------------------------------------------------------------
!
  IF( cmprsed .AND. filcmprs == 1 .AND. myproc == 0) THEN
    CALL cmprs( rstinf(1:lrstfn) )
  END IF

  GOTO 990

  999   CONTINUE
  WRITE(6,'(a)') ' Error reading restart data '//rstinf
  WRITE(6,'(a,i3,a)') ' Fortran channel ',rstiunt,' was used.'
  WRITE(6,'(a)') ' Job stopped in subroutine rstinsplit!'
  IF(myproc ==0)  CLOSE (UNIT=rstiunt)

  CALL arpsstop('arpsstop called from RSTINSPLIT error reading restart'//   &
                 'file ',1)

  990 CONTINUE

  DEALLOCATE(in1d,in2d,in3d)
  DEALLOCATE(in2di,in3di)
  DEALLOCATE(in2dns,in2dew)
  DEALLOCATE(in4d,in4dq)

  RETURN
END SUBROUTINE rstinsplit
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CPYARY3D                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE cpyary3d(nx,ny,nz, ain, aout)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Copy the contents of array 'ain' into 'aout'.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  11/10/1992.
!
!  MODIFICATION HISTORY:
!
!  2/11/93 (K. Droegemeier)
!  Added full documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    ain      Input array
!    nx       1st Dimension of input and output arrays.
!    ny       2nd Dimension of input and output arrays.
!    nz       3rd Dimension of input and output arrays.
!
!  OUTPUT:
!
!    aout     Output array.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: nx,ny,nz
  REAL,    INTENT(IN)  :: ain (nx,ny,nz) ! Input array to be copied in aout.
  REAL,    INTENT(OUT) :: aout(nx,ny,nz) ! Array whose value will be copied from ain.

  INTEGER :: i,j,k

!-----------------------------------------------------------------------
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        aout(i,j,k) = ain(i,j,k)
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE cpyary3d
