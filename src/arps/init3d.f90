!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE INITIAL                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE initial(mptr,nx,ny,nz,nzsoil,nxndg,nyndg,nzndg,nstyps,       &
           exbcbufsz,                                                   &
           u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,                       &
           udteb,udtwb,udtnb,udtsb,vdteb,vdtwb,vdtnb,vdtsb,             &
           wdteb,wdtwb,wdtnb,wdtsb,pdteb,pdtwb,pdtnb,pdtsb,             &
           sdteb,sdtwb,sdtnb,sdtsb,                                     &
           ubar,vbar,wbar,ptbar,pbar,ptbari,pbari,                      &
           rhostr,rhostri,qvbar,ppi,csndsq,                             &
           x,y,z,zp,zpsoil,hterain,mapfct,                              &
           j1,j2,j3,j3soil,aj3x,aj3y,aj3z,j3inv,j3soilinv,              &
           trigs1,trigs2,ifax1,ifax2,                                   &
           wsave1,wsave2,vwork1,vwork2,                                 &
           sinlat, kmh,kmv,rprntl,                                      &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,ptsfc,qvsfc,                    &
           ptcumsrc,qcumsrc,w0avg,nca,kfraincv,                         &
           cldefi,xland,bmjraincv,                                      &
           raing,rainc,prcrate,exbcbuf,bcrlx,radfrc,radsw,              &
           rnflx,radswnet,radlwin,usflx,vsflx,ptsflx,qvsflx,            &
           uincr,vincr,wincr,pincr,ptincr,qvincr,                       &
           qcincr,qrincr,qiincr,qsincr,qhincr,                          &
           tem1soil,tem2soil,tem3soil,tem4soil,tem5soil,                &
           temxy1,tem1,tem2,tem3,tem4,tem5,tem6,tem7,                   &
           tem8,tem9,tem10,tem11,tem12,tem13,                           &
           tem14,tem15,tem16,temscalar,                                 &
           tem1_0,tem2_0,tem3_0)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Initialize the model parameters and variables, including base state
!  variables, dependent variables and grid structure.
!
!  This is the main driver for model initializations.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  11/5/92.
!
!  MODIFICATION HISTORY:
!
!  5/02/92 (M. Xue)
!  Added full documentation.
!
!  5/03/92 (M. Xue)
!  Further documentation.
!
!  9/14/1992 (M. Xue)
!  Different surface drag coefficients defined for momentum, heat and
!  moisture.
!  Three options included for the Coriolis force calculations.
!
!  2/12/94 (Yuhe Liu)
!  Surface data and variables added for surface energy budget model.
!
!  6/10/94 (M. Xue &AS)
!  Added call to initpwr.
!
!  02/07/1995 (Yuhe Liu)
!  Added a new 2-D permanent array, veg(nx,ny), to the argument list
!
!  08/30/1995 (Yuhe Liu)
!  Moved the initialization of external boundary arrays from the
!  main program to this subroutine.
!
!  10/31/95   (D. Weber)
!  Added trigs1,trigs2,ifax1,ifax2 for use in the upper w-p
!  radiation condition.
!
!  1/22/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor to ARPS governing equations.
!
!  07/22/97 (D. Weber)
!  Added wsave1,wsave2,vwork1,vwork2 for use in the even fft version
!  of the upper w-p radiation condition (fftopt=2).
!
!  08/01/97 (Zonghui Huo)
!  Added Kain-fritsch cumulus parameterization scheme.
!
!  11/05/97 (D. Weber)
!  Added temxy5 array for use in the bottom boundary condition for
!  pertrubation pressure (hydrostatic).
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
!  9/15/1998 (D. Weber)
!  Added ptbari, pbari (inverse) for use in optimizing the code.
!
!  8/31/1998 (K. Brewster)
!  Added call to ININUDGE to initialize nudging terms.
!
!  11/18/1998 (K. Brewster)
!  Changed pibar to ppi (full pi) and moved initialization.
!
!  12/09/1998 (Donghai Wang)
!  Added the snow cover.
!
!  07/10/2001 (K. Brewster)
!  Added increment arrays to argument list and to call to ininudge.
!
!  03/13/2002 (Eric Kemp)
!  Added arrays for WRF BMJ cumulus scheme.
!
!  April 2002 (Fanyou Kong)
!  Added WRF new Kain-Fritsch (April 2002 version: KF_ETA) scheme
!  initialization (lookup table)
!
!  05/02/2002 (Dan Weber and Jerry Brotzge)
!  Added arrays for the new soil model.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    mptr     Grid identifier.
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of grid points in the soil model in the -z-direction
!
!    nxndg    Number of x grid points for nudging (1 or nx)
!    nyndg    Number of y grid points for nudging (1 or ny)
!    nzndg    Number of z grid points for nudging (1 or nz)
!
!  OUTPUT:
!
!    u        x-component of velocity at all time levels (m/s).
!    v        y-component of velocity at all time levels (m/s).
!    w        z-component of velocity at all time levels (m/s).
!    wcont    Contravariant vertical velocity (m/s)
!    ptprt    Perturbation potential temperature at all time levels (K)
!    pprt     Perturbation pressure at all time levels (Pascal)
!    qv       Water vapor specific humidity at all time levels (kg/kg)
!    qc       Cloud water mixing ratio at all time levels (kg/kg)
!    qr       Rainwater mixing ratio at all time levels (kg/kg)
!    qi       Cloud ice mixing ratio at all time levels (kg/kg)
!    qs       Snow mixing ratio at all time levels (kg/kg)
!    qh       Hail mixing ratio at all time levels (kg/kg)
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!
!    udteb    Time tendency of u field at east boundary (m/s**2)
!    udtwb    Time tendency of u field at west boundary (m/s**2)
!    udtnb    Time tendency of u field at north boundary (m/s**2)
!    udtsb    Time tendency of u field at south boundary (m/s**2)
!
!    vdteb    Time tendency of v field at east boundary (m/s**2)
!    vdtwb    Time tendency of v field at west boundary (m/s**2)
!    vdtnb    Time tendency of v field at north boundary (m/s**2)
!    vdtsb    Time tendency of v field at south boundary (m/s**2)
!
!    wdteb    Time tendency of w field at east boundary (m/s**2)
!    wdtwb    Time tendency of w field at west boundary (m/s**2)
!    wdtnb    Time tendency of w field at north boundary (m/s**2)
!    wdtsb    Time tendency of w field at south boundary (m/s**2)
!
!    pdteb    Time tendency of pprt field at east boundary (PASCAL/s)
!    pdtwb    Time tendency of pprt field at west boundary (PASCAL/s)
!    pdtnb    Time tendency of pprt field at north boundary (PASCAL/s)
!    pdtsb    Time tendency of pprt field at south boundary (PASCAL/s)
!
!    sdteb    Time tendency of a scalar field at east boundary (m/s**2)
!    sdtwb    Time tendency of a scalar field at west boundary (m/s**2)
!    sdtnb    Time tendency of a scalar field at north boundary (m/s**2)
!    sdtsb    Time tendency of a scalar field at south boundary (m/s**2)
!
!    ubar     Base state x-velocity component (m/s)
!    vbar     Base state y-velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    ptbari   Inverse Base state potential temperature (K)
!    pbari    Inverse Base state pressure (Pascal)
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    rhostri  Inverse base state density rhobar times j3 (kg/m**3)
!    qvbar    Base state water vapor specific humidity (kg/kg).
!    ppi      Exner function
!    csndsq   Sound wave speed squared.
!
!    x        x-coordinate of grid points in computational space (m)
!    y        y-coordinate of grid points in computational space (m)
!    z        z-coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!    zpsoil   Vertical coordinate of grid points in the soil model
!             in physical space (m).
!    hterain  Terrain height (m)
!    mapfct   Map factors at scalar, u and v points
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    j3soil   Soil coordinate transformation Jacobian  d(zpsoil)/dz
!    aj3x     Avgx of the coordinate transformation Jacobian  d(zp)/dz
!    aj3y     Avgy of the coordinate transformation Jacobian  d(zp)/dz
!    aj3z     Avgz of the coordinate transformation Jacobian  d(zp)/dz
!    j3inv    Inverse of the coordinate transformation j3
!    j3soilinv Inverse of the soil coordinate transformation j3soil
!
!    trigs1   Array containing pre-computed trig function for fftopt=1.
!    trigs2   Array containing pre-computed trig function for fftopt=1.
!    ifax1    Array containing the factors of nx for fftopt=1.
!    ifax2    Array containing the factors of ny for fftopt=1.
!
!    vwork1   2-D work array for fftopt=2 option.
!    vwork2   2-D work array for fftopt=2 option.
!    wsave1   Work array for fftopt=2 option.
!    wsave2   Work array for fftopt=2 option.
!
!    sinlat   Sin of latitude at each grid point
!
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!    rprntl   Reciprocal of Prandtl number
!
!    soiltyp  Soil type
!    stypfrct  Soil type fraction
!    vegtyp   Vegetation type
!    lai      Leaf Area Index
!    roufns   Surface roughness
!    veg      Vegetation fraction
!
!    tsoil    Soil temperature (K)
!    qsoil    Soil moisture (m**3/m**3)
!    wetcanp  Canopy water amount
!    ptsfc    Ground surface potential temperature (K)
!    qvsfc    Effective S.H. at sfc.
!
!    ptcumsrc Source term in pt-equation due to cumulus parameterization
!    qcumsrc  Source term in water equations due to cumulus parameterization
!
!    nca      Counter for CAPE release in the Kain-Fritsch scheme
!    kfraincv K-F convective precipitation rate
!    cldefi   BMJ cloud efficiency
!    xland    BMJ land/sea mask
!    bmjraincv BMJ convective precipitation amount
!
!    radfrc   Radiation forcing (K/s)
!    radsw    Solar radiation reaching the surface
!    rnflx    Net absorbed radiation by the surface
!    radswnet Net shortwave radiation
!    radlwin  Incoming longwave radiation
!
!    raing    Grid supersaturation rain
!    rainc    Cumulus convective rain
!
!    usflx    Surface flux of u-momentum (kg/(m*s**2))
!    vsflx    Surface flux of v-momentum (kg/(m*s**2))
!    ptsflx   Surface heat flux (K*kg/(m**2 * s ))
!    qvsflx   Surface moisture flux of (kg/(m**2 * s))
!
!    temxy1   Temporary work array
!
!    tem1soil Soil model temporary work array.
!    tem2soil Soil model temporary work array.
!    tem3soil Soil model temporary work array.
!    tem4soil Soil model temporary work array.
!    tem5soil Soil model temporary work array.
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
!    tem10    Temporary work array.
!    tem11    Temporary work array.
!    tem12    Temporary work array.
!    tem13    Temporary work array.
!    tem14    Temporary work array.
!    tem15    Temporary work array.
!    tem16    Temporary work array.
!    tem17    Temporary work array.
!    tem18    Temporary work array.
!    tem19    Temporary work array.
!    tem20    Temporary work array.
!    tem21    Temporary work array.
!    tem22    Temporary work array.
!    tem23    Temporary work array.
!    tem24    Temporary work array.
!    tem25    Temporary work array.
!    tem26    Temporary work array.
!
!    tem1_0   Temporary work array.
!    tem2_0   Temporary work array.
!    tem3_0   Temporary work array.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'timelvls.inc'
  INCLUDE 'globcst.inc'

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  INTEGER :: mptr              ! Grid identifier

  INTEGER :: nx,ny,nz          ! The number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of grid points in the -z-direction

  INTEGER :: nxndg,nyndg,nzndg ! The number of grid points in 3 directions

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s).

  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s).
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s).
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)
  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature
                               ! from that of base state atmosphere (Kelvin).
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure from that
                               ! of base state atmosphere (Pascal).
  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg).

  REAL :: qscalar(nx,ny,nz,nt,nscalar)

  REAL :: tke   (nx,ny,nz,nt)  ! Turbulent Kinetic Energy ((m/s)**2)

  REAL :: udteb (ny,nz)        ! T-tendency of u at e-boundary (m/s**2)
  REAL :: udtwb (ny,nz)        ! T-tendency of u at w-boundary (m/s**2)
  REAL :: udtnb (nx,nz)        ! T-tendency of u at n-boundary (m/s**2)
  REAL :: udtsb (nx,nz)        ! T-tendency of u at s-boundary (m/s**2)

  REAL :: vdteb (ny,nz)        ! T-tendency of v at e-boundary (m/s**2)
  REAL :: vdtwb (ny,nz)        ! T-tendency of v at w-boundary (m/s**2)
  REAL :: vdtnb (nx,nz)        ! T-tendency of v at n-boundary (m/s**2)
  REAL :: vdtsb (nx,nz)        ! T-tendency of v at s-boundary (m/s**2)

  REAL :: wdteb (ny,nz)        ! T-tendency of w at e-boundary (m/s**2)
  REAL :: wdtwb (ny,nz)        ! T-tendency of w at w-boundary (m/s**2)
  REAL :: wdtnb (nx,nz)        ! T-tendency of w at n-boundary (m/s**2)
  REAL :: wdtsb (nx,nz)        ! T-tendency of w at s-boundary (m/s**2)

  REAL :: pdteb (ny,nz)        ! T-tendency of pprt at e-boundary (PASCAL/s)
  REAL :: pdtwb (ny,nz)        ! T-tendency of pprt at w-boundary (PASCAL/s)
  REAL :: pdtnb (nx,nz)        ! T-tendency of pprt at n-boundary (PASCAL/s)
  REAL :: pdtsb (nx,nz)        ! T-tendency of pprt at s-boundary (PASCAL/s)

  REAL :: sdteb (ny,nz)        ! T-tendency of w at e-boundary (m/s**2)
  REAL :: sdtwb (ny,nz)        ! T-tendency of w at w-boundary (m/s**2)
  REAL :: sdtnb (nx,nz)        ! T-tendency of w at n-boundary (m/s**2)
  REAL :: sdtsb (nx,nz)        ! T-tendency of w at s-boundary (m/s**2)

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s).
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s).
  REAL :: wbar  (nx,ny,nz)     ! Base state w-velocity (m/s).
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: ptbari(nx,ny,nz)     ! Inverse Base state pot. temperature (K)
  REAL :: pbari (nx,ny,nz)     ! Inverse Base state pressure (Pascal).
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: rhostri(nx,ny,nz)    ! Inverse base state density rhobar times j3.
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                               ! (kg/kg).
  REAL :: ppi   (nx,ny,nz)     ! Exner function.
  REAL :: csndsq(nx,ny,nz)     ! Sound wave speed squared.

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point on the staggered grid.
  REAL :: zpsoil(nx,ny,nzsoil) ! The physical height coordinate defined
                               ! at the center of a soil layer(m).

  REAL :: hterain(nx,ny)       ! Terrain height (m).

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian -d(zp)/dx.
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian -d(zp)/dy.
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian  d(zp)/dz.
  REAL :: j3soil(nx,ny,nzsoil) ! Coordinate transformation Jacobian
                               ! defined as d( zpsoil )/d( zsoil ).
  REAL :: aj3x  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE X-DIR.
  REAL :: aj3y  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
  REAL :: aj3z  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Z-DIR.
  REAL :: j3inv (nx,ny,nz)     ! Inverse of J3
  REAL :: j3soilinv(nx,ny,nzsoil) ! Inverse of J3soil.

  REAL :: trigs1(3*(nx-1)/2+1) ! Array containing pre-computed trig
                               ! function for fftopt=1.
  REAL :: trigs2(3*(ny-1)/2+1) ! Array containing pre-computed trig
                               ! function for fftopt=1.
  INTEGER :: ifax1(13)         ! Array containing the factors of nx for
                               ! fftopt=1.
  INTEGER :: ifax2(13)         ! Array containing the factors of ny for
                               ! fftopt=1.

  REAL :: vwork1 (nx+1,ny+1)   ! 2-D work array for fftopt=2.
  REAL :: vwork2 (ny,nx+1)     ! 2-D work array for fftopt=2.
  REAL :: wsave1 (3*(ny-1)+15) ! Work array for fftopt=2.
  REAL :: wsave2 (3*(nx-1)+15) ! Work array for fftopt=2.

  REAL :: sinlat(nx,ny)        ! Sin of latitude at each grid point

  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: rprntl(nx,ny,nz)     ! Reciprocal of Prandtl number


  INTEGER :: nstyps                  ! Number of soil type
  INTEGER :: soiltyp(nx,ny,nstyps)   ! Soil types at grids
  REAL    :: stypfrct(nx,ny,nstyps)  ! Fraction of soil types
  INTEGER :: vegtyp (nx,ny)          ! Vegetation type
  REAL :: lai    (nx,ny)          ! Leaf Area Index
  REAL :: roufns (nx,ny)          ! Surface roughness
  REAL :: veg    (nx,ny)          ! Vegetation fraction

  REAL :: tsoil (nx,ny,nzsoil,0:nstyps) ! Soil layer temperature (K)
  REAL :: qsoil (nx,ny,nzsoil,0:nstyps) ! Soil layer moisture  (m**3/m**3)
  REAL :: wetcanp(nx,ny,0:nstyps)    ! Canopy water amount
  REAL :: snowdpth(nx,ny)            ! Snow depth (m)

  REAL :: qvsfc  (nx,ny,0:nstyps)    ! Effective qv at sfc.
  REAL :: ptsfc  (nx,ny)             ! Ground surface potential
                                     ! temperature (K)

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

                               !EMK BMJ
  REAL,INTENT(INOUT) :: cldefi(nx,ny)      ! BMJ cloud efficiency
  REAL,INTENT(INOUT) :: xland(nx,ny)       ! BMJ land mask
                                           !   (1.0 = land, 2.0 = sea)
  REAL,INTENT(INOUT) :: bmjraincv(nx,ny)   ! BMJ convective rainfall (cm)
                               !EMK END

  REAL :: radfrc(nx,ny,nz)     ! Radiation forcing (K/s)
  REAL :: radsw(nx,ny)         ! Solar radiation reacing the surface
  REAL :: rnflx(nx,ny)         ! Net absorbed radiation by the surface
  REAL :: radswnet(nx,ny)      ! Net shortwave radiation
  REAL :: radlwin(nx,ny)       ! Incoming longwave radiation

  REAL :: raing(nx,ny)         ! Grid supersaturation rain
  REAL :: rainc(nx,ny)         ! Cumulus convective rain
  REAL :: prcrate(nx,ny,4)     ! precipitation rates (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precipitation rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulus precip. rate
                               ! prcrate(1,1,4) = microphysics precip. rate

  REAL :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  INTEGER :: exbcbufsz         ! EXBC buffer size
  REAL :: exbcbuf( exbcbufsz ) ! EXBC buffer array
  REAL :: bcrlx(nx,ny)         ! EXBC relaxation coefficients

  REAL :: uincr(nxndg,nyndg,nzndg)      ! Analysis increment for u
  REAL :: vincr(nxndg,nyndg,nzndg)      ! Analysis increment for v
  REAL :: wincr(nxndg,nyndg,nzndg)      ! Analysis increment for w
  REAL :: pincr(nxndg,nyndg,nzndg)      ! Analysis increment for p
  REAL :: ptincr(nxndg,nyndg,nzndg)     ! Analysis increment for pt
  REAL :: qvincr(nxndg,nyndg,nzndg)     ! Analysis increment for qv
  REAL :: qcincr(nxndg,nyndg,nzndg)     ! Analysis increment for qc
  REAL :: qrincr(nxndg,nyndg,nzndg)     ! Analysis increment for qr
  REAL :: qiincr(nxndg,nyndg,nzndg)     ! Analysis increment for qi
  REAL :: qsincr(nxndg,nyndg,nzndg)     ! Analysis increment for qs
  REAL :: qhincr(nxndg,nyndg,nzndg)     ! Analysis increment for qh

  REAL :: temxy1(nx,ny)           ! Temporary work array

  REAL :: tem1soil(nx,ny,nzsoil)  ! Temporary soil model work array.
  REAL :: tem2soil(nx,ny,nzsoil)  ! Temporary soil model work array.
  REAL :: tem3soil(nx,ny,nzsoil)  ! Temporary soil model work array.
  REAL :: tem4soil(nx,ny,nzsoil)  ! Temporary soil model work array.
  REAL :: tem5soil(nx,ny,nzsoil)  ! Temporary soil model work array.


  REAL :: tem1  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem5  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem6  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem7  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem8  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem9  (nx,ny,nz)     ! Temporary work array.
  REAL :: tem10 (nx,ny,nz)     ! Temporary work array.
  REAL :: tem11 (nx,ny,nz)     ! Temporary work array.
  REAL :: tem12 (nx,ny,nz)     ! Temporary work array.
  REAL :: tem13 (nx,ny,nz)     ! Temporary work array.
  REAL :: tem14 (nx,ny,nz)     ! Temporary work array.
  REAL :: tem15 (nx,ny,nz)     ! Temporary work array.
  REAL :: tem16 (nx,ny,nz)     ! Temporary work array.

  REAL :: temscalar(nx,ny,ny,MAX(10,nscalar))

  REAL :: tem1_0(0:nx,0:ny,0:nz)     ! Temporary work array.
  REAL :: tem2_0(0:nx,0:ny,0:nz)     ! Temporary work array.
  REAL :: tem3_0(0:nx,0:ny,0:nz)     ! Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'bndry.inc'
  INCLUDE 'indtflg.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'exbc.inc'
  INCLUDE 'nudging.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j
  INTEGER :: ireturn
  REAL    :: tem

                              !EMK BMJ
  LOGICAL :: restart

  INTEGER,PARAMETER :: vegwaterflag = 14
  INTEGER,PARAMETER :: xland_waterflag = 2
  INTEGER,PARAMETER :: xland_landflag = 1
                              !EMK END

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  mgrid = mptr

  grdin = 0
  basin = 0
  varin = 0
  mstin = 0
  rainin= 0
  prcin = 0
  icein = 0
  trbin = 0
  sfcin = 0
  landin= 0
  radin = 0
  flxin = 0
!wdt update - init0 no longer necessary (arrays set to 0 after allocation)
!
!-----------------------------------------------------------------------
!
!  INITialize model array VARiables.
!
!-----------------------------------------------------------------------
!
  CALL initgrdvar(nx,ny,nz,nzsoil,nt,nstyps,exbcbufsz,                  &
                  x,y,z,zp,zpsoil,hterain,mapfct,                       &
                  j1,j2,j3,j3soil,aj3x,aj3y,aj3z,j3inv,j3soilinv,       &
                  u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,                &
                  udteb, udtwb, vdtnb, vdtsb,                           &
                  pdteb,pdtwb ,pdtnb ,pdtsb,                            &
                  trigs1,trigs2,ifax1,ifax2,                            &
                  wsave1,wsave2,vwork1,vwork2,                          &
                  ubar,vbar,wbar,ptbar,pbar,ptbari,pbari,               &
                  rhostr,rhostri,qvbar,ppi,csndsq,                      &
                  soiltyp,stypfrct,vegtyp,lai,roufns,veg,               &
                  tsoil,qsoil,wetcanp,snowdpth,qvsfc,                   &
                  ptcumsrc,qcumsrc,w0avg,nca,kfraincv,                  &
                  cldefi,xland,bmjraincv,                               &
                  raing,rainc,prcrate,exbcbuf,                          &
                  radfrc,radsw,rnflx,radswnet,radlwin,                  &
                  usflx,vsflx,ptsflx,qvsflx,                            &
                  tem1soil,tem2soil,tem3soil,tem4soil,tem5soil,         &
                  tem1,tem2,tem3,tem4,tem5,                             &
                  tem6,tem7,tem8,tem9)

  DO j=1,ny-1
    DO i=1,nx-1
      tem = 0.5 * ( pprt(i,j,1,2)+pbar(i,j,1)                           &
                  + pprt(i,j,2,2)+pbar(i,j,2) )
      ptsfc(i,j)=tsoil(i,j,1,0)*(p0/tem)**rddcp
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Calculate the sin of the lattitude of each grid point, to be used
!  in the calculation of latitude-dependent Coriolis parameters.
!
!-----------------------------------------------------------------------
!
  CALL gtsinlat(nx,ny,x,y, sinlat, tem1,tem2, tem3)
!
!-----------------------------------------------------------------------
!
!  Initialize arrays that store the lookup table data.
!
!-----------------------------------------------------------------------
!

  CALL initlktb

!
!-----------------------------------------------------------------------
!
!  Initialize the external boundary data array.
!
!-----------------------------------------------------------------------
!
  IF( lbcopt == 2 .AND. mptr == 1 ) THEN

    ireturn = 0

    !  DBW question why not soil model variables as well????

    CALL extbdtini(nx,ny,nz,                                            &
                   u,v,w,ptprt,pprt,                                    &
                   qv,qscalar(:,:,:,1,:),ptbar,pbar,                    &
                   exbcbuf(nu0exb),  exbcbuf(nv0exb),                   &
                   exbcbuf(nw0exb),  exbcbuf(npt0exb),                  &
                   exbcbuf(npr0exb), exbcbuf(nqv0exb),                  &
                   exbcbuf(nqscalar0exb(1)),                            &
                   exbcbuf(nudtexb), exbcbuf(nvdtexb),                  &
                   exbcbuf(nwdtexb), exbcbuf(nptdtexb),                 &
                   exbcbuf(nprdtexb),exbcbuf(nqvdtexb),                 &
                   exbcbuf(nqscalardtexb(1)),                           &
                   bcrlx,                                               &
                   tem1,tem2,tem3,tem4,tem5,tem6,                       &
                   temscalar,tem7,ireturn)

    IF( ireturn == 1 ) THEN
      WRITE (6,'(a/a)')                                                 &
          'Can not find the external boundary data. Dump the',          &
          'history file and restart file and then STOP the model.'
      CALL arpsstop('arpsstop called from initial with ext boundary file',1)
    ELSE IF( ireturn == 2 ) THEN
      WRITE (6,'(a/a)')                                                 &
          'Can not open the external boundary data. Dump the history',  &
          'file and restart file and then STOP the model.'
      CALL arpsstop('arpsstop called from initial with opening ext      &
            & boundary file ',1)
    ELSE IF( ireturn == 3 ) THEN
      WRITE (6,'(a/a)')                                                 &
          'Read errors in the external boundary data file. Dump the',   &
          'history file and restart file and then STOP the model.'
      CALL arpsstop('arpsstop called from initial with reading ext      &
            & boundary file ',1)
    ELSE IF( ireturn /= 0 ) THEN
      WRITE (6,'(a/a)')                                                 &
          'Other errors in getting the external boundary data. Dump the', &
          'history file and restart file and then STOP the model.'
      CALL arpsstop('arpsstop called from initial with the ext      &
            & boundary file ',1)
    END IF

  END IF

  IF( lbcopt == 1 .AND. mptr == 1 .AND. rbc_plbc == 2) THEN
    CALL extbdtini_rbc(nx,ny,nz, bcrlx)
  END IF

!
!-----------------------------------------------------------------------
!
!  Initialize nudging assimilation variables
!
!-----------------------------------------------------------------------
!
  IF( nudgopt > 0 )                                                     &
    CALL ininudge(nxndg,nyndg,nzndg,                                    &
                  uincr,vincr,wincr,pincr,ptincr,qvincr,                &
                  qcincr,qrincr,qiincr,qsincr,qhincr,ireturn)

!-----------------------------------------------------------------------
!
! Initialize KF_ETA arrays and look-up tables.
!
!-----------------------------------------------------------------------

  IF (cnvctopt == 5) THEN

    !...Now initialize KF_ETA look-up tables

    IF (initopt == 2 .or. initopt == 4) THEN
      restart = .TRUE.
    ELSE
      restart = .FALSE.
    END IF

    CALL interface_wrf_kfinit(nx,ny,nz,nca,restart)
  END IF ! IF (cnvctopt == 5) THEN

!-----------------------------------------------------------------------
!
! Initialize BMJ arrays and look-up tables.
!
!-----------------------------------------------------------------------

  ! Since xland is already declared, why do not make it general?
  DO j = 1,ny-1
    DO i = 1,nx-1
      IF (vegtyp(i,j) == vegwaterflag) THEN
        xland(i,j) = xland_waterflag
      ELSE
        xland(i,j) = xland_landflag
      END IF
    END DO ! DO i = 1,nx
  END DO ! DO j = 1,ny

  IF (cnvctopt == 4) THEN

    !...Now initialize BMJ look-up tables

    IF (initopt == 2 .or. initopt == 4) THEN
      restart = .TRUE.
    ELSE
      restart = .FALSE.
    END IF

    CALL interface_wrf_bmjinit(nx,ny,nz,cldefi,restart)
  END IF ! IF (cnvctopt == 4) THEN

!
!-----------------------------------------------------------------------
!
!  End of model initialization.
!
!-----------------------------------------------------------------------
!
  RETURN
END SUBROUTINE initial
