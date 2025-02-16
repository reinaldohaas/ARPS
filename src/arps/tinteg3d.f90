!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CORDINTG                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE cordintg(mptr,frstep, nx,ny,nz,nzsoil,nxndg,nyndg,nzndg,     &
           rbufsz,nstyps,exbcbufsz,                                     &
           u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,pbldpth,               &
           udteb,udtwb,udtnb,udtsb,vdteb,vdtwb,vdtnb,vdtsb,             &
           wdteb,wdtwb,wdtnb,wdtsb,pdteb,pdtwb,pdtnb,pdtsb,             &
           sdteb,sdtwb,sdtnb,sdtsb,                                     &
           ubar,vbar,ptbar,pbar,ptbari,pbari,                           &
           rhostr,rhostri,qvbar,ppi,csndsq,                             &
           x,y,z,zp,zpsoil,mapfct,j1,j2,j3,j3soil,aj3x,aj3y,aj3z,j3inv, &
           j3soilinv,trigs1,trigs2,ifax1,ifax2,                         &
           wsave1,wsave2,vwork1,vwork2,                                 &
           sinlat,kmh,kmv,rprntl,                                       &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,ptsfc,qvsfc,                    &
           ptcumsrc,qcumsrc,raing,rainc,prcrate,w0avg,nca,kfraincv,     &
           cldefi,xland,bmjraincv,                                      &
           radfrc,radsw,rnflx,radswnet,radlwin,rad2d,radbuf,            &
           exbcbuf,bcrlx,usflx,vsflx,ptsflx,qvsflx,                     &
           uincr,vincr,wincr,pincr,ptincr,qvincr,                       &
           qcincr,qrincr,qiincr,qsincr,qhincr,                          &
           tem1soil,tem2soil,tem3soil,tem4soil,tsdiffus,                &
           phydro,tem1,tem2,tem3,tem4,tem5,tem6,tem7,                   &
           tem8,tem9,tem10,tem11,tem12,tem13,                           &
           tem14,tem15,tem16,temscalar,                                 &
           tem1_0,tem2_0,tem3_0)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Coordinate the forward integration of all model time-dependent
!  variables.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/01/1991.
!
!  MODIFICATION HISTORY:
!
!  4/26/92 (K. Droegemeier)
!  Added full documentation.
!
!  4/10/93 (M. Xue & Hao Jin)
!  Add the terrain.
!
!  01/24/94 (Yuhe Liu)
!  Added array tsfc, qvsfc, tsoil, wetsfc, wetdp, wetcanp, and some of
!  temporary arrays to the arguments of subroutine SFCFLX, for
!  prediction of the surface temperature and specific humidity
!
!  9/10/94 (D. Weber & Y. Lu)
!  Cleaned up documentation.
!
!  02/07/1995 (Yuhe Liu)
!  Added a new 2-D array, veg(nx,ny), to the argument list
!
!  7/6/1995 (M. Xue)
!  Pressure detrending is applied to the base grid (as opposed to
!  nested grids) only.
!
!  10/31/95 (D. Weber)
!  Added trigs1,trigs2,ifax1,ifax2 for use in the upper w-p
!  radiation condition.
!
!  01/25/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor to ARPS governing equations.
!
!  01/31/1996 (V. Wong and X. Song)
!  Added a parameter, qpfgfrq, that controls the frequency of calling
!  the subroutine QPFGRID.
!
!  08/14/1996 (Yuhe Liu)
!  Moved the definition of radiation buffer to the ARPS main program
!  and passed it through argument list.
!
!  07/22/97 (D. Weber)
!  Added  wsave1,wsave2,vwork1,vwork2 for use in the upper w-p
!  radiation condition (fftopt=2).
!
!  08/01/97 (Zonghui Huo)
!  Added Kain-fritsch cumulus parameterization scheme.
!
!  10/21/97 (Donghai Wang)
!  Using total density (rho) in the calculation of the pressure
!  gradient force terms, and added the second order terms
!  in the linerized buoyancy terms.
!
!  11/05/97 (D. Weber)
!  Added phydro array for use in the bottom boundary condition for
!  perturbation pressure (hydrostatic).
!
!  11/06/97 (D. Weber)
!  Added three additional levels to the mapfct array.  The three
!  levels (4,5,6) represent the inverse of the first three in order.
!  The inverse map factors are computed to improve efficiency.
!
!  2/11/1998 (M. Xue)
!  Corrections made to ensure no moist process is activatived when
!  moist=0.
!
!  4/15/1998 (Donghai Wang)
!  Added the source terms to the right hand terms of the qc,qr,qi,qs
!  equations due to the K-F cumulus parameterization.
!
!  4/15/1998 (Donghai Wang)
!  Added the running average vertical velocity (array w0avg)
!  for the K-F cumulus parameterization scheme.
!
!  9/18/98 (D. Weber)
!  Added precomputed averages of j3 in the x,y, and z directions
!  to improve code efficiency.
!
!  8/31/1998 (K. Brewster)
!  Added call to NUDGEALL for nudging to observation increments.
!
!  11/18/98 (Keith Brewster)
!  Changed pibar to ppi (full pi).
!
!  12/09/1998 (Donghai Wang)
!  Added the snow cover.
!
!  12/05/2001 (Yunheng Wang)
!  Replaced sfcflx call with sfcphysics call which was a structure
!  change by M. Xue.
!
!  07/10/2001 (K. Brewster)
!  Added increment arrays to argument list and to call to nudgeall
!
!  07/23/2001 (K. Brewster)
!  Added mptr to call to tinteg, needed for printing of diagnostic noise
!  parameter.
!
!  03/13/2002 (Eric Kemp)
!  Added arrays for WRF BMJ cumulus scheme.
!
!  05/13/2002  (J. Brotzge)
!  Added soil scheme option; removed tsfc, wetsfc, wetdp, and
!  replaced with tsoil and qsoil variables
!
!  04/28/2005 (Nate Snook)
!  Added microphysics option #4 (Straka implementation of Lin et. al.
!  ice microphysics), updated input variables for that subroutine.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    mptr     Grid identifier.
!    frstep   Flag to determine if this is the initial time step.
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of grid points in the soil
!
!    u        x component of velocity at times tpast and tpresent (m/s)
!    v        y component of velocity at times tpast and tpresent (m/s)
!    w        Vertical component of Cartesian velocity at times
!             tpast and tpresent (m/s)
!    wcont    Contravariant vertical velocity (m/s)
!    ptprt    Perturbation potential temperature at times tpast and
!             tpresent (K)
!    pprt     Perturbation pressure at times tpast and tpresent (Pascal)
!    qv       Water vapor specific humidity at times tpast and tpresent (kg/kg)
!    qc       Cloud water mixing ratio at times tpast and tpresent (kg/kg)
!    qr       Rainwater mixing ratio at times tpast and tpresent (kg/kg)
!    qi       Cloud ice mixing ratio at times tpast and tpresent (kg/kg)
!    qs       Snow mixing ratio at times tpast and tpresent (kg/kg)
!    qh       Hail mixing ratio at times tpast and tpresent (kg/kg)
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!    pbldpth  Planetary boundary layer depth (m)
!
!    udteb    Time tendency of u field at east boundary (m/s**2)
!    udtwb    Time tendency of u field at west boundary (m/s**2)
!
!    vdtnb    Time tendency of v field at north boundary (m/s**2)
!    vdtsb    Time tendency of v field at south boundary (m/s**2)
!
!    pdteb    Time tendency of pprt field at east boundary (Pascal/s)
!    pdtwb    Time tendency of pprt field at west boundary (Pascal/s)
!    pdtnb    Time tendency of pprt field at north boundary (Pascal/s)
!    pdtsb    Time tendency of pprt field at south boundary (Pascal/s)
!
!    phydro   Big time step forcing term for use in computing the
!             hydrostatic pressure at k=1.
!
!    ubar     Base state zonal velocity component (m/s)
!    vbar     Base state meridional velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    ptbari   Inverse Base state potential temperature (K)
!    pbari    Inverse Base state pressure (Pascal)
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    rhostri  Inverse base state density rhobar times j3 (kg/m**3)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!    ppi      Exner function.
!    csndsq   Sound wave speed squared.
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!    zpsoil   Vertical coordinate of grid points in the soil (m)
!    mapfct   Map factors at scalar, u and v points
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    j3soil   Coordinate transformation Jacobian  d(zpsoil)/dz
!    aj3x     Avgx of the coordinate transformation Jacobian  d(zp)/dz
!    aj3y     Avgy of the coordinate transformation Jacobian  d(zp)/dz
!    aj3z     Avgz of the coordinate transformation Jacobian  d(zp)/dz
!    j3inv    Inverse of the coordinate transformation Jacobian  d(zp)/dz
!    j3soilinv Inverse of the coordinate transformation Jacobian  d(zpsoil)/dz
!
!    trigs1   Array containing pre-computed trig function for
!             fftopt=1.
!    trigs2   Array containing pre-computed trig function for
!             fftopt=1.
!    ifax1    Array containing the factors of nx for fftopt=1.
!    ifax2    Array containing the factors of ny for fftopt=1.
!
!    vwork1   2-D work array for fftopt=2.
!    vwork2   2-D work array for fftopt=2.
!    wsave1   Work array for fftopt=2.
!    wsave2   Work array for fftopt=2.
!
!    soiltyp  Soil type
!    vegtyp   Vegetation type
!    lai      Leaf Area Index
!    roufns   Surface roughness
!    veg      Vegetation fraction
!
!    tsoil    Soil temperature (K)
!    qvsfc    Effective qv at sfc.
!    qsoil    Soil moisture (m**3/m**3)
!    wetcanp  Canopy water amount
!
!    ptsfc    Ground surface potential temperature (K)
!    qvsfc    Effective S.H. at sfc.
!
!  OUTPUT:
!    u        Updated (by dtbig) x velocity component at times tpast
!             and tpresent (m/s).
!    v        Updated (by dtbig) y velocity component at times tpast
!             and tpresent (m/s).
!    w        Updated (by dtbig) z velocity component at times tpast
!             and tpresent (m/s).
!    ptprt    Updated (by dtbig) perturbation potential temperature
!             at times tpast and tpresent (K).
!    pprt     Updated (by dtbig) perturbation pressure at times tpast
!             and tpresent (Pascal).
!    qv       Updated (by dtbig) water vapor specific humidity at times
!             tpast and tpresent (kg/kg).
!    qc       Updated (by dtbig) cloud water mixing ratio at times
!             tpast and tpresent (kg/kg).
!    qr       Updated (by dtbig) rainwater mixing ratio at times tpast
!             and tpresent (kg/kg).
!    qi       Updated (by dtbig) cloud ice mixing ratio at times tpast
!             and tpresent (kg/kg).
!    qs       Updated (by dtbig) snow mixing ratio at times tpast and
!             tpresent (kg/kg).
!    qh       Updated (by dtbig) hail mixing ratio at times tpast and
!             tpresent (kg/kg).
!
!    udteb    Time tendency of u field at east boundary (m/s**2)
!    udtwb    Time tendency of u field at west boundary (m/s**2)
!
!    vdtnb    Time tendency of v field at north boundary (m/s**2)
!    vdtsb    Time tendency of v field at south boundary (m/s**2)
!
!    pdteb    Time tendency of pprt field at east boundary (Pascal/s)
!    pdtwb    Time tendency of pprt field at west boundary (Pascal/s)
!    pdtnb    Time tendency of pprt field at north boundary (Pascal/s)
!    pdtsb    Time tendency of pprt field at south boundary (Pascal/s)
!
!    sinlat   Sin of latitude at each grid point
!
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!    rprntl   Reciprocal of Prandtl number
!
!    tsfc     Soil temperature (K)
!    qvsfc    Effective qv at sfc.
!    qsoil    Soil moisture (m**3/m**3)
!    wetcanp  Canopy water amount
!
!    ptsfc    Ground surface potential temperature (K)
!    qvsfc    Effective S.H. at sfc.
!
!    ptcumsrc Source term in pt-equation due to cumulus parameterization
!    qcumsrc  Source term in water equations due to cumulus parameterization
!    raing    Grid supersaturation rain
!    rainc    Cumulus convective rain
!    kfraincv   K-F convective rainfall (cm)
!    nca      K-F counter for CAPE release
!    cldefi   BMJ cloud efficiency
!    xland    BMJ land/sea mask
!    bmjraincv   BMJ convective rainfall (cm)
!
!    radfrc   Radiation forcing (K/s)
!    radsw    Solar radiation reaching the ground surface (W/m**2)
!    rnflx    Net radiation flux at the ground surface (W/m**2)
!    radswnet Net shortwave radiation
!    radlwin  Incoming longwave radiation
!    rad2d    2-D arrays for radiation calculation
!    radbuf   Buffer array to carry temporary working arrays for
!             radiation calculation
!
!  WORK ARRAYS:
!
!    udtnb    Time tendency of u field at north boundary (m/s**2)
!    udtsb    Time tendency of u field at south boundary (m/s**2)
!
!    vdteb    Time tendency of v field at east boundary (m/s**2)
!    vdtwb    Time tendency of v field at west boundary (m/s**2)
!
!    wdteb    Time tendency of w field at east boundary (m/s**2)
!    wdtwb    Time tendency of w field at west boundary (m/s**2)
!    wdtnb    Time tendency of w field at north boundary (m/s**2)
!    wdtsb    Time tendency of w field at south boundary (m/s**2)
!
!    ptdteb   Time tendency of ptprt field at east boundary (K/s)
!    ptdtwb   Time tendency of ptprt field at west boundary (K/s)
!    ptdtnb   Time tendency of ptprt field at north boundary(K/s)
!    ptdtsb   Time tendency of ptprt field at south boundary(K/s)
!
!    sdteb    Time tendency of a scalar at e-boundary
!    sdtwb    Time tendency of a scalar at w-boundary
!    sdtnb    Time tendency of a scalar at n-boundary
!    sdtsb    Time tendency of a scalar at s-boundary
!
!    usflx    Surface flux of u-momentum (kg/(m*s**2))
!    vsflx    Surface flux of v-momentum (kg/(m*s**2))
!    ptsflx   Surface heat flux (K*kg/(m**2 * s ))
!    qvsflx   Surface moisture flux of (kg/(m**2 * s))
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

  IMPLICIT NONE             ! Force explicit declarations
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'     ! Global constants that control model execution
  INCLUDE 'bndry.inc'       ! Boundary condition control parameters
  INCLUDE 'phycst.inc'      ! Physical constants
  INCLUDE 'exbc.inc'        ! EXBC constants and parameters
  INCLUDE 'nudging.inc'     ! variables for nudging assimilation
  INCLUDE 'mp.inc'          ! Message passing parameters.
  INCLUDE 'timelvls.inc'
  INCLUDE 'radcst.inc'      ! Radiation constants and parameters


!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  INTEGER :: mptr              ! Grid identifier.

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of grid points in 3 directions

  INTEGER :: nxndg,nyndg,nzndg ! Number of grid points in 3 directions

  INTEGER :: frstep            ! Indicator of the first time step call

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)
  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature
                               ! from that of base state atmosphere (K)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure from that
                               ! of base state atmosphere (Pascal)
  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)

  REAL :: qscalar(nx,ny,nz,nt,nscalar)

  REAL :: tke   (nx,ny,nz,nt)  ! Turbulent Kinetic Energy ((m/s)**2)
  REAL :: pbldpth(nx,ny,nt)    ! Planetary boundary layer depth (m)

  REAL :: udteb (ny,nz)        ! Time tendency of u at e-boundary (m/s**2)
  REAL :: udtwb (ny,nz)        ! Time tendency of u at w-boundary (m/s**2)
  REAL :: udtnb (nx,nz)        ! Time tendency of u at n-boundary (m/s**2)
  REAL :: udtsb (nx,nz)        ! Time tendency of u at s-boundary (m/s**2)

  REAL :: vdteb (ny,nz)        ! Time tendency of v at e-boundary (m/s**2)
  REAL :: vdtwb (ny,nz)        ! Time tendency of v at w-boundary (m/s**2)
  REAL :: vdtnb (nx,nz)        ! Time tendency of v at n-boundary (m/s**2)
  REAL :: vdtsb (nx,nz)        ! Time tendency of v at s-boundary (m/s**2)

  REAL :: wdteb (ny,nz)        ! Time tendency of w at e-boundary (m/s**2)
  REAL :: wdtwb (ny,nz)        ! Time tendency of w at w-boundary (m/s**2)
  REAL :: wdtnb (nx,nz)        ! Time tendency of w at n-boundary (m/s**2)
  REAL :: wdtsb (nx,nz)        ! Time tendency of w at s-boundary (m/s**2)

  REAL :: pdteb (ny,nz)        ! Time tendency of pprt at e-boundary (Pascal/s)
  REAL :: pdtwb (ny,nz)        ! Time tendency of pprt at w-boundary (Pascal/s)
  REAL :: pdtnb (nx,nz)        ! Time tendency of pprt at n-boundary (Pascal/s)
  REAL :: pdtsb (nx,nz)        ! Time tendency of pprt at s-boundary (Pascal/s)

  REAL :: phydro(nx,ny)        ! Big time step forcing for computing
                               ! hydrostatic pprt at k=1.

  REAL :: sdteb (ny,nz)        ! Time tendency of a scalar at e-boundary
  REAL :: sdtwb (ny,nz)        ! Time tendency of a scalar at w-boundary
  REAL :: sdtnb (nx,nz)        ! Time tendency of a scalar at n-boundary
  REAL :: sdtsb (nx,nz)        ! Time tendency of a scalar at s-boundary

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: ptbari(nx,ny,nz)     ! Inverse Base state pot. temperature (K)
  REAL :: pbari (nx,ny,nz)     ! Inverse Base state pressure (Pascal).
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: rhostri(nx,ny,nz)    ! Inverse base state density rhobar times j3.
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific
                               ! humidity (kg/kg)
  REAL :: ppi   (nx,ny,nz)     ! Exner function.
  REAL :: csndsq(nx,ny,nz)     ! Sound wave speed squared.

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of staggered grid.
  REAL :: zpsoil (nx,ny,nzsoil)      ! The physical height coordinate defined at
                               ! w-point of the soil.

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian -d(zp)/d(x)
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian -d(zp)/d(y)
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian  d(zp)/d(z)
  REAL :: j3soil(nx,ny,nzsoil)       ! Coordinate transformation Jacobian  d(zp)/d(z)

  REAL :: aj3x  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE X-DIR.
  REAL :: aj3y  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
  REAL :: aj3z  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Z-DIR.

  REAL :: j3inv (nx,ny,nz)     ! Coordinate transformation Jacobian  d(zp)/d(z)
  REAL :: j3soilinv (nx,ny,nzsoil)   ! Coordinate transformation Jacobian  d(zpsoil)/d(z)

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

  INTEGER :: nstyps                 ! Number of soil types
  INTEGER :: soiltyp(nx,ny,nstyps)  ! Soil type at each point
  REAL    :: stypfrct(nx,ny,nstyps) ! Fraction of soil types
  INTEGER :: vegtyp (nx,ny)         ! Vegetation type at each point

  REAL :: lai    (nx,ny)          ! Leaf Area Index
  REAL :: roufns (nx,ny)          ! Surface roughness
  REAL :: veg    (nx,ny)          ! Vegetation fraction

  REAL :: tsoil (nx,ny,nzsoil,0:nstyps)! Soil temperature(K)
  REAL :: qsoil (nx,ny,nzsoil,0:nstyps) ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny,0:nstyps) ! Canopy water amount
  REAL :: snowdpth(nx,ny)         ! Snow depth (m)

  REAL :: ptsfc  (nx,ny)          ! Potential temperature at the ground level (K)
  REAL :: qvsfc  (nx,ny,0:nstyps) ! Effective qv at the surface (kg/kg)


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
  REAL :: kfraincv (nx,ny)     ! K-F convective rainfall (cm)
  INTEGER :: nca (nx,ny)       ! K-F counter for CAPE release

  REAL,INTENT(INOUT) :: cldefi(nx,ny) ! BMJ cloud efficiency
  REAL,INTENT(IN)    :: xland(nx,ny)  ! BMJ land mask
                                      ! (1.0 = land, 2.0 = sea)
  REAL,INTENT(INOUT) :: bmjraincv (nx,ny) ! BMJ convective rainfall
                                          ! (cm)

  REAL :: raing  (nx,ny)       ! Gridscale rainfall (mm)
  REAL :: rainc  (nx,ny)       ! Cumulus rainfall (mm)
  REAL :: prcrate(nx,ny,4)     ! precipitation rate (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precipitation rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulus precip. rate
                               ! prcrate(1,1,4) = microphysics precip. rate
!
!-----------------------------------------------------------------------
!
!  Arrays for radiation
!
!-----------------------------------------------------------------------
!
  INTEGER :: rbufsz            ! Radiation buffer size
  REAL :: radfrc(nx,ny,nz)     ! Radiation forcing (K/s)
  REAL :: radsw (nx,ny)        ! Solar radiation reaching the surface
  REAL :: radswnet(nx,ny)      ! Net solar radiation, (SWin - SWout)
  REAL :: radlwin(nx,ny)       ! Incoming longwave radiation
  REAL :: rnflx (nx,ny)        ! Net absorbed radiation flux at surface


  REAL :: rad2d(nx,ny,nrad2d)
    ! Buffur array to carry the variables calculated or used in
    ! radiation calculation. The last index defines the variables:
    ! 1  = nrsirbm,  Solar IR surface albedo for beam
    ! 2  = nrsirdf,  Solar IR surface albedo for diffuse
    ! 3  = nrsuvbm,  Solar UV surface albedo for beam
    ! 4  = nrsuvdf,  Solar UV surface albedo for diffuse
    ! 5  = ncosz,    Cosine of zenith
    ! 6  = ncosss,   Cosine of angle between sun light and slope
    ! 7  = nfdirir,  all-sky direct downward IR flux (0.7-10 micron)
    !                at the surface
    ! 8  = nfdifir,  all-sky diffuse downward IR flux
    !                at the surface
    ! 9  = nfdirpar, all-sky direct downward and par flux
    !                (0.4-0.7 micron) at the surface
    ! 10 = nfdifpar, all-sky diffuse downward and par flux
    !                at thediffuse downward and par flux
    !                at the surface

  REAL :: radbuf(rbufsz)       ! Buffer array storing temporary working
                               ! arrays for radiation computing

  INTEGER :: exbcbufsz         ! EXBC buffer size
  REAL :: exbcbuf( exbcbufsz ) ! EXBC buffer array
  REAL :: bcrlx(nx,ny)         ! EXBC relaxation coefficients

  REAL :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m**2*s))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

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
!
!-----------------------------------------------------------------------
!
!  Temporary arrays
!
!-----------------------------------------------------------------------
!
  REAL :: cdm   (nx,ny)        ! Drag coefficient for surface momentum flux
  REAL :: cdh   (nx,ny)        ! Drag coefficient for surface heat flux
  REAL :: cdq   (nx,ny)        ! Drag coefficient for surface moisture flux

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

  REAL, TARGET :: temscalar(nx,ny,nz,MAX(10,nscalar))

  REAL, POINTER :: tem17 (:,:,:)     ! Temporary work array.
  REAL, POINTER :: tem18 (:,:,:)     ! Temporary work array.
  REAL, POINTER :: tem19 (:,:,:)     ! Temporary work array.
  REAL, POINTER :: tem20 (:,:,:)     ! Temporary work array.
  REAL, POINTER :: tem21 (:,:,:)     ! Temporary work array.
  REAL, POINTER :: tem22 (:,:,:)     ! Temporary work array.
  REAL, POINTER :: tem23 (:,:,:)     ! Temporary work array.
  REAL, POINTER :: tem24 (:,:,:)     ! Temporary work array.
  REAL, POINTER :: tem25 (:,:,:)     ! Temporary work array.
  REAL, POINTER :: tem26 (:,:,:)     ! Temporary work array.

  REAL, POINTER :: temscalar1 (:,:,:)
  REAL, POINTER :: temscalar2 (:,:,:)
  REAL, POINTER :: temscalar3 (:,:,:)

  REAL :: tem1_0(0:nx,0:ny,0:nz)     ! Temporary work array.
  REAL :: tem2_0(0:nx,0:ny,0:nz)     ! Temporary work array.
  REAL :: tem3_0(0:nx,0:ny,0:nz)     ! Temporary work array.

  REAL :: tem1soil(nx,ny,nzsoil) ! Temporary work array.
  REAL :: tem2soil(nx,ny,nzsoil) ! Temporary work array.
  REAL :: tem3soil(nx,ny,nzsoil) ! Temporary work array.
  REAL :: tem4soil(nx,ny,nzsoil) ! Temporary work array.
  REAL :: tsdiffus(nx,ny,nzsoil) ! Temporary work array.

!
!-----------------------------------------------------------------------
!
!  Miscellaneous local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nocalls
  SAVE nocalls
  DATA nocalls / 0 /
  INTEGER :: i,j,k
  INTEGER :: nq
  INTEGER :: ireturn           ! Return status indicator
  REAL    :: dtbig1
  REAL    :: pisum
  INTEGER :: tstrtlvl

  CHARACTER :: ayear*4, amm*2, aday*2

  INTEGER :: mphyflg, mscheme

  REAL, ALLOCATABLE :: ptbar_z(:)
  REAL, ALLOCATABLE :: qvbar_z(:)
  REAL, ALLOCATABLE :: dzc_z(:)

  ! Added by DTD for additional microphysics array output

  CHARACTER (LEN=40 ) :: varunits
  CHARACTER (LEN=40 ) :: varname
  CHARACTER (LEN=6 ) :: varid

  CHARACTER (LEN=40) :: mpteqnunits(28)
  CHARACTER (LEN=40) :: mpteqnname(28)
  CHARACTER (LEN=40) :: mpteqnid(28)

  REAL, ALLOCATABLE :: mpteqnterms(:,:,:,:)
  REAL, ALLOCATABLE :: mptrate(:,:,:)
  REAL(KIND=8), ALLOCATABLE :: N0x(:,:,:,:)
  REAL, ALLOCATABLE :: N0x_sng(:,:,:,:)

  INTEGER :: nqscalar, istatus

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  tem17 => temscalar(:,:,:,1)
  tem18 => temscalar(:,:,:,2)
  tem19 => temscalar(:,:,:,3)
  tem20 => temscalar(:,:,:,4)
  tem21 => temscalar(:,:,:,5)
  tem22 => temscalar(:,:,:,6)
  tem23 => temscalar(:,:,:,7)
  tem24 => temscalar(:,:,:,8)
  tem25 => temscalar(:,:,:,9)
  tem26 => temscalar(:,:,:,10)

  mgrid = mptr
!
!-----------------------------------------------------------------------
!
!  Calculate the radiation forcing and surface fluxes if desired.
!
!-----------------------------------------------------------------------
!
  CALL set_acct(rad_acct)

  IF ( radopt /= 0 .AND. MOD(nstep-1, nradstp) == 0 ) THEN

    IF (myproc == 0)                &
      WRITE(6,'(a,i8,a,f10.2,a)')                                       &
                   ' Compute radiation at time step,', nstep,           &
                   ', model time=',curtim,' (s)'

    IF (radopt /= 3) THEN

      CALL radiation(nx,ny,nz,rbufsz,                                   &
                     ptprt(1,1,1,1),pprt(1,1,1,1),                      &
                     qv(1,1,1,1),qscalar(:,:,:,1,:),                    &
                     ptbar,pbar,ppi,rhostr,                             &
                     x,y,z,zp, j3inv,                                   &
                     soiltyp(1,1,1),tsoil(1,1,1,0),qsoil(1,1,1,0),      &
                     snowdpth,radfrc,radsw,rnflx,radswnet,radlwin,      &
                     rad2d(1,1,nrsirbm), rad2d(1,1,nrsirdf),            &
                     rad2d(1,1,nrsuvbm), rad2d(1,1,nrsuvdf),            &
                     rad2d(1,1,ncosz),   rad2d(1,1,ncosss),             &
                     rad2d(1,1,nfdirir), rad2d(1,1,nfdifir),            &
                     rad2d(1,1,nfdirpar),rad2d(1,1,nfdifpar),           &
                     tem1, tem2, tem3, tem4, tem5,                      &
                     tem6, tem7, tem8, tem9, tem10,                     &
                     tem11,tem12,tem13,tem14,tem15,tem16,               &
                     radbuf, tem17, tem18)

!   Note: Have not defined radswnet and radlwin explicitly for
!           radopt /= 1!

    ELSE IF (radopt == 3) THEN         !read Mesonet radiation data file (JAB)

      CALL initztime(ayear,amm,aday)

      CALL readjradd(nx,ny,nz,ayear,amm,aday,ztime,                    &
                     radsw,radswnet,radlwin,rnflx)

    ENDIF                              ! Radopt = 3 test block

    IF( lvldbg >= 2 ) THEN
      CALL checkshx(radfrc, nx,ny,nz,1,nx-1,1,ny-1,2,nz-2,              &
                    'radfrc', tem1)
      CALL checkshy(radfrc, nx,ny,nz,1,nx-1,1,ny-1,2,nz-2,              &
                    'radfrc', tem1)
      CALL checkshx(radsw,nx,ny,1,1,nx-1,1,ny-1,1,1,                    &
                    'radsw', tem1)
      CALL checkshy(radsw,nx,ny,1,1,nx-1,1,ny-1,1,1,                    &
                    'radsw', tem1)
      CALL checkshx(rnflx,nx,ny,1,1,nx-1,1,ny-1,1,1,                    &
                    'rnflx', tem1)
      CALL checkshy(rnflx,nx,ny,1,1,nx-1,1,ny-1,1,1,                    &
                    'rnflx', tem1)
    END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate the surface fluxes of momentum, heat, and water
!  vapor for use in the turbulent mixing subroutine
!
!-----------------------------------------------------------------------
!
  CALL set_acct(misc_acct)

  IF(sadvopt /= 4) THEN                   ! Leapfrop scheme
    tstrtlvl = tpast
  ELSE                                    ! Forward-based scheme
    tstrtlvl = tpresent
  END IF

  CALL sfcphysics(nx,ny,nz,nzsoil,nstyps,                               &
                  u(1,1,1,tstrtlvl),v(1,1,1,tstrtlvl),                  &
                  w(1,1,1,tstrtlvl),ptprt(1,1,1,tstrtlvl),              &
                  pprt(1,1,1,tstrtlvl),qv(1,1,1,tstrtlvl),              &
                  rhostr,ptbar,pbar,qvbar,                              &
                  x,y,z, zp,zpsoil, j1,j2,j3,j3soil,j3inv,              &
                  j3soilinv, prcrate(1,1,1),                            &
                  soiltyp,stypfrct,vegtyp,lai,roufns,veg,               &
                  tsoil,qsoil,wetcanp,snowdpth,qvsfc,                   &
                  radsw,rnflx,radswnet,radlwin,                         &
                  cdm,cdh,cdq,usflx,vsflx,ptsflx,qvsflx,pbldpth,        &
                  tsdiffus,tem1soil,tem2soil,tem3soil,tem4soil,         &
                  tem1,tem2,tem3,tem4,tem5,                             &
                  tem6,tem7,tem8,tem9,tem10,tem11,tem12,tem13,          &
                  tem14,tem15,tem16,tem17,tem18,tem19,tem20)

  CALL set_acct(bc_acct)
  IF( lbcopt == 2 .AND. mgrid == 1 ) THEN
!
!-----------------------------------------------------------------------
!
!  Read in variables from an external boundary file, and
!  calculate the linear time tendencies of the external data.
!
!-----------------------------------------------------------------------
!
    CALL extbdt(nx,ny,nz, ptbar,pbar, ireturn,                          &
                exbcbuf(nu0exb), exbcbuf(nv0exb),                       &
                exbcbuf(nw0exb), exbcbuf(npt0exb),                      &
                exbcbuf(npr0exb),exbcbuf(nqv0exb),                      &
                exbcbuf(nqscalar0exb(1)),                               &
                exbcbuf(nudtexb), exbcbuf(nvdtexb),                     &
                exbcbuf(nwdtexb), exbcbuf(nptdtexb),                    &
                exbcbuf(nprdtexb),exbcbuf(nqvdtexb),                    &
                exbcbuf(nqscalardtexb(1)),                              &
                tem1,tem2,tem3,tem4,tem5,tem6,                          &
                temscalar,tem7)

    CALL mpmaxi(ireturn)

    IF( ireturn == 0 ) THEN
      GO TO 112
    ELSE IF( ireturn == 1 ) THEN
      WRITE (6,'(a/a)')                                                 &
          'Can not find the external boundary data. Dump the',          &
          'history file and restart file and then STOP the model.'
      GO TO 111
    ELSE IF( ireturn == 2 ) THEN
      WRITE (6,'(a/a)')                                                 &
          'Can not open the external boundary data. Dump the history',  &
          'file and restart file and then STOP the model.'
      GO TO 111
    ELSE IF( ireturn == 3 ) THEN
      WRITE (6,'(a/a)')                                                 &
          'Read errors in the external boundary data file. Dump the',   &
          'history file and restart file and then STOP the model.'
      GO TO 111
    ELSE
      WRITE (6,'(a/a)')                                                 &
          'Other errors in getting the external boundary data. Dump the', &
          'history file and restart file and then STOP the model.'
      GO TO 111
    END IF

    111     CONTINUE

    CALL set_acct(output_acct)
!
!-----------------------------------------------------------------------
!
!  Write out external data arrays and stop the program if data read
!  failed.
!
!-----------------------------------------------------------------------
!
    IF (mp_opt > 0) THEN
      CALL arpsstop("arpstop called from CORDINTG  ",1)
    END IF

    CALL abortdmp(mptr,nx,ny,nz,nzsoil,nstyps,                          &
                  u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,                &
                  ubar,vbar,ptbar,pbar,rhostr,qvbar,kmh,kmv,            &
                  x,y,z,zp,zp(1,1,2),zpsoil,mapfct, j1,j2,j3,j3soil,    &
                  soiltyp,stypfrct,vegtyp,lai,roufns,veg,               &
                  tsoil,qsoil,wetcanp,snowdpth,                         &
                  raing,rainc,prcrate,                                  &
                  radfrc,radsw,rnflx,radswnet,radlwin,                  &
                  usflx,vsflx,ptsflx,qvsflx,                            &
                  tem1,tem2,tem3, tem4)

!  blocking inserted for ordering i/o for message passing
    DO i=0,nprocs-1,dumpstride
      IF(myproc >= i.AND.myproc <= i+dumpstride-1)THEN

        IF(mp_opt > 0 .AND. joindmp(FINDX_R) > 0) THEN

        CALL rstjoinout(nx,ny,nz,nzsoil,nstyps,exbcbufsz,               &
                u,v,w,ptprt,pprt,qv,qscalar,tke,                        &
                udteb, udtwb, vdtnb, vdtsb,                             &
                pdteb ,pdtwb ,pdtnb ,pdtsb,                             &
                ubar,vbar,ptbar,pbar,rhostr,qvbar,                      &
                x,y,z,zp,zpsoil,zp(1,1,2),mapfct,                       &
                soiltyp,stypfrct,vegtyp,lai,roufns,veg,                 &
                tsoil,qsoil,wetcanp,snowdpth,qvsfc,                     &
                ptcumsrc,qcumsrc,w0avg,nca,kfraincv,                    &
                cldefi,xland,bmjraincv,                                 &
                radfrc,radsw,rnflx,radswnet,radlwin,                    &
                raing,rainc,prcrate,exbcbuf,tem1)

        ELSE

        CALL rstout(nx,ny,nz,nzsoil,nstyps,exbcbufsz,                   &
                u,v,w,ptprt,pprt,qv,qscalar,tke,                        &
                udteb, udtwb, vdtnb, vdtsb,                             &
                pdteb ,pdtwb ,pdtnb ,pdtsb,                             &
                ubar,vbar,ptbar,pbar,rhostr,qvbar,                      &
                x,y,z,zp,zpsoil,zp(1,1,2),mapfct,                       &
                soiltyp,stypfrct,vegtyp,lai,roufns,veg,                 &
                tsoil,qsoil,wetcanp,snowdpth,qvsfc,                     &
                ptcumsrc,qcumsrc,w0avg,nca,kfraincv,                    &
                cldefi,xland,bmjraincv,                                 &
                radfrc,radsw,rnflx,radswnet,radlwin,                    &
                raing,rainc,prcrate,exbcbuf,tem1)
        END IF

      END IF
      IF (mp_opt > 0) CALL mpbarrier
    END DO

    CALL arpsstop("arpsstop called from CORDINTG run unstable, aborted",1)

    112     CONTINUE

  END IF
!
!-----------------------------------------------------------------------
!
!  Adjust potential temperature and Qv fields with grid condensation
!  or warm rain microphysics
!
!-----------------------------------------------------------------------
!
  CALL set_acct(qpfgrd_acct)

  IF ( (moist /= 0) .AND. (cnvctopt == 1) .AND.                         &
         (MOD(curtim+0.001,qpfgfrq) < (0.5*dtbig)) ) THEN

    CALL qpfgrid(nx,ny,nz,prcrate(1,1,2),                               &
         pprt(1,1,1,2),ptprt(1,1,1,2),qv(1,1,1,2),pbar,ptbar,           &
         rhostr,zp,j3,j3inv,raing, tem1(1,1,1),tem1(1,1,2),             &
         tem1(1,1,3),tem1(1,1,4))

    IF ( lbcopt == 2 .AND.  mgrid == 1 ) THEN
      CALL acct_interrupt(bc_acct)
      CALL exbcpt( nx,ny,nz, curtim, ptprt(1,1,1,2),                    &
                   exbcbuf(npt0exb),exbcbuf(nptdtexb) )
      CALL exbcq( nx,ny,nz, 0,curtim, qv(1,1,1,2),                      &
                   exbcbuf(nqv0exb),exbcbuf(nqvdtexb) )
      CALL acct_stop_inter
    END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  Integrate the conservation equations of momentum, heat, pressure
!  and water variables forward by one timestep
!
!-----------------------------------------------------------------------
!
  CALL set_acct(tinteg_acct)

  IF(tintegopt == 1) THEN

    CALL tinteg(mptr,frstep,nx,ny,nz,nzsoil,exbcbufsz,                  &
              u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,pbldpth,            &
              udteb,udtwb,udtnb,udtsb,vdteb,vdtwb,vdtnb,vdtsb,          &
              wdteb,wdtwb,wdtnb,wdtsb,pdteb,pdtwb,pdtnb,pdtsb,          &
              phydro,sdteb,sdtwb,sdtnb,sdtsb,                           &
              ubar,vbar,ptbar,pbar,ptbari,pbari,                        &
              rhostr,rhostri,qvbar,ppi,csndsq,                          &
              x,y,z,zp,zpsoil,mapfct,                                   &
              j1,j2,j3,j3soil,aj3x,aj3y,aj3z,j3inv,                     &
              trigs1,trigs2,ifax1,ifax2,                                &
              wsave1,wsave2,vwork1,vwork2,sinlat,                       &
              ptsfc,qvsfc(1,1,0),prcrate, radfrc,                       &
              usflx,vsflx,ptsflx,qvsflx,kmh,kmv,rprntl,                 &
              ptcumsrc,qcumsrc,raing,rainc,w0avg,nca,kfraincv,          &
              cldefi,xland,bmjraincv,roufns,                            &
              exbcbuf, bcrlx,                                           &
              tem1,tem2,tem3,tem4,tem5,tem6,tem7,                       &
              tem8,tem9,tem10,tem11,tem12,tem13,                        &
              tem14,tem15,tem16,tem17,tem18,tem19,                      &
              tem20,tem21,tem22,tem23,tem24,tem25,tem26,                &
              tem1_0,tem2_0,tem3_0)

  ELSE IF(tintegopt == 2 .OR. tintegopt == 3) THEN

    CALL tinteg_RK3(mptr,frstep,nx,ny,nz,nzsoil,exbcbufsz,              &
              u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,pbldpth,            &
              udteb,udtwb,udtnb,udtsb,vdteb,vdtwb,vdtnb,vdtsb,          &
              wdteb,wdtwb,wdtnb,wdtsb,pdteb,pdtwb,pdtnb,pdtsb,          &
              phydro,sdteb,sdtwb,sdtnb,sdtsb,                           &
              ubar,vbar,ptbar,pbar,ptbari,pbari,                        &
              rhostr,rhostri,qvbar,ppi,csndsq,                          &
              x,y,z,zp,zpsoil,mapfct,                                   &
              j1,j2,j3,j3soil,aj3x,aj3y,aj3z,j3inv,                     &
              trigs1,trigs2,ifax1,ifax2,                                &
              wsave1,wsave2,vwork1,vwork2,sinlat,                       &
              ptsfc,qvsfc(1,1,0),prcrate, radfrc,                       &
              usflx,vsflx,ptsflx,qvsflx,kmh,kmv,rprntl,                 &
              ptcumsrc,qcumsrc,raing,rainc,w0avg,nca,kfraincv,          &
              cldefi,xland,bmjraincv,roufns,                            &
              exbcbuf, bcrlx,                                           &
              tem1,tem2,tem3,tem4,tem5,tem6,tem7,                       &
              tem8,tem9,tem10,tem11,tem12,tem13,                        &
              tem14,tem15,tem16,tem17,tem18,tem19,                      &
              tem20,tem21,tem22,tem23,tem24,tem25,tem26,                &
              tem1_0,tem2_0,tem3_0)
  END IF

  CALL set_acct(misc_acct)
  IF( pdetrnd == 1 .AND. mptr == 1 ) THEN
!
!-----------------------------------------------------------------------
!
!  Detrend pressure. This should be used only when unwanted trend
!  develops in the pressure field. The precedure forces the
!  domain averaged Exner function to zero.
!
!  Pressure detrending is sometimes necessary when open lateral
!  boundary condition is used. Pressure detrending is applied
!  to the base grid (as opposed to nested grids) only.
!
!-----------------------------------------------------------------------
!
    pisum = 0.0
    DO k=2,nz-2
      DO j=2,ny-2
        DO i=2,nx-2
          pisum = pisum + ppi(i,j,k)/(pbar(i,j,k)+pprt(i,j,k,3))        &
                                    * pprt(i,j,k,3)
        END DO
      END DO
    END DO

    IF (mp_opt > 0) THEN
      CALL mptotal(pisum)
      pisum = pisum/(nproc_x*nproc_y)
    END IF

    pisum = pisum/((nx-3)*(ny-3)*(nz-3))
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          pprt(i,j,k,3) = pprt(i,j,k,3)                                 &
                        - pisum*(pbar(i,j,k)+pprt(i,j,k,3))/ppi(i,j,k)
        END DO
      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Call the microphysics routines which adjust the thermodynamic
!  and water variables.
!
!  Upon exiting this subroutine, supersaturation is removed.
!
!-----------------------------------------------------------------------
!
  IF( moist /= 0 ) THEN

    CALL set_acct(microph_acct)

    mphyflg = 0
    IF( MOD(nstep-1,nmphystp) == 0) mphyflg = 1

    IF(tintegopt == 1) THEN !Leapfrog time integration
      IF(frstep == 1) THEN
        dtbig1 = dtbig/2
      ELSE
        dtbig1 = dtbig
      END IF
    ELSE ! RK3 time integration
      dtbig1 = dtbig
    END IF

    !DTD: precalculate temperature and store in mptrate
    !This will be recalculated after the microphysics routine
    !and subtracted from the precalculated value to give
    !the temperature difference due to microphysics adjustment

    ALLOCATE(mptrate(nx,ny,nz))
    mptrate = 0.0

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
           mptrate(i,j,k) = (ptprt(i,j,k,3)+ptbar(i,j,k))*ppi(i,j,k)
!          IF (mptrate(i,j,k) /= 0.0) THEN
!            print*,'mptrate(',i,',',j,',',k,')',mptrate(i,j,k)
!          END IF
        END DO
      END DO
    END DO

    IF(mphyopt >= 0 .and. mphyopt <= 1) THEN
      nqscalar = 2
    ELSE IF(mphyopt >= 2 .and. mphyopt <= 7) THEN
      nqscalar = 5
    ELSE IF(mphyopt >= 8 .and. mphyopt <= 11) THEN
      nqscalar = 6
    ELSE
      nqscalar = 6
    END IF

    ALLOCATE(mpteqnterms(nx,ny,nz,28))
    mpteqnterms = 0.0
    ALLOCATE(N0x(nx,ny,nz,nqscalar))
    N0x = 0.0
    ALLOCATE(N0x_sng(nx,ny,nz,nqscalar))
    N0x_sng = 0.0

    IF( mphyopt == 1 .OR. (mphyopt == 0 .AND. moist == 1) ) THEN

      IF ( cnvctopt /= 1 ) THEN

!-----------------------------------------------------------------------
!
!  Kessler warm rain microphysics or saturation adjustment only.
!
!-----------------------------------------------------------------------

        CALL microph(nx,ny,nz,mphyflg,dtbig1,                           &
                     ptprt,pprt,qv,qscalar,                             &
                     raing,prcrate(1,1,4),                              &
                     rhostr,pbar,ptbar,ppi,j3,j3inv,                    &
                     tem1,tem2,tem3,tem4,tem5,tem6,                     &
                     tem20,tem21,tem22,tem23,tem24)

      END IF

    ELSE IF( mphyopt == 2 ) THEN    ! Ice microphysics.

!-----------------------------------------------------------------------
!
!  Lin ice microphysics scheme with modifications by W.G. Tao.
!  (Lin, Y.-L., R. D. Farley, and H. D. Orville, 1983:
!  Bulk parameterization of the snow field in a cloud model.
!  J. Climate Appl. Meteor., 22, 1065-1092.
!  Tao, W.-K., and J. Simpson, 1993: Goddard cumulus ensemble model.
!  Part I: Model description. Terres. Atmos. Ocean Sci., 4, 35-72.)
!
!-----------------------------------------------------------------------
!
      CALL microph_ice(nx,ny,nz,mphyflg,dtbig1,                         &
                       ptprt,pprt,qv,qscalar,                           &
                       raing,prcrate(1,1,4),                            &
                       rhostr,pbar,ptbar,qvbar,ppi,j3,j3inv,            &
                       tem1,tem2,tem3,tem4,tem5,tem6,tem7,              &
                       tem8,tem9,tem10,tem11,tem12,tem13,tem14,tem15,   &
                       mpteqnterms)

    ELSE IF ( mphyopt == 3 ) THEN

!-----------------------------------------------------------------------
!
!  Schultz ice microphysics scheme (Schultz, P., 1995: An explicit
!  cloud physics parameterization for operational numerical
!  weather prediction. Mon. Wea. Rev., 123, 3331-3343).
!
!-----------------------------------------------------------------------
!
      CALL microph_nem(nx,ny,nz,mphyflg,dtbig1,ptprt,pprt,qv,qscalar,   &
                       raing,prcrate(1,1,4),rhostr,                     &
                       pbar,ptbar,ppi,j3,j3inv,                         &
                       tem1,tem2,tem3,tem4,tem5,tem6,tem7,              &
                       tem8,tem9,tem10,tem11,tem12,tem13,               &
                       tem14,tem15,tem16,                               &
                       tem20,tem21,tem22,tem23,tem24)


    ELSE IF ( mphyopt == 4 ) THEN

!-----------------------------------------------------------------------
!
!  Straka implementation of Lin, Farley, Orville (1983) 3-ice scheme
!  Reference: Gilmore et al (2004) Mon. Wea. Rev.
!
!-----------------------------------------------------------------------
!
      ALLOCATE(qvbar_z(nz))
      ALLOCATE(ptbar_z(nz))
      ALLOCATE(dzc_z  (nz))

      DO k=1,nz
        ptbar_z(k)=ptbar(1,1,k)
        qvbar_z(k)=qvbar(1,1,k)
      END DO

      DO k=1,nz-1
         dzc_z(k)= 1.0/(zp(1,1,k+1)-zp(1,1,k))
      END DO
      dzc_z(nz)=dzc_z(nz-1)

      CALL lfo_ice_drive(ptbar_z,qvbar_z,                             &
              ptbar,ptprt(1,1,1,3),qv(1,1,1,3),qscalar(:,:,:,3,:),    &
              qvbar, pprt(1,1,1,3), pbar,                             &
              dtbig1,dzc_z,rhostr,j3inv,tem1,                         &
              nx, ny, nz , 0 )

      WHERE (qscalar(:,:,:,3,:) < 0.0) qscalar(:,:,:,3,:) = 0.0

      DEALLOCATE(ptbar_z,qvbar_z,dzc_z)

    ELSE IF ( mphyopt >= 5 .AND. mphyopt <= 7 ) THEN

      mscheme = mphyopt - 5

      CALL micro_wsm6_driver(mscheme,nx,ny,nz,dtbig1,zp,w(:,:,:,tfuture),&
              ptprt(:,:,:,tfuture),ptbar,pprt(:,:,:,tfuture),pbar,ppi,  &
              qv(:,:,:,tfuture),qvbar,qscalar(:,:,:,tfuture,:),         &
              raing,prcrate(:,:,4),                                     &
              tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,tem9,             &
              mpteqnterms,N0x_sng)


    ELSE IF ( mphyopt >= 8 .AND. mphyopt <= 11 ) THEN

      mscheme = mphyopt - 7

      IF (ny > MAX(nscalar,12)) THEN             ! we want to use continuous memeory
        temscalar1 => temscalar(:,:,:,1)
        temscalar2 => temscalar(:,:,:,2)
        temscalar3 => temscalar(:,:,:,3) ! only used when mscheme == 1
      ELSE IF (ny > 3) THEN              ! We know ny = 3 for 2-D simulations
        !temscalar1 => temscalar(:,1,:,:) ! pointers to discontinuous array sections
        !temscalar2 => temscalar(:,2,:,:)
        !temscalar3 => temscalar(:,3,:,:) ! only used when mscheme == 1

        ALLOCATE(temscalar1(nx,nz,nscalar))
        ALLOCATE(temscalar2(nx,nz,nscalar))
        ALLOCATE(temscalar3(nx,nz,12))

      ELSE                               ! impossible
        CALL arpsstop('No enough temporary space for temscalar before   &
                     & calling multimoment microphysics.',1)
      END IF

      CALL multimoment_driver(mscheme,nx,ny,nz,dtbig1,zp,w(:,:,:,tfuture),&
             ptprt,ptbar,pprt,pbar,ppi,qv,qvbar,qscalar,                &
             raing,prcrate(:,:,4),                                      &
             tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,tem9,tem10,        &
             temscalar1,temscalar2,temscalar3) ! ,mpteqnterms,N0x)

      IF (ny > 3 .and. ny <= MAX(nscalar,12)) THEN
        IF (ASSOCIATED(temscalar1)) DEALLOCATE(temscalar1, temscalar2, temscalar3)
      END IF

    ELSE IF ( mphyopt == 12 ) THEN

      CALL my2mom_driver(nx,ny,nz,dtbig1,w(:,:,:,tfuture),              &
              ptprt,ptbar,pprt,pbar,ppi,qv,qvbar,qscalar,               &
              raing,prcrate(:,:,4),istatus)
      IF (istatus /= 0) CALL arpsstop("arpsstop called from CORDINTG after mymom_driver.",1)

    ELSE                             ! Invalid option

      WRITE(6,'(1x,a,/1x,a,i3)')                                        &
           'Invalid option for microphysics parameterization.',         &
           'Job stopped here. MPHYOPT was ',mphyopt
      CALL arpsstop("arpsstop called from CORDINTG improper mphyopt",1)

    END IF

!
!-----------------------------------------------------------------------
!
!  Write out some microphysics fields to separate files for diagnostic
!  purposes (e.g. evaporation rates, intercept parameters).  Added by
!  DTD 01/07
!
!-----------------------------------------------------------------------
!
   IF (curtim >= tstrtdmp .and. MOD(curtim,thisdmp) == 0 .and. mpthermdiag == 1) THEN

   ! Write out the microphysical cooling and heating terms

   mpteqnid(1)='evapqc'
   mpteqnid(2)='evapqr'
   mpteqnid(3)='sublqi'
   mpteqnid(4)='sublqs'
   mpteqnid(5)='sublqg'
   mpteqnid(6)='sublqh'
   mpteqnid(7)='meltqi'
   mpteqnid(8)='meltqs'
   mpteqnid(9)='meltqg'
   mpteqnid(10)='meltqh'
   mpteqnid(11)='condqc'
   mpteqnid(12)='condqr'
   mpteqnid(13)='nuclqi'
   mpteqnid(14)='depoqi'
   mpteqnid(15)='depoqs'
   mpteqnid(16)='depoqg'
   mpteqnid(17)='depoqh'
   mpteqnid(18)='frzqci'
   mpteqnid(19)='colqci'
   mpteqnid(20)='colqcs'
   mpteqnid(21)='colqcg'
   mpteqnid(22)='colqch'
   mpteqnid(23)='frzqrh'
   mpteqnid(24)='colqri'
   mpteqnid(25)='colqrs'
   mpteqnid(26)='colqrg'
   mpteqnid(27)='colqrh'
   mpteqnid(28)='clipqv'
!   mpteqnid(29)='clipqc'
!   mpteqnid(30)='clipqr'
!   mpteqnid(31)='clipqi'
!   mpteqnid(32)='clipqs'
!   mpteqnid(33)='clipqg'
!   mpteqnid(34)='clipqh'

   mpteqnname(1)='Cloud evaporation rate'
   mpteqnname(2)='Rain evaporation rate'
   mpteqnname(3)='Ice sublimation rate'
   mpteqnname(4)='Snow sublimation rate'
   mpteqnname(5)='Graupel sublimation rate'
   mpteqnname(6)='Hail sublimation rate'
   mpteqnname(7)='Ice melting rate'
   mpteqnname(8)='Snow melting rate'
   mpteqnname(9)='Graupel melting rate'
   mpteqnname(10)='Hail melting rate'
   mpteqnname(11)='Cloud condensation rate'
   mpteqnname(12)='Rain condensation rate'
   mpteqnname(13)='Ice nucleation rate'
   mpteqnname(14)='Ice deposition rate'
   mpteqnname(15)='Snow deposition rate'
   mpteqnname(16)='Graupel deposition rate'
   mpteqnname(17)='Hail deposition rate'
   mpteqnname(18)='Cloud-to-ice freezing rate'
   mpteqnname(19)='Cloud-to-ice collection rate'
   mpteqnname(20)='Cloud-to-snow collection rate'
   mpteqnname(21)='Cloud-to-graupel collection rate'
   mpteqnname(22)='Cloud-to-hail collection rate'
   mpteqnname(23)='Rain-to-hail freezing rate'
   mpteqnname(24)='Rain-to-ice collection rate'
   mpteqnname(25)='Rain-to-snow collection rate'
   mpteqnname(26)='Rain-to-graupel collection rate'
   mpteqnname(27)='Rain-to-hail collection rate'
   mpteqnname(28)='qv clipping rate'
!   mpteqnname(29)='qr clipping rate'
!   mpteqnname(30)='qi clipping rate'
!   mpteqnname(31)='qs clipping rate'
!   mpteqnname(32)='qg clipping rate'
!   mpteqnname(33)='qh clipping rate'

   DO nq=1,28
     mpteqnunits(nq)='kg/kg/s'
   END DO

   DO nq=1,28

     varid = mpteqnid(nq)
     varname = mpteqnname(nq)
     varunits = mpteqnunits(nq)

     CALL wrtvar2(nx,ny,nz,mpteqnterms(:,:,:,nq), varid,varname,varunits, &
                  curtim,runname, dirname,3,2,1,istatus)
   END DO


   ! Write out total heating/cooling rate due to
   ! microphysics

   varid = 'mptrate'
   varname = 'Microphysical heating/cooling rate'
   varunits = 'C'

   DO k=1,nz-1
     DO j=1,ny-1
       DO i=1,nx-1
         mptrate(i,j,k) = ((ptprt(i,j,k,3)+ptbar(i,j,k))*ppi(i,j,k)-mptrate(i,j,k))/dtbig1
       END DO
     END DO
   END DO

   CALL wrtvar2(nx,ny,nz,mptrate, varid,varname,varunits,curtim,runname,    &
                dirname,3,2,1,istatus)

   ! Write out intercept parameter for the various species

!   IF(mphyopt >= 5) THEN
!     DO nq=1,nqscalar

!       varid = 'n0__'//qnames(nq)
!       varname = qdescp(nq)
!       varunits = 'm^-4'

!       IF(mphyopt >= 8) THEN
!       DO k=1,nz-1
!         DO j=1,ny-1
!           DO i=1,nx-1

!             N0x_sng(i,j,k,nq) = sngl(N0x(i,j,k,nq))

            ! IF(N0x(i,j,k,nq) /= 0.0) THEN
            !   print*,'nonzero N0x before writeout',N0x(i,j,k,nq)
            ! END IF
!           END DO
!         END DO
!       END DO
!       END IF

!       CALL wrtvar2(nx,ny,nz,N0x_sng(:,:,:,nq), varid,varname,varunits,curtim,runname,    &
!               dirname,3,2,1,istatus)

!     END DO
!   END IF

   END IF ! Writing additional mp variables

   DEALLOCATE(mpteqnterms)
   DEALLOCATE(N0x)
   DEALLOCATE(N0x_sng)
   DEALLOCATE(mptrate)
!
!-----------------------------------------------------------------------
!
!  The microphysics may change the boundary values of both the
!  potential temperature and specific humidity. For the choice of
!  external boundary conditions, we must adjust them to the external
!  BC.
!
!-----------------------------------------------------------------------

    IF ( lbcopt == 2 .AND. mgrid == 1 ) THEN

      CALL acct_interrupt(bc_acct)

      CALL exbcpt( nx,ny,nz, curtim+dtbig, ptprt(1,1,1,tfuture),      &
                   exbcbuf(npt0exb),exbcbuf(nptdtexb) )
      CALL exbcq( nx,ny,nz, 0,curtim+dtbig, qv(1,1,1,tfuture),        &
                   exbcbuf(nqv0exb),exbcbuf(nqvdtexb) )

      DO nq = 1, nscalar
        CALL exbcq( nx,ny,nz,nq,curtim+dtbig,qscalar(:,:,:,tfuture,nq),&
                  exbcbuf(nqscalar0exb(nq)),exbcbuf(nqscalardtexb(nq)))
      END DO

      CALL acct_stop_inter

    END IF

    DO j=1,ny-1
      DO i=1,nx-1
        prcrate(i,j,1)=prcrate(i,j,2)+prcrate(i,j,3)+prcrate(i,j,4)
      END DO
    END DO

  END IF     ! moist.ne.0
!
!-----------------------------------------------------------------------
!
!  Apply nudging
!
!-----------------------------------------------------------------------
!
  CALL set_acct(misc_acct)

  IF(nudgopt > 0 .AND. MOD(nstep,nudgstp) == 0 .AND.       &
     (curtim+dtbig) < ndstop .AND. (curtim+dtbig) > ndstart ) THEN

    IF (myproc == 0) WRITE(6,'(a,f10.1)')                               &
        ' Nudging forecast at ',(curtim+dtbig)

    CALL nudgeall(nx,ny,nz, nxndg,nyndg,nzndg,                          &
        u(1,1,1,tpresent),v(1,1,1,tpresent),                            &
        w(1,1,1,tpresent),pprt(1,1,1,tpresent),ptprt(1,1,1,tpresent),   &
        qv(1,1,1,tpresent),qscalar(:,:,:,tpresent,:),                   &
        uincr,vincr,wincr,pincr,ptincr,qvincr,                          &
        qcincr,qrincr,qiincr,qsincr,qhincr)


    CALL nudgeall(nx,ny,nz, nxndg,nyndg,nzndg,                          &
        u(1,1,1,tfuture),v(1,1,1,tfuture),                              &
        w(1,1,1,tfuture),pprt(1,1,1,tfuture),ptprt(1,1,1,tfuture),      &
        qv(1,1,1,tfuture),qscalar(:,:,:,tfuture,:),                     &
        uincr,vincr,wincr,pincr,ptincr,qvincr,                          &
        qcincr,qrincr,qiincr,qsincr,qhincr)

  END IF
!
!-----------------------------------------------------------------------
!
!  Apply the Asselin time filter to all variables
!
!-----------------------------------------------------------------------
!

  IF(tintegopt == 1) THEN
    IF( flteps /= 0 .AND. frstep /= 1 ) THEN

  !   DO nq = nscalarq+1, nscalar

        ! DTD: Test for Z variable issue: Multiply Z by 1.0e18 before passing
        ! to dynamics, and then reverting back afterwards.

  !      IF(nq >= 13) THEN
  !        DO k=1,nz
  !          DO j=1,ny
  !            DO i=1,nx
  !              qscalar(i,j,k,1,nq) = qscalar(i,j,k,1,nq)*1.0e18
  !              qscalar(i,j,k,2,nq) = qscalar(i,j,k,2,nq)*1.0e18
  !              qscalar(i,j,k,3,nq) = qscalar(i,j,k,3,nq)*1.0e18
  !            END DO
  !          END DO
  !        END DO
  !      END IF

  !   END DO

      CALL tfilt(nx,ny,nz, u,v,w,ptprt,pprt,qv,qscalar,tke)

  !   DO nq=nscalarq+1,nscalar
  !      IF(nq >= 13) THEN
  !        DO k=1,nz
  !          DO j=1,ny
  !            DO i=1,nx
  !              qscalar(i,j,k,1,nq) = qscalar(i,j,k,1,nq)*1.0e-18
  !              qscalar(i,j,k,2,nq) = qscalar(i,j,k,2,nq)*1.0e-18
  !              qscalar(i,j,k,3,nq) = qscalar(i,j,k,3,nq)*1.0e-18
  !            END DO
  !          END DO
  !        END DO
  !      END IF

  !   END DO


    END IF
  END IF

!-----------------------------------------------------------------------
!
!    Update Meso and flux variables              (JAB)
!
!------------------------------------------------------------------------

  IF (radopt == 3) THEN

    CALL initztime(ayear,amm,aday)              ! (JAB)


    CALL readjmeso(nx,ny,nz,ayear,amm,aday,ztime,                    &
       u(1,1,1,3),v(1,1,1,3),w(1,1,1,3),ptbar,pbar,ptprt(1,1,1,3),   &
       pprt(1,1,1,3),qv(1,1,1,3),prcrate)

  END IF

!
!-----------------------------------------------------------------------
!
!
!  Call the Klemp-Lilly (1978)/Durran (1983) version of the radiation
!
!
!-----------------------------------------------------------------------
!

  IF(rbcopt == 2 .OR. rbcopt == 3 .OR. rbcopt == 4)THEN

    CALL set_acct(bc_acct)

    IF(frstep == 1) THEN
      dtbig1 = dtbig/2
    ELSE
      dtbig1 = dtbig
    END IF

    CALL bdtu(nx,ny,nz,dtbig1, u,ubar,udteb,udtwb)

    CALL bdtv(nx,ny,nz,dtbig1, v,vbar,vdtnb,vdtsb)

  END IF
!
!-----------------------------------------------------------------------
!
!  To prepare for the next timestep integration, swap time levels
!  for all time-dependent variables.
!
!  The swap is as follows:
!
!  Fields at time tpresent are assigned to tpast,
!  Fields at time tfuture  are assigned to tpresent,
!  Fields at time tfuture  are not changed.
!
!-----------------------------------------------------------------------
!
  CALL set_acct(misc_acct)

  CALL tflip(nx,ny,nz, u,v,w,ptprt,pprt,qv,qscalar,tke)
!
!-----------------------------------------------------------------------
!
!  Print max./min. for the nested grid case.
!
!-----------------------------------------------------------------------

!  IF( mgrid.eq.1 .and. nestgrd.eq.1 ) THEN

  curtim = curtim + dtbig

!    CALL maxmin(mptr,nx,ny,nz,nzsoil,tlevel,rhobar,                     &
!                u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,kmh,kmv,          &
!                x,y,z,zp,zpsoil,mapfct,                                 &
!                tsoil(:,:,:,0),qsoil(:,:,:,0),wetcanp(:,:,0),           &
!                tem1,tem2,tem3,tem4)
!
!-----------------------------------------------------------------------
!
!  Check the stability of the integration.  If the model solution
!  appears to be exceeding the linear stability condition, then
!  perform a data dump for post-mortem.
!
!-----------------------------------------------------------------------
!
  CALL chkstab(mgrid,nx,ny,nz,nzsoil,nstyps,                            &
               u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,                   &
               ubar,vbar,ptbar,pbar,rhostr,qvbar,kmh,kmv,               &
               x,y,z,zp,zpsoil,zp(1,1,2),mapfct,j1,j2,j3,j3soil,        &
               soiltyp,stypfrct,vegtyp,lai,roufns,veg,                  &
               tsoil,qsoil,wetcanp,snowdpth,                            &
               raing,rainc,prcrate,                                     &
               radfrc,radsw,rnflx,radswnet,radlwin,                     &
               usflx,vsflx,ptsflx,qvsflx,                               &
               tem1,tem2,tem3, tem4)

  curtim = curtim - dtbig

!  ENDIF

  RETURN
END SUBROUTINE cordintg
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE TINTEG                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE tinteg(mptr, frstep, nx,ny,nz,nzsoil,exbcbufsz,              &
           u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,pbldpth,               &
           udteb,udtwb,udtnb,udtsb,vdteb,vdtwb,vdtnb,vdtsb,             &
           wdteb,wdtwb,wdtnb,wdtsb,pdteb,pdtwb,pdtnb,pdtsb,             &
           phydro,sdteb,sdtwb,sdtnb,sdtsb,                              &
           ubar,vbar,ptbar,pbar,ptbari,pbari,                           &
           rhostr,rhostri,qvbar,ppi,csndsq,                             &
           x,y,z,zp,zpsoil,mapfct,                                      &
           j1,j2,j3,j3soil,aj3x,aj3y,aj3z,j3inv,                        &
           trigs1,trigs2,ifax1,ifax2,                                   &
           wsave1,wsave2,vwork1,vwork2,sinlat,                          &
           ptsfc,qvsfc,prcrate, radfrc,                                 &
           usflx,vsflx,ptsflx,qvsflx,kmh,kmv,rprntl,                    &
           ptcumsrc,qcumsrc,raing,rainc,w0avg,nca,kfraincv,             &
           cldefi,xland,bmjraincv,roufns,                               &
           exbcbuf,bcrlx,                                               &
           rhofct,defsq,lenscl,                                         &
           uforce,vforce,wforce,pforce,ptforce,                         &
           tem1,tem2,tem3,tem4,tem5,tem6,                               &
           tem7,tem8,tem9,tem10,tem11,tem12,tem13,tem14,tem15,          &
           tem16,tem17,tem18,                                           &
           tem1_0,tem2_0,tem3_0)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Orchestrate the time integration of the dynamics of the basic governing
!  equations for a single time step. Note that the time tendencies of the
!  non-conservative processes (e.g., microphysics, radiation, etc.)
!  are not included in this routine but are applied in their
!  corresponding routines.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  3/21/92.
!
!  MODIFICATION HISTORY:
!
!  5/03/92 (K. Droegemeier)
!    Added full documentation.
!
!  5/20/92 (M. Xue)
!    Reformatted certain documentations.
!
!  4/10/93 (M. Xue & Hao Jin)
!    Add the terrain.
!
!  5/26/93 (M. Xue)
!    Added w into the solvpt argument list.
!
!  9/7/94 (M.Xue)
!    Call to ADVP modified.
!
!  9/10/94 (D. Weber & Y. Lu)
!  Cleaned up documentation.
!
!  8/20/95 (M. Xue)
!  Bug fix with the use of ptcumsrc.
!
!  10/31/95 (D. Weber)
!  Added trigs1,trigs2,ifax1,ifax2 for use in the upper w-p
!  radiation condition.
!
!  01/23/1996 (Donghai Wang, Ming Xue and Yuhe Liu)
!  Added the map projection factor to ARPS governing equations.
!
!  03/11/96 (M. Xue)
!  Modified the call to SOLVTKE, so that all three time levels of
!  ptprt, pprt etc are passed into the routine.
!
!  03/19/96 (Yuhe Liu)
!  Added radiation forcing
!
!  07/10/1997 (Fanyou Kong - CMRP)
!  Fixed a bug in 'solvtke' with corrected temporary arga
!
!  07/22/97 (D. Weber)
!  Added wsave1,wsave2,vwork1,vwork2 for use in the even fft version
!  of the upper w-p radiation condition (fftopt=2).
!
!  08/01/97 (Zonghui Huo)
!  Added Kain-fritsch cumulus parameterization scheme.
!
!  10/21/97 (Donghai Wang)
!  Using total density (rho) in the calculation of the pressure
!  gradient force terms, and added the second order terms
!  in the linerized buoyancy terms.
!
!  11/05/97 (D. Weber)
!  Added phydro array for use in the bottom boundary condition for
!  perturbation pressure (hydrostatic).
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
!  9/18/98 (D. Weber)
!  Added precomputed averages of j3 in the x, y, and z directions
!  to improve the code efficiency.
!
!  07/23/2001 (K. Brewster)
!  Added mptr to argument list, needed for printing of diagnostic
!  noise parameter.  Also added mptr to calls to smlstep, solvuw,
!  solvwpex, and solvwpim.
!
!  13 March 2002 (Eric Kemp)
!  Added arrays for WRF BMJ cumulus scheme.
!
!  April 2002 (Fanyou Kong)
!  Added cnvctopt=5 option for the new WRF K-F (KF_ETA) scheme
!
!  09/23/2011 (Dan Dawson)
!  Added code to convert the number concentration and reflectivity scalars (for the MY multimoment
!  scheme) to per kg units for the dynamics.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    frstep   Initial time step flag.
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of soil levels
!    u        x component of velocity at times tpast and tpresent (m/s)
!    v        y component of velocity at times tpast and tpresent (m/s)
!    w        Vertical component of Cartesian velocity at times
!             tpast and tpresent (m/s)
!    wcont    Contravariant vertical velocity (m/s)
!             computational coordinates (m/s)
!    ptprt    Perturbation potential temperature at times tpast and
!             tpresent (K)
!    pprt     Perturbation pressure at times tpast and tpresent (Pascal)
!    qv       Water vapor specific humidity at times tpast and tpresent (kg/kg)
!    qc       Cloud water mixing ratio at times tpast and tpresent (kg/kg)
!    qr       Rainwater mixing ratio at times tpast and tpresent (kg/kg)
!    qi       Cloud ice mixing ratio at times tpast and tpresent (kg/kg)
!    qs       Snow mixing ratio at times tpast and tpresent (kg/kg)
!    qh       Hail mixing ratio at times tpast and tpresent (kg/kg)
!
!    udteb    Time tendency of u field at east boundary (m/s**2)
!    udtwb    Time tendency of u field at west boundary (m/s**2)
!
!    vdtnb    Time tendency of v field at north boundary (m/s**2)
!    vdtsb    Time tendency of v field at south boundary (m/s**2)
!
!    pdteb    Time tendency of pprt field at east boundary (Pascal/s)
!    pdtwb    Time tendency of pprt field at west boundary (Pascal/s)
!    pdtnb    Time tendency of pprt field at north boundary (Pascal/s)
!    pdtsb    Time tendency of pprt field at south boundary (Pascal/s)
!
!    phydro   Big time step forcing term for use in computing the
!             hydrostatic pressure at k=1.
!
!    ubar     Base state zonal velocity component (m/s)
!    vbar     Base state meridional velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    ptbari   Inverse Base state potential temperature (K)
!    pbari    Inverse Base state pressure (Pascal)
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    rhostri  Inverse base state density rhobar times j3 (kg/m**3)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!    cnsd     Sound wave speed.
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!    zpsoil   Height of soil levels (m)
!    mapfct   Map factors at scalar, u and v points
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    j3soil   Coordinate transformation Jacobian  d(zpsoil)/dz
!    aj3x     Avgx of the coordinate transformation Jacobian  d(zp)/dz
!    aj3y     Avgy of the coordinate transformation Jacobian  d(zp)/dz
!    aj3z     Avgz of the coordinate transformation Jacobian  d(zp)/dz
!    j3inv    Inverse of the coordinate transformation j3
!    trigs1   Array containing pre-computed trig function for
!             fftopt=1.
!    trigs2   Array containing pre-computed trig function for
!             fftopt=1.
!    ifax1    Array containing the factors of nx for fftopt=1.
!    ifax2    Array containing the factors of ny for fftopt=1.
!
!    vwork1   2-D work array for fftopt=2.
!    vwork2   2-D work array for fftopt=2.
!    wsave1   Work array for fftopt=2.
!    wsave2   Work array for fftopt=2.
!
!    sinlat   Sin of latitude at each grid point
!
!    usflx    Surface flux of u-momentum (kg/(m*s**2))
!    vsflx    Surface flux of v-momentum (kg/(m*s**2))
!    ptsflx   Surface heat flux (K*kg/(m**2 * s ))
!    qvsflx   Surface moisture flux (kg/(m**2 * s))
!
!    ptcumsrc Source term in pt-equation due to cumulus parameterization
!    qcumsrc  Source term in water equations due to cumulus parameterization
!    kfraincv   K-F convective rainfall (cm)
!    nca      K-F counter for CAPE release
!    cldefi   BMJ cloud efficiency
!    xland    BMJ land/sea mask
!    bmjraincv   BMJ convective rainfall (cm)
!
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!    rprntl   Reciprocal of Prandtl number
!
!    radfrc   Radiation forcing (K/s)
!
!  OUTPUT:
!
!    u        x component of velocity at time tfuture (m/s)
!    v        y component of velocity at time tfuture (m/s)
!    w        Vertical component of Cartesian velocity at tfuture (m/s)
!    ptprt    Perturbation potential temperature at time tfuture (K)
!    pprt     Perturbation pressure at time tfuture (Pascal)
!    qv       Water vapor specific humidity at time tfuture (kg/kg)
!    qc       Cloud water mixing ratio at time tfuture (kg/kg)
!    qr       Rainwater mixing ratio at time tfuture (kg/kg)
!    qi       Cloud ice mixing ratio at time tfuture (kg/kg)
!    qs       Snow mixing ratio at time tfuture (kg/kg)
!    qh       Hail mixing ratio at time tfuture (kg/kg)
!
!    udteb    Time tendency of u field at east boundary (m/s**2)
!    udtwb    Time tendency of u field at west boundary (m/s**2)
!
!    vdtnb    Time tendency of v field at north boundary (m/s**2)
!    vdtsb    Time tendency of v field at south boundary (m/s**2)
!
!    pdteb    Time tendency of pprt field at east boundary (Pascal/s)
!    pdtwb    Time tendency of pprt field at west boundary (Pascal/s)
!    pdtnb    Time tendency of pprt field at north boundary (Pascal/s)
!    pdtsb    Time tendency of pprt field at south boundary (Pascal/s)
!
!    phydro   Big time step forcing term for use in computing the
!             hydrostatic pressure at k=1.
!
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!    pbldpth  Planetary boundary layer depth (m)
!
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!    rprntl   Reciprocal of Prandtl number
!
!  WORK ARRAYS:
!
!    defsq    Deformation squared (1/s**2)
!    lenscl   Turbulent mixing length scale (m)
!    uforce   Acoustically inactive forcing terms in u-Eq. ((kg/(m*s)**2)
!    vforce   Acoustically inactive forcing terms in v-Eq. ((kg/(m*s)**2)
!    wforce   Acoustically inactive forcing terms in w-Eq. ((kg/(m*s)**2)
!    pforce   Acoustically inactive forcing terms in p-eq. (Pascal/s)
!    ptforce  Gravity wave inactive forcing terms in pt-eq. (K*kg/(m**3*s))
!
!    rhofct   rho-factor: rhobar/rho
!
!    udtnb    Time tendency of u field at north boundary (m/s**2)
!    udtsb    Time tendency of u field at south boundary (m/s**2)
!
!    vdteb    Time tendency of v field at east boundary (m/s**2)
!    vdtwb    Time tendency of v field at west boundary (m/s**2)
!
!    wdteb    Time tendency of w field at east boundary (m/s**2)
!    wdtwb    Time tendency of w field at west boundary (m/s**2)
!    wdtnb    Time tendency of w field at north boundary (m/s**2)
!    wdtsb    Time tendency of w field at south boundary (m/s**2)
!
!    ptdteb   Time tendency of ptprt field at east boundary (K/s)
!    ptdtwb   Time tendency of ptprt field at west boundary (K/s)
!    ptdtnb   Time tendency of ptprt field at north boundary(K/s)
!    ptdtsb   Time tendency of ptprt field at south boundary(K/s)
!
!    sdteb    Time tendency of a scalar at e-boundary
!    sdtwb    Time tendency of a scalar at w-boundary
!    sdtnb    Time tendency of a scalar at n-boundary
!    sdtsb    Time tendency of a scalar at s-boundary
!
!    tem1     Temporary work array
!    tem2     Temporary work array
!    tem3     Temporary work array
!    tem4     Temporary work array
!    tem5     Temporary work array
!    tem6     Temporary work array
!    tem7     Temporary work array
!    tem8     Temporary work array
!    tem9     Temporary work array
!    tem10    Temporary work array
!    tem11    Temporary work array
!    tem12    Temporary work array
!    tem13    Temporary work array
!    tem14    Temporary work array
!    tem15    Temporary work array
!    tem16    Temporary work array
!    tem17    Temporary work array
!    tem18    Temporary work array
!
!    tem1_0   Temporary work array.
!    tem2_0   Temporary work array.
!    tem3_0   Temporary work array.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE        ! Force explicit declarations
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'bndry.inc'
  INCLUDE 'exbc.inc'
  INCLUDE 'timelvls.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  INTEGER :: mptr

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of soil levels

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)
  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)
  REAL :: tke   (nx,ny,nz,nt)  ! Turbulent Kinetic Energy ((m/s)**2)
  REAL :: pbldpth(nx,ny,nt)    ! Planetary boundary layer depth (m)

  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)
  REAL :: qscalar(nx,ny,nz,nt,nscalar)

  REAL :: udteb (ny,nz)        ! Time tendency of u at e-boundary (m/s**2)
  REAL :: udtwb (ny,nz)        ! Time tendency of u at w-boundary (m/s**2)
  REAL :: udtnb (nx,nz)        ! Time tendency of u at n-boundary (m/s**2)
  REAL :: udtsb (nx,nz)        ! Time tendency of u at s-boundary (m/s**2)

  REAL :: vdteb (ny,nz)        ! Time tendency of v at e-boundary (m/s**2)
  REAL :: vdtwb (ny,nz)        ! Time tendency of v at w-boundary (m/s**2)
  REAL :: vdtnb (nx,nz)        ! Time tendency of v at n-boundary (m/s**2)
  REAL :: vdtsb (nx,nz)        ! Time tendency of v at s-boundary (m/s**2)

  REAL :: wdteb (ny,nz)        ! Time tendency of w at e-boundary (m/s**2)
  REAL :: wdtwb (ny,nz)        ! Time tendency of w at w-boundary (m/s**2)
  REAL :: wdtnb (nx,nz)        ! Time tendency of w at n-boundary (m/s**2)
  REAL :: wdtsb (nx,nz)        ! Time tendency of w at s-boundary (m/s**2)

  REAL :: pdteb (ny,nz)        ! Time tendency of pprt at e-boundary (Pascal/s)
  REAL :: pdtwb (ny,nz)        ! Time tendency of pprt at w-boundary (Pascal/s)
  REAL :: pdtnb (nx,nz)        ! Time tendency of pprt at n-boundary (Pascal/s)
  REAL :: pdtsb (nx,nz)        ! Time tendency of pprt at s-boundary (Pascal/s)

  REAL :: phydro(nx,ny)        ! Pressure at k=1 computed using pert.
                               ! hydrostatic relation.

  REAL :: sdteb (ny,nz)        ! Time tendency of a variable at e-boundary
  REAL :: sdtwb (ny,nz)        ! Time tendency of a variable at w-boundary
  REAL :: sdtnb (nx,nz)        ! Time tendency of a variable at n-boundary
  REAL :: sdtsb (nx,nz)        ! Time tendency of a variable at s-boundary

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: ptbari(nx,ny,nz)     ! Inverse Base state pot. temperature (K)
  REAL :: pbari (nx,ny,nz)     ! Inverse Base state pressure (Pascal).
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: rhostri(nx,ny,nz)    ! Inverse base state density rhobar times j3.
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                               ! (kg/kg)
  REAL :: ppi   (nx,ny,nz)     ! Exner function
  REAL :: csndsq(nx,ny,nz)     ! Sound wave speed squared.

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
  REAL :: zpsoil(nz,ny,nzsoil) ! Height of soil levels.

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! - d( zp )/d( x ).
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! - d( zp )/d( y ).
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! d( zp )/d( z ).
  REAL :: j3soil(nx,ny,nzsoil) ! Coordinate transformation Jacobian defined as
                               ! d( zpsoil )/d( z ).
  REAL :: aj3x  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE X-DIR.
  REAL :: aj3y  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
  REAL :: aj3z  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Z-DIR.
  REAL :: j3inv (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! d( zp )/d( z ).
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
  REAL :: wsave1 (3*(ny-1)+15) ! Work array for fftopt =2.
  REAL :: wsave2 (3*(nx-1)+15) ! Work array for fftopt =2.

  REAL :: sinlat(nx,ny)        ! Sin of latitude at each grid point

  REAL :: ptsfc  (nx,ny)       ! Potential temperature at the ground level (K)
  REAL :: qvsfc  (nx,ny)       ! Effective qv at the surface (kg/kg)

  REAL :: prcrate(nx,ny,4)     ! precipitation rate (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precipitation rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulus precip. rate
                               ! prcrate(1,1,4) = microphysics precip. rate

  REAL :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m**2*s))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: rprntl(nx,ny,nz)     ! Reciprocal of Prandtl number

  REAL :: defsq (nx,ny,nz)     ! Deformation squared (1/s**2)
  REAL :: lenscl(nx,ny,nz)     ! Turbulent mixing length scale (m)

  REAL :: uforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in u-momentum equation (kg/(m*s)**2)
  REAL :: vforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in v-momentum equation (kg/(m*s)**2)
  REAL :: wforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in w-momentum equation (kg/(m*s)**2)
  REAL :: pforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in pressure equation (Pascal/s)
  REAL :: ptforce(nx,ny,nz)    ! Gravity wave inactive forcing terms
                               ! in pressure equation (Pascal/s)

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

  REAL, INTENT(INOUT) :: cldefi(nx,ny) ! BMJ cloud efficiency
  REAL, INTENT(IN)    :: xland(nx,ny)  ! BMJ land mask
                                      ! (1.0 = land, 2.0 = sea)
  REAL, INTENT(INOUT) :: bmjraincv(nx,ny) ! BMJ convective rainfall
                                         ! (cm)

  REAL, INTENT(IN)   :: roufns(nx,ny)

  REAL :: raing(nx,ny)         ! Grid supersaturation rain
  REAL :: rainc(nx,ny)         ! Cumulus convective rain

  REAL :: radfrc(nx,ny,nz)     ! Radiation forcing (K/s)

  INTEGER :: exbcbufsz         ! EXBC buffer size
  REAL    :: exbcbuf( exbcbufsz ) ! EXBC buffer array
  REAL    :: bcrlx (nx,ny)        ! EXBC relaxation coefficients

  REAL :: rhofct(nx,ny,nz)     ! rho-factor: rhobar/rho

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array
  REAL :: tem5  (nx,ny,nz)     ! Temporary work array
  REAL :: tem6  (nx,ny,nz)     ! Temporary work array
  REAL :: tem7  (nx,ny,nz)     ! Temporary work array
  REAL :: tem8  (nx,ny,nz)     ! Temporary work array
  REAL :: tem9  (nx,ny,nz)     ! Temporary work array
  REAL :: tem10 (nx,ny,nz)     ! Temporary work array
  REAL :: tem11 (nx,ny,nz)     ! Temporary work array
  REAL :: tem12 (nx,ny,nz)     ! Temporary work array
  REAL :: tem13 (nx,ny,nz)     ! Temporary work array
  REAL :: tem14 (nx,ny,nz)     ! Temporary work array
  REAL :: tem15 (nx,ny,nz)     ! Temporary work array
  REAL :: tem16 (nx,ny,nz)     ! Temporary work array
  REAL :: tem17 (nx,ny,nz)     ! Temporary work array
  REAL :: tem18 (nx,ny,nz)     ! Temporary work array

  REAL :: tem1_0(0:nx,0:ny,0:nz)     ! Temporary work array.
  REAL :: tem2_0(0:nx,0:ny,0:nz)     ! Temporary work array.
  REAL :: tem3_0(0:nx,0:ny,0:nz)     ! Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL    :: dtbig1            ! Local value of big time step size
  REAL    :: dtsml1            ! Local value of small time step size
  INTEGER :: frstep            ! Flag for the initial time step
  INTEGER :: i,j,k
  INTEGER :: qflag             ! Indicator for the water/ice type
                               ! when calling SOLVQ.
  INTEGER :: ntst
  REAL    :: tst
  INTEGER :: RK3step

  LOGICAL :: convertunits
  CHARACTER(LEN=40) :: qvarname

  INTEGER :: istatus
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  RK3step = 0    ! Not needed in this version of tinteg
!
!-----------------------------------------------------------------------
!
!  To initialize the three-time level scheme, we assume for
!  the first step that the values of all variables at time
!  past equal the values at time present.  We then perform a
!  forward-in-time integration for the first step only and
!  then, using the same algorithm, we perform a centered-in-time
!  integration for all subsequent steps. The size of the
!  timestep is adjusted accordingly, such that the leapfrog step
!  size is twice that of the first forward-in-time step.
!
!-----------------------------------------------------------------------
!
  IF(frstep == 1) THEN      ! frstep=1 on the first time step
                            ! only - indicating a forward-in-
                            ! time integration.
    dtbig1 = dtbig/2
    dtsml1 = dtsml/2

  ELSE

    dtbig1 = dtbig
    dtsml1 = dtsml

  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate wcont at time tpresent. wcont will be used in FRCUVW
!  and FRCP. rhostr averaged to u, v, w points. They are stored
!  in tem1, tem2 and tem3.
!
!-----------------------------------------------------------------------
!
  CALL rhouvw(nx,ny,nz,rhostr,tem1,tem2,tem3)

  CALL wcontra(nx,ny,nz,                                                &
               u(1,1,1,tpresent),v(1,1,1,tpresent),                     &
               w(1,1,1,tpresent),mapfct,j1,j2,j3,aj3z,                  &
               rhostr,tem1,tem2,tem3,wcont,tem4,tem5)

!-----------------------------------------------------------------------
!
! Calculate Planetary boundary layer tendency
!
!-----------------------------------------------------------------------

  !IF (pblopt > 0) THEN
  IF (tmixopt > 4) THEN
    !ALLOCATE(rublten (nx,ny,nz), STAT = istatus)
    !ALLOCATE(rvblten (nx,ny,nz), STAT = istatus)
    !ALLOCATE(rthblten(nx,ny,nz), STAT = istatus)
    !ALLOCATE(rqvblten(nx,ny,nz), STAT = istatus)
    !ALLOCATE(rqcblten(nx,ny,nz), STAT = istatus)
    !ALLOCATE(rqiblten(nx,ny,nz), STAT = istatus)
    !rublten  = 0.0
    !rvblten  = 0.0
    !rthblten = 0.0
    !rqvblten = 0.0
    !rqcblten = 0.0
    !rqiblten = 0.0

    CALL pbl_driver(tintegopt,tmixopt-4,nx,ny,nz,                       &
                  u,v,w,ptprt,pprt,ppi,qv,qscalar,                      &
                  ptbar,pbar,qvbar,rhostr,tem1,tem2,dtbig1,             &
                  !rublten,rvblten,rthblten,rqvblten,rqcblten,rqiblten,  &
                  zp,roufns,ptsflx,qvsflx,xland,ptsfc,                  &
                  pbldpth,kmv,rprntl,                                   &
                  tem3,tem4,tem5,tem6,tem7,tem8,tem9,                   &
                  tem10,tem11,tem12,tem13,tem14,istatus)
!WRITE(*,*) '==1==',MAXVAL(kmv)
  ELSE

    !ALLOCATE(rublten (1,1,1), STAT = istatus)
    !ALLOCATE(rvblten (1,1,1), STAT = istatus)
    !ALLOCATE(rthblten(1,1,1), STAT = istatus)
    !ALLOCATE(rqvblten(1,1,1), STAT = istatus)
    !ALLOCATE(rqcblten(1,1,1), STAT = istatus)
    !ALLOCATE(rqiblten(1,1,1), STAT = istatus)
    !rublten  = 0.0
    !rvblten  = 0.0
    !rthblten = 0.0
    !rqvblten = 0.0
    !rqcblten = 0.0
    !rqiblten = 0.0

  END IF
!
!-----------------------------------------------------------------------
!
!  Compute the acoustically inactive terms in the momentum and
!  pressure equations that are held fixed during the small time
!  step computations.  This includes advection, buoyancy, mixing
!  (both physical and computational), and the Coriolis terms.
!  These forcing terms are accumulated into arrays for each
!  of the momentum equations, e.g., uforce for the u-equation,
!  vforce for the v-equation, wforce for the w-equation and
!  pforce for the p-equation.
!
!  Note: Arrays, pforce and ptforce, are used as work arrays in
!        subroutine frcuvw.
!
!-----------------------------------------------------------------------
!
  CALL frcuvw(nx,ny,nz,nzsoil,exbcbufsz,dtbig1,                         &
              u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,pbldpth,            &
              ubar,vbar,ptbar,pbar,ptbari,pbari,rhostr,qvbar,           &
              usflx,vsflx, x,y,z,zp,zpsoil,mapfct,                      &
              j1,j2,j3,j3soil,aj3x,aj3y,aj3z,j3inv, sinlat, ptsfc,      &
              uforce,vforce,wforce,kmh,kmv,rprntl,lenscl,defsq,         &
              exbcbuf, bcrlx,rhofct,phydro,                             &
              pforce,ptforce,                                           &
              tem1,tem2,tem3,tem4,tem5,tem6,tem7,                       &
              tem8,tem9,tem10,tem11,tem12,tem13,tem14)

!WRITE(*,*) '==2==',MAXVAL(kmv)
  IF (tmixopt == 4 .OR. tmixopt == 5) THEN
!
!-----------------------------------------------------------------------
!
!  Integrate the TKE equation.
!  NOTE: Arrays lenscl and defsq should NOT be changed between
!  calls to frcuvw and solvtke.
!  Note: pforce and ptforce are used as work arrays in the following call.
!
!-----------------------------------------------------------------------
!
    CALL solvtke(RK3step,nx,ny,nz,dtbig1,u,v,wcont,ptprt,pprt,          &
                 qv,qscalar,tke,                                        &
                 ubar,vbar,ptbar,pbar,ptbari,rhostr,rhostri,qvbar,      &
                 x,y,z,zp, mapfct, j1,j2,j3,aj3x,aj3y,j3inv,            &
                 kmh,kmv,rprntl,lenscl,defsq,                           &
                 ptsflx,qvsflx,                                         &
                 tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,               &
                 tem9,tem10,tem11,ptforce, pforce, & ! pforce and ptforce used as tem arrays
                 tem1_0,tem2_0,tem3_0)

  END IF

  CALL frcp(nx,ny,nz,exbcbufsz, dtbig1,                                 &
            u,v,w,wcont,ptprt,pprt,qv,qscalar,                          &
            ptbar,pbar,rhostr,qvbar,mapfct,j1,j2,j3,aj3x,aj3y,aj3z,     &
            pforce,                                                     &
            exbcbuf, bcrlx,                                             &
            tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,tem9)

!
!-----------------------------------------------------------------------
!
!  Calculate gravity wave or acoustic wave inactive terms in the
!  potential temperature equation.
!
!  The force terms are stored in ptforce.
!
!-----------------------------------------------------------------------
!
  CALL frcpt(nx,ny,nz, exbcbufsz, dtbig1,ptprt,u,v,w,wcont,             &
             ptbar,rhostr,rhostri,kmh,kmv,rprntl,                       &
             usflx,vsflx,ptsflx,pbldpth,                                &
             x,y,z,zp,mapfct, j1,j2,j3,aj3x,aj3y,j3inv,ptsfc,           &
             ptforce,                                                   &
             exbcbuf, bcrlx,                                            &
             tem1,tem2,tem3,tem4,tem5,tem6,tem7,                        &
             tem8,tem9,tem10,tem11,                                     &
             tem1_0,tem2_0,tem3_0,tem12)

  CALL set_acct(tinteg_acct)
  IF ( radopt == 2 ) THEN
    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          ptforce(i,j,k) = ptforce(i,j,k)                               &
                         + rhostr(i,j,k)*radfrc(i,j,k)
        END DO
      END DO
    END DO
  END IF

!
!-----------------------------------------------------------------------
!
!  Calculate source and sink terms in temperature (ptprt) and
!  moisture (qv) equations that are due to subgrid scale cumulus
!  convection.
!
!-----------------------------------------------------------------------
!
  IF (cnvctopt == 1 .OR. cnvctopt == 2 .OR. cnvctopt == 3 .OR. &
      cnvctopt == 4 .OR. cnvctopt == 5) THEN

    CALL set_acct(cum_acct)

!-----------------------------------------------------------------------
!
! Calculate w0avg, is close to a running mean vertical velocity,
! tst is the number of time steps in 10min for K-F scheme
!
!-----------------------------------------------------------------------
!
    IF (cnvctopt == 3 .OR. cnvctopt == 5) THEN

      ntst=nint( 600.0/dtbig )
      tst=FLOAT(ntst)

      IF( (curtim-tstart) <= 300.0 .AND. initopt /= 2 ) THEN
        DO k=1,nz
          DO j=1,ny-1
            DO i=1,nx-1
              w0avg(i,j,k)= (w(i,j,k,1)+w(i,j,k,2))*0.5
            END DO
          END DO
        END DO

      ELSE IF ( (curtim-tstart) > 300.0 .OR. initopt == 2) THEN

        DO k=1,nz
          DO j=1,ny-1
            DO i=1,nx-1
              w0avg(i,j,k)=(w0avg(i,j,k)*(tst-1.)+w(i,j,k,2))/tst
            END DO
          END DO
        END DO
      END IF
    END IF
!
!-----------------------------------------------------------------------
!
!  Calculate vertical velocity for KUO scheme, or running-average
!  vertical velocity (m/s) for Kain Fritsch scheme
!
!-----------------------------------------------------------------------
!
!    IF( cnvctopt == 3 .AND. MOD(curtim+0.001,confrq) <= (0.5*dtbig)     &
!             .AND. ((curtim-tstart) > 300.0                             &
    IF( (cnvctopt == 3 .OR. cnvctopt == 5) .AND.                        &
        MOD(curtim+0.001,confrq) <= (0.5*dtbig) .AND.                   &
        ((curtim-tstart) > 300.0 .OR. initopt == 2)) THEN

      DO k=1,nz
        DO j=1,ny-1
          DO i=1,nx-1
!        tem1(i,j,k) = (w(i,j,k,1)+w(i,j,k,2))*0.5
            tem1(i,j,k) = w0avg(i,j,k)
          END DO
        END DO
      END DO


    ELSE IF( MOD(curtim+0.001,confrq) <= (0.5*dtbig) .OR. nstep == 1 ) THEN
      DO k=1,nz
        DO j=1,ny-1
          DO i=1,nx-1
            tem1(i,j,k) = wcont(i,j,k)*j3(i,j,k)
          END DO
        END DO
      END DO
    END IF
!
!-----------------------------------------------------------------------
!
!  Call cumulus parameterization schemes
!  Make sure to reset kfraincv, ptcumsrc and qcumsrc to 0 once nca<=0
!
!-----------------------------------------------------------------------
!
    IF( (moist /= 0) .AND. ((curtim-tstart) > 300.0 .OR. initopt == 2) .AND. &
          (MOD(curtim+0.001,confrq) <= (0.5*dtbig))) THEN

      DO j=1,ny-1
        DO i=1,nx-1
          IF (nca(i,j) <= 0) THEN
            kfraincv(i,j) = 0.0
            prcrate(i,j,3) = 0.0
          END IF
        END DO
      END DO

      DO k=2,nz-2
        DO j=1,ny-1
          DO i=1,nx-1
            IF (nca(i,j) <= 0) THEN
              ptcumsrc(i,j,k) = 0.0
              qcumsrc(i,j,k,1) = 0.0
              qcumsrc(i,j,k,2) = 0.0
              qcumsrc(i,j,k,3) = 0.0
              qcumsrc(i,j,k,4) = 0.0
              qcumsrc(i,j,k,5) = 0.0
            END IF
          END DO
        END DO
      END DO

      CALL qpfcums(nx,ny,nz,prcrate(1,1,3),qvsflx,                      &
                   u(1,1,1,2),v(1,1,1,2),tem1,                          &
                   pprt(1,1,1,2),ptprt(1,1,1,2),qv(1,1,1,2),            &
                   pbar,ptbar,qvbar,rhostr,zp,j3,                       &
                   ptcumsrc,qcumsrc,rainc,nca,kfraincv,                 &
                   cldefi,xland,bmjraincv,                              &
                   tem2,tem3,tem4,tem5,tem6,                            &
                   tem7,tem8,tem9,tem10,tem11)
    END IF
!
!  Accumulate rainc and update nca. Note: raincv is in cm.
!
    DO j=1,ny-1
      DO i=1,nx-1
        rainc(i,j) = rainc(i,j) + 10.0*kfraincv(i,j) + 10.0*bmjraincv(i,j)
        IF ( nca(i,j) >= 1 ) nca(i,j) = nca(i,j) - 1
      END DO
    END DO

    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          ptforce(i,j,k) = ptforce(i,j,k)                               &
                         + rhostr(i,j,k)*ptcumsrc(i,j,k)
        END DO
      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  End of cumulus parameterization
!
!-----------------------------------------------------------------------
!
  CALL set_acct(misc_acct)

  IF( lvldbg >= 2 ) THEN
    CALL checkuhx(uforce, nx,ny,nz,2,nx-1,1,ny-1,2,nz-2,                &
                  'ufrcx', tem1)
    CALL checkuhy(uforce, nx,ny,nz,2,nx-1,1,ny-1,2,nz-2,                &
                  'ufrcy', tem1)
    CALL checkvhx(vforce, nx,ny,nz,1,nx-1,2,ny-1,2,nz-2,                &
                  'vfrcx', tem1)
    CALL checkvhy(vforce, nx,ny,nz,1,nx-1,2,ny-1,2,nz-2,                &
                  'vfrcy', tem1)
    CALL checkwhx(wforce, nx,ny,nz,1,nx-1,1,ny-1,2,nz-1,                &
                  'wfrcx', tem1)
    CALL checkwhy(wforce, nx,ny,nz,1,nx-1,1,ny-1,2,nz-1,                &
                  'wfrcy', tem1)
    CALL checkshx(pforce, nx,ny,nz,1,nx-1,1,ny-1,2,nz-2,                &
                  'pfrcx', tem1)
    CALL checkshy(pforce, nx,ny,nz,1,nx-1,1,ny-1,2,nz-2,                &
                  'pfrcy', tem1)
    CALL checkshx(ptforce, nx,ny,nz,1,nx-1,1,ny-1,2,nz-2,               &
                  'ptfrcx', tem1)
    CALL checkshy(ptforce, nx,ny,nz,1,nx-1,1,ny-1,2,nz-2,               &
                  'ptfrcy', tem1)
  END IF
!
!-----------------------------------------------------------------------
!
!  Integrate the momentum and pressure equations (i.e., the
!  acoustically-active equations) in time using a mode-splitting
!  approach.  The momentum components are advanced in time using
!  a forward scheme (relative to the pressure gradient force terms)
!  and the pressure is advanced using a backward scheme (relative to
!  the divergence term, which is evaluated using the newly updated
!  velocities hence the backward scheme). These small time steps
!  bring the momentum and pressure from time tpast to tfuture.
!  During these steps, the acoustically-inactive forcing terms (e.g.,
!  uforce, pforce, etc.) are held fixed at time tpresent (i.e. they
!  are evaluated only once, at time tpresent, for each large time
!  step integration), hence the time integration is leapfrog-
!  in-time relative to these forcing terms.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  When this solver is called for a nested grid (not the base grid),
!  the information stored in u, v, w and pprt is used to set the
!  time tendency of these variables at the lateral boundaries.
!
!-----------------------------------------------------------------------
!
  IF( mgrid /= 1 .AND. nestgrd == 1 ) THEN

    CALL nestbdt(nx,ny,nz,u, 1,nx,1,ny-1,1,nz-1,dtbig,                  &
                 udteb,udtwb,udtnb,udtsb)

    CALL nestbdt(nx,ny,nz,v, 1,nx-1,1,ny,1,nz-1,dtbig,                  &
                 vdteb,vdtwb,vdtnb,vdtsb)

    CALL nestbdt(nx,ny,nz,w, 1,nx-1,1,ny-1,1,nz,dtbig,                  &
                 wdteb,wdtwb,wdtnb,wdtsb)

    CALL nestbdt(nx,ny,nz,pprt, 1,nx-1,1,ny-1,1,nz-1,dtbig,             &
                 pdteb,pdtwb,pdtnb,pdtsb)

    CALL nestbdt(nx,ny,nz,ptprt,1,nx-1,1,ny-1,1,nz-1,dtbig,             &
                 sdteb,sdtwb,sdtnb,sdtsb)

!
!-----------------------------------------------------------------------
!
!    For the first forward time step, make sure that the fields at
!    tpast equal those at tpresent. This can only be done after
!    the time tendencies on the boundaries are calculated and stored
!    in the time tendency arrays.
!
!-----------------------------------------------------------------------
!
    IF( frstep == 1 ) THEN

      DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            u    (i,j,k,tpast)=u    (i,j,k,tpresent)
            v    (i,j,k,tpast)=v    (i,j,k,tpresent)
            w    (i,j,k,tpast)=w    (i,j,k,tpresent)
            pprt (i,j,k,tpast)=pprt (i,j,k,tpresent)
            ptprt(i,j,k,tpast)=ptprt(i,j,k,tpresent)
          END DO
        END DO
      END DO

    END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  Note here defsq is used as a work array by ACOUST.
!
!-----------------------------------------------------------------------
!
  CALL set_acct(smlstp_acct)

  CALL smlstep(mptr, nx,ny,nz, exbcbufsz, dtbig1,dtsml1,                &
               u,v,w,wcont,pprt,ptprt,                                  &
               udteb,udtwb,udtnb,udtsb,vdteb,vdtwb,vdtnb,vdtsb,         &
               wdteb,wdtwb,wdtnb,wdtsb,pdteb,pdtwb,pdtnb,pdtsb,         &
               phydro,sdteb,sdtwb,sdtnb,sdtsb,                          &
               ubar,vbar,ptbar,pbar,ptbari,pbari,rhostr,rhostri,        &
               csndsq,                                                  &
               uforce,vforce,wforce,pforce,ptforce,                     &
               x,y,z,zp,mapfct, j1,j2,j3,aj3x,aj3y,aj3z,j3inv,          &
               trigs1,trigs2,ifax1,ifax2,                               &
               wsave1,wsave2,vwork1,vwork2,                             &
               exbcbuf,rhofct,                                          &
               tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,tem9,defsq,      &
               tem10,tem11,tem12,tem13,tem14,tem15,tem16,tem17,         &
               tem18,tem1_0,tem2_0,tem3_0)

  CALL set_acct(misc_acct)

  IF( lvldbg >= 1 ) THEN
    CALL checkuhx(u(1,1,1,tfuture),                                     &
                  nx,ny,nz,1,nx,1,ny-1,1,nz-1, 'uxbig', tem1)
    CALL checkuhy(u(1,1,1,tfuture),                                     &
                  nx,ny,nz,1,nx,1,ny-1,1,nz-1, 'uybig', tem1)
    CALL checkvhx(v(1,1,1,tfuture),                                     &
                  nx,ny,nz,1,nx-1,1,ny,1,nz-1, 'vxbig', tem1)
    CALL checkvhy(v(1,1,1,tfuture),                                     &
                  nx,ny,nz,1,nx-1,1,ny,1,nz-1, 'vybig', tem1)
    CALL checkwhx(w(1,1,1,tfuture),                                     &
                  nx,ny,nz,1,nx-1,1,ny-1,1,nz, 'wxbig', tem1)
    CALL checkwhy(w(1,1,1,tfuture),                                     &
                  nx,ny,nz,1,nx-1,1,ny-1,1,nz, 'wybig', tem1)
    CALL checkshx(pprt(1,1,1,tfuture),                                  &
                  nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, 'pxbig', tem1)
    CALL checkshy(pprt(1,1,1,tfuture),                                  &
                  nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, 'pybig', tem1)
    CALL checkshx(ptprt(1,1,1,tfuture),                                 &
                  nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, 'ptxbig', tem1)
    CALL checkshy(ptprt(1,1,1,tfuture),                                 &
                  nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, 'ptybig', tem1)
  END IF
!
!-----------------------------------------------------------------------
!
!  Since wcont was reset in SMLSTP. Now need to re-calculate its
!  value at tpresent. Which will be used in SOLVQV and SOLVQ.
!
!-----------------------------------------------------------------------
!
  CALL set_acct(tinteg_acct)

  CALL rhouvw(nx,ny,nz,rhostr,tem1,tem2,tem3)

  CALL wcontra(nx,ny,nz,                                                &
               u(1,1,1,tpresent),v(1,1,1,tpresent),                     &
               w(1,1,1,tpresent),mapfct,j1,j2,j3,aj3z,                  &
               rhostr,tem1,tem2,tem3,wcont,tem4,tem5)

!
!-----------------------------------------------------------------------
!
!  Since pressure was reset in SMLSTP. Now need to re-calculate its
!  value at tfuture which will be used in microphysics.
!
!-----------------------------------------------------------------------
!
  CALL setppi(nx,ny,nz,nt,tfuture,pprt,pbar,ppi)

  IF( lvldbg >= 1 ) THEN
    CALL acct_interrupt(misc_acct)

    CALL checkwhx(wcont, nx,ny,nz,1,nx-1,1,ny-1,2,nz-1,                 &
                  'wcxbig', tem1)
    CALL checkwhy(wcont, nx,ny,nz,1,nx-1,1,ny-1,2,nz-1,                 &
                  'wcybig', tem1)
    CALL acct_stop_inter
  END IF
!
!-----------------------------------------------------------------------
!
!  Integrate the liquid water substance continuity equations
!  forward one timestep.
!
!  Sources and sinks due to phase changes and radiative processes
!  are handled separately.
!
!  If the simulation is for dry dynamics, then we skip over
!  the moisture computations to save computer resources.
!  The flag for this choice is "moist", i.e.,
!
!  Moist = 0 for dry run
!  Moist = 1 for moist run
!
!-----------------------------------------------------------------------
!

  IF( moist == 1) THEN   ! Determine if this run is dry or moist.

!
!-----------------------------------------------------------------------
!
!  Water vapor equation. Array uforce and vforce are used as work
!  arrays.
!  Before we call the solve routine for the water varaibles, we
!  need to compute ustr,vstr,wstr and store them into tem9,10,and 11
!  These arrays cannot be altered and will be used in the
!  following six moisture solve calls.
!
!
!-----------------------------------------------------------------------
!
    IF (sadvopt == 1 .OR. sadvopt == 2 .OR. sadvopt == 3 ) THEN
                             ! 2nd or 4th order advection

      CALL rhouvw(nx,ny,nz,rhostr,tem1,tem2,tem3)

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx
            tem9(i,j,k)=u(i,j,k,2)*tem1(i,j,k)
          END DO
        END DO
      END DO

      DO k=1,nz-1
        DO j=1,ny
          DO i=1,nx-1
            tem10(i,j,k)=v(i,j,k,2)*tem2(i,j,k)
          END DO
        END DO
      END DO

      DO k=1,nz
        DO j=1,ny-1
          DO i=1,nx-1
            tem11(i,j,k)=wcont(i,j,k)*tem3(i,j,k)
          END DO
        END DO
      END DO

    ELSE IF( sadvopt == 4 .OR. sadvopt == 5) THEN  ! FCT advection

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem3_0(i,j,k)=rhostr(i,j,k)
          END DO
        END DO
      END DO

      CALL acct_interrupt(bc_acct)

      CALL extndsbc(tem3_0,nx,ny,nz,0,ebc,wbc,nbc,sbc,tbc,bbc)

      CALL acct_stop_inter

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx
            tem9(i,j,k)=u(i,j,k,2)*(tem3_0(i-1,j,k)+tem3_0(i,j,k))      &
                      *mapfct(i,j,5)*0.5
          END DO
        END DO
      END DO

      DO k=1,nz-1
        DO j=1,ny
          DO i=1,nx-1
            tem10(i,j,k)=v(i,j,k,2)*(tem3_0(i,j-1,k)+tem3_0(i,j,k))     &
                      *mapfct(i,j,6)*0.5
          END DO
        END DO
      END DO

      DO k=1,nz
        DO j=1,ny-1
          DO i=1,nx-1
            tem11(i,j,k)=wcont(i,j,k)                                   &
                      *(tem3_0(i,j,k-1)+tem3_0(i,j,k))*0.5
          END DO
        END DO
      END DO

    END IF

!-----------------------------------------------------------------------
!
!  Solve transport-diffusion equation for water vapor mixing ratio.
!  Arrays uforce and vforce are used as work arrays.
!  Tem9, tem10, tem11 contain ustr,vstr, and wstr calculated
!  earlier, and will be needed gain by next call to solvq.
!
!-----------------------------------------------------------------------

    CALL solvqv(RK3step,nx,ny,nz, exbcbufsz, dtbig1,                    &
                 qv,u,v,wcont, tem9,tem10,tem11,                        &
                 sdteb,sdtwb,sdtnb,sdtsb,                               &
                 rhostr,rhostri,qvbar,kmh,kmv,rprntl,qvsflx,pbldpth,    &
                 x,y,z,zp,mapfct, j1,j2,j3,aj3x,aj3y,j3inv,             &
                 qcumsrc(1,1,1,1),                                      &
                 usflx,vsflx,ptsflx,ptsfc,qvsfc,ptbar,ptprt,            &
                 exbcbuf, bcrlx,                                        &
                 uforce,vforce,                                         &
                 tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,               &
                 tem1_0,tem2_0,tem3_0)

!-----------------------------------------------------------------------
!
!  Solve transport-diffusion equations for moments of cloud and
!  hydrometeors.  Array uforce and vforce are used as work
!  arrays. Tem9, tem10, tem11 contain ustr,vstr, and wstr calculated
!  earlier.  Note that the special treatment of number concentration and
!  reflectivity factor is no longer needed, since the new version of the
!  MY scheme passes out per mass variables instead of per volume.
!
!-----------------------------------------------------------------------

    DO qflag = 1,nscalar

      ! DTD: For the number concentration and reflectivity scalars,
      ! first convert from /m^3 to /kg by dividing through by air density
      convertunits = .FALSE.
      qvarname = qnames(qflag)

      ! Assume that scalars names contains 2 characters and the first
      ! one is 'n' then it is number concentration for MY scheme.
      ! The first character is 'z', then it is reflectivity scalars for MY scheme.

      IF ( LEN_TRIM(qvarname) == 2 .AND.                                &
           (qvarname(1:1) == 'n' .OR. qvarname(1:1) == 'z') ) THEN

        DO k=1,nz
          DO j=1,ny
            DO i=1,nx
              qscalar(i,j,k,1,qflag) = qscalar(i,j,k,1,qflag)*rhostri(i,j,k)*j3(i,j,k)
              qscalar(i,j,k,2,qflag) = qscalar(i,j,k,2,qflag)*rhostri(i,j,k)*j3(i,j,k)
              qscalar(i,j,k,3,qflag) = qscalar(i,j,k,3,qflag)*rhostri(i,j,k)*j3(i,j,k)
            END DO
          END DO
        END DO

        convertunits = .TRUE.

      END IF

      CALL solvq(RK3step,nx,ny,nz, exbcbufsz, dtbig1, qflag,            &
                 qscalar(:,:,:,:,qflag),u,v,wcont, tem9,tem10,tem11,    &
                 sdteb,sdtwb,sdtnb,sdtsb, rhostr,rhostri,               &
                 kmh,kmv,rprntl,x,y,z,zp,mapfct,                        &
                 j1,j2,j3,aj3x,aj3y,j3inv,                              &
                 qcumsrc(1,1,1,2),qcumsrc(1,1,1,3),                     &
                 qcumsrc(1,1,1,4),qcumsrc(1,1,1,5),                     &
                 exbcbuf, bcrlx,                                        &
                 uforce,vforce,                                         &
                 tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,               &
                 tem1_0,tem2_0,tem3_0)

      ! Convert back to /m^3 from /kg for N and Z
      IF ( convertunits ) THEN

        DO k=1,nz
          DO j=1,ny
            DO i=1,nx
              qscalar(i,j,k,1,qflag) = qscalar(i,j,k,1,qflag)*rhostr(i,j,k)*j3inv(i,j,k)
              qscalar(i,j,k,2,qflag) = qscalar(i,j,k,2,qflag)*rhostr(i,j,k)*j3inv(i,j,k)
              qscalar(i,j,k,3,qflag) = qscalar(i,j,k,3,qflag)*rhostr(i,j,k)*j3inv(i,j,k)
            END DO
          END DO
        END DO

      END IF

    END DO

    IF( lvldbg >= 1 ) THEN
      CALL acct_interrupt(misc_acct)
      CALL checkshx(qv(1,1,1,tfuture),                                  &
                    nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, 'qvx', tem1)
      CALL checkshy(qv(1,1,1,tfuture),                                  &
                    nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, 'qvy', tem1)

      IF (P_QC > 0) THEN
        CALL checkshx(qscalar(1,1,1,tfuture,P_QC),                      &
                      nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, 'qcx', tem1)
        CALL checkshy(qscalar(1,1,1,tfuture,P_QC),                      &
                      nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, 'qcy', tem1)
      END IF

      IF (P_QR > 0) THEN
        CALL checkshx(qscalar(1,1,1,tfuture,P_QR),                      &
                      nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, 'qrx', tem1)
        CALL checkshy(qscalar(1,1,1,tfuture,P_QR),                      &
                      nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, 'qry', tem1)
      END IF

      CALL acct_stop_inter
    END IF

  END IF

  RETURN
END SUBROUTINE tinteg
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE TINTEG_RK3                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE tinteg_RK3(mptr, frstep, nx,ny,nz,nzsoil,exbcbufsz,          &
           u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,pbldpth,               &
           udteb,udtwb,udtnb,udtsb,vdteb,vdtwb,vdtnb,vdtsb,             &
           wdteb,wdtwb,wdtnb,wdtsb,pdteb,pdtwb,pdtnb,pdtsb,             &
           phydro,sdteb,sdtwb,sdtnb,sdtsb,                              &
           ubar,vbar,ptbar,pbar,ptbari,pbari,                           &
           rhostr,rhostri,qvbar,ppi,csndsq,                             &
           x,y,z,zp,zpsoil,mapfct,                                      &
           j1,j2,j3,j3soil,aj3x,aj3y,aj3z,j3inv,                        &
           trigs1,trigs2,ifax1,ifax2,                                   &
           wsave1,wsave2,vwork1,vwork2,sinlat,                          &
           ptsfc,qvsfc,prcrate, radfrc,                                 &
           usflx,vsflx,ptsflx,qvsflx,kmh,kmv,rprntl,                    &
           ptcumsrc,qcumsrc,raing,rainc,w0avg,nca,kfraincv,             &
           cldefi,xland,bmjraincv,roufns,                               &
           exbcbuf,bcrlx,                                               &
           rhofct,defsq,lenscl,                                         &
           uforce,vforce,wforce,pforce,ptforce,                         &
           tem1,tem2,tem3,tem4,tem5,tem6,                               &
           tem7,tem8,tem9,tem10,tem11,tem12,tem13,tem14,tem15,          &
           tem16,tem17,tem18,                                           &
           tem1_0,tem2_0,tem3_0)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Orchestrate the time integration of the dynamics of the basic governing
!  equations for a single time step. Note that the time tendencies of the
!  non-conservative processes (e.g., microphysics, radiation, etc.)
!  are not included in this routine but are applied in their
!  corresponding routines.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  3/21/92.
!
!  MODIFICATION HISTORY:
!
!  5/03/92 (K. Droegemeier)
!    Added full documentation.
!
!  5/20/92 (M. Xue)
!    Reformatted certain documentations.
!
!  4/10/93 (M. Xue & Hao Jin)
!    Add the terrain.
!
!  5/26/93 (M. Xue)
!    Added w into the solvpt argument list.
!
!  9/7/94 (M.Xue)
!    Call to ADVP modified.
!
!  9/10/94 (D. Weber & Y. Lu)
!  Cleaned up documentation.
!
!  8/20/95 (M. Xue)
!  Bug fix with the use of ptcumsrc.
!
!  10/31/95 (D. Weber)
!  Added trigs1,trigs2,ifax1,ifax2 for use in the upper w-p
!  radiation condition.
!
!  01/23/1996 (Donghai Wang, Ming Xue and Yuhe Liu)
!  Added the map projection factor to ARPS governing equations.
!
!  03/11/96 (M. Xue)
!  Modified the call to SOLVTKE, so that all three time levels of
!  ptprt, pprt etc are passed into the routine.
!
!  03/19/96 (Yuhe Liu)
!  Added radiation forcing
!
!  07/10/1997 (Fanyou Kong - CMRP)
!  Fixed a bug in 'solvtke' with corrected temporary arga
!
!  07/22/97 (D. Weber)
!  Added wsave1,wsave2,vwork1,vwork2 for use in the even fft version
!  of the upper w-p radiation condition (fftopt=2).
!
!  08/01/97 (Zonghui Huo)
!  Added Kain-fritsch cumulus parameterization scheme.
!
!  10/21/97 (Donghai Wang)
!  Using total density (rho) in the calculation of the pressure
!  gradient force terms, and added the second order terms
!  in the linerized buoyancy terms.
!
!  11/05/97 (D. Weber)
!  Added phydro array for use in the bottom boundary condition for
!  perturbation pressure (hydrostatic).
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
!  9/18/98 (D. Weber)
!  Added precomputed averages of j3 in the x, y, and z directions
!  to improve the code efficiency.
!
!  07/23/2001 (K. Brewster)
!  Added mptr to argument list, needed for printing of diagnostic
!  noise parameter.  Also added mptr to calls to smlstep, solvuw,
!  solvwpex, and solvwpim.
!
!  13 March 2002 (Eric Kemp)
!  Added arrays for WRF BMJ cumulus scheme.
!
!  April 2002 (Fanyou Kong)
!  Added cnvctopt=5 option for the new WRF K-F (KF_ETA) scheme
!
!  10/2009 (Dan Dawson and Ming Xue)
!  Split off from original subroutine.  This version orchestrates
!  the split RK3 time integration scheme (Wicker and Skamarock 2002)
!
!  09/23/2011 (Dan Dawson)
!  Added code to convert the number concentration and reflectivity scalars (for the MY multimoment
!  scheme) to per kg units for the dynamics.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    frstep   Initial time step flag.
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of soil levels
!    u        x component of velocity at times tpast and tpresent (m/s)
!    v        y component of velocity at times tpast and tpresent (m/s)
!    w        Vertical component of Cartesian velocity at times
!             tpast and tpresent (m/s)
!    wcont    Contravariant vertical velocity (m/s)
!             computational coordinates (m/s)
!    ptprt    Perturbation potential temperature at times tpast and
!             tpresent (K)
!    pprt     Perturbation pressure at times tpast and tpresent (Pascal)
!    qv       Water vapor specific humidity at times tpast and tpresent (kg/kg)
!    qc       Cloud water mixing ratio at times tpast and tpresent (kg/kg)
!    qr       Rainwater mixing ratio at times tpast and tpresent (kg/kg)
!    qi       Cloud ice mixing ratio at times tpast and tpresent (kg/kg)
!    qs       Snow mixing ratio at times tpast and tpresent (kg/kg)
!    qh       Hail mixing ratio at times tpast and tpresent (kg/kg)
!
!    udteb    Time tendency of u field at east boundary (m/s**2)
!    udtwb    Time tendency of u field at west boundary (m/s**2)
!
!    vdtnb    Time tendency of v field at north boundary (m/s**2)
!    vdtsb    Time tendency of v field at south boundary (m/s**2)
!
!    pdteb    Time tendency of pprt field at east boundary (Pascal/s)
!    pdtwb    Time tendency of pprt field at west boundary (Pascal/s)
!    pdtnb    Time tendency of pprt field at north boundary (Pascal/s)
!    pdtsb    Time tendency of pprt field at south boundary (Pascal/s)
!
!    phydro   Big time step forcing term for use in computing the
!             hydrostatic pressure at k=1.
!
!    ubar     Base state zonal velocity component (m/s)
!    vbar     Base state meridional velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    ptbari   Inverse Base state potential temperature (K)
!    pbari    Inverse Base state pressure (Pascal)
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    rhostri  Inverse base state density rhobar times j3 (kg/m**3)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!    cnsd     Sound wave speed.
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!    zpsoil   Height of soil levels (m)
!    mapfct   Map factors at scalar, u and v points
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    j3soil   Coordinate transformation Jacobian  d(zpsoil)/dz
!    aj3x     Avgx of the coordinate transformation Jacobian  d(zp)/dz
!    aj3y     Avgy of the coordinate transformation Jacobian  d(zp)/dz
!    aj3z     Avgz of the coordinate transformation Jacobian  d(zp)/dz
!    j3inv    Inverse of the coordinate transformation j3
!    trigs1   Array containing pre-computed trig function for
!             fftopt=1.
!    trigs2   Array containing pre-computed trig function for
!             fftopt=1.
!    ifax1    Array containing the factors of nx for fftopt=1.
!    ifax2    Array containing the factors of ny for fftopt=1.
!
!    vwork1   2-D work array for fftopt=2.
!    vwork2   2-D work array for fftopt=2.
!    wsave1   Work array for fftopt=2.
!    wsave2   Work array for fftopt=2.
!
!    sinlat   Sin of latitude at each grid point
!
!    usflx    Surface flux of u-momentum (kg/(m*s**2))
!    vsflx    Surface flux of v-momentum (kg/(m*s**2))
!    ptsflx   Surface heat flux (K*kg/(m**2 * s ))
!    qvsflx   Surface moisture flux (kg/(m**2 * s))
!
!    ptcumsrc Source term in pt-equation due to cumulus parameterization
!    qcumsrc  Source term in water equations due to cumulus parameterization
!    kfraincv   K-F convective rainfall (cm)
!    nca      K-F counter for CAPE release
!    cldefi   BMJ cloud efficiency
!    xland    BMJ land/sea mask
!    bmjraincv   BMJ convective rainfall (cm)
!
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!    rprntl   Reciprocal of Prandtl number
!
!    radfrc   Radiation forcing (K/s)
!
!  OUTPUT:
!
!    u        x component of velocity at time tfuture (m/s)
!    v        y component of velocity at time tfuture (m/s)
!    w        Vertical component of Cartesian velocity at tfuture (m/s)
!    ptprt    Perturbation potential temperature at time tfuture (K)
!    pprt     Perturbation pressure at time tfuture (Pascal)
!    qv       Water vapor specific humidity at time tfuture (kg/kg)
!    qc       Cloud water mixing ratio at time tfuture (kg/kg)
!    qr       Rainwater mixing ratio at time tfuture (kg/kg)
!    qi       Cloud ice mixing ratio at time tfuture (kg/kg)
!    qs       Snow mixing ratio at time tfuture (kg/kg)
!    qh       Hail mixing ratio at time tfuture (kg/kg)
!
!    udteb    Time tendency of u field at east boundary (m/s**2)
!    udtwb    Time tendency of u field at west boundary (m/s**2)
!
!    vdtnb    Time tendency of v field at north boundary (m/s**2)
!    vdtsb    Time tendency of v field at south boundary (m/s**2)
!
!    pdteb    Time tendency of pprt field at east boundary (Pascal/s)
!    pdtwb    Time tendency of pprt field at west boundary (Pascal/s)
!    pdtnb    Time tendency of pprt field at north boundary (Pascal/s)
!    pdtsb    Time tendency of pprt field at south boundary (Pascal/s)
!
!    phydro   Big time step forcing term for use in computing the
!             hydrostatic pressure at k=1.
!
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!    pbldpth  Planetary boundary layer depth (m)
!
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!    rprntl   Reciprocal of Prandtl number
!
!  WORK ARRAYS:
!
!    defsq    Deformation squared (1/s**2)
!    lenscl   Turbulent mixing length scale (m)
!    uforce   Acoustically inactive forcing terms in u-Eq. ((kg/(m*s)**2)
!    vforce   Acoustically inactive forcing terms in v-Eq. ((kg/(m*s)**2)
!    wforce   Acoustically inactive forcing terms in w-Eq. ((kg/(m*s)**2)
!    pforce   Acoustically inactive forcing terms in p-eq. (Pascal/s)
!    ptforce  Gravity wave inactive forcing terms in pt-eq. (K*kg/(m**3*s))
!
!    rhofct   rho-factor: rhobar/rho
!
!    udtnb    Time tendency of u field at north boundary (m/s**2)
!    udtsb    Time tendency of u field at south boundary (m/s**2)
!
!    vdteb    Time tendency of v field at east boundary (m/s**2)
!    vdtwb    Time tendency of v field at west boundary (m/s**2)
!
!    wdteb    Time tendency of w field at east boundary (m/s**2)
!    wdtwb    Time tendency of w field at west boundary (m/s**2)
!    wdtnb    Time tendency of w field at north boundary (m/s**2)
!    wdtsb    Time tendency of w field at south boundary (m/s**2)
!
!    ptdteb   Time tendency of ptprt field at east boundary (K/s)
!    ptdtwb   Time tendency of ptprt field at west boundary (K/s)
!    ptdtnb   Time tendency of ptprt field at north boundary(K/s)
!    ptdtsb   Time tendency of ptprt field at south boundary(K/s)
!
!    sdteb    Time tendency of a scalar at e-boundary
!    sdtwb    Time tendency of a scalar at w-boundary
!    sdtnb    Time tendency of a scalar at n-boundary
!    sdtsb    Time tendency of a scalar at s-boundary
!
!    tem1     Temporary work array
!    tem2     Temporary work array
!    tem3     Temporary work array
!    tem4     Temporary work array
!    tem5     Temporary work array
!    tem6     Temporary work array
!    tem7     Temporary work array
!    tem8     Temporary work array
!    tem9     Temporary work array
!    tem10    Temporary work array
!    tem11    Temporary work array
!    tem12    Temporary work array
!    tem13    Temporary work array
!    tem14    Temporary work array
!    tem15    Temporary work array
!    tem16    Temporary work array
!    tem17    Temporary work array
!    tem18    Temporary work array
!
!    tem1_0   Temporary work array.
!    tem2_0   Temporary work array.
!    tem3_0   Temporary work array.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE        ! Force explicit declarations
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'bndry.inc'
  INCLUDE 'exbc.inc'
  INCLUDE 'timelvls.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  INTEGER :: mptr

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of soil levels

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)
  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)
  REAL :: tke   (nx,ny,nz,nt)  ! Turbulent Kinetic Energy ((m/s)**2)
  REAL :: pbldpth(nx,ny,nt)    ! Planetary boundary layer depth (m)

  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)
  REAL :: qscalar(nx,ny,nz,nt,nscalar)

  REAL :: udteb (ny,nz)        ! Time tendency of u at e-boundary (m/s**2)
  REAL :: udtwb (ny,nz)        ! Time tendency of u at w-boundary (m/s**2)
  REAL :: udtnb (nx,nz)        ! Time tendency of u at n-boundary (m/s**2)
  REAL :: udtsb (nx,nz)        ! Time tendency of u at s-boundary (m/s**2)

  REAL :: vdteb (ny,nz)        ! Time tendency of v at e-boundary (m/s**2)
  REAL :: vdtwb (ny,nz)        ! Time tendency of v at w-boundary (m/s**2)
  REAL :: vdtnb (nx,nz)        ! Time tendency of v at n-boundary (m/s**2)
  REAL :: vdtsb (nx,nz)        ! Time tendency of v at s-boundary (m/s**2)

  REAL :: wdteb (ny,nz)        ! Time tendency of w at e-boundary (m/s**2)
  REAL :: wdtwb (ny,nz)        ! Time tendency of w at w-boundary (m/s**2)
  REAL :: wdtnb (nx,nz)        ! Time tendency of w at n-boundary (m/s**2)
  REAL :: wdtsb (nx,nz)        ! Time tendency of w at s-boundary (m/s**2)

  REAL :: pdteb (ny,nz)        ! Time tendency of pprt at e-boundary (Pascal/s)
  REAL :: pdtwb (ny,nz)        ! Time tendency of pprt at w-boundary (Pascal/s)
  REAL :: pdtnb (nx,nz)        ! Time tendency of pprt at n-boundary (Pascal/s)
  REAL :: pdtsb (nx,nz)        ! Time tendency of pprt at s-boundary (Pascal/s)

  REAL :: phydro(nx,ny)        ! Pressure at k=1 computed using pert.
                               ! hydrostatic relation.

  REAL :: sdteb (ny,nz)        ! Time tendency of a variable at e-boundary
  REAL :: sdtwb (ny,nz)        ! Time tendency of a variable at w-boundary
  REAL :: sdtnb (nx,nz)        ! Time tendency of a variable at n-boundary
  REAL :: sdtsb (nx,nz)        ! Time tendency of a variable at s-boundary

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: ptbari(nx,ny,nz)     ! Inverse Base state pot. temperature (K)
  REAL :: pbari (nx,ny,nz)     ! Inverse Base state pressure (Pascal).
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: rhostri(nx,ny,nz)    ! Inverse base state density rhobar times j3.
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                               ! (kg/kg)
  REAL :: ppi   (nx,ny,nz)     ! Exner function
  REAL :: csndsq(nx,ny,nz)     ! Sound wave speed squared.

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.
  REAL :: zpsoil(nz,ny,nzsoil) ! Height of soil levels.

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! - d( zp )/d( x ).
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! - d( zp )/d( y ).
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! d( zp )/d( z ).
  REAL :: j3soil(nx,ny,nzsoil) ! Coordinate transformation Jacobian defined as
                               ! d( zpsoil )/d( z ).
  REAL :: aj3x  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE X-DIR.
  REAL :: aj3y  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
  REAL :: aj3z  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Z-DIR.
  REAL :: j3inv (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! d( zp )/d( z ).
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
  REAL :: wsave1 (3*(ny-1)+15) ! Work array for fftopt =2.
  REAL :: wsave2 (3*(nx-1)+15) ! Work array for fftopt =2.

  REAL :: sinlat(nx,ny)        ! Sin of latitude at each grid point

  REAL :: ptsfc  (nx,ny)       ! Potential temperature at the ground level (K)
  REAL :: qvsfc  (nx,ny)       ! Effective qv at the surface (kg/kg)

  REAL :: prcrate(nx,ny,4)     ! precipitation rate (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precipitation rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulus precip. rate
                               ! prcrate(1,1,4) = microphysics precip. rate

  REAL :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m**2*s))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: rprntl(nx,ny,nz)     ! Reciprocal of Prandtl number

  REAL :: defsq (nx,ny,nz)     ! Deformation squared (1/s**2)
  REAL :: lenscl(nx,ny,nz)     ! Turbulent mixing length scale (m)

  REAL :: uforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in u-momentum equation (kg/(m*s)**2)
  REAL :: vforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in v-momentum equation (kg/(m*s)**2)
  REAL :: wforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in w-momentum equation (kg/(m*s)**2)
  REAL :: pforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in pressure equation (Pascal/s)
  REAL :: ptforce(nx,ny,nz)    ! Gravity wave inactive forcing terms
                               ! in pressure equation (Pascal/s)

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

  REAL, INTENT(INOUT) :: cldefi(nx,ny) ! BMJ cloud efficiency
  REAL, INTENT(IN)    :: xland(nx,ny)  ! BMJ land mask
                                       ! (1.0 = land, 2.0 = sea)
  REAL, INTENT(INOUT) :: bmjraincv(nx,ny) ! BMJ convective rainfall
                                          ! (cm)
  REAL, INTENT(IN)    :: roufns(nx,ny)

  REAL :: raing(nx,ny)         ! Grid supersaturation rain
  REAL :: rainc(nx,ny)         ! Cumulus convective rain

  REAL :: radfrc(nx,ny,nz)     ! Radiation forcing (K/s)

  INTEGER :: exbcbufsz         ! EXBC buffer size
  REAL :: exbcbuf( exbcbufsz ) ! EXBC buffer array
  REAL :: bcrlx (nx,ny)        ! EXBC relaxation coefficients

  REAL :: rhofct(nx,ny,nz)     ! rho-factor: rhobar/rho

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array
  REAL :: tem5  (nx,ny,nz)     ! Temporary work array
  REAL :: tem6  (nx,ny,nz)     ! Temporary work array
  REAL :: tem7  (nx,ny,nz)     ! Temporary work array
  REAL :: tem8  (nx,ny,nz)     ! Temporary work array
  REAL :: tem9  (nx,ny,nz)     ! Temporary work array
  REAL :: tem10 (nx,ny,nz)     ! Temporary work array
  REAL :: tem11 (nx,ny,nz)     ! Temporary work array
  REAL :: tem12 (nx,ny,nz)     ! Temporary work array
  REAL :: tem13 (nx,ny,nz)     ! Temporary work array
  REAL :: tem14 (nx,ny,nz)     ! Temporary work array
  REAL :: tem15 (nx,ny,nz)     ! Temporary work array
  REAL :: tem16 (nx,ny,nz)     ! Temporary work array
  REAL :: tem17 (nx,ny,nz)     ! Temporary work array
  REAL :: tem18 (nx,ny,nz)     ! Temporary work array

  REAL :: tem1_0(0:nx,0:ny,0:nz)     ! Temporary work array.
  REAL :: tem2_0(0:nx,0:ny,0:nz)     ! Temporary work array.
  REAL :: tem3_0(0:nx,0:ny,0:nz)     ! Temporary work array.

  ! Allocatable temporary arrays

  REAL,ALLOCATABLE :: utem  (:,:,:,:)
  REAL,ALLOCATABLE :: vtem  (:,:,:,:)
  REAL,ALLOCATABLE :: wtem  (:,:,:,:)
  REAL,ALLOCATABLE :: pprttem  (:,:,:,:)
  REAL,ALLOCATABLE :: ptprttem  (:,:,:,:)
  REAL,ALLOCATABLE :: tketem   (:,:,:,:)
  REAL,ALLOCATABLE :: qvtem (:,:,:,:)
  REAL,ALLOCATABLE :: qscalartem (:,:,:,:,:)

  REAL,ALLOCATABLE :: uadv(:,:,:)       ! Advection forcing in u-momentum
                                        ! equation
  REAL,ALLOCATABLE :: vadv(:,:,:)       ! Advection forcing in v-momentum
                                        ! equation
  REAL,ALLOCATABLE :: wadv(:,:,:)       ! Advection forcing in w-momentum
                                        ! equation
  REAL,ALLOCATABLE :: wbuoy(:,:,:)      ! Buoyancy forcing in w-momentum equation
  REAL,ALLOCATABLE :: ptadv(:,:,:)      ! Advection forcing in pt equation
  REAL,ALLOCATABLE :: padv(:,:,:)       ! Advection forcing in p equation
  REAL,ALLOCATABLE :: tkeadv(:,:,:)     ! Advection forcing in tke equation
  REAL,ALLOCATABLE :: tkeforce(:,:,:)   ! Other forcing in the tke equation
  REAL,ALLOCATABLE :: qvadv(:,:,:)      ! Advection forcing in the qv equation
  REAL,ALLOCATABLE :: qvforce(:,:,:)    ! Other forcing in the qv equation
  REAL,ALLOCATABLE :: qadv(:,:,:,:)    ! Advection forcing in the scalar equation(s)
  REAL,ALLOCATABLE :: qforce(:,:,:,:)  ! Other forcing in the scalar equation(s)

  REAL,ALLOCATABLE :: utot(:,:,:)
  REAL,ALLOCATABLE :: vtot(:,:,:)
  REAL,ALLOCATABLE :: wtot(:,:,:)
  REAL,ALLOCATABLE :: pttot(:,:,:)
  REAL,ALLOCATABLE :: ptot(:,:,:)

  REAL :: mp_tem(MAX(nx+1,ny+1)*(nz+1)) ! Temporary message passing array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL    :: dtbig1            ! Local value of big time step size
  REAL    :: dtsml1            ! Local value of small time step size
  INTEGER :: frstep            ! Flag for the initial time step
  INTEGER :: i,j,k,t
  INTEGER :: qflag             ! Indicator for the water/ice type
                               ! when calling SOLVQ.
  INTEGER :: ntst
  REAL    :: tst
  INTEGER :: RK3step           ! Integer indicating RK3 substep

  INTEGER :: istatus
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!  Added 11/15/2003 - Dan Dawson.  Modified subroutine tinteg for time
!  integration
!  using Runge-Kutta 3rd order (RK3) time differencing after Wicker and
!  Skamarock (2002).
!
!-----------------------------------------------------------------------
!
! Allocate temporary arrays
  ALLOCATE(utem  (nx,ny,nz,nt))
  ALLOCATE(vtem  (nx,ny,nz,nt))
  ALLOCATE(wtem  (nx,ny,nz,nt))
  ALLOCATE(pprttem  (nx,ny,nz,nt))
  ALLOCATE(ptprttem  (nx,ny,nz,nt))
  ALLOCATE(tketem  (nx,ny,nz,nt))
  ALLOCATE(qvtem  (nx,ny,nz,nt))
  ALLOCATE(qscalartem  (nx,ny,nz,nt,nscalar))

  ALLOCATE(uadv  (nx,ny,nz))
  ALLOCATE(vadv  (nx,ny,nz))
  ALLOCATE(wadv  (nx,ny,nz))
  ALLOCATE(wbuoy  (nx,ny,nz))
  ALLOCATE(ptadv  (nx,ny,nz))
  ALLOCATE(padv  (nx,ny,nz))
  IF (tmixopt == 4 .OR. tmixopt == 5) THEN
    ALLOCATE(tkeadv  (nx,ny,nz))
    ALLOCATE(tkeforce  (nx,ny,nz))
  END IF
  IF (moist == 1) THEN
    ALLOCATE(qvadv  (nx,ny,nz))
    ALLOCATE(qvforce  (nx,ny,nz))
    ALLOCATE(qadv  (nx,ny,nz,nscalar))
    ALLOCATE(qforce  (nx,ny,nz,nscalar))
  END IF

  ALLOCATE(utot (nx,ny,nz))
  ALLOCATE(vtot (nx,ny,nz))
  ALLOCATE(wtot (nx,ny,nz))
  ALLOCATE(pttot (nx,ny,nz))
  ALLOCATE(ptot (nx,ny,nz))
!
!-----------------------------------------------------------------------
!
! Initialize the large and small time steps
!
!-----------------------------------------------------------------------
!
  dtbig1 = dtbig
  dtsml1 = dtsml

! There are two options for the RK3 integration.  In the first option
! (tintegopt = 2), all forcing terms except for advection are calculated
! ahead of time and held constant during the RK3 loop.  In the second option
! (tintegopt = 3), all forcing terms (except for radiation, cumulus parameterization,
! surface physics, microphysics, etc.) are updated with each substep of the
! RK3 loop.  The latter is more expensive, but also more accurate and provides
! for better numerical convergence of the temporal truncation error.

  IF(tintegopt == 2) THEN
    RK3step = 0
!
!-----------------------------------------------------------------------
!
!  Calculate wcont at time tpresent. wcont will be used in FRCUVW
!  and FRCP. rhostr averaged to u, v, w points. They are stored
!  in tem1, tem2 and tem3.
!
!-----------------------------------------------------------------------
!
  CALL rhouvw(nx,ny,nz,rhostr,tem1,tem2,tem3)
  CALL wcontra(nx,ny,nz,                                                &
               u(1,1,1,tpresent),v(1,1,1,tpresent),                     &
               w(1,1,1,tpresent),mapfct,j1,j2,j3,aj3z,                  &
               rhostr,tem1,tem2,tem3,wcont,tem4,tem5)

!-----------------------------------------------------------------------
!
! Calculate Planetary boundary layer tendency
!
!-----------------------------------------------------------------------

    !IF (pblopt > 0) THEN
    IF (tmixopt > 4) THEN
      !ALLOCATE(rublten (nx,ny,nz), STAT = istatus)
      !ALLOCATE(rvblten (nx,ny,nz), STAT = istatus)
      !ALLOCATE(rthblten(nx,ny,nz), STAT = istatus)
      !ALLOCATE(rqvblten(nx,ny,nz), STAT = istatus)
      !ALLOCATE(rqcblten(nx,ny,nz), STAT = istatus)
      !ALLOCATE(rqiblten(nx,ny,nz), STAT = istatus)
      !rublten  = 0.0
      !rvblten  = 0.0
      !rthblten = 0.0
      !rqvblten = 0.0
      !rqcblten = 0.0
      !rqiblten = 0.0

      CALL pbl_driver(tintegopt,tmixopt-4,nx,ny,nz,                       &
                    u,v,w,ptprt,pprt,ppi,qv,qscalar,                      &
                    ptbar,pbar,qvbar,rhostr,dtbig1,                       &
                    !rublten,rvblten,rthblten,rqvblten,rqcblten,rqiblten,  &
                    zp,roufns,ptsflx,qvsflx,xland,ptsfc,                  &
                    pbldpth,kmv,rprntl,                                   &
                    tem1,tem2,tem3,tem4,tem5,tem6,tem7,                   &
                    tem8,tem9,tem10,tem11,tem12,istatus)
    ELSE

      !ALLOCATE(rublten (1,1,1), STAT = istatus)
      !ALLOCATE(rvblten (1,1,1), STAT = istatus)
      !ALLOCATE(rthblten(1,1,1), STAT = istatus)
      !ALLOCATE(rqvblten(1,1,1), STAT = istatus)
      !ALLOCATE(rqcblten(1,1,1), STAT = istatus)
      !ALLOCATE(rqiblten(1,1,1), STAT = istatus)
      !rublten  = 0.0
      !rvblten  = 0.0
      !rthblten = 0.0
      !rqvblten = 0.0
      !rqcblten = 0.0
      !rqiblten = 0.0

    END IF
!
!-----------------------------------------------------------------------
!
!  Compute the acoustically inactive terms in the momentum and
!  pressure equations that are held fixed during the small time
!  step computations.  This includes buoyancy and mixing
!  (both physical and computational), and the Coriolis terms.
!  (The advection terms will be updated within the RK3 loop).
!  Each of the forcing terms are calculated at time tpresent.
!  (In the RK3 scheme, time level tpast is not used).
!  These forcing terms are accumulated into arrays for each
!  of the momentum equations, e.g., uforce for the u-equation,
!  vforce for the v-equation, wforce for the w-equation and
!  pforce for the p-equation.
!
!  Note: Arrays, pforce and ptforce, are used as work arrays in
!        subroutine frcuvw.
!
!-----------------------------------------------------------------------
!
  CALL frcuvw(nx,ny,nz,nzsoil,exbcbufsz,dtbig1,                         &
              u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,pbldpth,            &
              ubar,vbar,ptbar,pbar,ptbari,pbari,rhostr,qvbar,           &
              usflx,vsflx, x,y,z,zp,zpsoil,mapfct,                      &
              j1,j2,j3,j3soil,aj3x,aj3y,aj3z,j3inv, sinlat, ptsfc,      &
              uforce,vforce,wforce,kmh,kmv,rprntl,lenscl,defsq,         &
              exbcbuf, bcrlx,rhofct,phydro,                             &
              pforce,ptforce,                                           &
              tem1,tem2,tem3,tem4,tem5,tem6,tem7,                       &
              tem8,tem9,tem10,tem11,tem12,tem13,tem14)
!
!-----------------------------------------------------------------------
!
!  Calculate gravity wave or acoustic wave inactive terms in the
!  potential temperature equation.
!
!  The force terms are stored in ptforce.
!
!-----------------------------------------------------------------------
!
  CALL frcpt(nx,ny,nz, exbcbufsz, dtbig1,ptprt,u,v,w,wcont,             &
             ptbar,rhostr,rhostri,kmh,kmv,rprntl,                       &
             usflx,vsflx,ptsflx,pbldpth,                                &
             x,y,z,zp,mapfct, j1,j2,j3,aj3x,aj3y,j3inv,ptsfc,           &
             ptforce,                                                   &
             exbcbuf, bcrlx,                                            &
             tem1,tem2,tem3,tem4,tem5,tem6,tem7,                        &
             tem8,tem9,tem10,tem11,                                     &
             tem1_0,tem2_0,tem3_0,tem12)

  CALL set_acct(tinteg_acct)

  ELSE
    ptforce = 0.0  ! Make sure this is initialized in the case of tintegopt == 3
  END IF

  IF ( radopt == 2 ) THEN
    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          ptforce(i,j,k) = ptforce(i,j,k)                               &
                         + rhostr(i,j,k)*radfrc(i,j,k)
        END DO
      END DO
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate source and sink terms in temperature (ptprt) and
!  moisture (qv) equations that are due to subgrid scale cumulus
!  convection.
!
!-----------------------------------------------------------------------
!
  IF (cnvctopt == 1 .OR. cnvctopt == 2 .OR. cnvctopt == 3 .OR. &
      cnvctopt == 4 .OR. cnvctopt == 5) THEN

    CALL set_acct(cum_acct)

!-----------------------------------------------------------------------
!
! Calculate w0avg, is close to a running mean vertical velocity,
! tst is the number of time steps in 10min for K-F scheme
!
!-----------------------------------------------------------------------
!
    IF (cnvctopt == 3 .OR. cnvctopt == 5) THEN

      ntst=nint( 600.0/dtbig )
      tst=FLOAT(ntst)

      IF( (curtim-tstart) <= 300.0 .AND. initopt /= 2 ) THEN
        DO k=1,nz
          DO j=1,ny-1
            DO i=1,nx-1
              w0avg(i,j,k)= (w(i,j,k,1)+w(i,j,k,2))*0.5
            END DO
          END DO
        END DO

      ELSE IF ( (curtim-tstart) > 300.0 .OR. initopt == 2) THEN

        DO k=1,nz
          DO j=1,ny-1
            DO i=1,nx-1
              w0avg(i,j,k)=(w0avg(i,j,k)*(tst-1.)+w(i,j,k,2))/tst
            END DO
          END DO
        END DO
      END IF
    END IF
!
!-----------------------------------------------------------------------
!
!  Calculate vertical velocity for KUO scheme, or running-average
!  vertical velocity (m/s) for Kain Fritsch scheme
!
!-----------------------------------------------------------------------
!
!    IF( cnvctopt == 3 .AND. MOD(curtim+0.001,confrq) <= (0.5*dtbig)     &
!             .AND. ((curtim-tstart) > 300.0                             &
    IF( (cnvctopt == 3 .OR. cnvctopt == 5) .AND. MOD(curtim+0.001,confrq) &
             <= (0.5*dtbig) .AND. ((curtim-tstart) > 300.0                             &
                    .OR. initopt == 2)) THEN

      DO k=1,nz
        DO j=1,ny-1
          DO i=1,nx-1
!        tem1(i,j,k) = (w(i,j,k,1)+w(i,j,k,2))*0.5
            tem1(i,j,k) = w0avg(i,j,k)
          END DO
        END DO
      END DO


    ELSE IF( MOD(curtim+0.001,confrq) <= (0.5*dtbig).OR.                &
              nstep == 1 ) THEN
      DO k=1,nz
        DO j=1,ny-1
          DO i=1,nx-1
            tem1(i,j,k) = wcont(i,j,k)*j3(i,j,k)
          END DO
        END DO
      END DO
    END IF
!
!-----------------------------------------------------------------------
!
!  Call cumulus parameterization schemes
!  Make sure to reset kfraincv, ptcumsrc and qcumsrc to 0 once nca<=0
!
!-----------------------------------------------------------------------
!
    IF( (moist /= 0) .AND. ((curtim-tstart) > 300.0 .OR. initopt == 2) .AND. &
          (MOD(curtim+0.001,confrq) <= (0.5*dtbig))) THEN

      DO j=1,ny-1
        DO i=1,nx-1
          IF (nca(i,j) <= 0) THEN
            kfraincv(i,j) = 0.0
            prcrate(i,j,3) = 0.0
          END IF
        END DO
      END DO

      DO k=2,nz-2
        DO j=1,ny-1
          DO i=1,nx-1
            IF (nca(i,j) <= 0) THEN
              ptcumsrc(i,j,k) = 0.0
              qcumsrc(i,j,k,1) = 0.0
              qcumsrc(i,j,k,2) = 0.0
              qcumsrc(i,j,k,3) = 0.0
              qcumsrc(i,j,k,4) = 0.0
              qcumsrc(i,j,k,5) = 0.0
            END IF
          END DO
        END DO
      END DO

      CALL qpfcums(nx,ny,nz,prcrate(1,1,3),qvsflx,                      &
                   u(1,1,1,2),v(1,1,1,2),tem1,                          &
                   pprt(1,1,1,2),ptprt(1,1,1,2),qv(1,1,1,2),            &
                   pbar,ptbar,qvbar,rhostr,zp,j3,                       &
                   ptcumsrc,qcumsrc,rainc,nca,kfraincv,                 &
                   cldefi,xland,bmjraincv,                              &
                   tem2,tem3,tem4,tem5,tem6,                            &
                   tem7,tem8,tem9,tem10,tem11)
    END IF
!
!  Accumulate rainc and update nca. Note: raincv is in cm.
!
    DO j=1,ny-1
      DO i=1,nx-1
        rainc(i,j) = rainc(i,j) + 10.0*kfraincv(i,j) + 10.0*bmjraincv(i,j)
        IF ( nca(i,j) >= 1 ) nca(i,j) = nca(i,j) - 1
      END DO
    END DO

    DO k=2,nz-2
      DO j=1,ny-1
        DO i=1,nx-1
          ptforce(i,j,k) = ptforce(i,j,k)                               &
                         + rhostr(i,j,k)*ptcumsrc(i,j,k)
        END DO
      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  End of cumulus parameterization
!
!-----------------------------------------------------------------------
!
  CALL set_acct(misc_acct)

  IF( lvldbg >= 2 ) THEN
    CALL checkuhx(uforce, nx,ny,nz,2,nx-1,1,ny-1,2,nz-2,                &
                  'ufrcx', tem1)
    CALL checkuhy(uforce, nx,ny,nz,2,nx-1,1,ny-1,2,nz-2,                &
                  'ufrcy', tem1)
    CALL checkvhx(vforce, nx,ny,nz,1,nx-1,2,ny-1,2,nz-2,                &
                  'vfrcx', tem1)
    CALL checkvhy(vforce, nx,ny,nz,1,nx-1,2,ny-1,2,nz-2,                &
                  'vfrcy', tem1)
    CALL checkwhx(wforce, nx,ny,nz,1,nx-1,1,ny-1,2,nz-1,                &
                  'wfrcx', tem1)
    CALL checkwhy(wforce, nx,ny,nz,1,nx-1,1,ny-1,2,nz-1,                &
                  'wfrcy', tem1)
    CALL checkshx(pforce, nx,ny,nz,1,nx-1,1,ny-1,2,nz-2,                &
                  'pfrcx', tem1)
    CALL checkshy(pforce, nx,ny,nz,1,nx-1,1,ny-1,2,nz-2,                &
                  'pfrcy', tem1)
    CALL checkshx(ptforce, nx,ny,nz,1,nx-1,1,ny-1,2,nz-2,               &
                  'ptfrcx', tem1)
    CALL checkshy(ptforce, nx,ny,nz,1,nx-1,1,ny-1,2,nz-2,               &
                  'ptfrcy', tem1)
  END IF
!
!-----------------------------------------------------------------------
!
!  Integrate the momentum and pressure equations (i.e., the
!  acoustically-active equations) in time using a mode-splitting
!  approach.  The momentum components are advanced in time using
!  a forward scheme (relative to the pressure gradient force terms)
!  and the pressure is advanced using a backward scheme (relative to
!  the divergence term, which is evaluated using the newly updated
!  velocities hence the backward scheme). These small time steps
!  bring the momentum and pressure from time tpresent to tfuture.
!  During these steps, the acoustically-inactive forcing terms (e.g.,
!  uforce, pforce, etc.) are held fixed (i.e. they
!  are evaluated only once, for each RK3 large time
!  step integration).  For the first iteration, they are evaluated at
!  time tpresent, for the second, at time tpresent + dtbig/3, and for
!  the third, time tpresent + dtbig/2.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  When this solver is called for a nested grid (not the base grid),
!  the information stored in u, v, w and pprt is used to set the
!  time tendency of these variables at the lateral boundaries.
!
!-----------------------------------------------------------------------
!
  IF( mgrid /= 1 .AND. nestgrd == 1 ) THEN

    CALL nestbdt(nx,ny,nz,u, 1,nx,1,ny-1,1,nz-1,dtbig,                  &
                 udteb,udtwb,udtnb,udtsb)

    CALL nestbdt(nx,ny,nz,v, 1,nx-1,1,ny,1,nz-1,dtbig,                  &
                 vdteb,vdtwb,vdtnb,vdtsb)

    CALL nestbdt(nx,ny,nz,w, 1,nx-1,1,ny-1,1,nz,dtbig,                  &
                 wdteb,wdtwb,wdtnb,wdtsb)

    CALL nestbdt(nx,ny,nz,pprt, 1,nx-1,1,ny-1,1,nz-1,dtbig,             &
                 pdteb,pdtwb,pdtnb,pdtsb)

    CALL nestbdt(nx,ny,nz,ptprt,1,nx-1,1,ny-1,1,nz-1,dtbig,             &
                 sdteb,sdtwb,sdtnb,sdtsb)

!
!-----------------------------------------------------------------------
!
!    For the first forward time step, make sure that the fields at
!    tpast equal those at tpresent. This can only be done after
!    the time tendencies on the boundaries are calculated and stored
!    in the time tendency arrays. (Not needed for RK3)
!
!-----------------------------------------------------------------------
!
!    IF( frstep == 1 ) THEN

!      DO k=1,nz
!        DO j=1,ny
!          DO i=1,nx
!            u    (i,j,k,tpast)=u    (i,j,k,tpresent)
!            v    (i,j,k,tpast)=v    (i,j,k,tpresent)
!            w    (i,j,k,tpast)=w    (i,j,k,tpresent)
!            pprt (i,j,k,tpast)=pprt (i,j,k,tpresent)
!            ptprt(i,j,k,tpast)=ptprt(i,j,k,tpresent)
!          END DO
!        END DO
!      END DO

!    END IF

  END IF
!
!-----------------------------------------------------------------------
!
!  Begin the Runge-Kutta 3rd order time integration loop
!
!-----------------------------------------------------------------------
!

  CALL set_acct(tinteg_acct)

  ! Assign the values of the dependent variables to temporary arrays which will be used in the Runge-Kutta procedure

  DO t=1,3
   DO k=1,nz
    DO j=1,ny
     DO i=1,nx
      utem(i,j,k,t) = u(i,j,k,t)
      vtem(i,j,k,t) = v(i,j,k,t)
      wtem(i,j,k,t) = w(i,j,k,t)
      pprttem(i,j,k,t) = pprt(i,j,k,t)
      ptprttem(i,j,k,t) = ptprt(i,j,k,t)
      IF(tmixopt == 4 .OR. tmixopt == 5) THEN
        tketem(i,j,k,t) = tke(i,j,k,t)
      END IF
      IF(moist == 1) THEN
        qvtem(i,j,k,t) = qv(i,j,k,t)
        DO qflag=1,nscalar
          qscalartem(i,j,k,t,qflag) = qscalar(i,j,k,t,qflag)
        END DO
      END IF
     END DO
    END DO
   END DO
  END DO

  DO t=1,3

    ! Zero out the arrays containing the advection forcing terms, as they
    ! are recalculated each iteration.  Also store the dependent variables at the proper time level into
    ! temporary arrays to be passed to the advection subroutines.
    ! In the case of tintegopt == 3, these arrays will also hold the mixing terms

    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          uadv(i,j,k) = 0.0
          vadv(i,j,k) = 0.0
          wadv(i,j,k) = 0.0
          ptadv(i,j,k) = 0.0
          pforce(i,j,k) = 0.0
          IF(tmixopt == 4 .OR. tmixopt == 5) THEN
            tkeadv(i,j,k) = 0.0
          END IF
          IF(moist == 1) THEN
            qvadv(i,j,k) = 0.0
            DO qflag=1,nscalar
              qadv(i,j,k,qflag) = 0.0
            END DO
          END IF
        END DO
      END DO
    END DO

    ! Note, the below adjustments to dtsml within the given RK3 interval
    ! ensure that the number of small time steps fits evenly within the
    ! intervals.
    IF(t == 1) THEN
      dtbig1 = dtbig/3.0
      dtsml1 = MIN(dtsml,dtbig1)
      ! Check to see if new dtsml1 fits evenly within the dtbig/3 interval
      ! and adjust it accordingly
      i = MAX(1,nint(dtbig1/dtsml1))
      dtsml1 = dtbig1/i
      !dtsml1=dtsml/3.0
      RK3step = 1
    ELSE IF(t == 2) THEN
      dtbig1 = dtbig/2.0
      dtsml1 = MIN(dtsml,dtbig1)
      !dtsml1 = dtsml/2.0
      RK3step = 2
    ELSE IF(t == 3) THEN
      dtbig1 = dtbig
      dtsml1 = dtsml
      RK3step = 3
    END IF

    CALL rhouvw(nx,ny,nz,rhostr,tem1,tem2,tem3)

    CALL wcontra(nx,ny,nz,                                              &
               utem(1,1,1,tpresent),vtem(1,1,1,tpresent),               &
               wtem(1,1,1,tpresent),mapfct,j1,j2,j3,aj3z,               &
               rhostr,tem1,tem2,tem3,wcont,tem4,tem5)


    IF(tintegopt == 2) THEN

    ! Calculate the advection terms in the momentum equations

      CALL set_acct(advuvw_acct)

      CALL advuvw(nx,ny,nz,utem,vtem,wtem,wcont,rhostr,ubar,vbar, mapfct,  &
              uadv,vadv,wadv,                                           &
              tem1,tem2,tem3,                                           &
              tem4,tem5,tem6,tem7,tem8,tem9)

    ELSE IF(tintegopt == 3) THEN

!-----------------------------------------------------------------------
!
! Calculate Planetary boundary layer tendency
!
!-----------------------------------------------------------------------

      !IF (pblopt > 0) THEN
      IF (tmixopt > 4) THEN
        !ALLOCATE(rublten (nx,ny,nz), STAT = istatus)
        !ALLOCATE(rvblten (nx,ny,nz), STAT = istatus)
        !ALLOCATE(rthblten(nx,ny,nz), STAT = istatus)
        !ALLOCATE(rqvblten(nx,ny,nz), STAT = istatus)
        !ALLOCATE(rqcblten(nx,ny,nz), STAT = istatus)
        !ALLOCATE(rqiblten(nx,ny,nz), STAT = istatus)
        !rublten  = 0.0
        !rvblten  = 0.0
        !rthblten = 0.0
        !rqvblten = 0.0
        !rqcblten = 0.0
        !rqiblten = 0.0

        CALL pbl_driver(tintegopt,tmixopt-4,nx,ny,nz,                       &
                      u,v,w,ptprt,pprt,ppi,qv,qscalar,                      &
                      ptbar,pbar,qvbar,rhostr,dtbig1,                       &
                      !rublten,rvblten,rthblten,rqvblten,rqcblten,rqiblten,  &
                      zp,roufns,ptsflx,qvsflx,xland,ptsfc,                  &
                      pbldpth,kmv,rprntl,                                   &
                      tem1,tem2,tem3,tem4,tem5,tem6,tem7,                   &
                      tem8,tem9,tem10,tem11,tem12,istatus)
      ELSE

        !ALLOCATE(rublten (1,1,1), STAT = istatus)
        !ALLOCATE(rvblten (1,1,1), STAT = istatus)
        !ALLOCATE(rthblten(1,1,1), STAT = istatus)
        !ALLOCATE(rqvblten(1,1,1), STAT = istatus)
        !ALLOCATE(rqcblten(1,1,1), STAT = istatus)
        !ALLOCATE(rqiblten(1,1,1), STAT = istatus)
        !rublten  = 0.0
        !rvblten  = 0.0
        !rthblten = 0.0
        !rqvblten = 0.0
        !rqcblten = 0.0
        !rqiblten = 0.0

      END IF
    !
    !-----------------------------------------------------------------------
    !
    !  Compute the acoustically inactive terms in the momentum and
    !  pressure equations that are held fixed during the small time
    !  step computations.  This includes advection, buoyancy, mixing
    !  (both physical and computational), and the Coriolis terms.
    !  These forcing terms are accumulated into arrays for each
    !  of the momentum equations, e.g., uforce for the u-equation,
    !  vforce for the v-equation, wforce for the w-equation and
    !  pforce for the p-equation.
    !
    !  Note: Arrays, ptot and pttot, are used as work arrays in
    !        subroutine frcuvw.
    !
    !-----------------------------------------------------------------------
    !
      CALL frcuvw(nx,ny,nz,nzsoil,exbcbufsz,dtbig1,                         &
              utem,vtem,wtem,wcont,ptprttem,pprttem,qvtem,qscalartem,tketem,pbldpth,            &
              ubar,vbar,ptbar,pbar,ptbari,pbari,rhostr,qvbar,           &
              usflx,vsflx, x,y,z,zp,zpsoil,mapfct,                      &
              j1,j2,j3,j3soil,aj3x,aj3y,aj3z,j3inv, sinlat, ptsfc,      &
              uadv,vadv,wadv,kmh,kmv,rprntl,lenscl,defsq,               &
              exbcbuf, bcrlx,rhofct,phydro,                             &
              ptot,pttot,                                           &
              tem1,tem2,tem3,tem4,tem5,tem6,tem7,                       &
              tem8,tem9,tem10,tem11,tem12,tem13,tem14)

    END IF

    CALL set_acct(misc_acct)

    IF( lvldbg >= 4 ) THEN
      CALL checkuhx(uforce, nx,ny,nz,2,nx-1,1,ny-1,2,nz-2,              &
                  'uforce after advu', tem9)
      CALL checkuhy(uforce, nx,ny,nz,2,nx-1,1,ny-1,2,nz-2,              &
                  'uforce after advu', tem9)
      CALL checkvhx(vforce, nx,ny,nz,1,nx-1,2,ny-1,2,nz-2,              &
                  'vforce after advv', tem9)
      CALL checkvhy(vforce, nx,ny,nz,1,nx-1,2,ny-1,2,nz-2,              &
                  'vforce after advv', tem9)
      CALL checkwhx(wforce, nx,ny,nz,1,nx-1,1,ny-1,2,nz-1,              &
                  'wforce after advw', tem9)
      CALL checkwhy(wforce, nx,ny,nz,1,nx-1,1,ny-1,2,nz-1,              &
                  'wforce after advw', tem9)
    END IF

!
!-----------------------------------------------------------------------
!
!  Integrate the TKE equation.
!  NOTE: Arrays lenscl and defsq should NOT be changed between
!  calls to frcuvw and solvtke.
!
!  Upon exit of this subroutine, tke will be updated by one RK3
!  substep.  In the case of tintegopt == 2, and if this is the first RK3
!  step, tkeforce will contain all the forcing terms except advection
!  calculated at the beginning of the RK3 interval.  These will be held
!  constant for the rest of the RK3 steps, so tkeforce MUST NOT BE
!  CHANGED until after the RK3 loop.  In the case of tintegopt = 3,
!  all the forcing terms are updated for each RK3 subinterval.
!-----------------------------------------------------------------------
!
    IF (tmixopt == 4 .OR. tmixopt == 5) THEN
      CALL solvtke(RK3step,nx,ny,nz,dtbig1,utem,vtem,wcont,ptprttem,pprttem, &
                 qvtem,qscalartem,tke,                                  &
                 ubar,vbar,ptbar,pbar,ptbari,rhostr,rhostri,qvbar,      &
                 x,y,z,zp, mapfct, j1,j2,j3,aj3x,aj3y,j3inv,            &
                 kmh,kmv,rprntl,lenscl,defsq,                           &
                 ptsflx,qvsflx,                                         &
                 tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,               &
                 tem9,tem10,tem11,tkeadv,tkeforce,                      &
                 tem1_0,tem2_0,tem3_0)
    END IF
    ! Calculate the forcing terms for the pressure equation.  Since this
    ! includes only advection, we don't need to make a separate call to
    ! advp here, and it holds for both tinteg options.

    CALL frcp(nx,ny,nz,exbcbufsz, dtbig1,                               &
            utem,vtem,wtem,wcont,ptprttem,pprttem,qvtem,qscalartem,     &
            ptbar,pbar,rhostr,qvbar,mapfct,j1,j2,j3,aj3x,aj3y,aj3z,     &
            pforce,                                                     &
            exbcbuf, bcrlx,                                             &
            tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,tem9)

    IF(tintegopt == 3) THEN

!
!-----------------------------------------------------------------------
!
!  Calculate gravity wave or acoustic wave inactive terms in the
!  potential temperature equation.
!
!  The force terms are stored in ptadv.
!
!-----------------------------------------------------------------------
!
    CALL frcpt(nx,ny,nz, exbcbufsz, dtbig1,ptprttem,utem,vtem,wtem,wcont,             &
            ptbar,rhostr,rhostri,kmh,kmv,rprntl,                       &
            usflx,vsflx,ptsflx,pbldpth,                                &
            x,y,z,zp,mapfct, j1,j2,j3,aj3x,aj3y,j3inv,ptsfc,           &
            ptadv,                                                     &
            exbcbuf, bcrlx,                                            &
            tem1,tem2,tem3,tem4,tem5,tem6,tem7,                        &
            tem8,tem9,tem10,tem11,                                     &
            tem1_0,tem2_0,tem3_0,tem12)

    ELSE

    ! Calculate the advection terms in the potential temperature equation

    CALL set_acct(advs_acct)

    ! For sadvopt 4, to save computational expense, for now we will only
    ! Call the FCT advection subroutine for the final RK3 step.
    ! For the first 2 steps, the normal advection subroutine is called

    IF(sadvopt == 4 .and. t == 3) THEN

      CALL advptfct(nx,ny,nz,dtbig1,ptprttem,utem,vtem,wtem,wcont,      &
                  rhostr,rhostri,ptbar,mapfct,j3,j3inv,                 &
                  ptadv,                                                &
                  tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,              &
                  tem9,tem10,tem11,tem1_0,tem2_0,tem3_0,mp_tem)
    ELSE
      CALL advpt(nx,ny,nz,ptprttem,utem,vtem,wtem,wcont, rhostr,ptbar,mapfct,         &
                  j3,j3inv,ptadv,                                          &
                  tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8)
    END IF

    END IF

    CALL set_acct(tinteg_acct)
    ! Sum the advection (and mixing for tintegopt = 3) terms with the other forcing terms that were
    ! previously calculated outside of the RK3 loop.  Store the result
    ! in the total forcing arrays.  This is to ensure that the total forcings
    ! are recalculated each iteration of the loop

    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          IF(tintegopt == 2) THEN
            utot(i,j,k) = uadv(i,j,k) + uforce(i,j,k)
            vtot(i,j,k) = vadv(i,j,k) + vforce(i,j,k)
            wtot(i,j,k) = wadv(i,j,k) + wforce(i,j,k)
            pttot(i,j,k) = -ptadv(i,j,k) + ptforce(i,j,k)
          ELSE
            utot(i,j,k) = uadv(i,j,k)
            vtot(i,j,k) = vadv(i,j,k)
            wtot(i,j,k) = wadv(i,j,k)
            pttot(i,j,k) = ptadv(i,j,k) + ptforce(i,j,k)
          END IF
          ptot(i,j,k) = pforce(i,j,k)
        END DO
      END DO
    END DO

    ! Treat the vertically-implicit mixing.  This is done here because we need
    ! the updated advective forcings.

    IF(tintegopt == 2 .and. trbvimp == 1) THEN

      ! Vertically-implicit mixing for u,v,w

      CALL set_acct(tmix_acct)

      CALL vmiximpuvw(nx,ny,nz,dtbig1,utem(1,1,1,1),vtem(1,1,1,1),wtem(1,1,1,1),   &
                    rhostr,kmv,j1,j2,j3inv,                             &
                    utot,vtot,wtot,                               &
                    tem1,tem2,tem3,tem4,                                &
                    tem5,tem6,tem7,tem8,tem9,tem10,tem11)

      ! Vertically-implicit mixing for pt

      CALL set_acct(tmix_acct)

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem1(i,j,k)=ptbar(i,j,k)+ptprttem(i,j,k,tpresent)
          END DO
        END DO
      END DO

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem2(i,j,k)=rhostr(i,j,k)*kmv(i,j,k)*rprntl(i,j,k)            &
                        *j3inv(i,j,k)*j3inv(i,j,k)
          END DO
        END DO
      END DO

      CALL vmiximps(nx,ny,nz,dtbig1*0.5,tem1,rhostr,tem2,                 &
                    pttot,tem3,tem4,tem5,tem6)

    END IF

    CALL set_acct(smlstp_acct)

    !Call the forward-backward time integration scheme for the small
    !time steps, using the forcing terms calculated previously.  After
    !smlstep exits, all dependent variables should be updated to the
    !current future time in the RK3 scheme.  The next iteration will
    !then recalculate the forcings based on these updated variables

    CALL smlstep(mptr, nx,ny,nz, exbcbufsz, dtbig1,dtsml1,              &
               u,v,w,wcont,pprt,ptprt,                                  &
               udteb,udtwb,udtnb,udtsb,vdteb,vdtwb,vdtnb,vdtsb,         &
               wdteb,wdtwb,wdtnb,wdtsb,pdteb,pdtwb,pdtnb,pdtsb,         &
               phydro,sdteb,sdtwb,sdtnb,sdtsb,                          &
               ubar,vbar,ptbar,pbar,ptbari,pbari,rhostr,rhostri,        &
               csndsq,                                                  &
               utot,vtot,wtot,ptot,pttot,                         &
               x,y,z,zp,mapfct, j1,j2,j3,aj3x,aj3y,aj3z,j3inv,          &
               trigs1,trigs2,ifax1,ifax2,                               &
               wsave1,wsave2,vwork1,vwork2,                             &
               exbcbuf,rhofct,                                          &
               tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,tem9,defsq,      &
               tem10,tem11,tem12,tem13,tem14,tem15,tem16,tem17,         &
               tem18,tem1_0,tem2_0,tem3_0)

    CALL set_acct(misc_acct)

  IF( lvldbg >= 1 ) THEN
    CALL checkuhx(u(1,1,1,tfuture),                                     &
                  nx,ny,nz,1,nx,1,ny-1,1,nz-1, 'uxbig', tem1)
    CALL checkuhy(u(1,1,1,tfuture),                                     &
                  nx,ny,nz,1,nx,1,ny-1,1,nz-1, 'uybig', tem1)
    CALL checkvhx(v(1,1,1,tfuture),                                     &
                  nx,ny,nz,1,nx-1,1,ny,1,nz-1, 'vxbig', tem1)
    CALL checkvhy(v(1,1,1,tfuture),                                     &
                  nx,ny,nz,1,nx-1,1,ny,1,nz-1, 'vybig', tem1)
    CALL checkwhx(w(1,1,1,tfuture),                                     &
                  nx,ny,nz,1,nx-1,1,ny-1,1,nz, 'wxbig', tem1)
    CALL checkwhy(w(1,1,1,tfuture),                                     &
                  nx,ny,nz,1,nx-1,1,ny-1,1,nz, 'wybig', tem1)
    CALL checkshx(pprt(1,1,1,tfuture),                                  &
                  nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, 'pxbig', tem1)
    CALL checkshy(pprt(1,1,1,tfuture),                                  &
                  nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, 'pybig', tem1)
    CALL checkshx(ptprt(1,1,1,tfuture),                                 &
                  nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, 'ptxbig', tem1)
    CALL checkshy(ptprt(1,1,1,tfuture),                                 &
                  nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, 'ptybig', tem1)
  END IF

!
!-----------------------------------------------------------------------
!
!  Since wcont was reset in SMLSTP. Now need to re-calculate its
!  value at tpresent. Which will be used in SOLVQV and SOLVQ.
!
!-----------------------------------------------------------------------
!
  CALL set_acct(tinteg_acct)

  CALL rhouvw(nx,ny,nz,rhostr,tem1,tem2,tem3)

  CALL wcontra(nx,ny,nz,                                                &
               utem(1,1,1,tpresent),vtem(1,1,1,tpresent),               &
               wtem(1,1,1,tpresent),mapfct,j1,j2,j3,aj3z,               &
               rhostr,tem1,tem2,tem3,wcont,tem4,tem5)
!
!-----------------------------------------------------------------------
!
!  Integrate the liquid water substance continuity equations
!  forward one timestep.
!
!  Sources and sinks due to phase changes and radiative processes
!  are handled separately.
!
!  If the simulation is for dry dynamics, then we skip over
!  the moisture computations to save computer resources.
!  The flag for this choice is "moist", i.e.,
!
!  Moist = 0 for dry run
!  Moist = 1 for moist run
!
!-----------------------------------------------------------------------
!

  IF( moist == 1) THEN   ! Determine if this run is dry or moist.

!
!-----------------------------------------------------------------------
!
!  Water vapor equation. Array uforce and vforce are used as work
!  arrays.
!  Before we call the solve routine for the water varaibles, we
!  need to compute ustr,vstr,wstr and store them into tem9,10,and 11
!  These arrays cannot be altered and will be used in the
!  following six moisture solve calls.
!
!  In the case of sadvopt == 4, we will only use the FCT advection
!  subroutine for the final RK3 step, and the normal advection subroutine
!  for the first 2 steps.
!
!-----------------------------------------------------------------------
!
    IF (sadvopt == 1 .OR. sadvopt == 2 .OR. sadvopt == 3 .or. (sadvopt == 4 .and. t /= 3)) THEN
                             ! 2nd or 4th order advection
!      CALL rhouvw(nx,ny,nz,rhostr,tem1,tem2,tem3)

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx
            tem9(i,j,k)=utem(i,j,k,2)*tem1(i,j,k)
          END DO
        END DO
      END DO

      DO k=1,nz-1
        DO j=1,ny
          DO i=1,nx-1
            tem10(i,j,k)=vtem(i,j,k,2)*tem2(i,j,k)
          END DO
        END DO
      END DO

      DO k=1,nz
        DO j=1,ny-1
          DO i=1,nx-1
            tem11(i,j,k)=wcont(i,j,k)*tem3(i,j,k)
          END DO
        END DO
      END DO

    ELSE IF( (sadvopt == 4 .OR. sadvopt == 5) .AND. t == 3) THEN  ! FCT advection

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem3_0(i,j,k)=rhostr(i,j,k)
          END DO
        END DO
      END DO

      CALL acct_interrupt(bc_acct)

      CALL extndsbc(tem3_0,nx,ny,nz,0,ebc,wbc,nbc,sbc,tbc,bbc)

      CALL acct_stop_inter

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx
            tem9(i,j,k)=utem(i,j,k,2)*(tem3_0(i-1,j,k)+tem3_0(i,j,k))      &
                      *mapfct(i,j,5)*0.5
          END DO
        END DO
      END DO

      DO k=1,nz-1
        DO j=1,ny
          DO i=1,nx-1
            tem10(i,j,k)=vtem(i,j,k,2)*(tem3_0(i,j-1,k)+tem3_0(i,j,k))     &
                      *mapfct(i,j,6)*0.5
          END DO
        END DO
      END DO

      DO k=1,nz
        DO j=1,ny-1
          DO i=1,nx-1
            tem11(i,j,k)=wcont(i,j,k)                                   &
                      *(tem3_0(i,j,k-1)+tem3_0(i,j,k))*0.5
          END DO
        END DO
      END DO

    END IF

!-----------------------------------------------------------------------
!
!  Solve transport-diffusion equation for water vapor mixing ratio.
!  Arrays uforce and vforce are used as work arrays.
!  Tem9, tem10, tem11 contain ustr,vstr, and wstr calculated
!  earlier, and will be needed gain by next call to solvq.
!
!  At the end of this subroutine, qv will be updated by one RK3 substep
!-----------------------------------------------------------------------

    CALL solvqv(t,nx,ny,nz, exbcbufsz, dtbig1,                          &
                 qv,utem,vtem,wcont, tem9,tem10,tem11,                  &
                 sdteb,sdtwb,sdtnb,sdtsb,                               &
                 rhostr,rhostri,qvbar,kmh,kmv,rprntl,qvsflx,pbldpth,    &
                 x,y,z,zp,mapfct, j1,j2,j3,aj3x,aj3y,j3inv,             &
                 qcumsrc(1,1,1,1),                                      &
                 usflx,vsflx,ptsflx,ptsfc,qvsfc,ptbar,ptprttem,         &
                 exbcbuf, bcrlx,                                        &
                 qvadv,qvforce,                                         &
                 tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,               &
                 tem1_0,tem2_0,tem3_0)

!-----------------------------------------------------------------------
!
!  Solve transport-diffusion equations for the moments of cloud and
!  hydrometeors.  Array uforce and vforce are used as work
!  arrays. Tem9, tem10, tem11 contain ustr,vstr, and wstr calculated
!  earlier.Note that the special treatment of number concentration and
!  reflectivity factor is no longer needed, since the new version of the
!  MY scheme passes out per mass variables instead of per volume.
!
!  At the end of this subroutine, the scalars in qscalar will be updated
!  by one RK3 substep
!-----------------------------------------------------------------------

    DO qflag = 1,nscalar

      ! DTD: For the number concentration and reflectivity scalars,
      ! first convert from /m^3 to /kg by dividing through by air density

      IF(qflag >= 7) THEN ! Assume that scalars from 7 and up are Nt or Z; this may need to be made
                          ! more general in the future, but is correct for MY scheme which is the only
                          ! scheme that has number concentration and reflectivity scalars at this time.
      DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            qscalar(i,j,k,1,qflag) = qscalar(i,j,k,1,qflag)*rhostri(i,j,k)*j3(i,j,k)
            qscalar(i,j,k,2,qflag) = qscalar(i,j,k,2,qflag)*rhostri(i,j,k)*j3(i,j,k)
            qscalar(i,j,k,3,qflag) = qscalar(i,j,k,3,qflag)*rhostri(i,j,k)*j3(i,j,k)
          END DO
        END DO
      END DO

      END IF

      CALL solvq(t,nx,ny,nz, exbcbufsz, dtbig1, qflag,            &
                 qscalar(:,:,:,:,qflag),utem,vtem,wcont, tem9,tem10,tem11,    &
                 sdteb,sdtwb,sdtnb,sdtsb, rhostr,rhostri,               &
                 kmh,kmv,rprntl,x,y,z,zp,mapfct,                        &
                 j1,j2,j3,aj3x,aj3y,j3inv,                              &
                 qcumsrc(1,1,1,2),qcumsrc(1,1,1,3),                     &
                 qcumsrc(1,1,1,4),qcumsrc(1,1,1,5),                     &
                 exbcbuf, bcrlx,                                        &
                 qadv(1,1,1,qflag),qforce(1,1,1,qflag),                 &
                 tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,               &
                 tem1_0,tem2_0,tem3_0)

      ! Convert back to /m^3 from /kg for N and Z

      IF(qflag >= 7) THEN ! Assume that scalars from 7 and up are Nt or Z; this may need to be made
                          ! more general in the future, but is correct for MY scheme which is the only
                          ! scheme that has number concentration and reflectivity scalars at this time.
        DO k=1,nz
          DO j=1,ny
            DO i=1,nx
              qscalar(i,j,k,1,qflag) = qscalar(i,j,k,1,qflag)*rhostr(i,j,k)*j3inv(i,j,k)
              qscalar(i,j,k,2,qflag) = qscalar(i,j,k,2,qflag)*rhostr(i,j,k)*j3inv(i,j,k)
              qscalar(i,j,k,3,qflag) = qscalar(i,j,k,3,qflag)*rhostr(i,j,k)*j3inv(i,j,k)
            END DO
          END DO
        END DO

      END IF

    END DO

    IF( lvldbg >= 1 ) THEN
      CALL acct_interrupt(misc_acct)
      CALL checkshx(qv(1,1,1,tfuture),                                  &
                    nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, 'qvx', tem1)
      CALL checkshy(qv(1,1,1,tfuture),                                  &
                    nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, 'qvy', tem1)

      IF (P_QC > 0) THEN
        CALL checkshx(qscalar(1,1,1,tfuture,P_QC),                      &
                      nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, 'qcx', tem1)
        CALL checkshy(qscalar(1,1,1,tfuture,P_QC),                      &
                      nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, 'qcy', tem1)
      END IF

      IF (P_QR > 0) THEN
        CALL checkshx(qscalar(1,1,1,tfuture,P_QR),                      &
                      nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, 'qrx', tem1)
        CALL checkshy(qscalar(1,1,1,tfuture,P_QR),                      &
                      nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, 'qry', tem1)
      END IF

      CALL acct_stop_inter
    END IF

  END IF ! moist = 1
    !Swap the dependent variable's time levels in the temporary arrays to prepare for the next Runge-Kutta integration

    DO k=1,nz
     DO j=1,ny
      DO i=1,nx
        utem(i,j,k,tpresent) = u(i,j,k,tfuture)
        vtem(i,j,k,tpresent) = v(i,j,k,tfuture)
        wtem(i,j,k,tpresent) = w(i,j,k,tfuture)
        pprttem(i,j,k,tpresent) = pprt(i,j,k,tfuture)
        ptprttem(i,j,k,tpresent) = ptprt(i,j,k,tfuture)
        IF(tmixopt == 4 .OR. tmixopt == 5 ) THEN
          tketem(i,j,k,tpresent) = tke(i,j,k,tfuture)
        END IF
        IF(moist == 1) THEN
          qvtem(i,j,k,tpresent) = qv(i,j,k,tfuture)
          DO qflag = 1,nscalar
            qscalartem(i,j,k,tpresent,qflag) = qscalar(i,j,k,tfuture,qflag)
          END DO
        END IF
      END DO
     END DO
    END DO

    CALL set_acct(tinteg_acct)

  END DO ! RK3 loop
!
!-----------------------------------------------------------------------
!
!  Since pressure was reset in SMLSTP. Now need to re-calculate its
!  value at tfuture which will be used in microphysics.
!
!-----------------------------------------------------------------------
!
  CALL setppi(nx,ny,nz,nt,tfuture,pprt,pbar,ppi)

  IF( lvldbg >= 1 ) THEN
    CALL acct_interrupt(misc_acct)

    CALL checkwhx(wcont, nx,ny,nz,1,nx-1,1,ny-1,2,nz-1,                 &
                  'wcxbig', tem1)
    CALL checkwhy(wcont, nx,ny,nz,1,nx-1,1,ny-1,2,nz-1,                 &
                  'wcybig', tem1)
    CALL acct_stop_inter
  END IF

! Deallocate temporary arrays

  DEALLOCATE(utem)
  DEALLOCATE(vtem)
  DEALLOCATE(wtem)
  DEALLOCATE(pprttem)
  DEALLOCATE(ptprttem)
  DEALLOCATE(tketem)
  DEALLOCATE(qvtem)
  DEALLOCATE(qscalartem)

  DEALLOCATE(uadv)
  DEALLOCATE(vadv)
  DEALLOCATE(wadv)
  DEALLOCATE(wbuoy)
  DEALLOCATE(ptadv)
  DEALLOCATE(padv)
  IF (tmixopt == 4 .OR. tmixopt == 5)  THEN
    DEALLOCATE(tkeadv)
    DEALLOCATE(tkeforce)
  END IF
  IF (moist == 1) THEN
    DEALLOCATE(qvadv)
    DEALLOCATE(qvforce)
    DEALLOCATE(qadv)
    DEALLOCATE(qforce)
  END IF

  DEALLOCATE(utot)
  DEALLOCATE(vtot)
  DEALLOCATE(wtot)
  DEALLOCATE(pttot)
  DEALLOCATE(ptot)

  RETURN
END SUBROUTINE tinteg_RK3
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SMLSTEP                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE smlstep(mptr, nx,ny,nz, exbcbufsz, dtbig1,dtsml1,            &
           u,v,w,wcont,pprt,ptprt,                                      &
           udteb,udtwb,udtnb,udtsb,vdteb,vdtwb,vdtnb,vdtsb,             &
           wdteb,wdtwb,wdtnb,wdtsb,pdteb,pdtwb,pdtnb,pdtsb,             &
           phydro,ptdteb,ptdtwb,ptdtnb,ptdtsb,                          &
           ubar,vbar,ptbar,pbar,ptbari,pbari,rhostr,rhostri,            &
           csndsq,                                                      &
           uforce,vforce,wforce,pforce,ptforce,                         &
           x,y,z,zp,mapfct, j1,j2,j3,aj3x,aj3y,aj3z,j3inv,              &
           trigs1,trigs2,ifax1,ifax2,                                   &
           wsave1,wsave2,vwork1,vwork2,                                 &
           exbcbuf,rhofct,                                              &
           rhostru,rhostrv,rhostrw,                                     &
           tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,tem9,tem10,          &
           tem11,tem12,tem13,tem14,tem15,tem16,tem17,tem18,tem19)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Coordinate the integration of the acoustically active parts of
!  the dynamic equations.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  3/21/92.
!
!  MODIFICATION HISTORY:
!
!  5/20/92 (K. Droegemeier and M. Xue)
!  Added full documentation.
!
!  4/10/93 (M. Xue & Hao Jin)
!  Add the terrain.
!
!  9/10/94 (D. Weber & Y. Lu)
!  Cleaned up documentation.
!
!  10/31/95 (D. Weber)
!  Added trigs1,trigs2,ifax1,ifax2 for use in the upper w-p
!  radiation condition.
!
!  1/25/96 (Donghai Wang & Yuhe Liu)
!  Added the map projection factor to ARPS governing equations.
!
!  07/22/97 (D. Weber)
!  Added wsave1,wsave2,vwork1,vwork2 for use in the even fft version
!  of the upper w-p radiation condition (fftopt=2).
!
!  10/21/97 (Donghai Wang)
!  Using total density (rho) in the calculation of the pressure
!  gradient force terms, and added the second order terms
!  in the linerized buoyancy terms.
!
!  11/05/97 (D. Weber)
!  Added phydro array for use in the bottom boundary condition for
!  perturbation pressure (hydrostatic).
!
!  11/06/97 (D. Weber)
!  Added three additional levels to the mapfct array.  The three
!  levels (4,5,6) represent the inverse of the first three in order.
!  The inverse map factors are computed to improve efficiency.
!
!  9/18/98 (D. Weber)
!  Added arrays aj3x,y,z, rhostri,ptbari,pbari, and tem9-12 for
!  use in optimizing this routine.
!  Tem9-19 is used to store commonly used groups of constants.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    dtbig1   Large time step size for this call.
!    dtsml1   Small time step size for this call.
!
!    u        x component of velocity at times tpast and tpresent (m/s)
!    v        y component of velocity at times tpast and tpresent (m/s)
!    w        Vertical component of Cartesian velocity at times
!             tpast and tpresent (m/s)
!    wcont    Contravariant vertical velocity (m/s)
!    pprt     Perturbation pressure at times tpast and tpresent (Pascal)
!    ptprt    Perturbation potential temperature at times tpast and
!             tpresent (K)
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
!    pdteb    Time tendency of pprt field at east boundary (Pascal/s)
!    pdtwb    Time tendency of pprt field at west boundary (Pascal/s)
!    pdtnb    Time tendency of pprt field at north boundary (Pascal/s)
!    pdtsb    Time tendency of pprt field at south boundary (Pascal/s)
!
!    ptdteb   Time tendency of ptprt field at east boundary (K/s)
!    ptdtwb   Time tendency of ptprt field at west boundary (K/s)
!    ptdtnb   Time tendency of ptprt field at north boundary(K/s)
!    ptdtsb   Time tendency of ptprt field at south boundary(K/s)
!
!    phydro   Big time step forcing term for use in computing the
!             hydrostatic pressure at k=1.
!
!    ubar     Base state zonal velocity component (m/s)
!    vbar     Base state meridional velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    ptbari   Inverse base state potential temperature (K)
!    pbari    Inverse base state pressure (Pascal)
!    rhostr   Base state density rhobar times j3 (kg/m**3)
!    rhostri  Inverse base state density rhobar times j3 (kg/m**3)
!    csndsq   Sound wave speed squared.
!
!    uforce   Acoustically inactive forcing terms in u-eq. (kg/(m*s)**2)
!    vforce   Acoustically inactive forcing terms in v-eq. (kg/(m*s)**2)
!    wforce   Acoustically inactive forcing terms in w-eq. (kg/(m*s)**2)
!    pforce   Acoustically inactive forcing terms in p-eq. (Pascal/s)
!    ptforce  Gravity wave inactive forcing terms in pt-eq. (K*kg/(m**3*s))
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space (m)
!
!    mapfct   Map factors at scalar, u and v points
!
!    j1       Coordinate transformation Jacobian -d(zp)/dx
!    j2       Coordinate transformation Jacobian -d(zp)/dy
!    j3       Coordinate transformation Jacobian  d(zp)/dz
!    aj3x     Avgx of the coordinate transformation Jacobian  d(zp)/dz
!    aj3y     Avgy of the coordinate transformation Jacobian  d(zp)/dz
!    aj3z     Avgz of the coordinate transformation Jacobian  d(zp)/dz
!    j3inv    Inverse of the coordinate transformation j3
!    trigs1   Array containing pre-computed trig function for
!             fftopt=1.
!    trigs2   Array containing pre-computed trig function for
!             fftopt=1.
!    ifax1    Array containing the factors of nx for fftopt=1.
!    ifax2    Array containing the factors of ny for fftopt=1.
!
!    vwork1   2-D work array for fftopt = 2.
!    vwork2   2-D work array for fftopt = 2.
!    wsave1   Work array for fftopt = 2.
!    wsave2   Work array for fftopt = 2.
!
!    rhostru  Averaged rhostr at u points (kg/m**3).
!    rhostrv  Averaged rhostr at v points (kg/m**3).
!    rhostrw  Averaged rhostr at w points (kg/m**3).
!
!  OUTPUT:
!
!    u        x component of velocity at time tfuture (m/s)
!    v        y component of velocity at time tfuture (m/s)
!    w        Vertical component of Cartesian velocity at tfuture (m/s)
!    pprt     Perturbation pressure at time tfuture (Pascal)
!    ptprt    Perturbation potential temperature at time tfuture (K)
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array
!    tem2     Temporary work array
!    tem3     Temporary work array
!    tem4     Temporary work array
!    tem5     Temporary work array
!    tem6     Temporary work array
!    tem7     Temporary work array
!    tem8     Temporary work array
!    tem9     Temporary work array
!    tem10    Temporary work array
!    tem11    Temporary work array
!    tem12    Temporary work array
!    tem13    Temporary work array
!    tem14    Temporary work array
!    tem15    Temporary work array
!    tem16    Temporary work array
!    tem17    Temporary work array
!    tem18    Temporary work array
!    tem19    Temporary work array
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE        ! Force explicit declarations

  INCLUDE 'timelvls.inc'

  INTEGER :: mptr

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  REAL :: dtbig1               ! Large time step size for this call.
  REAL :: dtsml1               ! Small time step size for this call.

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: wcont (nx,ny,nz)     ! Contravariant vertical velocity (m/s)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)
  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)

  REAL :: udteb (ny,nz)        ! Time tendency of u at e-boundary (m/s**2)
  REAL :: udtwb (ny,nz)        ! Time tendency of u at w-boundary (m/s**2)
  REAL :: udtnb (nx,nz)        ! Time tendency of u at n-boundary (m/s**2)
  REAL :: udtsb (nx,nz)        ! Time tendency of u at s-boundary (m/s**2)

  REAL :: vdteb (ny,nz)        ! Time tendency of v at e-boundary (m/s**2)
  REAL :: vdtwb (ny,nz)        ! Time tendency of v at w-boundary (m/s**2)
  REAL :: vdtnb (nx,nz)        ! Time tendency of v at n-boundary (m/s**2)
  REAL :: vdtsb (nx,nz)        ! Time tendency of v at s-boundary (m/s**2)

  REAL :: wdteb (ny,nz)        ! Time tendency of w at e-boundary (m/s**2)
  REAL :: wdtwb (ny,nz)        ! Time tendency of w at w-boundary (m/s**2)
  REAL :: wdtnb (nx,nz)        ! Time tendency of w at n-boundary (m/s**2)
  REAL :: wdtsb (nx,nz)        ! Time tendency of w at s-boundary (m/s**2)

  REAL :: pdteb (ny,nz)        ! Time tendency of pprt at e-boundary (Pascal/s)
  REAL :: pdtwb (ny,nz)        ! Time tendency of pprt at w-boundary (Pascal/s)
  REAL :: pdtnb (nx,nz)        ! Time tendency of pprt at n-boundary (Pascal/s)
  REAL :: pdtsb (nx,nz)        ! Time tendency of pprt at s-boundary (Pascal/s)

  REAL :: phydro(nx,ny)        ! Big time step forcing for computing
                               ! hydrostatic pprt at k=1.

  REAL :: ptdteb(ny,nz)        ! Time tendency of ptprt at e-boundary (K/s)
  REAL :: ptdtwb(ny,nz)        ! Time tendency of ptprt at w-boundary (K/s)
  REAL :: ptdtnb(nx,nz)        ! Time tendency of ptprt at n-boundary (K/s)
  REAL :: ptdtsb(nx,nz)        ! Time tendency of ptprt at s-boundary (K/s)

  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: ptbari(nx,ny,nz)     ! Inverse base state pot. temperature (K)
  REAL :: pbari (nx,ny,nz)     ! Inverse base state pressure (Pascal).
  REAL :: rhostr(nx,ny,nz)     ! Base state density rhobar times j3.
  REAL :: rhostri(nx,ny,nz)    ! Inverse base state density rhobar times j3.
  REAL :: csndsq(nx,ny,nz)     ! Sound wave speed squared.

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined at
                               ! w-point of the staggered grid.

  REAL :: mapfct(nx,ny,8)      ! Map factors at scalar, u and v points

  REAL :: j1    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! - d( zp )/d( x ).
  REAL :: j2    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! - d( zp )/d( y ).
  REAL :: j3    (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! d( zp )/d( z ).
  REAL :: aj3x  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE X-DIR.
  REAL :: aj3y  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Y-DIR.
  REAL :: aj3z  (nx,ny,nz)     ! Coordinate transformation Jacobian defined
                               ! as d( zp )/d( z ) AVERAGED IN THE Z-DIR.
  REAL :: j3inv (nx,ny,nz)     ! Coordinate transformation Jacobian defined as
                               ! d( zp )/d( z ).

  REAL :: trigs1(3*(nx-1)/2+1) ! Array containing pre-computed trig
                               ! function for fftopt=1.
  REAL :: trigs2(3*(ny-1)/2+1) ! Array containing pre-computed trig
                               ! function for fftopt=1.
  INTEGER :: ifax1(13)         ! Array containing the factors of nx for
                               ! fftopt=1.
  INTEGER :: ifax2(13)         ! Array containing the factors of ny for
                               ! fftopt=1.

  REAL :: vwork1 (nx+1,ny+1)   ! 2-D work array for fftopt = 2.
  REAL :: vwork2 (ny,nx+1)     ! 2-D work array for fftopt = 2.
  REAL :: wsave1 (3*(ny-1)+15) ! Work array for fftopt = 2.
  REAL :: wsave2 (3*(nx-1)+15) ! Work array for fftopt = 2.

  REAL :: uforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in u-momentum equation (kg/(m*s)**2)
  REAL :: vforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in v-momentum equation (kg/(m*s)**2)
  REAL :: wforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in w-momentum equation (kg/(m*s)**2)
  REAL :: pforce(nx,ny,nz)     ! Acoustically inactive forcing terms
                               ! in pressure equation (Pascal/s)
  REAL :: ptforce(nx,ny,nz)    ! Gravity wave inactive forcing terms
                               ! in potential temperature eq. (K*kg/(m**3*s))

  REAL :: rhostru(nx,ny,nz)    ! Averaged rhostr at u points (kg/m**3).
  REAL :: rhostrv(nx,ny,nz)    ! Averaged rhostr at v points (kg/m**3).
  REAL :: rhostrw(nx,ny,nz)    ! Averaged rhostr at w points (kg/m**3).

  INTEGER :: exbcbufsz         ! EXBC buffer size
  REAL :: exbcbuf( exbcbufsz ) ! EXBC buffer array

  REAL :: rhofct(nx,ny,nz)     ! rho-factor: rhobar/rho

  REAL :: tem1  (nx,ny,nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
  REAL :: tem3  (nx,ny,nz)     ! Temporary work array
  REAL :: tem4  (nx,ny,nz)     ! Temporary work array
  REAL :: tem5  (nx,ny,nz)     ! Temporary work array
  REAL :: tem6  (nx,ny,nz)     ! Temporary work array
  REAL :: tem7  (nx,ny,nz)     ! Temporary work array
  REAL :: tem8  (nx,ny,nz)     ! Temporary work array
  REAL :: tem9  (nx,ny,nz)     ! Temporary work array
  REAL :: tem10 (nx,ny,nz)     ! Temporary work array
  REAL :: tem11 (nx,ny,nz)     ! Temporary work array
  REAL :: tem12 (nx,ny,nz)     ! Temporary work array
  REAL :: tem13 (nx,ny,nz)     ! Temporary work array
  REAL :: tem14 (nx,ny,nz)     ! Temporary work array
  REAL :: tem15 (nx,ny,nz)     ! Temporary work array
  REAL :: tem16 (nx,ny,nz)     ! Temporary work array
  REAL :: tem17 (nx,ny,nz)     ! Temporary work array
  REAL :: tem18 (nx,ny,nz)     ! Temporary work array
  REAL :: tem19 (nx,ny,nz)     ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k
  INTEGER :: ismstp

  REAL :: curtsml               ! Current time during small time step
                                ! integration.
  REAL :: tem,tema,temb
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'bndry.inc'
  INCLUDE 'exbc.inc'
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
!  Integrate the momentum equations using a forward-in-time (relative
!  to the pressure gradient terms) scheme and the pressure equation
!  using a backward-in-time (relative to the mass divergence term)
!  scheme.
!
!  dtsml = (2*dtbig)/nsmstp is used so that after nsmstp small time
!  step integrations, the fields are brought from time t-dtbig
!  to t+dtbig, i.e. from time tpast to time future.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  The first time through the small timestep loop, set the values
!  at time tfuture equal to those at tpast.  By doing this, we
!  simplify the forward-in-time integration because the equations
!  take the form:
!
!  u(i,j,k,future) = u(i,j,k,future)
!    :                + dtsml*( acoustic +  non-acoustic forcing terms )
!
!  i.e., we only use a single time level in this two-time level
!  scheme and automatically finish with u at time tfuture.
!
!-----------------------------------------------------------------------
!
!call test_dump(w(1,1,tfuture),"XXXsml_wfut",nx,ny,nz,3,1)
!call test_dump(w(1,1,tpast),"XXXsml_wpast",nx,ny,nz,3,1)
!call test_dump(w(1,1,tpresent),"XXXsml_wpres",nx,ny,nz,3,1)

 IF(tintegopt == 1) THEN
  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        u    (i,j,k,tfuture)=u    (i,j,k,tpast)
        v    (i,j,k,tfuture)=v    (i,j,k,tpast)
        w    (i,j,k,tfuture)=w    (i,j,k,tpast)
        pprt (i,j,k,tfuture)=pprt (i,j,k,tpast)
      END DO
    END DO
  END DO
 ELSE IF(tintegopt == 2 .or. tintegopt == 3) THEN
  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        u    (i,j,k,tfuture)=u    (i,j,k,tpresent)
        v    (i,j,k,tfuture)=v    (i,j,k,tpresent)
        w    (i,j,k,tfuture)=w    (i,j,k,tpresent)
        pprt (i,j,k,tfuture)=pprt (i,j,k,tpresent)
      END DO
    END DO
  END DO
 END IF


  IF(ptsmlstp == 1) THEN

    IF( sadvopt == 4 ) THEN

!-----------------------------------------------------------------------
!
!  When FCT is used to advect ptprt, all (e.g., advection,
!  mixing) except for ptbar advection terms were evaluated at tpresent,
!  and the large time step integration is forward in time.
!  We have to do somthing special here, so that the same ptprt(tfuture)
!  results had the ptbar advection not been included when ptprt is
!  is updated in small steps. Therefore, we require
!
!         (ptprt(tfuture)-ptprt(tpast))/(2*dtbig1)
!       = (ptprt(tfuture)-ptprt(tpresent))/dtbig
!
!  TODO: Update this section of code for RK3 integration
!-----------------------------------------------------------------------

      tem = 2*dtbig1/dtbig

      IF(nestgrd == 1.AND.mgrid /= 1) THEN
                ! Don't touch the boundary values at tpast for nested grids.

        DO k=1,nz-1
          DO j=2,ny-2
            DO i=2,nx-2
              ptprt(i,j,k,tpast)=ptprt(i,j,k,tfuture)                   &
                   -(ptprt(i,j,k,tfuture)-ptprt(i,j,k,tpresent))*tem
            END DO
          END DO
        END DO

      ELSE

        DO k=1,nz-1
          DO j=1,ny-1
            DO i=1,nx-1
              ptprt(i,j,k,tpast)=ptprt(i,j,k,tfuture)                   &
                   -(ptprt(i,j,k,tfuture)-ptprt(i,j,k,tpresent))*tem
            END DO
          END DO
        END DO

      END IF

    END IF

    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
           IF(tintegopt == 1) THEN
             ptprt(i,j,k,tfuture)=ptprt(i,j,k,tpast)
           ELSE IF(tintegopt == 2 .or. tintegopt == 3) THEN
             ptprt(i,j,k,tfuture)=ptprt(i,j,k,tpresent)
           END IF
        END DO
      END DO
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
!  Calculate wcont at future time level.
!
!-----------------------------------------------------------------------
!
  CALL rhouvw(nx,ny,nz,rhostr,rhostru,rhostrv,rhostrw)

  CALL wcontra(nx,ny,nz,                                                &
             u(1,1,1,tfuture),v(1,1,1,tfuture),                         &
             w(1,1,1,tfuture),mapfct,j1,j2,j3,aj3z,                     &
             rhostr,rhostru,rhostrv,rhostrw,wcont,tem1,tem2)

  IF( lvldbg >= 3 ) THEN
    CALL checkwhx(wcont, nx,ny,nz,1,nx-1,1,ny-1,2,nz-1,                 &
                  'wcxsml', tem1)
    CALL checkwhy(wcont, nx,ny,nz,1,nx-1,1,ny-1,2,nz-1,                 &
                  'wcysml', tem1)
  END IF

!-----------------------------------------------------------------------
!
!  Compute a number of static variables (wrt small time step) and
!  pass into the appropriate subroutines.
!  NOTE:  ALL OF THE COMPUTATIONS ASSUME THAT THE VARIABLES USED
!         ARE CONSTANT DURING THE SMALL TIME STEP!!!!!!!
!         THEY CANNOT BE OVERWRITTEN UNTIL AFTER THE W-P SOLVER.
!
!-----------------------------------------------------------------------

  temb = 0.5*dtsml1
  DO k=2,nz-2
    DO j=1,ny-1
      DO i=2,nx-1
        tema = 1.0/rhostru(i,j,k)
        tem13(i,j,k) =temb*(rhofct(i,j,k)+rhofct(i-1,j,k))              &
                          *mapfct(i,j,2)*tema
        uforce(i,j,k) = dtsml1*uforce(i,j,k)*tema
      END DO
    END DO
  END DO

  DO k=2,nz-2
    DO j=2,ny-1
      DO i=1,nx-1
        tema = 1.0/rhostrv(i,j,k)
        tem14(i,j,k) =temb*(rhofct(i,j,k)+rhofct(i,j-1,k))              &
                           *mapfct(i,j,3)*tema
        vforce(i,j,k) = dtsml1*vforce(i,j,k)*tema
      END DO
    END DO
  END DO

  DO k=2,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem17(i,j,k)=1.0/(pbar(i,j,k)+pbar(i,j,k-1)) ! used in solvwp
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem12(i,j,k)=1.0/rhostrw(i,j,k)              ! used in solvwp
        tem9 (i,j,k)=rhostr(i,j,k)*j3inv(i,j,k)      ! used in solvwp
        tem10(i,j,k)=csndsq(i,j,k)*j3inv(i,j,k)*tem9(i,j,k)  ! solvwp
        tem11(i,j,k)=rhostr(i,j,k)*pbari(i,j,k)      ! used in solvwp
        tem18(i,j,k)=rhostr(i,j,k)*ptbari(i,j,k)     ! used in solvwp
      END DO
    END DO
  END DO

  IF(vimplct == 0)THEN  ! pre-compute wforce, etc.. for wp-ex.
    DO k=2,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem19(i,j,k) =temb*(rhofct(i,j,k)+rhofct(i,j,k-1))            &
                             *tem12(i,j,k)
          wforce(i,j,k) = dtsml1*wforce(i,j,k)*tem12(i,j,k)
        END DO
      END DO
    END DO
  END IF

  tema = dtsml1*dxinv
  DO k=2,nz-2
    DO j=1,ny-1
      DO i=1,nx
        tem15(i,j,k)=tema*aj3x(i,j,k)*mapfct(i,j,5) ! used in solvwp's
      END DO
    END DO
  END DO

  temb = dtsml1*dyinv
  DO k=2,nz-2
    DO j=1,ny
      DO i=1,nx-1
        tem16(i,j,k)=temb*aj3y(i,j,k)*mapfct(i,j,6) ! used in solvwp's
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!    For peqopt=2, set
!
!    pforce=pforce (containging exbc damping only)
!      + rhostr*cs**2/ptbar*d(pt)/dt
!       =pforce (containging exbc damping only)
!      + cs**2/ptbar * (ptforce + buoyancy if not already included).
!
!    Note that adiabet heating/colling is not included in ptforce,
!    therefore is not accounted for in p equation.
!
!-----------------------------------------------------------------------
!
  IF(peqopt == 2) THEN

    IF( ptsmlstp == 0 ) THEN
      DO k=2,nz-2
        DO j=1,ny-1
          DO i=1,nx-1
            pforce(i,j,k)=pforce(i,j,k)                                 &
                         +csndsq(i,j,k)*ptforce(i,j,k)*ptbari(i,j,k)
          END DO
        END DO
      END DO
    ELSE
      DO k=2,nz-2
        DO j=1,ny-1
          DO i=1,nx-1
            pforce(i,j,k)=pforce(i,j,k)+csndsq(i,j,k)*ptbari(i,j,k)     &
                     *(ptforce(i,j,k)-tem9(i,j,k)*                      &
                       0.25*(ptbar(i,j,k+1)-ptbar(i,j,k-1))*dzinv*      &
                       (w(i,j,k+1,2)+w(i,j,k,2)) )
          END DO
        END DO
      END DO
    END IF

  END IF

  DO k=2,nz-2   ! preparing pforce...for all solvwp codes...
    DO j=1,ny-1
      DO i=1,nx-1
        pforce(i,j,k) = dtsml1*pforce(i,j,k)*j3inv(i,j,k)
      END DO
    END DO
  END DO


  DO ismstp =1, nsmstp
!
!-----------------------------------------------------------------------
!
!  Calculate the current time within small time step. It will be
!  used to calculate the external boundary fields.
!
!-----------------------------------------------------------------------
!
    curtsml = curtim + dtbig - 2*dtbig1 + ismstp*dtsml1
!
!-----------------------------------------------------------------------
!
!  Advance the three momentum equations one small timestep.
!  (Arguments such as u(1,1,1,tfuture) indicate the passing of
!  array values at time tfuture into the subroutine. Inside the
!  subroutine, these arrays are defined at a single time level.
!  u, v, w at tfuture and wcont are updated in time exiting
!  this routine.
!
!  Advance the pressure (p at tfuture) one small timestep.
!
!  tem7 is used to transfer wpgrad between subroutine SOLVUV and
!  subroutine SOLVWPIM or SOLVWPEX.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Calculate the boundary time tendencies of u and v using Klemp &
!  Wilhelmson type open boundary condition (rbcopt=1).
!
!-----------------------------------------------------------------------
!
    IF( rbcopt == 1 .AND. mgrid == 1 ) THEN

      CALL bkwsmldt(nx,ny,nz, u(1,1,1,tfuture),v(1,1,1,tfuture),        &
                    ubar,vbar,udteb,udtwb,vdtnb,vdtsb)

    END IF

    CALL solvuv(mptr, nx,ny,nz, exbcbufsz, dtsml1, curtsml,             &
                u(1,1,1,tfuture),v(1,1,1,tfuture),                      &
                wcont,pprt(1,1,1,tfuture),                              &
                udteb,udtwb,udtnb,udtsb,                                &
                vdteb,vdtwb,vdtnb,vdtsb,                                &
                rhostr, uforce,vforce,ubar,vbar,                        &
                x,y,z,zp,mapfct, j1,j2,j3,j3inv,                        &
                exbcbuf,rhofct,                                         &
                rhostru,rhostrv,rhostrw,tem13,tem14,tem7,               &
                tem1,tem2,tem3,tem4,tem5)
!
!-----------------------------------------------------------------------
!
!  Warning: tem7 is carrying wpgrad can not be overwritten before
!  the subroutine for w and p-equations are called.
!
!-----------------------------------------------------------------------
!

    IF( lvldbg >= 3 ) THEN
      CALL checkuhx(u(1,1,1,tfuture),                                   &
                    nx,ny,nz,1,nx,  1,ny-1,1,nz-1, 'uxsml', tem1)
      CALL checkuhy(u(1,1,1,tfuture),                                   &
                    nx,ny,nz,1,nx,  1,ny-1,1,nz-1, 'uysml', tem1)
      CALL checkvhx(v(1,1,1,tfuture),                                   &
                    nx,ny,nz,1,nx-1,1,ny,  1,nz-1, 'vxsml', tem1)
      CALL checkvhy(v(1,1,1,tfuture),                                   &
                    nx,ny,nz,1,nx-1,1,ny,  1,nz-1, 'vysml', tem1)
    END IF
!
!-----------------------------------------------------------------------
!
!  Advance the w and pressure (at tfuture) one small timestep.
!
!-----------------------------------------------------------------------
!
    IF( vimplct == 0 ) THEN   ! Not implicit
!
!-----------------------------------------------------------------------
!
!  Integrate the w and pressure equations one small time step
!  using explicit scheme.
!
!  NOTE: tem11,tem18,tem12,tem15,tem16,tem9 are assumed to be
!        constant during the small time step!!!
!
!-----------------------------------------------------------------------
!
      CALL solvwpex(mptr, nx,ny,nz, exbcbufsz, dtsml1, curtsml,         &
                    u(1,1,1,tfuture),v(1,1,1,tfuture),                  &
                    w(1,1,1,tfuture),wcont,                             &
                    ptprt(1,1,1,tfuture),pprt(1,1,1,tfuture),phydro,    &
                    wdteb,wdtwb,wdtnb,wdtsb,pdteb,pdtwb,pdtnb,pdtsb,    &
                    rhostr,ptbar,ptbari,pbari,csndsq,                   &
                    wforce,tem7,pforce,                                 &
                    x,y,z,zp, mapfct, j1,j2,j3,aj3x,aj3y,aj3z,j3inv,    &
                    rhostru,rhostrv,rhostrw,                            &
                    exbcbuf,rhofct,                                     &
                    tem1,tem2,tem3,tem4,tem5,                           &
                    tem11,tem18,tem19,tem15,tem16,tem9)
    ELSE   ! Implicit
!
!-----------------------------------------------------------------------
!
!  Integrate the w and pressure equations one small time step
!  using vertically implicit scheme.
!
!  NOTE: tem9,tem10,tem11,tem12,tem15,tem16,tem17,tem18
!        are assumed to be constant during the small time step!!!
!
!-----------------------------------------------------------------------
!
      IF(peqopt == 1) THEN

        CALL solvwpim(mptr,nx,ny,nz,exbcbufsz, dtsml1, curtsml,         &
                    u(1,1,1,tfuture),v(1,1,1,tfuture),                  &
                    w(1,1,1,tfuture),wcont,                             &
                    ptprt(1,1,1,tfuture),pprt(1,1,1,tfuture),phydro,    &
                    wdteb,wdtwb,wdtnb,wdtsb,pdteb,pdtwb,pdtnb,pdtsb,    &
                    rhostr,ptbar,ptbari,pbari,csndsq,                   &
                    wforce,tem7,pforce,                                 &
                    x,y,z,zp,mapfct,                                    &
                    j1,j2,j3,aj3x,aj3y,aj3z,j3inv,                      &
                    trigs1,trigs2,ifax1,ifax2,                          &
                    wsave1,wsave2,vwork1,vwork2,                        &
                    rhostru,rhostrv,rhostrw,                            &
                    exbcbuf,rhofct,                                     &
                    tem1,tem2,tem3,tem4,tem5,tem6,tem8,                 &
                    tem9,tem10,tem11,tem12,tem15,tem16,tem17)

      ELSE

        CALL solvwpim1(nx,ny,nz,exbcbufsz, dtsml1, curtsml,             &
                    u(1,1,1,tfuture),v(1,1,1,tfuture),                  &
                    w(1,1,1,tfuture),wcont,                             &
                    ptprt(1,1,1,tfuture),pprt(1,1,1,tfuture),phydro,    &
                    wdteb,wdtwb,wdtnb,wdtsb,pdteb,pdtwb,pdtnb,pdtsb,    &
                    rhostr,ptbar,ptbari,pbari,csndsq,                   &
                    wforce,tem7,pforce,                                 &
                    x,y,z,zp,mapfct,                                    &
                    j1,j2,j3,aj3z,j3inv,                                &
                    trigs1,trigs2,ifax1,ifax2,                          &
                    wsave1,wsave2,vwork1,vwork2,                        &
                    rhostru,rhostrv,rhostrw,                            &
                    exbcbuf,rhofct,                                     &
                    tem1,tem2,tem3,tem4,tem5,tem6,tem8,                 &
                    tem11,tem18,tem12,tem17,tem9)
      END IF

    END IF

    IF( lvldbg >= 3 ) THEN
      CALL checkwhx(w(1,1,1,tfuture),                                   &
                    nx,ny,nz,1,nx-1,1,ny-1,1,nz,   'wxsml', tem1)
      CALL checkwhy(w(1,1,1,tfuture),                                   &
                    nx,ny,nz,1,nx-1,1,ny-1,1,nz,   'wysml', tem1)
      CALL checkshx(pprt(1,1,1,tfuture),                                &
                    nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, 'pxsml', tem1)
      CALL checkshy(pprt(1,1,1,tfuture),                                &
                    nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, 'pysml', tem1)
    END IF

!-----------------------------------------------------------------------
!
!  Integrate the potential temperature equation one small time step
!  when ptsmlstp=1.
!
!-----------------------------------------------------------------------
!

    IF( ptsmlstp == 1 ) THEN

      CALL solvpt_sml(nx,ny,nz, exbcbufsz, dtbig1,dtsml1,curtsml,       &
                      ptprt,w, ptdteb,ptdtwb,ptdtnb,ptdtsb,             &
                      ptbar,rhostr,rhostri,rhostrw,j3,aj3z,             &
                      ptforce, exbcbuf,                                 &
                      tem1,tem2)

      IF( lvldbg >= 3) THEN
        CALL checkshx(ptprt(1,1,1,tfuture),                             &
                      nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, 'ptsmlx',tem1)
        CALL checkshy(ptprt(1,1,1,tfuture),                             &
                      nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, 'ptsmly',tem1)
      END IF

    END IF


  END DO

!
!-----------------------------------------------------------------------
!
!  Integrate the potential temperature equation one large time step
!  when ptsmlstp=0.
!
!-----------------------------------------------------------------------
!
  IF( ptsmlstp == 0 ) THEN

    CALL solvpt_lrg(nx,ny,nz, exbcbufsz, dtbig1, ptprt,                 &
                    ptdteb,ptdtwb,ptdtnb,ptdtsb,rhostr,rhostri,         &
                    ptforce,exbcbuf,tem1)

  END IF

  RETURN
END SUBROUTINE smlstep
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE TFILT                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE tfilt(nx,ny,nz,                                              &
           u,v,w,ptprt,pprt,qv,qscalar,tke)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Apply the Asselin time filter to all time-dependent variables.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91.
!
!  MODIFICATION HISTORY:
!
!  5/20/92 (K. Droegemeier and M. Xue)
!  Added full documentation.
!
!  9/10/94 (D. Weber & Y. Lu)
!  Cleaned up documentation.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    u        x component of velocity at all time levels (m/s)
!    v        y component of velocity at all time levels (m/s)
!    w        Vertical component of Cartesian velocity at all time
!             levels (m/s)
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
!  OUTPUT:
!
!    u        x component of velocity at time tpresent (m/s)
!    v        y component of velocity at time tpresent(m/s)
!    w        Vertical component of Cartesian velocity at tpresent (m/s)
!    ptprt    Perturbation potential temperature at time tpresent (K)
!    pprt     Perturbation pressure at time tpresent (Pascal)
!    qv       Water vapor specific humidity at time tpresent (kg/kg)
!    qc       Cloud water mixing ratio at time tpresent (kg/kg)
!    qr       Rainwater mixing ratio at time tpresent (kg/kg)
!    qi       Cloud ice mixing ratio at time tpresent (kg/kg)
!    qs       Snow mixing ratio at time tpresent (kg/kg)
!    qh       Hail mixing ratio at time tpresent (kg/kg)
!    tke      Turbulent Kinetic Energy ((m/s)**2)
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
  INCLUDE 'timelvls.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)
  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)

  REAL :: qscalar(nx,ny,nz,nt,nscalar)

  REAL :: tke   (nx,ny,nz,nt)  ! Turbulent kinetic energy

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: nq
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  CALL aselin(u, nx,ny,nz, 1,nx,1,ny-1,1,nz-1, flteps)

  CALL aselin(v, nx,ny,nz, 1,nx-1,1,ny,1,nz-1, flteps)

  CALL aselin(w, nx,ny,nz, 1,nx-1,1,ny-1,1,nz, flteps)

  IF(sadvopt /= 4) CALL aselin(ptprt, nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, flteps)

  CALL aselin(pprt, nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, flteps)

  IF( tmixopt == 4 .OR. tmixopt == 5) THEN
    IF(sadvopt /= 4) CALL aselin(tke, nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, flteps)
  END IF

  IF( moist /= 0 .AND. sadvopt /= 4 ) THEN

    CALL aselin(qv, nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, flteps)

    DO nq = 1, nscalar
      CALL aselin(qscalar(:,:,:,:,nq),nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, flteps)
    END DO

  ELSE IF (nscalar > 0) THEN   ! Added cc for no moisture process

    DO nq = 1, nscalar
      CALL aselin(qscalar(:,:,:,:,nq),nx,ny,nz,1,nx-1,1,ny-1,1,nz-1, flteps)
    END DO

  END IF

  RETURN
END SUBROUTINE tfilt

!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ASELIN                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE aselin(a, nx,ny,nz, nx1,nx2,ny1,ny2,nz1,nz2, flteps)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Apply Asselin time filter.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91.
!
!  MODIFICATION HISTORY:
!
!  5/20/92 (K. Droegemeier and M. Xue)
!  Added full documentation.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    a        A variable at three time levels.
!
!    nx1,nx2  i-index array bounds where the time filter is applied
!    ny1,ny2  j-index array bounds where the time filter is applied
!    nz1,nz2  k-index array bounds where the time filter is applied
!
!    flteps   The asselin time filter coefficient.
!             Defined in globcst.inc
!
!  OUTPUT:
!
!    a        The new array values at time tpresent after the time
!             filter is applied.
!
!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'timelvls.inc'

  INTEGER :: nx, ny, nz   ! Number of grid points in 3 directions

  REAL :: a(nx,ny,nz,nt)

  INTEGER :: nx1, nx2
  INTEGER :: ny1, ny2
  INTEGER :: nz1, nz2
  REAL    :: flteps       ! Coefficient of the Asselin time filter
!
!-----------------------------------------------------------------------
!
!  Misc. local variable:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=nz1,nz2
    DO j=ny1,ny2
      DO i=nx1,nx2
        a(i,j,k,tpresent)=a(i,j,k,tpresent)+flteps*                     &
            (a(i,j,k,tfuture)-2*a(i,j,k,tpresent)+a(i,j,k,tpast))
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE aselin

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE TFLIP                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE tflip(nx,ny,nz, u,v,w,ptprt,pprt,qv,qscalar,tke)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Update the fields u, v, w, ptprt ,pprt ,qv ,qc ,qr ,qi ,qs ,qh and
!  tke in time. The fields at tpast are replaced by those at tpresent.
!  The fields at tpresent are replaced by those at tfuture. The fields
!  at tfuture are not changed.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91.
!
!  MODIFICATION HISTORY:
!
!  5/20/92 (K. Droegemeier and M. Xue)
!  Added full documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    u        x component of velocity at all time levels (m/s)
!    v        y component of velocity at all time levels (m/s)
!    w        Vertical component of Cartesian velocity at all time
!             levels (m/s)
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
!  OUTPUT:
!
!    u        x-velocity at tpast and tpresent updated in time (m/s)
!    v        y-velocity at tpast and tpresent updated in time (m/s)
!    w        Vertical component of Cartesian velocity at tpast
!             and tpresent updated in time (m/s)
!    ptprt    Perturbation potential temperature at tpast and tpresent
!             updated in time (K)
!    pprt     Perturbation pressure at tpast and tpresent
!             updated in time (Pascal)
!    qv       Water vapor specific humidity at tpast and tpresent
!             updated in time (kg/kg)
!    qc       Cloud water mixing ratio at tpast and tpresent
!             updated in time (kg/kg)
!    qr       Rainwater mixing ratio at tpast and tpresent
!             updated in time (kg/kg)
!    qi       Cloud ice mixing ratio at tpast and tpresent
!             updated in time (kg/kg)
!    qs       Snow mixing ratio at tpast and tpresent
!             updated in time (kg/kg)
!    qh       Hail mixing ratio at tpast and tpresent
!             updated in time (kg/kg)
!    tke      Turbulent Kinetic Energy ((m/s)**2)
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
  INCLUDE 'timelvls.inc'
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  REAL :: u     (nx,ny,nz,nt)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz,nt)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz,nt)  ! Total w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz,nt)  ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz,nt)  ! Perturbation pressure (Pascal)

  REAL :: qv    (nx,ny,nz,nt)  ! Water vapor specific humidity (kg/kg)

  REAL :: qscalar(nx,ny,nz,nt,nscalar)

  REAL :: tke   (nx,ny,nz,nt)  ! Turbulent kinetic energy

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: nq
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  CALL tswap(nx,ny,nz, u)
  CALL tswap(nx,ny,nz, v)
  CALL tswap(nx,ny,nz, w)
  CALL tswap(nx,ny,nz, ptprt)
  CALL tswap(nx,ny,nz, pprt)

  IF( tmixopt == 4 .OR. tmixopt == 5) THEN
    CALL tswap(nx,ny,nz,tke )
  END IF

  IF( moist /= 0 ) THEN
    CALL tswap(nx,ny,nz, qv)
  END IF

  DO nq = 1,nscalar
    CALL tswap(nx,ny,nz,qscalar(:,:,:,:,nq))
  END DO

  RETURN
END SUBROUTINE tflip

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE TSWAP                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE tswap(nx,ny,nz, a)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Update the time levels of a time-dependent array.
!
!  a(*,*,*,tpast   ) = a(*,*,*,tpresent)
!  a(*,*,*,tpresent) = a(*,*,*,tfuture )
!  a(*,*,*,tfuture ) is not changed.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/10/91.
!
!  MODIFICATION HISTORY:
!
!  5/20/92 (K. Droegemeier and M. Xue)
!  Added full documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    a        A 3-D array at three time levels that will be updated
!             in time.
!
!  OUPTUT:
!
!    a        Its new values at time tpast and tfuture.
!
!-----------------------------------------------------------------------

!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'timelvls.inc'

  INTEGER :: nx, ny, nz   ! Number of grid points in 3 directions

  REAL :: a(nx,ny,nz,nt)  ! Array whose values will be updated in time
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, k
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
  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        a(i,j,k,tpast   )=a(i,j,k,tpresent)
        a(i,j,k,tpresent)=a(i,j,k,tfuture )
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE tswap


!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE READJMESO                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE readjmeso(nx,ny,nz,ayear,amm,aday,arpstime,  &
           u,v,w,ptbar,pbar,ptprt,pprt,qv,prcrate)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in the mesonet data qc'd and processed by Jerry Brotzge.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Weber and Jerry Brotzge
!  8/20/2001.
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
!
!  OUTPUT:
!
!    u        Updated (by dtbig) x velocity component at times tpast
!             and tpresent (m/s).
!    v        Updated (by dtbig) y velocity component at times tpast
!             and tpresent (m/s).
!    w        Updated (by dtbig) z velocity component at times tpast
!             and tpresent (m/s).
!    ptprt    Updated (by dtbig) perturbation potential temperature
!             at times tpast and tpresent (K).
!    pprt     Updated (by dtbig) perturbation pressure at times tpast
!             and tpresent (Pascal).
!    qv       Updated (by dtbig) water vapor specific humidity at times
!             tpast and tpresent (kg/kg).
!    qc       Updated (by dtbig) cloud water mixing ratio at times
!             tpast and tpresent (kg/kg).

!    raing    Grid supersaturation rain
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!-----------------------------------------------------------------------
!

  IMPLICIT NONE             ! Force explicit declarations

!
!  INPUT
!

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  INTEGER, PARAMETER :: nmeso=288

  REAL :: uair  (nmeso)     ! Mesonet u-velocity (m/s)
  REAL :: vair  (nmeso)     ! Mesonet v-velocity (m/s)
  REAL :: wair  (nmeso)     ! Mesonet w-velocity (m/s)
  REAL :: pres  (nmeso)     ! Mesonet pressure (Pa)
  REAL :: tair  (nmeso)     ! Mesonet potential temperature (K)
  REAL :: rain  (nmeso)     ! Mesonet rainfall rate (kg/m^2/s)
  REAL :: relh  (nmeso)     ! Mesonet mixing ratio (kg/kg)
  REAL :: ws2m  (nmeso)     ! Mesonet 2 m wind speed (m/s)
  CHARACTER(LEN=4) :: temp1 (nmeso) ! Mesonet station name
  INTEGER :: stnm1 (nmeso)  ! Mesonet station number
  INTEGER :: mesotime(nmeso)  ! Mesonet time (minutes then seconds)
  REAL :: mesot(nmeso)   ! Mesonet time (seconds)

!
!  OUTPUT
!

  REAL :: u     (nx,ny,nz)  ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz)  ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz)  ! Total w-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)  ! Base state potential temperature (K)
  REAL :: pbar (nx,ny,nz)   ! Base state pressure (Pa)
  REAL :: ptprt (nx,ny,nz)  ! Perturbation potential temperature
                            ! from that of base state atmosphere (K)
  REAL :: pprt  (nx,ny,nz)  ! Perturbation pressure from that
                            ! of base state atmosphere (Pascal)
  REAL :: qv    (nx,ny,nz)  ! Water vapor specific humidity (kg/kg)
!  REAL :: raing (nx,ny)     ! Gridscale rainfall (mm)
  REAL :: prcrate(nx,ny,4)  ! Precipitation rate

!
!-----------------------------------------------------------------------
!
!  Miscellaneous local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,n
  CHARACTER :: ayear*4, amm*2, aday*2

  CHARACTER :: atemp*4, btemp*2, ctemp*2, dtemp*3
  CHARACTER :: astid*4, astnm*4, aime*4, atair*4
  CHARACTER :: arelh*4, aws2m*4, arain*4, apres*4
  CHARACTER :: auair*4, avair*4, awair*4

  INTEGER :: cmm, cdd, filelength

  REAL :: tema, temb
  REAL :: arpstime
  REAL :: temarpstime
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'     ! Global constants that control model execution
!-----------------------------------------------------------------------


!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  CALL initztime(ayear,amm,aday)

  filelength = len_trim(sitemeso)

  OPEN(unit=60,file=sitemeso(1:filelength)//     &
      ayear//amm//aday//'.mts',status='old',form='formatted')

!    loop that reads Jerry's mesonet data file (qc'd)
!    note that the ARPS time must match the mesonet datafile time.

  read(60,411) dtemp
  read(60,413) btemp, atemp, cmm, cdd, ctemp, ctemp, ctemp
  read(60,401) astid, astnm, aime, atair, arelh, arain,       &
               aws2m, apres, auair, avair, awair

  401   FORMAT (a4,7(1x,a4),3(1x,a1))
  402   FORMAT (1x,a4,4x,i3,5x,i4,2(1x,f9.4),1x,f10.6,1x,f9.4,       &
            1x,f11.4,3(1x,f9.4))
  411   FORMAT (6x,a3)
  413   FORMAT (7x,a2,1x,a4,1x,a2,1x,a2,1x,a2,1x,a2,1x,a2)


  DO n=1,288   ! perform the file read....

    READ(60,402) temp1(n),stnm1(n),mesotime(n),                         &
                 tair(n), relh(n), rain(n),                                 &
                 ws2m(n), pres(n),uair(n),vair(n),wair(n)


    relh(n) = relh(n) / (relh(n) + 1.0)   !Converts w to q

!  test to see if we have a time match...
    mesotime(n) = mesotime(n) * 60
    mesot(n) = real(mesotime(n))

    temarpstime = arpstime + 60.0

    IF (temarpstime.le.mesot(n)) THEN
      IF(temarpstime.eq.mesot(n)) THEN   ! we have an exact match

        DO j = 1,ny
          DO i = 1,nx
            ptprt(i,j,2) = tair (n) - ptbar(i,j,2)
            pprt(i,j,2) = pres(n) - pbar(i,j,2)
            qv(i,j,2) = relh(n)
            prcrate(i,j,1) = rain(n)
            u(i,j,2) = uair(n)
            v(i,j,2) = vair(n)
            w(i,j,2) = wair(n)

            ptprt(i,j,1) = tair (n) - ptbar(i,j,2)
            pprt(i,j,1) = pres(n) - pbar(i,j,2)
            qv(i,j,1) = relh(n)
            prcrate(i,j,1) = rain(n)
            u(i,j,1) = uair(n)
            v(i,j,1) = vair(n)
            w(i,j,1) = wair(n)

          END DO
        END DO

      ELSE IF(temarpstime.gt.mesot(n-1).and.temarpstime.lt.mesot(n))THEN


!    load the mesonet data into the arps arrays.
!    note we need to convert the rh to qv.....using ARPS routines!!!!!

        tema = (arpstime-mesot(n-1))/300.0
        temb = 1.0-tema

        DO j = 1,ny
          DO i = 1,nx
            ptprt(i,j,2) = temb*tair(n-1)+tema*tair(n) - ptbar(i,j,2)
            pprt(i,j,2) = temb*pres(n-1)+tema*pres(n) - pbar(i,j,2)
            qv(i,j,2) = temb*relh(n-1)+tema*relh(n)
            prcrate(i,j,1) =  temb*rain(n-1)+tema*rain(n)
            u(i,j,2) = temb*uair(n-1) + tema*uair(n)
            v(i,j,2) = temb*vair(n-1) + tema*vair(n)
            w(i,j,2) = temb*wair(n-1) + tema*wair(n)

            ptprt(i,j,1) = temb*tair(n-1)+tema*tair(n) - ptbar(i,j,2)
            pprt(i,j,1) = temb*pres(n-1)+tema*pres(n) - pbar(i,j,2)
            qv(i,j,1) = temb*relh(n-1)+tema*relh(n)
            u(i,j,1) = temb*uair(n-1) + tema*uair(n)
            v(i,j,1) = temb*vair(n-1) + tema*vair(n)
            w(i,j,1) = temb*wair(n-1) + tema*wair(n)


          END DO
        END DO

      END IF   !  end of the time if loop
    END IF      ! end of the time if loop

  END DO                                   !Time loop

  CLOSE (60)

  RETURN
END SUBROUTINE readjmeso

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE READJFLUX                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE readjflux(nx,ny,nz,ayear,amm,aday,arpstime,  &
           radfrc,ptsflx,qvsflx)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in the mesonet data qc'd and processed by Jerry Brotzge.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Weber and Jerry Brotzge
!  8/20/2001.
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
!
!  OUTPUT:
!
!    radfrc   Radiation forcing (K/s)
!
!    ptsflx   Surface heat flux (K*kg/ (m**2*s))
!
!    qvsflx   Surface moisture flux (kg/ (m**2*s))
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!-----------------------------------------------------------------------
!

  IMPLICIT NONE             ! Force explicit declarations

!
!  INPUT
!

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  INTEGER, PARAMETER :: nmeso=288

  REAL :: rnet  (nmeso)     ! Mesonet CNR1 rnet (W/m**2)
  REAL :: hflx  (nmeso)     ! Mesonet sonic sensible heat flux (W/m**2)
  REAL :: lflx  (nmeso)     ! Mesonet sonic latent heat flux (W/m**2)
  REAL :: gflx  (nmeso)     ! Mesonet ground heat flux (W/m**2)
  REAL :: l2fx  (nmeso)     ! Mesonet residual latent heat flux (W/m**2)
  REAL :: irtx  (nmeso)     ! Mesonet rainfall rate (kg/m^2/s)
  REAL :: prof  (nmeso)     ! Mesonet profile sensible flux (W/m**2)
  CHARACTER(LEN=4) :: temp1 (nmeso) ! Mesonet station name
  INTEGER :: stnm1 (nmeso)  ! Mesonet station number
  INTEGER :: mesotime(nmeso)   ! Mesonet time (minutes then seconds)
  REAL :: mesot(nmeso)      ! Mesonet time (seconds)

!
!  OUTPUT
!

  REAL :: radfrc (nx,ny,nz)  ! Total u-velocity (K/s)
  REAL :: ptsflx (nx,ny)     ! Total v-velocity (K*kg/(m**2*s))
  REAL :: qvsflx (nx,ny)     ! Total w-velocity (kg/(m**2*s))

!
!-----------------------------------------------------------------------
!
!  Miscellaneous local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,n
  CHARACTER :: ayear*4, amm*2, aday*2
  CHARACTER :: atemp*4, btemp*2, ctemp*2, dtemp*3
  CHARACTER :: astid*4, astnm*4, aime*4, airtx*4
  CHARACTER :: arnet*4, ahflx*4, alflx*4, al2fx*4
  CHARACTER :: agflx*4

  INTEGER :: cmm, cdd, filelength

  REAL :: tema, temb
  REAL :: arpstime

!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'     ! Global constants that control model execution
!-----------------------------------------------------------------------



!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!


  filelength = len_trim(siteflux)

  OPEN(unit=70,file=siteflux(1:filelength)//          &
         ayear//amm//aday//'.mts',status='old',       &
         form='formatted')

!    loop that reads Jerry's mesonet data file (qc'd)
!    note that the ARPS time must match the mesonet datafile time.


  read(70,511) dtemp
  read(70,513) btemp, atemp, cmm, cdd, ctemp, ctemp, ctemp
  read(70,501) astid, astnm, aime, arnet, ahflx, alflx,       &
               agflx, al2fx, airtx

  DO n=1,288   ! perform the file read....

    read(70,502) temp1(n),stnm1(n),mesotime(n),        &
                 rnet(n), hflx(n), lflx(n),            &
                 gflx(n), l2fx(n),irtx(n), prof(n)

    hflx(n) = hflx(n) / 1004.0     !Converts W/m**2 to (K*kg/(s*m**2))
    lflx(n) = lflx(n) / 2250000.0  !Converts W/m**2 to (kg/(s*m**2))
    l2fx(n) = l2fx(n) / 2250000.0  !Converts W/m**2 to (kg/(s*m**2))

    501   FORMAT (a4,8(1x,a4))
    511   FORMAT (6x,a3)
    513   FORMAT (7x,a2,1x,a4,1x,a2,1x,a2,1x,a2,1x,a2,1x,a2)
    502   format (1x,a4,4x,i3,5x,i4,7(1x,f9.4))

!  test to see if we have a time match...
    mesotime(n) = mesotime(n) * 60
    mesot(n) = real(mesotime(n))

    IF (arpstime.le.mesot(n)) THEN

      IF(arpstime.eq.mesot(n)) THEN   ! we have an exact match

        DO j = 1,ny
          DO i = 1,nx
            ptsflx(i,j) = hflx(n)
            qvsflx(i,j) = l2fx(n)
            radfrc(i,j,2) = rnet(n)
            radfrc(i,j,1) = rnet(n)
          END DO
        END DO

      ELSE IF(arpstime.gt.mesot(n-1).and.arpstime.le.mesot(n))THEN

!    load the mesonet data into the arps arrays.
!    note we need to convert the rh to qv.....using ARPS routines!!!!!

        tema = (arpstime-mesot(n-1))/300.0
        temb = 1.0-tema

        DO j = 1,ny
          DO i = 1,nx
             ptsflx(i,j) = temb*hflx(n-1)+tema*hflx(n)
             qvsflx(i,j) = temb*l2fx(n-1)+tema*l2fx(n)
             radfrc(i,j,2) = temb*rnet(n-1)+tema*rnet(n)
             radfrc(i,j,1) = temb*rnet(n-1)+tema*rnet(n)
           END DO
        END DO

      END IF   !  end of the time if loop
    END IF       ! end of the time if loop

  END DO                                   !Time loop

  CLOSE (70)

  RETURN
END SUBROUTINE readjflux

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE READJRADD                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE readjradd(nx,ny,nz,ayear,amm,aday,arpstime,  &
           radsw,radswnet,radlwin,rnflx)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in the mesonet radt'n data qc'd and processed by Jerry Brotzge.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Weber and Jerry Brotzge
!  8/28/2001.
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
!
!  OUTPUT:
!
!    radsw   Surface solar radiation flux (W/m**2)
!    rnflx   Surface net radiation flux (W/m**2)
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!-----------------------------------------------------------------------
!

  IMPLICIT NONE             ! Force explicit declarations

!
!  INPUT
!

  INTEGER :: nx, ny, nz        ! Number of grid points in 3 directions

  INTEGER, PARAMETER :: nmeso=288

  REAL    :: rnet (nmeso)      !Mesonet CNR1 rnet (W/m**2)
  REAL    :: swou (nmeso)      !Mesonet sw outgoing (W/m**2)
  REAL    :: swin (nmeso)      !Mesonet sw incoming (W/m**2)
  REAL    :: lwin (nmeso)      !Mesonet lw incoming (W/m**2)
  REAL    :: lwou (nmeso)      !Mesonet lw outgoing (W/m**2)

  CHARACTER*4 :: temp1 (nmeso) ! Mesonet station name
  INTEGER     :: stnm1 (nmeso)  ! Mesonet station number
  INTEGER     :: mesotime(nmeso) !Mesonet time (minutes then seconds)
  REAL        :: mesot(nmeso)   ! Mesonet time (seconds)

  !
  !  OUTPUT
  !

  REAL :: radsw(nx,ny)   !Surface solar radiation flux (W/m**2)
  REAL :: radswou(nx,ny) !Reflected solar radiation flux (W/m**2)
  REAL :: radlwin(nx,ny)  !Incoming longwave radiation flux (W/m**2)
  REAL :: rnflx(nx,ny)   !Surface net radiation flux (W/m**2)
  REAL :: radswnet(nx,ny) !Net solar radiation flux (W/m**2)

!
!-----------------------------------------------------------------------
!
!  Miscellaneous local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,n
  CHARACTER :: ayear*4, amm*2, aday*2

  CHARACTER :: atemp*4, btemp*2, ctemp*2, dtemp*3
  CHARACTER :: astid*4, astnm*4, aime*4
  CHARACTER :: arnet*4, aswin*4, aswou*4, alwin*4
  CHARACTER :: alwou*4

  INTEGER :: cmm, cdd, filelength

  REAL :: tema, temb
  REAL :: arpstime

!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  CALL initztime(ayear,amm,aday)

  filelength = len_trim(siternet)

  OPEN(unit=20,file=siternet(1:filelength)//                            &
      ayear//amm//aday//'.mts',status='old',form='formatted')

  !    loop that reads Jerry's mesonet data file (qc'd)
  !    note that the ARPS time must match the mesonet datafile time.

  read(20,311) dtemp
  read(20,313) btemp, atemp, cmm, cdd, ctemp, ctemp, ctemp
  read(20,301) astid, astnm, aime, aswin, aswou, alwin,                 &
               alwou, arnet

  301   FORMAT (a4,7(1x,a4))
  302   FORMAT (1x,a4,4x,i3,5x,i4,5(1x,f9.4))
  311   FORMAT (6x,a3)
  313   FORMAT (7x,a2,1x,a4,1x,a2,1x,a2,1x,a2,1x,a2,1x,a2)

  DO n=1,288   ! perform the file read....

    READ(20,302) temp1(n),stnm1(n),mesotime(n),        &
                 swin(n), swou(n), lwin(n),            &
                 lwou(n), rnet(n)

!  test to see if we have a time match...
    mesotime(n) = mesotime(n) * 60
    mesot(n) = real(mesotime(n))

    IF (arpstime.le.mesot(n)) THEN
       IF(arpstime.eq.mesot(n)) THEN   ! we have an exact match

         DO j = 1,ny
           DO i = 1,nx
             radsw(i,j) = swin(n)
             radswou(i,j) = swou(n)
             radlwin(i,j) = lwin(n)
             rnflx(i,j) = rnet(n)
             radswnet(i,j) = radsw(i,j) - radswou(i,j)
           END DO
         END DO

       ELSE IF(arpstime.gt.mesot(n-1).and.arpstime.lt.mesot(n))THEN

!    load the mesonet data into the arps arrays.

         tema = (arpstime-mesot(n-1))/300.0
         temb = 1.0-tema

         DO j = 1,ny
           DO i = 1,nx
             radsw(i,j) = temb*swin(n-1)+tema*swin(n)
             radswou(i,j) = temb*swou(n-1)+tema*swou(n)
             radlwin(i,j) = temb*lwin(n-1)+tema*lwin(n)
             rnflx(i,j) = temb*rnet(n-1)+tema*rnet(n)
             radswnet(i,j) = radsw(i,j) - radswou(i,j)
           END DO
         END DO

       END IF   !  end of the time if loop
    END IF    !  end of the time if loop for matching arpstime and mesot

  END DO    !Time loop

  CLOSE (20)

  RETURN
END SUBROUTINE readjradd
