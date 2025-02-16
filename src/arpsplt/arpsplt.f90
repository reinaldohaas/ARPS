PROGRAM arpsplt
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                PROGRAM ARPSPLT                       ######
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
!  This is a graphic analysis/plotting program for display ARPS
!  history format data set.
!
!  It is based on Ming Xue's graphics package ZXPLOT, which is an
!  independent package from the ARPS code. The object code library
!  of ZXPLOT for most platforms are freely available from
!  ftp://ftp.caps.ou.edu/pub/ZXPLOT3. Documentation and other info
!  on ZXPLOT can be found at http://www.caps.ou.edu/ZXPLOT.
!
!  ZXPLOT can be interfaced with NCAR graphics to produce metafile
!  output consistent with NCAR graphics viewing facilities
!  or linked to the Postscript driver (pure fortran program) to
!  produce Postscript output directly (i.e., no NCAR graphics needed).
!
!  Documentation on the control parameters for plotting can be found
!  in input/arpsplt.input.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue, CAPS/OU.
!    2/19/1992.
!
!  MODIFICATION HISTORY:
!    See HISTORY file.
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
!    qc       Cloud water mixing ratio (kg/kg)
!    qr       Rainwater mixing ratio (kg/kg)
!    qi       Cloud ice mixing ratio (kg/kg)
!    qs       Snow mixing ratio (kg/kg)
!    qh       Hail mixing ratio (kg/kg)
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
!    psl      Sea level pressure (mb)
!
!  CALCULATED DATA ARRAYS:
!
!    u        x-component of velocity (m/s)
!    v        y-component of velocity (m/s)
!    w        z-component of velocity (m/s)
!    pt       potential temperature (K)
!    qv       water vapor mixing ratio (kg/kg)
!    td       dew-point temperature (C)
!    cape     CAPE  (J/kg)
!    cin      CIN   (J/kg)
!    thet     theta_E (K)
!    heli     helicity (m2/s2)
!    uh       updraft helicity (m2/s2)
!    srlfl    storm-relative low-level flow (0-2km AGL)
!    srmfl    storm-relative mid-level flow (2-9km AGL)
!    shr37    7km - 3km wind shear
!    ustrm    Estimated storm motion (Bob Johns)
!    vstrm    Estimated storm motion (Bob Johns)
!    capst    CAPE strength
!    blcon    boundary layer convergence
!    ct       convective temperature
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
!  Arrays for plots on constant pressure levels
!
!-----------------------------------------------------------------------
!
!   tz       Temperature (K) on computational grids
!   t700     Temperature (K) on 700mb pressure grids
!   zps3d    negative log pressure(Pascal) at ARPS grid points
!   algpzc   -log(pressure) at scalar grid points
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  USE module_arbitrary_vario

  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Grid dimensions.
  INTEGER :: nzprofc ! the maximum vertical index in height (zpc/zpsoilc)
                 ! and variables to be profiled when calling vprofil
                 ! subroutine. 06/10/2002, Zuwen He
                 ! In the atmosphere model, the vertical index is
                 ! typically nz-1, while in the soil model, it's nzsoil.

  INTEGER :: nzsoil            ! levels of soil model
  INTEGER :: nstyps            ! Maximum number of soil types.

  INTEGER :: hinfmt
  INTEGER :: nhisfile_max,nhisfile
  PARAMETER (nhisfile_max=200)
  CHARACTER (LEN=256) :: grdbasfn
  CHARACTER (LEN=256) :: hisfile(nhisfile_max)
  CHARACTER (LEN=256) :: hdmpftrailer

  COMMON /init2_hisf/ hinfmt,nhisfile, grdbasfn, hisfile, hdmpftrailer
  INTEGER :: first_frame
  COMMON /frstfrm/ first_frame

  INTEGER :: lengbf,nf,lenfil

  INTEGER, PARAMETER :: max_dim=200
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
  INCLUDE 'arpsplt.inc'
  INCLUDE 'alloc.inc'

  INCLUDE 'mp.inc'
  INCLUDE 'arpstrajc.inc'
!
!-----------------------------------------------------------------------
!
!  Arrays to be read in:
!
!-----------------------------------------------------------------------
!
  REAL, DIMENSION(:), POINTER :: x   ! The x-coord. of the physical and
                                     ! computational grid. Defined at u-point.
  REAL, DIMENSION(:), POINTER :: y   ! The y-coord. of the physical and
                                     ! computational grid. Defined at v-point.
  REAL, DIMENSION(:), POINTER :: z   ! The z-coord. of the computational grid.
                                     ! Defined at w-point on the staggered grid.
  REAL, DIMENSION(:,:,:), POINTER :: zp     ! The physical height coordinate defined at
                                            ! w-point of the staggered grid.
  REAL, DIMENSION(:,:,:), POINTER :: zpsoil ! The physical height coordinate defined at
                                            ! w-point of the staggered grid for soil model.

  REAL, DIMENSION(:,:,:), POINTER :: uprt   ! Perturbation u-velocity (m/s)
  REAL, DIMENSION(:,:,:), POINTER :: vprt   ! Perturbation v-velocity (m/s)
  REAL, DIMENSION(:,:,:), POINTER :: wprt   ! Perturbation w-velocity (m/s)
  REAL, DIMENSION(:,:,:), POINTER :: ptprt  ! Perturbation potential temperature
                                            ! from that of base state atmosphere (K)
  REAL, DIMENSION(:,:,:), POINTER :: pprt   ! Perturbation pressure from that
                                            ! of base state atmosphere (Pascal)
  REAL, DIMENSION(:,:,:), POINTER :: qvprt
  REAL, DIMENSION(:,:,:), POINTER :: qv     ! Water vapor specific humidity (kg/kg)
  REAL, DIMENSION(:,:,:,:), POINTER :: qscalar
  REAL, DIMENSION(:,:,:), POINTER :: tke    ! Turbulent Kinetic Energy ((m/s)**2)
  REAL, DIMENSION(:,:,:), POINTER :: kmh    ! Horizontal turb. mixing coef. for
                                            ! momentum. ( m**2/s )
  REAL, DIMENSION(:,:,:), POINTER :: kmv    ! Vertical turb. mixing coef. for
                                            ! momentum. ( m**2/s )
  REAL, DIMENSION(:,:,:), POINTER :: ubar   ! Base state u-velocity (m/s)
  REAL, DIMENSION(:,:,:), POINTER :: vbar   ! Base state v-velocity (m/s)
  REAL, DIMENSION(:,:,:), POINTER :: wbar   ! Base state w-velocity (m/s)
  REAL, DIMENSION(:,:,:), POINTER :: ptbar  ! Base state potential temperature (K)
  REAL, DIMENSION(:,:,:), POINTER :: pbar   ! Base state pressure (Pascal)
  REAL, DIMENSION(:,:,:), POINTER :: rhobar ! Base state density rhobar
  REAL, DIMENSION(:,:,:), POINTER :: qvbar  ! Base state water vapor specific humidity
                                            ! (kg/kg)

  INTEGER, DIMENSION(:,:,:), POINTER :: soiltyp ! Soil type
  INTEGER, DIMENSION(:,:),   POINTER :: vegtyp  ! Vegetation type
  REAL, DIMENSION(:,:,:),  POINTER :: stypfrct  ! Soil type fraction
  REAL, DIMENSION(:,:),    POINTER :: lai       ! Leaf Area Index
  REAL, DIMENSION(:,:),    POINTER :: roufns   ! Surface roughness
  REAL, DIMENSION(:,:),    POINTER :: veg      ! Vegetation fraction

  REAL, DIMENSION(:,:,:,:), POINTER :: tsoil   ! soil temperature (K)
  REAL, DIMENSION(:,:,:,:), POINTER :: qsoil   ! soil moisture
  REAL, DIMENSION(:,:,:), POINTER :: wetcanp   ! Canopy water amount
  REAL, DIMENSION(:,:), POINTER :: snowdpth  ! Snow depth (m)
  REAL, DIMENSION(:,:), POINTER :: raing     ! Grid supersaturation rain
  REAL, DIMENSION(:,:), POINTER :: rainc     ! Cumulus convective rain
  REAL, DIMENSION(:,:), POINTER :: raint     ! Total rain (rainc+raing)
  REAL, DIMENSION(:,:,:), POINTER :: prcrate ! precipitation rate (kg/(m**2*s))
                                             ! prcrate(1,1,1) = total precip. rate
                                             ! prcrate(1,1,2) = grid scale precip. rate
                                             ! prcrate(1,1,3) = cumulus precip. rate
                                             ! prcrate(1,1,4) = microphysics precip. rate

  REAL, DIMENSION(:,:,:), POINTER :: radfrc  ! Radiation forcing (K/s)
  REAL, DIMENSION(:,:), POINTER :: radsw     ! Solar radiation reaching the surface
  REAL, DIMENSION(:,:), POINTER :: rnflx     ! Net radiation flux absorbed by surface
  REAL, DIMENSION(:,:), POINTER :: radswnet  ! Net shortwave radiation
  REAL, DIMENSION(:,:), POINTER :: radlwin   ! Incoming longwave radiation


  REAL, DIMENSION(:,:), POINTER :: usflx     ! Surface flux of u-momentum (kg/(m*s**2))
  REAL, DIMENSION(:,:), POINTER :: vsflx     ! Surface flux of v-momentum (kg/(m*s**2))
  REAL, DIMENSION(:,:), POINTER :: ptsflx    ! Surface heat flux (K*kg/(m*s**2))
  REAL, DIMENSION(:,:), POINTER :: qvsflx    ! Surface moisture flux (kg/(m**2*s))
!
!-----------------------------------------------------------------------
!
!  Arrays derived from the read-in arrays
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: u     (:,:,:)        ! Total u-velocity (m/s)
  REAL, ALLOCATABLE :: v     (:,:,:)        ! Total v-velocity (m/s)
  REAL, ALLOCATABLE :: w     (:,:,:)        ! Total w-velocity (m/s)
  REAL, ALLOCATABLE :: pt    (:,:,:)  ! Total poten

  REAL, ALLOCATABLE :: xc   (:,:,:)   ! x-coor of scalar point (km)
  REAL, ALLOCATABLE :: yc   (:,:,:)   ! y-coor of scalar point (km)
  REAL, ALLOCATABLE :: zc   (:,:,:)   ! z-coor of scalar point in computational
                                      ! space (km)
  REAL, ALLOCATABLE :: zpc  (:,:,:)   ! z-coor of scalar point in physical
                                      ! space (km)
  REAL, ALLOCATABLE :: zpsoilc (:,:,:) ! zsoil-coor of scalar point in physical
                                       ! space (m)

  REAL, ALLOCATABLE :: acc_raing(:,:,:) ! Accumulated grid supersaturation rain
                                        ! between 2 successive time levels
  REAL, ALLOCATABLE :: acc_rainc(:,:,:) ! Accumulated cumulus convective rain
                                        ! between 2 successive time levels
  REAL, ALLOCATABLE :: acc_raint(:,:,:) ! Accumulated rain (rainc+raing)
                                        ! between 2 successive time levels

  REAL, ALLOCATABLE :: hterain(:,:)   ! The height of the terrain.

  REAL, ALLOCATABLE :: psl  (:,:)     ! Sea level pressure (mb)
  REAL, ALLOCATABLE :: td   (:,:,:)   ! dew-point temperature (C)
!
!-----------------------------------------------------------------------
!
!  Arrays for plots on constant pressure levels
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: tz  (:,:,:)    ! Temperature (K) on computational grids
  REAL, ALLOCATABLE :: t700  (:,:)    ! Temperature (K) on 700mb pressure grids

  REAL, ALLOCATABLE :: algpzc(:,:,:)  ! -log(pressure) at scalar grid points
  REAL, ALLOCATABLE :: zps3d(:,:,:)   !
  REAL, ALLOCATABLE :: zpsoils3d(:,:,:)   !

!-----------------------------------------------------------------------
!
!  Array for CAPE , CIN
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: li(:,:)     ! Lifted Index (K)
  REAL, ALLOCATABLE :: cape(:,:)   ! CAPE  (J/kg)
  REAL, ALLOCATABLE :: cin(:,:)    ! CIN   (J/kg)
  REAL, ALLOCATABLE :: thet(:,:)   ! theta_E (K)
  REAL, ALLOCATABLE :: heli(:,:)   ! helicity
  REAL, ALLOCATABLE :: uh(:,:)     ! updraft helicity
  REAL, ALLOCATABLE :: brn(:,:)    ! Bulk Richardson Number (Weisman and Klemp)
  REAL, ALLOCATABLE :: brnu(:,:)   ! Shear parameter of BRN, "U-squared"
  REAL, ALLOCATABLE :: srlfl(:,:)  ! storm-relative low-level flow (0-2km AGL)
  REAL, ALLOCATABLE :: srmfl(:,:)  ! storm-relative mid-level flow (2-9km AGL)
  REAL, ALLOCATABLE :: ustrm(:,:)  ! Estimated storm motion (Bob Johns)
  REAL, ALLOCATABLE :: vstrm(:,:)  ! Estimated storm motion (Bob Johns)

  REAL, ALLOCATABLE :: capst(:,:)  ! cap strength
  REAL, ALLOCATABLE :: blcon(:,:)  ! boundary llayer convergence
  REAL, ALLOCATABLE :: ct(:,:)     ! convective temperature

  REAL, ALLOCATABLE :: sinlat(:,:) ! Sin of latitude at each grid point

  REAL, ALLOCATABLE :: xs(:),ys(:)

  REAL :: dxkm, dykm, dzkm
  REAL :: dzsoilcm         ! Zuwen He, 05/31/2002 in cm.
  REAL :: alttostpr
!
!-----------------------------------------------------------------------
!
!  Temporary work arrays for general use
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: tem1(:,:,:)
  REAL, ALLOCATABLE :: tem2(:,:,:)
  REAL, ALLOCATABLE :: tem3(:,:,:)
  REAL, ALLOCATABLE :: tem4(:,:,:)
  REAL, ALLOCATABLE :: tem5(:,:,:)
  REAL, ALLOCATABLE :: tem6(:,:,:)
  REAL, ALLOCATABLE :: tem7(:,:,:)
  REAL, ALLOCATABLE :: tem8(:,:,:)
  REAL, ALLOCATABLE :: tem9(:,:,:)
!
!-----------------------------------------------------------------------
!
!  Work arrays used by profile plotting.
!
!-----------------------------------------------------------------------
!
  REAL :: xprof(max_dim)
  REAL :: yprof(max_dim)
!
!-----------------------------------------------------------------------
!
!  Work arrays to be used in interpolation subroutines
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: b1(:,:),b2(:,:)
  REAL, ALLOCATABLE :: u1(:,:),v1(:,:)
  REAL, ALLOCATABLE :: u2(:,:),v2(:,:),w2(:,:),zs2(:,:)
  REAL, ALLOCATABLE :: xp(:),yp(:)

!
!-----------------------------------------------------------------------
!
!  Some universal constants
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: kappa=287.053/cp,       &
                     gamma=6.5,              & ! 6.5 K/km
                     ex1=0.1903643,          & ! R*gamma/g
                     ex2=5.2558774,          & ! g/R/gamma
                     mbtopa=100.

  INTEGER, PARAMETER :: nscalarmax=30     ! max number of scalars for plotting
  REAL :: factor

!-----------------------------------------------------------------------
!
!  Plotting control parameters entered by user
!
!-----------------------------------------------------------------------
!
  INTEGER :: nchin  ! input file unit number
  INTEGER :: n, iorig

  INTEGER :: vtrplt, vtpplt, vagplt,vtrstrm,vtpstrm, xuvplt, strmplt
  INTEGER :: uplot, vplot, wplot, ptplot, pplot, hplot, tplot, vhplot,  &
             vsplot

  INTEGER :: qvplot,kmhplt,          &
             kmvplt,tkeplt,rhplot,rfplot,pteplt,tdplot,qwplot,          &
             rfcplt,qtplot, rhiplot, zdrplt, kdpplt, zdpplt, rhvplt

  INTEGER :: qscalarplot(nscalarmax)

  INTEGER :: dualpol
  CHARACTER (LEN=256) :: rsadir
  REAL :: wavelen
  INTEGER :: rfopt, wrfopt
  INTEGER :: upplot,vpplot,wpplot,ptpplt,ppplot,qvpplt,vorpplt,divpplt
  INTEGER :: divqplt,gricplt, avorplt
  INTEGER :: viqcplt,viqrplt,viqiplt,viqhplt,viqsplt,vilplt,viiplt,     &
             vicplt, vitplt, pwplt,tprplt, gprplt, cprplt

!
!-----------------------------------------------------------------------
!
!  Contour intervals
!
!-----------------------------------------------------------------------
!
  REAL :: vtrunit, vtpunit, vagunit,xuvunit, strmunit
  REAL :: uinc, vinc, winc, ptinc, pinc, hinc, tinc, vhinc,vsinc
  REAL :: qvinc,qwinc,qtinc
  REAL :: qscalarinc(nscalarmax)
  REAL :: kmhinc,kmvinc,tkeinc,rhinc,rfinc,pteinc,tdinc,rfcinc, rhiinc,&
          zdrinc,kdpinc,zdpinc,rhvinc
  REAL :: upinc,vpinc,wpinc,ptpinc,ppinc,qvpinc,vorpinc,divpinc,divqinc
  REAL :: gricinc,avorinc
  REAL :: viqcinc,viqrinc,viqiinc,viqhinc,viqsinc,vilinc,viiinc,vicinc
  REAL :: vitinc, pwinc,tprinc, gprinc, cprinc

!
!-----------------------------------------------------------------------
!
!  Limited variable minimum and maximum for color contour shade
!
!-----------------------------------------------------------------------
!
  REAL :: wcpminc, wcpmaxc, trnminc, trnmaxc
  REAL :: raincminc, raincmaxc, raingminc, raingmaxc,                   &
          raintminc,raintmaxc
  REAL :: rainicminc, rainicmaxc, rainigminc, rainigmaxc,               &
          rainitminc,rainitmaxc

  REAL :: capeminc, capemaxc, cinminc, cinmaxc, thetminc, thetmaxc,     &
          heliminc, helimaxc, uhminc, uhmaxc, capsminc, capsmaxc,       &
          blcominc, blcomaxc, ctcminc, ctcmaxc
  REAL :: uhmnhgt,uhmxhgt
  REAL :: brnminc, brnmaxc, bruminc, brumaxc, srlminc, srlmaxc,         &
          srmminc, srmmaxc, liminc, limaxc
  REAL :: uminc, umaxc, vminc, vmaxc, wminc, wmaxc,ptminc, ptmaxc
  REAL :: vhminc, vhmaxc,vsminc, vsmaxc
  REAL :: pminc, pmaxc, qvminc, qvmaxc
  REAL :: qscalarminc(nscalarmax),qscalarmaxc(nscalarmax)
  REAL :: qwminc, qwmaxc, qtminc, qtmaxc
  REAL :: rhminc, rhmaxc, rhiminc, rhimaxc
  REAL :: kmhminc, kmhmaxc,kmvminc, kmvmaxc,tkeminc, tkemaxc
  REAL :: rfminc, rfmaxc, pteminc, ptemaxc, upminc, upmaxc
  REAL :: rfcminc, rfcmaxc
  REAL :: zdrminc, kdpminc, zdpminc, rhvminc
  REAL :: zdrmaxc, kdpmaxc, zdpmaxc, rhvmaxc
  REAL :: vpminc, vpmaxc, wpminc, wpmaxc, ptpminc, ptpmaxc
  REAL :: ppminc, ppmaxc, qvpminc, qvpmaxc, vorpminc, vorpmaxc
  REAL :: divpminc, divpmaxc, divqminc, divqmaxc
  REAL :: gricminc, gricmaxc,avorminc, avormaxc
  REAL :: hminc, hmaxc, tminc, tmaxc,pslmaxc,pslminc
  REAL :: tdminc, tdmaxc
  REAL :: soiltpminc,soiltpmaxc,vegtpminc,vegtpmaxc,laiminc,laimaxc,    &
          rouminc,roumaxc,vegminc,vegmaxc,snowdminc,snowdmaxc
  REAL :: viqcminc,viqrminc,viqiminc,viqhminc,viqsminc,vilminc,viiminc, &
          vicminc, vitminc, pwminc, tprminc, gprminc, cprminc
  REAL :: viqcmaxc,viqrmaxc,viqimaxc,viqhmaxc,viqsmaxc,vilmaxc,viimaxc, &
          vicmaxc, vitmaxc, pwmaxc, tprmaxc, gprmaxc, cprmaxc

!-----------------------------------------------------------------------
!
!  Overlay control parameters
!
!-----------------------------------------------------------------------
!
  INTEGER :: vtrovr,vtpovr,vagovr,vtrstmovr,vtpstmovr,xuvovr,strmovr
  INTEGER :: uovr , vovr , wovr , ptovr , povr, hovr, tovr,vhovr,vsovr
  INTEGER :: qvovr ,kmhovr ,         &
             kmvovr,tkeovr,rhovr ,rfovr ,pteovr,tdovr, qwovr, qtovr,    &
             rfcovr, rhiovr, zdrovr, kdpovr, zdpovr, rhvovr
  INTEGER :: qscalarovr(nscalarmax)
  INTEGER :: upovr ,vpovr ,wpovr ,ptpovr,ppovr ,qvpovr,vorpovr,divpovr
  INTEGER :: trnovr,wcpovr,racovr,ragovr,ratovr,pslovr,                 &
             capovr,cinovr,theovr,helovr,uhovr,brnovr,brnuovr,          &
             srlfovr,srmfovr,liovr,capsovr,blcoovr,ctcovr
  INTEGER :: raicovr,raigovr,raitovr
  INTEGER :: divqovr,gricovr,avorovr
  INTEGER :: styovr,vtyovr,laiovr,rouovr,vegovr,snowdovr
  INTEGER :: viqcovr,viqrovr,viqiovr,viqhovr,viqsovr,vilovr,viiovr,     &
             vicovr, vitovr, pwovr, tprovr, gprovr, cprovr
!
!-----------------------------------------------------------------------
!
!  highlighting frequency for contour parameters
!
!-----------------------------------------------------------------------
!
  INTEGER :: uhlf , vhlf , whlf , pthlf , phlf, hhlf, thlf,vhhlf,vshlf
  INTEGER :: qvhlf ,kmhhlf ,         &
             kmvhlf,tkehlf,rhhlf ,rfhlf ,ptehlf,tdhlf, qwhlf, qthlf,    &
             rfchlf, rhihlf, zdrhlf, kdphlf, zdphlf, rhvhlf
  INTEGER :: qscalarhlf(nscalarmax)
  INTEGER :: uphlf ,vphlf ,wphlf ,ptphlf,pphlf ,qvphlf,vorphlf,divphlf
  INTEGER :: trnhlf,wcphlf,rachlf,raghlf,rathlf,pslhlf,                 &
             caphlf,cinhlf,thehlf,helhlf,uhhlf,brnhlf,brnuhlf,          &
             srlfhlf,srmfhlf,lihlf,capshlf,blcohlf,ctchlf
  INTEGER :: divqhlf,grichlf,avorhlf
  INTEGER :: styhlf,vtyhlf,laihlf,rouhlf,veghlf,snowdhlf
  INTEGER :: viqchlf,viqrhlf,viqihlf,viqhhlf,viqshlf,vilhlf,viihlf,     &
             vichlf, vithlf, pwhlf,tprhlf, gprhlf, cprhlf
  INTEGER :: raichlf,raighlf,raithlf

!
!-----------------------------------------------------------------------
!
!  define the attributes of zero contour to be plotted parameters
!
!-----------------------------------------------------------------------
!
  INTEGER :: uzro , vzro , wzro , ptzro , pzro, hzro, tzro,vhzro,vszro
  INTEGER :: qvzro ,kmhzro ,         &
             kmvzro,tkezro,rhzro ,rfzro ,ptezro,tdzro, qwzro, qtzro,    &
             rfczro, rhizro, zdrzro, kdpzro, zdpzro, rhvzro
  INTEGER :: qscalarzro(nscalarmax)
  INTEGER :: upzro ,vpzro ,wpzro ,ptpzro,ppzro ,qvpzro,vorpzro,divpzro
  INTEGER :: trnzro,wcpzro,raczro,ragzro,ratzro,pslzro,                 &
             capzro,cinzro,thezro,helzro,uhzro,brnzro,brnuzro,          &
             srlfzro,srmfzro,lizro,capszro,blcozro,ctczro
  INTEGER :: divqzro,griczro,avorzro
  INTEGER :: styzro,vtyzro,laizro,rouzro,vegzro,snowdzro
  INTEGER :: viqczro,viqrzro,viqizro,viqhzro,viqszro,vilzro,viizro,     &
             viczro, vitzro, pwzro, tprzro, gprzro, cprzro
  INTEGER :: raiczro,raigzro,raitzro

!
!-----------------------------------------------------------------------
!
!  Define the option for contour line stypes.
!
!-----------------------------------------------------------------------
!
  INTEGER :: usty , vsty , wsty , ptsty , psty, hsty, tsty,vhsty,vssty
  INTEGER :: qvsty ,kmhsty ,         &
             kmvsty,tkesty,rhsty ,rfsty ,ptesty,tdsty, qwsty, qtsty,    &
             rfcsty, rhisty, zdrsty, kdpsty, zdpsty, rhvsty
  INTEGER :: qscalarsty(nscalarmax)
  INTEGER :: upsty ,vpsty ,wpsty ,ptpsty,ppsty ,qvpsty,vorpsty,divpsty
  INTEGER :: trnsty,wcpsty,racsty,ragsty,ratsty,pslsty,                 &
             capsty,cinsty,thesty,helsty,uhsty,brnsty,brnusty,          &
             srlfsty,srmfsty,listy,capssty,blcosty,ctcsty
  INTEGER :: divqsty,gricsty,avorsty
  INTEGER :: stysty,vtysty,laisty,rousty,vegsty,snowdsty
  INTEGER :: viqcsty,viqrsty,viqisty,viqhsty,viqssty,vilsty,viisty,     &
             vicsty, vitsty, pwsty, tprsty, gprsty, cprsty
  INTEGER :: raicsty,raigsty,raitsty

  INTEGER :: msfplt,msfovr,msfcol1,msfcol2,msfprio,msfhlf,msfzro,msfsty
  REAL    :: msfinc, msfminc, msfmaxc

  INTEGER :: thkplt,thkovr,thkcol1,thkcol2,thkprio,thkhlf,thkzro,thksty
  REAL    :: thkinc, thkminc, thkmaxc

  INTEGER :: ipvplt,ipvovr,ipvcol1,ipvcol2,ipvprio,ipvhlf,ipvzro,ipvsty
  REAL    :: ipvinc, ipvminc, ipvmaxc
!
!-----------------------------------------------------------------------
!
!  Profile Plotting control parameters entered by user
!
!-----------------------------------------------------------------------
!
  INTEGER :: uprof, vprof, wprof, ptprof, pprof
  INTEGER :: qvprof,qcprof,qrprof,qiprof,qsprof,qhprof,kmhprof,         &
             kmvprof,tkeprof,rhprof,rfprof,pteprf
  INTEGER :: upprof,vpprof,wpprof,ptpprf,ppprof,qvpprf,vorpprf,divpprf
  INTEGER :: tsoilprof,qsoilprof  ! Zuwen He 05/31/2002
!
!-----------------------------------------------------------------------
!
!  Profile plot lower bound
!
!-----------------------------------------------------------------------
!
  REAL :: uprmin, vprmin, wprmin, ptprmin, pprmin
  REAL :: qvprmin,qcpmin,qrpmin,qipmin,qspmin,qhpmin
  REAL :: kmhpmin,kmvpmin,tkepmin
  REAL :: rhpmin,rfpmin,ptepmin,vorppmin,divppmin,avormin
  REAL :: uppmin,vppmin,wppmin,ptppmin,pppmin,qvppmin
  REAL :: tsoilprofmin,qsoilprofmin  ! Zuwen He 06/03/2002
!
!-----------------------------------------------------------------------
!
!  Profile plot upper bound
!
!-----------------------------------------------------------------------
!
  REAL :: uprmax, vprmax, wprmax, ptprmax, pprmax
  REAL :: qvprmax,qcpmax,qrpmax,qipmax,qspmax,qhpmax
  REAL :: kmhpmax, kmvpmax, tkepmax
  REAL :: rhpmax,rfpmax,ptepmax,vorppmax,divppmax, avormax
  REAL :: uppmax,vppmax,wppmax,ptppmax,pppmax,qvppmax
  REAL :: tsoilprofmax,qsoilprofmax  ! Zuwen He 06/03/2002
!
!-----------------------------------------------------------------------
!
!  3-D wireframe plotting
!
!-----------------------------------------------------------------------
!
  INTEGER :: w3dplt, q3dplt
  REAL    :: wisosf,qisosf

  INTEGER :: idisplay
  INTEGER :: imove, inwfrm
!
!-----------------------------------------------------------------------
!
!  Common blocks for plotting control parameters
!
!-----------------------------------------------------------------------
!
  INTEGER :: layover               ! set by OVERLAY, used in ctr2d,vtr2d
  REAL    :: ctinc,ctmin,ctmax,vtunt  ! contour interval and vector unit
  INTEGER :: icolor,icolor1,lbcolor,trcolor    ! required color
  CHARACTER (LEN=12) :: varname                ! save the plot variable

  REAL :: x01,y01                  ! first interpolation point for a
                                   ! vertical slice specified by two pts
  REAL :: x02,y02                  ! second interpolation point for a
                                   ! vertical slice specified by two pts
  REAL :: zlevel                   ! altitude (meters) of a horizontal
                                   ! slice so specified
  REAL :: sinaf,cosaf,dist,sqrtdxy

  COMMON /laypar/  layover
  COMMON /incunt/  ctinc,ctmin,ctmax,vtunt
  COMMON /recolor/ icolor,icolor1,lbcolor,trcolor
  COMMON /varplt1/ varname
  COMMON /slicev/  x01,y01,x02,y02,sinaf,cosaf,dist,sqrtdxy
  COMMON /sliceh/  zlevel

  REAL :: yxstrch    ! Stretching factor for x-y plots
  REAL :: zxstrch    ! Stretching factor for x-z plots
  REAL :: zystrch    ! Stretching factor for y-z plots
  REAL :: zhstrch    ! Stretching factor for arbitrary vertical slices
  REAL :: zsoilxstrch    ! Stretching factor for x-z plots for the soil model
  REAL :: zsoilystrch    ! Stretching factor for y-z plots for the soil model

  REAL :: aspratio
  COMMON /yratio/ aspratio      ! scaling factor: the y/x ratio.
!
!  Pass the plotting window parameters into subroutine encodwd.
!
  INTEGER :: iskip, jskip
  COMMON /pltwdw/ xbgn,xend,ybgn,yend,iskip,jskip

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!

  CHARACTER (LEN=80)  :: label

  CHARACTER (LEN=256) :: outfilename
  INTEGER :: iwtype,exchange

  !
  !  Variables defined by Youngsun Jung for dual-pol calculation
  !
  REAL :: pi
  REAL :: cx(6)     ! (PI/6)*rho_qx
  REAL, ALLOCATABLE :: alpha(:,:,:,:)   ! Shape parameter
!
!-----------------------------------------------------------------------
!
  INTEGER :: length
  CHARACTER (LEN=6) :: stem2
  CHARACTER (LEN=1) :: stem1
  INTEGER :: flagsin
  REAL    :: f_cputime,cpu0,cpu1,cpu2
  DOUBLE PRECISION :: f_walltime,second1, second2

  INTEGER :: ireturn
  INTEGER :: slicopt
  INTEGER :: i,j,k,ibgn,iend,jbgn,jend,kbgn,kend,ist,jst,kst,iob
  INTEGER :: ksoilbgn,ksoilend ! Zuwen He, 05/31/2002
  INTEGER :: ibgn1, iend1
  REAL    :: x11, y11
  INTEGER :: nxpic,nypic,islice,jslice,kslice,layout
  INTEGER :: idxfmt,filetm
  INTEGER :: isoilslice,jsoilslice,ksoilslice   ! Zuwen 05/31/2002
  REAL :: time,angl,aspect
  REAL :: time1c, time2c, time1g, time2g, time1t, time2t, timediff

  REAL :: xr,yr,zr,x1,x2,y1,y2,z1,z2,xlimit,ylimit
  REAL :: zmax,zmin
  REAL :: zsoilr,zsoil1,zsoil2
  REAL :: zsoilmax,zsoilmin  ! Zuwen He, 05/31/2002
  REAL :: pmb,drot
  INTEGER :: ifile
  INTEGER :: imax,jmax,kmax,imin,jmin,kmin
  REAL :: wmin, wmax, tk, tdk, psta
  REAL :: xbgn,xend,ybgn,yend,zbgn,zend
  REAL :: zsoilbgn,zsoilend ! Zuwen He, 05/31/2002  <0
  LOGICAL :: fexist
  LOGICAL :: rdallhist
  INTEGER :: onvf
  LOGICAL :: need_uv
  INTEGER :: ii, jj, nzmax
  CHARACTER(LEN=3), PARAMETER :: fmtstr(11) = (/'bin','asc','hdf','pak',&
                              'svi','bn2','net','net','gad','grb','v5d'/)

  REAL :: utmax,utmin,vtmax,vtmin,wtmax,wtmin,vsmax,vsmin
  REAL :: qvmax,qvmin
  REAL :: qscalarmax(nscalarmax),qscalarmin(nscalarmax)
  REAL :: qwmax,qwmin
  REAL :: rhmax,rhmin,rfmax,rfmin,ptemax,ptemin,upmax,upmin
  REAL :: vpmax,vpmin,wpmax,wpmin,ptpmax,ptpmin,ppmax,ppmin
  REAL :: qvpmax,qvpmin,vormax,vormin,divmax,divmin
  REAL :: hmin, hmax, tdmin,tdmax, rhimin, rhimax

  INTEGER :: nprof, profopt, nxprpic, nyprpic, npicprof
  REAL :: zprofbgn, zprofend
  REAL :: zsoilprofbgn, zsoilprofend  ! Zuwen He, 05/31/2002
  REAL :: time0

  INTEGER :: trnplt,wetcanplt
  INTEGER :: soiltpplt,vegtpplt,laiplt,rouplt,vegplt,snowdplt
  INTEGER :: soiltpn       ! number of soil type 1 to 4
  INTEGER :: pslplt
  INTEGER :: raincplt,raingplt,raintplt
  INTEGER :: capeplt, cinplt, thetplt, heliplt, uhplt
  INTEGER :: brnplt, brnuplt, srlfplt, srmfplt, liplt
  INTEGER :: capsplt, blcoplt, ctcplt
  INTEGER :: rainicplt,rainigplt,rainitplt

  INTEGER :: ip,ipp,ipriority(nprio),iptemp(nprio),sigplt(nprio)

  REAL :: trninc,wcpinc
  REAL :: soiltpinc,vegtpinc,laiinc,rouinc,veginc,snowdinc
  REAL :: pslinc
  REAL :: raincinc, rainginc,raintinc
  REAL :: capeinc, cininc, thetinc, heliinc, uhinc
  REAL :: brninc, brnuinc, srlfinc, srmfinc, liinc
  REAL :: capsinc, blcoinc, ctcinc
  REAL :: rainicinc,rainiginc,rainitinc

!
! 05/30/2002 Zuwen He
!
! soil plot options and parameters
!
  INTEGER :: tsoilplt
  REAL    :: tsoilinc,tsoilminc,tsoilmaxc
  INTEGER :: tsoilovr,tsoilhlf,tsoilzro,tsoilcol1,tsoilcol2,tsoilprio

  INTEGER :: qsoilplt
  REAL    :: qsoilinc,qsoilminc,qsoilmaxc
  INTEGER :: qsoilovr,qsoilhlf,qsoilzro,qsoilcol1,qsoilcol2,qsoilprio

!
! END Zuwen He
!

  INTEGER :: ovrmap,mapgrid,mapcol(maxmap),mapline_style(maxmap),       &
             mapgridcol,nmapfile
  REAL    :: latgrid,longrid
  CHARACTER (LEN=256) :: mapfile(maxmap)
  INTEGER :: lmapfile
  COMMON /mappar / ovrmap
  COMMON /mappar1/ nmapfile,mapcol,mapline_style,mapfile
  COMMON /mappar2/ mapgrid,mapgridcol, latgrid,longrid

  REAL    :: ztmin,ztmax
  INTEGER :: ovrtrn
  COMMON /trnpar/ trnplt,ovrtrn,trninc,trnminc,trnmaxc,                 &
                  ztmin,ztmax

  INTEGER :: ovrobs,obsset,obscol,obs_marktyp
  REAL :: obs_marksz
  COMMON /obspar/ ovrobs,obsset,obscol,obs_marktyp, obs_marksz
  CHARACTER (LEN=256) :: sfcobfl(mxsfcobfl)
  INTEGER :: nsfcobfl,obunit,lsfcobfl

  INTEGER :: ovrstaopt
  INTEGER :: ovrstam,staset,ovrstan,ovrstav,stacol,markprio,wrtstax
  INTEGER :: nsta_typ,sta_typ(30),sta_marktyp(30),sta_markcol(30)
  REAL :: sta_marksz(30)
  REAL :: wrtstad
  CHARACTER (LEN=256) :: stalofl
  COMMON /sta_par/ ovrstaopt,ovrstam,staset,ovrstan,ovrstav,stacol,     &
         markprio,nsta_typ,sta_typ,sta_marktyp,                         &
         sta_markcol,sta_marksz,stalofl,wrtstax,wrtstad
  INTEGER :: stunit,lstalofl
!
!-----------------------------------------------------------------------
!
!  Surface (single-level) read-in observation variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=24) :: atime
  CHARACTER (LEN=5)  :: stn(mxsfcob)
  INTEGER :: n_meso_g,                                                  &
             n_meso_pos,n_sao_g,n_sao_pos_g,n_sao_b,n_sao_pos_b,        &
             n_obs_g,n_obs_pos_g,n_obs_b,n_obs_pos_b,nobs
  CHARACTER (LEN=8) :: obstype(mxsfcob)
  CHARACTER (LEN=8) :: wx(mxsfcob)
  CHARACTER (LEN=1) :: store_emv(mxsfcob,5)
  CHARACTER (LEN=4) :: store_amt(mxsfcob,5)
  REAL :: latob(mxsfcob),lonob(mxsfcob),elevob(mxsfcob)
  REAL :: tob(mxsfcob),tdob(mxsfcob),dd(mxsfcob),ff(mxsfcob)
  REAL :: ddg(mxsfcob),ffg(mxsfcob)
  REAL :: pstn(mxsfcob),pmsl(mxsfcob),alt(mxsfcob)
  REAL :: ceil(mxsfcob),lowcld(mxsfcob),cover(mxsfcob)
  REAL :: vis(mxsfcob),rad(mxsfcob)
  REAL :: store_hgt(mxsfcob,5)
  INTEGER :: kloud(mxsfcob),idp3(mxsfcob)
  INTEGER :: obstime(mxsfcob)
  REAL :: obs1(mxsfcob),obs2(mxsfcob)

  COMMON /sfc_obs1/ nobs
  COMMON /sfc_obs2/ latob,lonob,obs1,obs2

  INTEGER :: nsta,nstapro(mxstalo),nstatyp(mxstalo)
  REAL :: latsta(mxstalo), lonsta(mxstalo)
  CHARACTER (LEN=2 ) :: s_state(mxstalo)
  CHARACTER (LEN=5 ) :: s_name(mxstalo)
  CHARACTER (LEN=20) :: s_site(mxstalo)
  INTEGER :: s_elev(mxstalo)
  COMMON /sta_loc/latsta,lonsta,nstatyp,nstapro,nsta
  COMMON /sta_loc1/s_name

  REAL :: lblmag   ! A global magnification factor for labels.
  REAL :: ctrlbsiz, axlbsiz
  COMMON /labmag/ lblmag, ctrlbsiz, axlbsiz

  REAL :: winsiz   ! A global factor for window size
  COMMON /windows/ winsiz
  REAL :: margnx, margny   ! margin
  INTEGER :: pcolbar ! position of color bar
  INTEGER :: tickopt
  INTEGER :: axlbfmt   ! format for the axis's label
  INTEGER :: fontopt   ! the font of character
  INTEGER :: ctrlbopt  ! Contour labelling option
  INTEGER :: ctrstyle
  INTEGER :: ctrlbfrq
  COMMON /clb_frq/ ctrlbfrq

  INTEGER :: lbmaskopt
  INTEGER :: haxisu, vaxisu
  INTEGER :: lbaxis
  INTEGER :: lnmag
  REAL :: hmintick,vmajtick,vmintick,hmajtick
  COMMON /var_par/ fontopt,haxisu, vaxisu,lbaxis,tickopt,               &
          hmintick,vmajtick,vmintick,hmajtick,axlbfmt
  INTEGER :: presaxis_no
  REAL :: pres_val(20), pres_z(20)
  COMMON /pressbar_par/presaxis_no,pres_val,pres_z

  INTEGER :: ntitle,titcol, wpltime
  REAL :: titsiz
  CHARACTER (LEN=256) :: title(3), footer_l, footer_c, footer_r

  COMMON /titpar1/title, footer_l, footer_c, footer_r
  COMMON /titpar2/ntitle,titcol,wpltime, nxpic, nypic
  COMMON /titpar3/titsiz

  COMMON /tmphc1/x1,x2,y1,y2,z1,z2
!
!  Work arrays used by contouring and contour color fill routines
!  m*n=<50000, 8*m=<m if the array to be contours is demensioned m x n.
!

  INTEGER :: sovrlay
  LOGICAL :: plotovr
!
!-----------------------------------------------------------------------
!
!  Common block used by the Ncar Graphic streamline routine
!  INITA  Used to precondition grid boxes to be eligible to
!         start a streamline. For example, a value of 4
!         means that every fourth grid box is eligible ;
!         a value of 2 means that every other grid box is eligible.
!         (see INITB)
!
!  INITB  Used to precondition grid boxes to be eligible for
!         direction arrows. If the user changes the default values
!         of INITA and/or INITB, it should be done such that
!         MOD(INITA,INITB) = 0. For a dense grid try INITA=4 and
!         INITB=2 to reduce the CPU time.
!
!-----------------------------------------------------------------------
!
!  INTEGER :: inita , initb , iterp , iterc , igflg, imsg , icyc
!  REAL :: arowl , uvmsg , displ , dispc , cstop
!  COMMON /str03/  inita , initb , arowl , iterp , iterc , igflg,        &
!                  imsg , uvmsg , icyc , displ , dispc , cstop

!-----------------------------------------------------------------------
!
!  Namelist for plot input
!
!-----------------------------------------------------------------------
!
  INTEGER :: maxslicopt

  PARAMETER (maxslicopt=12)
!
! slicopt = 1: xy slice
! slicopt = 2: xz slice
! slicopt = 3: yz slice
! slicopt = 4: constant height
! slicopt = 5: vertical slice through two points
! slicopt = 6: constant pressure
! slicopt = 7: isentropic level
! slicopt = 8: surface xy slice
! slicopt = 9: xy-soil-slice for the soil model
! slicopt = 10: xz-soil-slice for the soil model
! slicopt = 11: yz-soil-slice for the soil model
! slicopt = 12: Vertical slices along trajectories


  INTEGER :: nslice(maxslicopt), indxslic
  INTEGER :: iplot(maxslicopt)

  INTEGER :: nslice_xy,nslice_xz, nslice_yz, nslice_h, nslice_v,        &
             nslice_p,nslice_pt

  INTEGER :: slice_xy(max_dim),slice_xz(max_dim), slice_yz(max_dim)

!
! 05/30/2002 Zuwen He
! Slice for soil variables
!

  INTEGER :: nslice_xy_soil
  INTEGER :: slice_xy_soil(max_dim)

  INTEGER :: nslice_xz_soil
  INTEGER :: slice_xz_soil(max_dim)

  INTEGER :: nslice_yz_soil
  INTEGER :: slice_yz_soil(max_dim)
!
! END Zuwen
!
  REAL :: xi1(max_dim),yi1(max_dim),xi2(max_dim),yi2(max_dim)
  REAL :: slice_h(max_dim), slice_p(max_dim),  slice_pt(max_dim)
  REAL :: xpnt1(max_dim),ypnt1(max_dim), xpnt2(max_dim),ypnt2(max_dim)

  INTEGER :: ucol1,vcol1,wcol1,ptcol1,pcol1,qvcol1,qscalarcol1(nscalarmax),              &
          kmhcol1,kmvcol1,tkecol1,          &
          rhcol1,rfcol1,rfccol1,ptecol1,upcol1,vpcol1,wpcol1,           &
          ptpcol1,ppcol1,qvpcol1,vorpcol1,divpcol1,vtrcol1,             &
          vtpcol1,vagcol1,vtrstmcol1,vtpstmcol1,hcol1,tcol1,            &
          divqcol1,zdrcol1,kdpcol1,zdpcol1,rhvcol1,                     &
          tdcol1,qwcol1,qtcol1,trncol1,                                 &
          wcpcol1,raccol1,ragcol1,ratcol1,pslcol1,stycol1,              &
          vtycol1,laicol1,roucol1,vegcol1,snowdcol1,                    &
          capcol1,cincol1,thecol1,helcol1,uhcol1,                       &
          vhcol1, vscol1,brncol1,brnucol1,srlfcol1,                     &
          srmfcol1,licol1, capscol1, blcocol1,griccol1,avorcol1,        &
          viqccol1,viqrcol1,viqicol1,viqhcol1,viqscol1,vilcol1,         &
          viicol1,viccol1, xuvcol1, strmcol1, ctccol1, rhicol1,         &
          vitcol1, pwcol1, tprcol1, gprcol1, cprcol1,                   &
          raiccol1,raigcol1,raitcol1

  INTEGER :: ucol2,vcol2,wcol2,ptcol2,pcol2,qvcol2,qscalarcol2(nscalarmax), &
          kmhcol2,kmvcol2,tkecol2,                                      &
          rhcol2,rfcol2,rfccol2,ptecol2,upcol2,vpcol2,wpcol2,           &
          ptpcol2,ppcol2,qvpcol2,vorpcol2,divpcol2,vtrcol2,             &
          vtpcol2,vagcol2,vtrstmcol2,vtpstmcol2,hcol2,tcol2,            &
          divqcol2,zdrcol2,kdpcol2,zdpcol2,rhvcol2,                     &
          tdcol2,qwcol2,qtcol2,trncol2,                                 &
          wcpcol2,raccol2,ragcol2,ratcol2,pslcol2,stycol2,              &
          vtycol2,laicol2,roucol2,vegcol2,snowdcol2,                    &
          capcol2,cincol2,thecol2,helcol2,                              &
          uhcol2,vhcol2,vscol2,brncol2,brnucol2,srlfcol2,               &
          srmfcol2,licol2, capscol2, blcocol2,griccol2,avorcol2,        &
          viqccol2,viqrcol2,viqicol2,viqhcol2,viqscol2,vilcol2,         &
          viicol2,viccol2, xuvcol2, strmcol2, ctccol2, rhicol2,         &
          vitcol2, pwcol2, tprcol2, gprcol2, cprcol2,                   &
          raiccol2,raigcol2,raitcol2

  INTEGER :: uprio,vprio,wprio,ptprio,pprio,qvprio,qscalarprio(nscalarmax), &
          kmhprio,kmvprio,tkeprio,                                      &
          rhprio,rfprio,rfcprio,rhiprio,                                &
          zdrprio,kdpprio,zdpprio,rhvprio,                              &
          pteprio,upprio,vpprio,wpprio,ptpprio,ppprio,qvpprio,          &
          vorpprio,divpprio,vtrprio,vtpprio,vagprio,vtrstmprio,         &
          vtpstmprio,hprio,tprio,divqprio,tdprio,qwprio,qtprio,         &
          trnprio,wcpprio,racprio, ragprio,ratprio,pslprio,             &
          styprio,vtyprio,laiprio,rouprio,vegprio,snowdprio,            &
          capprio,cinprio,theprio,helprio,uhprio,vhprio,vsprio,         &
          brnprio,bruprio,srlprio,srmprio,liprio,                       &
          capsprio,blcoprio,gricprio,avorprio,                          &
          viqcprio,viqrprio,viqiprio,viqhprio,viqsprio,vilprio,         &
          viiprio,vicprio, xuvprio, strmprio, ctcprio, vitprio,         &
          pwprio, tprprio, gprprio, cprprio,                            &
          raicprio,raigprio,raitprio

  INTEGER :: col_table,smooth
  CHARACTER (LEN=256) :: color_map
  INTEGER :: number_of_boxes, boxcol
  REAL    :: bx1(10), bx2(10),by1(10),by2(10)

  INTEGER :: number_of_polys, polycol
  REAL    :: vertx(max_verts,max_polys), verty(max_verts,max_polys)

  INTEGER :: istride,jstride,kstride

  COMMON /coltable/col_table,pcolbar
  COMMON /smoothopt/smooth
  COMMON /boxesopt/number_of_boxes,boxcol,bx1,bx2,by1,by2
  COMMON /polysopt/number_of_polys,polycol,vertx,verty

  CHARACTER (LEN=1) :: tunits  ! units for temperature F or C
  CHARACTER (LEN=1) :: tdunits ! units for dew-point temp F or C

  INTEGER :: vhunits,vtrunits,vtpunits,vagunits,xuvunits,strmunits
  INTEGER :: vtrtype,vtptype,vagtype,xuvtype, strmtype

  INTEGER :: tprunits, gprunits, cprunits

  INTEGER :: iunits, itype
  COMMON /windvtr/iunits, itype

  INTEGER :: racunit, ragunit, ratunit
  INTEGER :: raicunit,raigunit,raitunit

  INTEGER :: ovrlaymulopt, ovrmul_num
  CHARACTER (LEN=12) :: ovrmulname(50)
  CHARACTER (LEN=12) :: ovrname

  INTEGER :: setcontopt ,setcontnum
  CHARACTER (LEN=12) :: setcontvar(maxuneva)
  REAL :: setconts(maxunevm,maxuneva)
  COMMON /setcont_var/setcontvar
  COMMON /setcon_par/setcontopt,setcontnum,setconts

  COMMON /init5_coltab/ color_map

  INTEGER :: arbvaropt   ! plot arbitrary variable

  INTEGER :: istatus

  INTEGER :: finfmt3d(maxarbvar), finfmt2d(maxarbvar)

  CHARACTER (LEN=256) :: dirname3d(maxarbvar),dirname2d(maxarbvar)
  CHARACTER (LEN=256) :: filename3d(maxarbvar),filename2d(maxarbvar)
  CHARACTER (LEN=40) :: varunits, var_name
  CHARACTER (LEN=6) :: var3d(maxarbvar)
  INTEGER :: var3dnum
  INTEGER :: var3dplot(maxarbvar)
  REAL    :: var3dinc(maxarbvar), var3dminc(maxarbvar),                 &
             var3dmaxc(maxarbvar)
  INTEGER :: var3dovr(maxarbvar),var3dcol1(maxarbvar),                  &
          var3dcol2(maxarbvar),var3dprio(maxarbvar),                    &
          var3dhlf(maxarbvar),var3dzro(maxarbvar),                      &
          var3dsty(maxarbvar)

  REAL, ALLOCATABLE :: var3dv(:,:,:)

  CHARACTER (LEN=6) :: var2d(maxarbvar)
  INTEGER :: var2dnum
  INTEGER :: var2dplot(maxarbvar)
  REAL :: var2dinc(maxarbvar), var2dminc(maxarbvar),                    &
       var2dmaxc(maxarbvar)
  INTEGER :: var2dovr(maxarbvar),var2dcol1(maxarbvar),                  &
          var2dcol2(maxarbvar), var2dprio(maxarbvar),                   &
          var2dhlf(maxarbvar),var2dzro(maxarbvar),                      &
          var2dsty(maxarbvar)

  INTEGER :: iast, jast
  INTEGER :: vtr2dnum
  CHARACTER (LEN=256) :: diruv2d(maxarbvar)
  CHARACTER (LEN=256) :: filenameu2d(maxarbvar),filenamev2d(maxarbvar)
  INTEGER :: finfmtuv2d(maxarbvar)
  CHARACTER (LEN=6)  :: vtru2d(maxarbvar),vtrv2d(maxarbvar)
  INTEGER :: iastride(maxarbvar),jastride(maxarbvar)
  INTEGER :: vtraplt(maxarbvar), magaplt(maxarbvar)
  REAL :: vtraunit(maxarbvar)
  REAL :: magainc(maxarbvar),magaminc(maxarbvar),magamaxc(maxarbvar)
  INTEGER :: vtraovr(maxarbvar), magaovr(maxarbvar)
  INTEGER :: magahlf(maxarbvar), magazro(maxarbvar)
  INTEGER :: vtracol1(maxarbvar), magacol1(maxarbvar)
  INTEGER :: vtracol2(maxarbvar), magacol2(maxarbvar)
  INTEGER :: vtraprio(maxarbvar), magaprio(maxarbvar)
  INTEGER :: vtraunits(maxarbvar), magaunits(maxarbvar)
  INTEGER :: vtratype(maxarbvar), magasty(maxarbvar)

  REAL, ALLOCATABLE :: var2du(:,:)
  REAL, ALLOCATABLE :: var2dv(:,:)

  INTEGER :: missfill_opt,missval_colind    ! miss value color index
  COMMON /multi_value/ missfill_opt, missval_colind

  INTEGER :: xnwpic_called
  COMMON /callnwpic/xnwpic_called

  COMMON /init8_cntl1/       &
      hplot,  hinc,   hminc,   hmaxc,  hovr,   hcol1,hcol2, hprio,      &
            hhlf, hzro, hsty,                                           &
      msfplt,msfinc,msfminc, msfmaxc,msfovr,msfcol1,msfcol2,msfprio,    &
            msfhlf, msfzro,  msfsty,                                    &
      thkplt,thkinc,thkminc, thkmaxc,thkovr,thkcol1,thkcol2,thkprio,    &
            thkhlf, thkzro,  thksty,                                    &
      tplot,  tinc,   tminc,   tmaxc,  tovr,   tcol1,tcol2, tprio,      &
            thlf, tzro, tsty,                                   &
      uplot,  uinc,   uminc,   umaxc,  uovr,   ucol1,ucol2, uprio,      &
            uhlf, uzro, usty,                                           &
      vplot,  vinc,   vminc,   vmaxc,  vovr,   vcol1,vcol2, vprio,      &
            vhlf, vzro, vsty,                                           &
      vhplot, vhinc,  vhminc,  vhmaxc, vhovr,  vhcol1,vhcol2,vhprio,    &
            vhunits, vhhlf, vhzro, vhsty,                               &
      vsplot, vsinc,  vsminc,  vsmaxc, vsovr,  vscol1,vscol2,vsprio,    &
            vshlf, vszro, vssty,                                        &
      wplot,  winc,   wminc,   wmaxc,  wovr,   wcol1,wcol2, wprio,      &
            whlf, wzro, wsty,                                           &
      ptplot, ptinc,  ptminc,  ptmaxc, ptovr,  ptcol1,ptcol2,ptprio,    &
            pthlf, ptzro, ptsty,                                        &
      pplot , pinc,   pminc,   pmaxc,  povr,   pcol1,pcol2, pprio,      &
            phlf, pzro, psty,                                           &
      ipvplt,ipvinc,ipvminc, ipvmaxc,ipvovr,ipvcol1,ipvcol2,ipvprio,    &
            ipvhlf, ipvzro,  ipvsty

  COMMON /init9_cntl2/                                              &
      qvplot, qvinc,  qvminc,  qvmaxc, qvovr,  qvcol1,qvcol2,qvprio,    &
            qvhlf, qvzro, qvsty,                                        &
      qscalarplot,qscalarinc,qscalarminc,qscalarmaxc,qscalarovr,        &
            qscalarcol1,qscalarcol2,qscalarprio,qscalarhlf,             &
            qscalarzro,qscalarsty,                                      &
      qwplot, qwinc,  qwminc,  qwmaxc, qwovr,  qwcol1,qwcol2,qwprio,    &
            qwhlf, qwzro, qwsty,                                        &
      qtplot, qtinc,  qtminc,  qtmaxc, qtovr,  qtcol1,qtcol2,qtprio,    &
            qthlf, qtzro, qtsty
  COMMON /init10_cntl3/                                              &
      kmhplt, kmhinc, kmhminc,kmhmaxc,kmhovr,kmhcol1,kmhcol2,kmhprio,   &
            kmhhlf, kmhzro, kmhsty,                                     &
      kmvplt, kmvinc, kmvminc,kmvmaxc,kmvovr,kmvcol1,kmvcol2,kmvprio,   &
            kmvhlf, kmvzro, kmvsty,                                     &
      tkeplt, tkeinc, tkeminc, tkemaxc,tkeovr,tkecol1,tkecol2,tkeprio,  &
            tkehlf, tkezro, tkesty,                                     &
      rhplot, rhinc,  rhminc,  rhmaxc,  rhovr,   rhcol1,rhcol2,rhprio,  &
            rhhlf, rhzro,   rhsty,                                      &
      tdplot, tdinc,  tdminc,  tdmaxc,  tdovr,   tdcol1,tdcol2,tdprio,  &
            tdhlf, tdzro, tdsty,                                        &
      dualpol, rfopt, rsadir, wavelen,                                  &
      rfplot, rfinc,  rfminc,  rfmaxc,  rfovr,   rfcol1,rfcol2,rfprio,  &
            rfhlf, rfzro, rfsty,                                        &
      rfcplt, rfcinc, rfcminc, rfcmaxc,rfcovr,rfccol1,rfccol2,rfcprio,  &
            rfchlf, rfczro, rfcsty,                                     &
      pteplt, pteinc, pteminc, ptemaxc,pteovr,ptecol1,ptecol2,pteprio,  &
            ptehlf, ptezro, ptesty,                                     &
      zdrplt, zdrinc, zdrminc, zdrmaxc,zdrovr,zdrcol1,zdrcol2,zdrprio,  &
            zdrhlf, zdrzro, zdrsty,                                     &
      kdpplt, kdpinc, kdpminc, kdpmaxc,kdpovr,kdpcol1,kdpcol2,kdpprio,  &
            kdphlf, kdpzro, kdpsty,                                     &
      zdpplt, zdpinc, zdpminc, zdpmaxc,zdpovr,zdpcol1,zdpcol2,zdpprio,  &
            zdphlf, zdpzro, zdpsty,                                     &
      rhvplt, rhvinc, rhvminc, rhvmaxc,rhvovr,rhvcol1,rhvcol2,rhvprio,  &
            rhvhlf, rhvzro, rhvsty

  COMMON /init810_char_units/ tunits, tdunits

  COMMON /init11_cntl_prt1/                                             &
      upplot, upinc,  upminc,   upmaxc,   upovr,upcol1,upcol2,upprio,   &
            uphlf, upzro, upsty,                                        &
      vpplot, vpinc,  vpminc,   vpmaxc,   vpovr,vpcol1,vpcol2,vpprio,   &
            vphlf, vpzro, vpsty,                                        &
      wpplot, wpinc,  wpminc,   wpmaxc,   wpovr,wpcol1,wpcol2,wpprio,   &
            wphlf, wpzro, wpsty,                                        &
      ptpplt, ptpinc, ptpminc,ptpmaxc,ptpovr,ptpcol1,ptpcol2,ptpprio,   &
            ptphlf, ptpzro, ptpsty,                                     &
      ppplot, ppinc,  ppminc, ppmaxc,  ppovr,   ppcol1,ppcol2,ppprio,   &
            pphlf, ppzro, ppsty,                                        &
      qvpplt, qvpinc, qvpminc,qvpmaxc,qvpovr,qvpcol1,qvpcol2,qvpprio,   &
            qvphlf, qvpzro, qvpsty,                                     &
      vorpplt,vorpinc,vorpminc, vorpmaxc, vorpovr, vorpcol1,vorpcol2,   &
            vorphlf,  vorpprio, vorpzro, vorpsty,                       &
      divpplt,divpinc,divpminc, divpmaxc, divpovr, divpcol1,divpcol2,   &
            divphlf,  divpprio, divpzro, divpsty,                       &
      divqplt,divqinc,divqminc, divqmaxc, divqovr, divqcol1,divqcol2,   &
            divqhlf,divqprio, divqzro,divqsty
  COMMON /init12_cntl_prt2/                                             &
      gricplt,gricinc,gricminc, gricmaxc, gricovr, griccol1,griccol2,   &
            grichlf,gricprio, griczro, gricsty,                         &
      avorplt,avorinc,avorminc, avormaxc, avorovr, avorcol1,avorcol2,   &
            avorhlf,avorprio, avorzro, avorsty,                         &
      rhiplot, rhiinc,  rhiminc,  rhimaxc,  rhiovr,   rhicol1,rhicol2,  &
            rhiprio,  rhihlf, rhizro , rhisty

  COMMON /init13_cntl_vctr/ istride,jstride,kstride,                    &
      vtrplt, vtrunit,vtrovr,vtrcol1,vtrcol2,vtrprio,vtrunits,vtrtype,  &
      vtpplt, vtpunit, vtpovr,vtpcol1,vtpcol2,vtpprio,vtpunits,vtptype, &
      xuvplt, xuvunit,xuvovr,xuvcol1,xuvcol2,xuvprio,xuvunits,xuvtype,  &
      strmplt,strmunit,strmovr,strmcol1,strmcol2,strmprio,strmunits,    &
      strmtype,                                                         &
      vagplt, vagunit,vagovr,vagcol1,vagcol2,vagprio,vagunits,vagtype

  COMMON /init14_cntl_strm/                                             &
      vtrstrm, vtrstmovr, vtrstmcol1, vtrstmcol2, vtrstmprio,           &
      vtpstrm, vtpstmovr, vtpstmcol1, vtpstmcol2, vtpstmprio

  COMMON /init4_plotset/ iorig,zbgn,zend,zsoilbgn,zsoilend,             &
                      yxstrch,zxstrch,zystrch,                          &
                      zhstrch,zsoilxstrch,zsoilystrch,                  &
                      margnx,margny
  COMMON /init6_style/ lnmag, ctrlbopt, ctrstyle

  COMMON /init7_slice/ nslice_xy, nslice_xz, nslice_yz, nslice_h,        &
                       nslice_v,  nslice_p,  nslice_pt,                  &
                       nslice_xy_soil, nslice_xz_soil, nslice_yz_soil,   &
                       slice_xy, slice_xz, slice_yz, slice_h,            &
                       slice_p,  slice_pt, xpnt1, ypnt1, xpnt2, ypnt2,   &
                       slice_xy_soil, slice_xz_soil, slice_yz_soil,      &
                       imove

  COMMON /init15_sfc/                                                   &
      trnovr,trncol1,trncol2,trnprio,trnhlf, trnzro, trnsty,            &
      wetcanplt,wcpinc,wcpminc,wcpmaxc,wcpovr,wcpcol1,wcpcol2,wcpprio,  &
            wcphlf, wcpzro, wcpsty,                                     &
      raincplt,raincinc,raincminc,raincmaxc,racovr,raccol1,raccol2,     &
            rachlf, racprio, raczro, racsty, racunit,                   &
      raingplt,rainginc,raingminc,raingmaxc,ragovr,ragcol1,ragcol2,     &
            raghlf,  ragprio, ragzro, ragsty, ragunit,                  &
      raintplt,raintinc,raintminc,raintmaxc,ratovr,ratcol1,ratcol2,     &
            rathlf,  ratprio, ratzro, ratsty, ratunit,                  &
      rainicplt,rainicinc,rainicminc,rainicmaxc,raicovr,raiccol1,       &
            raiccol2,raichlf,raicprio,raiczro,raicsty,raicunit,         &
      rainigplt,rainiginc,rainigminc,rainigmaxc,raigovr,raigcol1,       &
            raigcol2,raighlf,raigprio,raigzro,raigsty,raigunit,         &
      rainitplt,rainitinc,rainitminc,rainitmaxc,raitovr,raitcol1,       &
            raitcol2,raithlf,raitprio,raitzro,raitsty,raitunit

  COMMON /init19_soil/                                                  &
      tsoilplt,tsoilinc,tsoilminc,tsoilmaxc,tsoilovr,                   &
            tsoilcol1,tsoilcol2,tsoilhlf,tsoilprio,tsoilzro,            &
      qsoilplt,qsoilinc,qsoilminc,qsoilmaxc,qsoilovr,                   &
            qsoilcol1,qsoilcol2,qsoilhlf,qsoilprio,qsoilzro

  COMMON /init16_sfc/                                                  &
      pslplt ,pslinc, pslminc, pslmaxc,pslovr,pslcol1,pslcol2,pslprio,  &
            pslhlf, pslzro, pslsty,                                     &
      capeplt,capeinc,capeminc,capemaxc,capovr,capcol1,capcol2,capprio, &
            caphlf, capzro, capsty,                                     &
      cinplt, cininc, cinminc, cinmaxc, cinovr,cincol1,cincol2,cinprio, &
            cinhlf, cinzro, cinsty,                                     &
      thetplt,thetinc,thetminc,thetmaxc,theovr,thecol1,thecol2,theprio, &
            thehlf, thezro, thesty,                                     &
      heliplt,heliinc,heliminc,helimaxc,helovr,helcol1,helcol2,helprio, &
            helhlf, helzro, helsty,                                     &
      uhplt,uhinc,uhminc,uhmaxc,uhovr,uhcol1,uhcol2,uhprio,             &
            uhhlf, uhzro, uhsty,uhmnhgt,uhmxhgt,                        &
      brnplt, brninc, brnminc, brnmaxc, brnovr,brncol1,brncol2,brnprio, &
            brnhlf, brnzro, brnsty,                                     &
      brnuplt, brnuinc, bruminc, brumaxc, brnuovr, brnucol1,brnucol2,   &
            brnuhlf,  brnuzro, brnusty, bruprio,                        &
      srlfplt, srlfinc, srlminc, srlmaxc, srlfovr, srlfcol1,srlfcol2,   &
            srlfhlf,  srlfzro, srlfsty, srlprio,                        &
      srmfplt, srmfinc, srmminc, srmmaxc, srmfovr, srmfcol1,srmfcol2,   &
            srmfhlf, srmfzro, srmfsty, srmprio
  COMMON /init17_sfc/                                                  &
      liplt, liinc, liminc, limaxc, liovr, licol1,licol2,liprio,        &
            lihlf, lizro, listy,                                        &
      capsplt, capsinc, capsminc, capsmaxc, capsovr, capscol1,capscol2, &
            capshlf, capszro, capssty, capsprio,                        &
      blcoplt, blcoinc, blcominc, blcomaxc, blcoovr, blcocol1,blcocol2, &
            blcohlf, blcozro, blcosty, blcoprio,                        &
      viqcplt, viqcinc, viqcminc, viqcmaxc, viqcovr, viqccol1,viqccol2, &
            viqchlf, viqczro, viqcsty, viqcprio,                        &
      viqiplt, viqiinc, viqiminc, viqimaxc, viqiovr, viqicol1,viqicol2, &
            viqihlf, viqizro, viqisty, viqiprio,                        &
      viqrplt, viqrinc, viqrminc, viqrmaxc, viqrovr, viqrcol1,viqrcol2, &
            viqrhlf, viqrzro, viqrsty, viqrprio,                        &
      viqsplt, viqsinc, viqsminc, viqsmaxc, viqsovr, viqscol1,viqscol2, &
            viqshlf, viqszro, viqssty,viqsprio,                         &
      viqhplt, viqhinc, viqhminc, viqhmaxc, viqhovr, viqhcol1,viqhcol2, &
            viqhhlf, viqhzro, viqhsty,viqhprio,                         &
      vilplt, vilinc, vilminc, vilmaxc, vilovr, vilcol1,vilcol2,        &
            vilhlf, vilzro, vilsty, vilprio
  COMMON /init18_sfc/                                                  &
      viiplt, viiinc, viiminc, viimaxc, viiovr, viicol1,viicol2,        &
            viihlf, viizro, viisty,  viiprio,                           &
      vicplt, vicinc, vicminc, vicmaxc, vicovr, viccol1,viccol2,        &
            vichlf, viczro, vicsty, vicprio,                            &
      ctcplt, ctcinc, ctcminc, ctcmaxc, ctcovr, ctccol1,ctccol2,        &
            ctchlf, ctczro, ctcsty, ctcprio,                            &
      vitplt, vitinc, vitminc, vitmaxc, vitovr, vitcol1,vitcol2,        &
            vithlf, vitzro, vitsty, vitprio,                            &
      pwplt, pwinc, pwminc, pwmaxc, pwovr, pwcol1,pwcol2,               &
            pwhlf, pwzro, pwsty, pwprio,                                &
      tprplt, tprinc, tprminc, tprmaxc, tprovr, tprcol1,tprcol2,        &
            tprhlf, tprzro, tprsty, tprprio, tprunits,                  &
      gprplt, gprinc, gprminc, gprmaxc, gprovr, gprcol1,gprcol2,        &
            gprhlf, gprzro, gprsty, gprprio, gprunits,                  &
      cprplt, cprinc, cprminc, cprmaxc, cprovr, cprcol1,cprcol2,        &
            cprhlf, cprzro, cprsty, cprprio, cprunits

  COMMON /init20_sfccha/                                               &
      soiltpplt,soiltpinc,soiltpminc,soiltpmaxc,styovr,stycol1,stycol2, &
            styhlf, styzro, stysty,styprio,soiltpn,                     &
      vegtpplt,vegtpinc,vegtpminc,vegtpmaxc,vtyovr,vtycol1,vtycol2,     &
            vtyhlf, vtyzro, vtysty,vtyprio,                             &
      laiplt,laiinc,laiminc,laimaxc,laiovr,laicol1,laicol2,laiprio,     &
            laihlf, laizro, laisty,                                     &
      rouplt,rouinc,rouminc,roumaxc,rouovr,roucol1,roucol2,rouprio,     &
            rouhlf, rouzro, rousty,                                     &
      vegplt,veginc,vegminc,vegmaxc,vegovr,vegcol1,vegcol2,vegprio,     &
            veghlf, vegzro, vegsty,                                     &
      snowdplt,snowdinc,snowdminc,snowdmaxc,snowdovr,snowdcol1,         &
            snowdcol2, snowdprio,snowdhlf, snowdzro, snowdsty
  COMMON /init23_wirfrm/ w3dplt, wisosf, q3dplt, qisosf

  REAL :: paprlnth
  COMMON /init3_page/ layout, inwfrm, paprlnth

  COMMON /init25_prof/ profopt, nprof, xprof, yprof,                 &
      npicprof, uprof, uprmin, uprmax, vprof, vprmin, vprmax,           &
      wprof,wprmin,wprmax,  ptprof,ptprmin,ptprmax,                     &
      pprof,pprmin,pprmax,  qvprof,qvprmin,qvprmax,                     &
      qcprof,qcpmin,qcpmax, qrprof,qrpmin,qrpmax,                       &
      qiprof,qipmin,qipmax, qsprof,qspmin,qspmax,                       &
      qhprof,qhpmin,qhpmax, rhprof,rhpmin,rhpmax,                       &
      kmhprof,kmhpmin,kmhpmax, kmvprof,kmvpmin,kmvpmax,                 &
      tkeprof,tkepmin,tkepmax,                                          &
      rfprof,rfpmin,rfpmax, pteprf,ptepmin,ptepmax,                     &
      upprof,uppmin,uppmax, vpprof,vppmin,vppmax,                       &
      wpprof,wppmin,wppmax, ptpprf,ptppmin,ptppmax,                     &
      ppprof,pppmin,pppmax, qvpprf,qvppmin,qvppmax,                     &
      vorpprf, vorppmin, vorppmax, divpprf, divppmin, divppmax,         &
      zprofbgn,zprofend,                                                &
      tsoilprof,tsoilprofmin,tsoilprofmax,                              &
      qsoilprof,qsoilprofmin,qsoilprofmax,                              &
      zsoilprofbgn,zsoilprofend,                                        &
      nxprpic, nyprpic

  COMMON /init24_obs/ nsfcobfl,sfcobfl

  COMMON /init22_ovrlay/ovrlaymulopt,ovrname,ovrmul_num,ovrmulname

  COMMON /init21_cntl_arbvar/arbvaropt,                                 &
      var3dnum,dirname3d,finfmt3d,filename3d,                           &
      var3d,var3dplot, var3dinc, var3dminc,var3dmaxc,                   &
      var3dovr, var3dhlf, var3dzro,var3dsty,var3dcol1, var3dcol2,       &
      var3dprio, var2dnum,dirname2d,finfmt2d,filename2d,                &
      var2d,var2dplot, var2dinc, var2dminc,var2dmaxc,                   &
      var2dovr, var2dhlf, var2dzro, var2dsty, var2dcol1, var2dcol2,     &
      var2dprio
  COMMON /init26_cntl_vctra/ vtr2dnum,diruv2d,vtru2d,vtrv2d,            &
      finfmtuv2d,filenameu2d,filenamev2d,iastride,jastride,             &
      vtraplt, vtraunit,vtraovr,vtracol1,vtracol2,                      &
      vtraprio,vtraunits,vtratype,                                      &
      magaplt,magainc,magaminc,magamaxc,magaovr,magahlf,magazro,        &
      magasty,magacol1,magacol2,magaunits

  INTEGER :: zxplot_called,zxout_unit
  REAL :: xorig_0, yorig_0

  REAL :: omega2, r2d
  REAL :: relvort

  REAL, ALLOCATABLE :: fcorio(:,:), mapfct(:,:)
  REAL :: umove_readin, vmove_readin,uadd,vadd
!
!-----------------------------------------------------------------------
!
! Variables for mpi jobs
!
!-----------------------------------------------------------------------

  INTEGER :: ncompressx, ncompressy ! compression in x and y direction:
                                    ! ncompressx=nprocx_in/nproc_x
                                    ! ncompressy=nprocy_in/nproc_y
  INTEGER :: nproc_node
  INTEGER :: nprocx_lw, nprocy_lw
  COMMON /init1_mpi/ ncompressx, ncompressy, nproc_node,nprocx_lw,nprocy_lw

  INTEGER :: nxlg, nylg             ! global domain

  INTEGER :: xpbgn, xpend, ypbgn, ypend
  INTEGER :: xpbgn1,xpend1,ypbgn1,ypend1
  INTEGER :: ibgnl,iendl,jbgnl,jendl
  INTEGER :: ibgnl1, iendl1
  INTEGER :: ibgnl2, jbgnl2      ! hold the value of ibgnl and jbgnl
                                 ! for vector plot
  INTEGER :: ibgnl21
  INTEGER :: ibgnla, jbgnla      ! hold the value of ibgnl and jbgnl
                                 ! for arbitrary 2D vector plot

  COMMON /processors/ xpbgn,xpend,ypbgn,ypend

  INTEGER :: agl
  COMMON /agl_or_asl/ agl

!
!-----------------------------------------------------------------------
!
!  Function f_qvsat and inline directives for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_qvsat

  REAL :: amin, amax, tem, pref
  REAL :: xr2, yr2
  REAL :: xor_current, yor_current

  REAL :: xshift, yshift

  INTEGER :: nch
  COMMON /XOUTCH/ nch

  INTEGER :: itrajc_start,itrajc_end
  INTEGER :: jpoint, jpoint_final
  INTEGER :: nq, posdefflag, id, jd, kd

  INTEGER :: P_GH

!fpp$ expand (f_qvsat)
!!dir$ inline always f_qvsat
!*$*  inline routine (f_qvsat)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  cpu0 = f_cputime()

  xi1= 50.0
  yi1= 50.0
  xi2= 50.0
  yi2= 50.0
  trcolor = 183   ! xz slice terrain color
  plotovr = .false.
  sovrlay = 0
  lbcolor = 1     ! laber color 1-black
!  imsg = 1        ! Activate missing value checking for streamline routine
!  uvmsg = -9999.0 ! Set the missing data value for streamline routine
!  inita = 8
!  initb = 8

  layover = 0
  ctinc = 0.0
  vtunt = 0.0
  aspratio = 1.0

  xnwpic_called =0

  first_frame = 1
!
!-----------------------------------------------------------------------
!
!  Set the max. and min. values for producing HDF images.
!
!-----------------------------------------------------------------------
!
  utmax =  40.0
  utmin = -40.0
  vtmax =  40.0
  vtmin = -40.0
  wtmax =  40.0
  wtmin = -40.0
  vsmax =  40.0
  vsmin = -40.0
  qvmax =  15.0
  qvmin =   0.0
!  qcmax =   5.0
!  qcmin =   0.0
!  qrmax =  10.0
!  qrmin =   0.0
!  qimax =  10.0
!  qimin =   0.0
!  qsmax =  10.0
!  qsmin =   0.0
!  qhmax =  10.0
!  qhmin =   0.0
  qwmax =  10.0
  qwmin =   0.0
  rhmax =   1.0
  rhmin =   0.0
  rhimax =   1.0
  rhimin =   0.0
  rfmax =  75.0
  rfmin =   0.0
  ptemax= 400.0
  ptemin= 250.0
  upmax =  40.0
  upmin = -40.0
  vpmax =  40.0
  vpmin = -40.0
  wpmax =  5.0
  wpmin = -5.0
  ptpmax=  0.0
  ptpmin= -8.0
  ppmax = 200.0
  ppmin =-200.0
  qvpmax=  10.0
  qvpmin=   0.0
  vormax=  0.02
  vormin= -0.02
  avormax=  0.02
  avormin= -0.02
  divmax=  0.02
  divmin= -0.02
  hmax =  2000.0
  hmin =    0.0
  tdmin = 200.
  tdmax = 400.

  time1c = 0.0
  time2c = 0.0
  time1g = 0.0
  time2g = 0.0
  time1t = 0.0
  time2t = 0.0
  timediff = 0.0

  iwtype = 1

  CALL initpltpara(nx,ny,nz,nzsoil,nstyps,outfilename,iwtype,istatus)

!-----------------------------------------------------------------------
!
!  Set some parameters
!
!-----------------------------------------------------------------------
  trcolor = trncol1

  nxlg = (nx-3)*nproc_x + 3
  nylg = (ny-3)*nproc_y + 3

  ipriority(1)=hprio
  ipriority(2)=tprio
  ipriority(3)=uprio
  ipriority(4)=vprio
  ipriority(5)=vhprio
  ipriority(6)=wprio
  ipriority(7)=ptprio
  ipriority(8)=pprio
  ipriority(9)=qvprio

  !DTD: priority numbers 10-40 reserved for qscalar array


  DO nq=1,nscalar
    ipriority(nq+9)=qscalarprio(nq)
  END DO

!  ipriority(10)=qcprio
!  ipriority(11)=qrprio
!  ipriority(12)=qiprio
!  ipriority(13)=qsprio
!  ipriority(14)=qhprio


  ipriority(41)=qwprio
  ipriority(42)=kmhprio
  ipriority(43)=kmvprio
  ipriority(44)=tkeprio
  ipriority(45)=rhprio
  ipriority(46)=tdprio
  ipriority(47)=rfprio
  ipriority(48)=pteprio
  ipriority(49)=upprio
  ipriority(50)=vpprio
  ipriority(51)=wpprio
  ipriority(52)=ptpprio
  ipriority(53)=ppprio
  ipriority(54)=qvpprio
  ipriority(55)=vorpprio
  ipriority(56)=divpprio
  ipriority(57)=divqprio

  ipriority(58)=vtrprio
  ipriority(59)=vtpprio

  ipriority(60)=vtrstmprio
  ipriority(61)=vtpstmprio

  ipriority(62)=rfcprio
  ipriority(63)=vsprio
  ipriority(64)=gricprio
  ipriority(65)=avorprio
  ipriority(66)=xuvprio
  ipriority(67)=qtprio
  ipriority(68)=rhiprio
  ipriority(69)=tsoilprio
  ipriority(70)=qsoilprio
  ipriority(71)=zdrprio
  ipriority(72)=kdpprio
  ipriority(73)=zdpprio
  ipriority(74)=rhvprio

  ipriority(81)=trnprio
  ipriority(82)=wcpprio
  ipriority(83)=racprio
  ipriority(84)=ragprio
  ipriority(85)=ratprio
  ipriority(86)=pslprio
  ipriority(87)=liprio
  ipriority(88)=capprio
  ipriority(89)=cinprio
  ipriority(90)=theprio
  ipriority(91)=helprio
  ipriority(92)=uhprio
  ipriority(93)=ctcprio
  ipriority(94)=brnprio
  ipriority(95)=bruprio
  ipriority(96)=srlprio
  ipriority(97)=srmprio
  ipriority(98)=capsprio
  ipriority(99)=blcoprio
  ipriority(100)=viqcprio
  ipriority(101)=viqrprio
  ipriority(102)=viqiprio
  ipriority(103)=viqsprio
  ipriority(104)=viqhprio
  ipriority(105)=vilprio
  ipriority(106)=viiprio
  ipriority(107)=vicprio
  ipriority(108)=strmprio
  ipriority(109)=vitprio
  ipriority(110)=pwprio
  ipriority(111)=tprprio
  ipriority(112)=gprprio
  ipriority(113)=cprprio

  ipriority(121)=styprio
  ipriority(122)=vtyprio
  ipriority(123)=laiprio
  ipriority(124)=rouprio
  ipriority(125)=vegprio
  ipriority(126)=snowdprio
  ipriority(127)=msfprio
  ipriority(128)=ipvprio
  ipriority(129)=vagprio
  ipriority(130)=thkprio

  ipriority(131)=raicprio
  ipriority(132)=raigprio
  ipriority(133)=raitprio

  DO i=1, var3dnum
    ipriority(150+i) = var3dprio(i)
  END DO

  DO i=1, var2dnum
    ipriority(170+i) = var2dprio(i)
  END DO

  nslice(1)=nslice_xy
  nslice(2)=nslice_xz
  nslice(3)=nslice_yz

  IF( nslice_h >= 0 ) THEN
    nslice(4)=nslice_h
  ELSE
    nslice(4)=0
  ENDIF

  nslice(5)=nslice_v
  nslice(6)=nslice_p
  nslice(7)=nslice_pt
  nslice(8)= 1             ! surface plots
  iplot(8) = 1


  iplot(1:7)=1
  IF(nslice(1) == 0) iplot(1)=0
  IF(nslice(2) == 0) iplot(2)=0
  IF(nslice(3) == 0) iplot(3)=0
  IF(nslice(4) == 0) iplot(4)=0
  IF(nslice(5) == 0) iplot(5)=0
  IF(nslice(6) == 0) iplot(6)=0
  IF(nslice(7) == 0) iplot(7)=0

  nslice(9) =nslice_xy_soil
  nslice(10)=nslice_xz_soil
  nslice(11)=nslice_yz_soil

  iplot(9:11)=1
  IF(nslice(9)  == 0) iplot(9)=0
  IF(nslice(10) == 0) iplot(10)=0
  IF(nslice(11) == 0) iplot(11)=0

  DO indxslic = 1,nslice(5)
    xi1(indxslic)=xpnt1(indxslic)
    yi1(indxslic)=ypnt1(indxslic)
    xi2(indxslic)=xpnt2(indxslic)
    yi2(indxslic)=ypnt2(indxslic)
  END DO

  IF(arbvaropt < 10) THEN
    rdallhist=.true.
  ELSE
    arbvaropt=mod(arbvaropt,10)
    rdallhist=.false.
  END IF

  umove_readin = umove
  vmove_readin = vmove
!
!-----------------------------------------------------------------------
!
!  When priority equal to zero , the order of the plotting is default.
!
!-----------------------------------------------------------------------
!
  DO i=1,nprio
    iptemp(i)=i
    sigplt(i)=0
  END DO
  iptemp(nprio)=0

!
!-----------------------------------------------------------------------
!
!  Allocate arrays.
!
!-----------------------------------------------------------------------
!
  nzmax = max(nzsoil,nz)

  ALLOCATE(x     (nx),STAT=istatus)
  ALLOCATE(y     (ny),STAT=istatus)
  ALLOCATE(z     (nz),STAT=istatus)
  ALLOCATE(zp    (nx,ny,nz),STAT=istatus)
  ALLOCATE(zpsoil(nx,ny,nzsoil),STAT=istatus)

  ALLOCATE(u     (nx,ny,nz),STAT=istatus)
  ALLOCATE(v     (nx,ny,nz),STAT=istatus)
  ALLOCATE(w     (nx,ny,nz),STAT=istatus)
  ALLOCATE(ptprt (nx,ny,nz),STAT=istatus)
  ALLOCATE(pprt  (nx,ny,nz),STAT=istatus)
  ALLOCATE(qv    (nx,ny,nz),STAT=istatus)
  ALLOCATE(qscalar(nx,ny,nz,nscalar),STAT=istatus)
  ALLOCATE(tke   (nx,ny,nz),STAT=istatus)
  ALLOCATE(kmh   (nx,ny,nz),STAT=istatus)
  ALLOCATE(kmv   (nx,ny,nz),STAT=istatus)

  ALLOCATE(ubar  (nx,ny,nz),STAT=istatus)
  ALLOCATE(vbar  (nx,ny,nz),STAT=istatus)
  ALLOCATE(wbar  (nx,ny,nz),STAT=istatus)
  ALLOCATE(ptbar (nx,ny,nz),STAT=istatus)
  ALLOCATE(pbar  (nx,ny,nz),STAT=istatus)
  ALLOCATE(rhobar(nx,ny,nz),STAT=istatus)
  ALLOCATE(qvbar (nx,ny,nz),STAT=istatus)

  ALLOCATE(soiltyp (nx,ny,nstyps),STAT=istatus)
  ALLOCATE(stypfrct(nx,ny,nstyps),STAT=istatus)
  ALLOCATE(vegtyp  (nx,ny),STAT=istatus)
  ALLOCATE(lai     (nx,ny),STAT=istatus)
  ALLOCATE(roufns  (nx,ny),STAT=istatus)
  ALLOCATE(veg     (nx,ny),STAT=istatus)

  ALLOCATE(tsoil   (nx,ny,nzsoil,0:nstyps),STAT=istatus)
  ALLOCATE(qsoil   (nx,ny,nzsoil,0:nstyps),STAT=istatus)
  ALLOCATE(wetcanp (nx,ny,0:nstyps),STAT=istatus)
  ALLOCATE(snowdpth(nx,ny),STAT=istatus)

  ALLOCATE(raing   (nx,ny),STAT=istatus)
  ALLOCATE(rainc   (nx,ny),STAT=istatus)
  ALLOCATE(raint   (nx,ny),STAT=istatus)
  ALLOCATE(prcrate (nx,ny,4),STAT=istatus)

  ALLOCATE(radfrc(nx,ny,nz),STAT=istatus)

  ALLOCATE(radsw (nx,ny),   STAT=istatus)
  ALLOCATE(rnflx (nx,ny),   STAT=istatus)
  ALLOCATE(radswnet(nx,ny), STAT=istatus)
  ALLOCATE(radlwin (nx,ny), STAT=istatus)
  ALLOCATE(usflx (nx,ny),   STAT=istatus)
  ALLOCATE(vsflx (nx,ny),   STAT=istatus)
  ALLOCATE(ptsflx(nx,ny),   STAT=istatus)
  ALLOCATE(qvsflx(nx,ny),   STAT=istatus)

  ALLOCATE(uprt  (nx,ny,nz),STAT=istatus)
  ALLOCATE(vprt  (nx,ny,nz),STAT=istatus)
  ALLOCATE(wprt  (nx,ny,nz),STAT=istatus)
  ALLOCATE(pt    (nx,ny,nz),STAT=istatus)
  ALLOCATE(qvprt (nx,ny,nz),STAT=istatus)

  ALLOCATE(xc    (nx,ny,nzmax),STAT=istatus)   ! shared with soil
  ALLOCATE(yc    (nx,ny,nzmax),STAT=istatus)   ! shared with soil
  ALLOCATE(zc    (nx,ny,nz),            STAT=istatus)
  ALLOCATE(zpc   (nx,ny,nz),            STAT=istatus)
  ALLOCATE(zpsoilc(nx,ny,nzsoil),       STAT=istatus)   !Zuwen

  ALLOCATE(acc_rainc(nx,ny,2),STAT=istatus)
  ALLOCATE(acc_raing(nx,ny,2),STAT=istatus)
  ALLOCATE(acc_raint(nx,ny,2),STAT=istatus)

  ALLOCATE(hterain(nx,ny),    STAT=istatus)

  ALLOCATE(psl   (nx,ny),     STAT=istatus)
  ALLOCATE(td    (nx,ny,nz),  STAT=istatus)

  ALLOCATE(tz    (nx,ny,nz),  STAT=istatus)
  ALLOCATE(t700  (nx,ny),     STAT=istatus)

  ALLOCATE(algpzc(nx,ny,nz),  STAT=istatus)
  ALLOCATE(zps3d (nx,ny,nz),  STAT=istatus)
  ALLOCATE(zpsoils3d (nx,ny,nzsoil),STAT=istatus)

  ALLOCATE(li   (nx,ny),STAT=istatus)
  ALLOCATE(cape (nx,ny),STAT=istatus)
  ALLOCATE(cin  (nx,ny),STAT=istatus)
  ALLOCATE(thet (nx,ny),STAT=istatus)
  ALLOCATE(heli (nx,ny),STAT=istatus)
  ALLOCATE(uh   (nx,ny),STAT=istatus)
  ALLOCATE(brn  (nx,ny),STAT=istatus)
  ALLOCATE(brnu (nx,ny),STAT=istatus)
  ALLOCATE(srlfl(nx,ny),STAT=istatus)
  ALLOCATE(srmfl(nx,ny),STAT=istatus)
  ALLOCATE(ustrm(nx,ny),STAT=istatus)
  ALLOCATE(vstrm(nx,ny),STAT=istatus)

  ALLOCATE(capst(nx,ny),STAT=istatus)
  ALLOCATE(blcon(nx,ny),STAT=istatus)
  ALLOCATE(ct   (nx,ny),STAT=istatus)
  ALLOCATE(sinlat(nx,ny),STAT=istatus)

  ALLOCATE(xs(nx),STAT=istatus)
  ALLOCATE(ys(ny),STAT=istatus)

  ALLOCATE(tem1(nx,ny,nzmax),STAT=istatus)
  ALLOCATE(tem2(nx,ny,nzmax),STAT=istatus)
  ALLOCATE(tem3(nx,ny,nzmax),STAT=istatus)
  ALLOCATE(tem4(nx,ny,nzmax),STAT=istatus)
  ALLOCATE(tem5(nx,ny,nzmax),STAT=istatus)
  ALLOCATE(tem6(nx,ny,nzmax),STAT=istatus)
  ALLOCATE(tem7(nx,ny,nzmax),STAT=istatus)
  ALLOCATE(tem8(nx,ny,nzmax),STAT=istatus)
  ALLOCATE(tem9(nx,ny,nzmax),STAT=istatus)

  ALLOCATE(b1(nx,ny),STAT=istatus)
  ALLOCATE(b2(nx,ny),STAT=istatus)
  ALLOCATE(u1(nx,ny),STAT=istatus)
  ALLOCATE(v1(nx,ny),STAT=istatus)

  ALLOCATE(u2 (nx+ny,nz),STAT=istatus)
  ALLOCATE(v2 (nx+ny,nz),STAT=istatus)
  ALLOCATE(w2 (nx+ny,nz),STAT=istatus)
  ALLOCATE(zs2(nx+ny,nz),STAT=istatus)
  ALLOCATE(xp (nx+ny),STAT=istatus)
  ALLOCATE(yp (nx+ny),STAT=istatus)

  ALLOCATE(var3dv(nx,ny,nz),STAT=istatus)
  ALLOCATE(var2du(nx,ny),STAT=istatus)
  ALLOCATE(var2dv(nx,ny),STAT=istatus)
  ALLOCATE(fcorio(nx,ny),STAT=istatus)
  ALLOCATE(mapfct(nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "arpsplt:mapfct")

  x     =0.0
  y     =0.0
  z     =0.0
  zp    =0.0
  zpsoil=0.0

  u     =0.0
  v     =0.0
  w     =0.0
  ptprt =0.0
  pprt  =0.0
  qv    =0.0
  qscalar=0.0
  tke   =0.0
  kmh   =0.0
  kmv   =0.0

  ubar  =0.0
  vbar  =0.0
  wbar  =0.0
  ptbar =0.0
  pbar  =0.0
  rhobar=0.0
  qvbar =0.0

  soiltyp =0.0
  stypfrct=0.0
  vegtyp  =0.0
  lai     =0.0
  roufns  =0.0
  veg     =0.0

  tsoil   =0.0
  qsoil   =0.0
  wetcanp =0.0
  snowdpth=0.0

  raing   =0.0
  rainc   =0.0
  raint   =0.0
  prcrate =0.0

  radfrc=0.0

  radsw =0.0
  rnflx =0.0
  usflx =0.0
  vsflx =0.0
  ptsflx=0.0
  qvsflx=0.0

  uprt  =0.0
  vprt  =0.0
  wprt  =0.0
  pt    =0.0
  qvprt =0.0

  xc    =0.0
  yc    =0.0
  zc    =0.0
  zpc   =0.0
  zpsoilc   =0.0

  acc_rainc = 0.0
  acc_raing = 0.0
  acc_raint = 0.0

  hterain=0.0

  psl   =0.0
  td    =0.0

  tz    =0.0
  t700  =0.0

  algpzc=0.0
  zps3d =0.0
  zpsoils3d =0.0

  li   =0.0
  cape =0.0
  cin  =0.0
  thet =0.0
  heli =0.0
  uh   =0.0
  brn  =0.0
  brnu =0.0
  srlfl=0.0
  srmfl=0.0
  ustrm=0.0
  vstrm=0.0

  capst=0.0
  blcon=0.0
  ct   =0.0
  sinlat=0.0

  xs=0.0
  ys=0.0

  tem1=0.0
  tem2=0.0
  tem3=0.0
  tem4=0.0
  tem5=0.0
  tem6=0.0
  tem7=0.0
  tem8=0.0
  tem9=0.0

  b1=0.0
  b2=0.0
  u1=0.0
  v1=0.0

  u2 =0.0
  v2 =0.0
  w2 =0.0
  zs2=0.0
  xp =0.0
  yp =0.0

  var3dv=0.0
  var2du=0.0
  var2dv=0.0
  fcorio=0.0

!-----------------------------------------------------------------------
!  Calculate alpha for DSD to calculte polarimetric variables
!  Added by Youngsun Jung
!-----------------------------------------------------------------------
  ALLOCATE(alpha(nx,ny,nz,6))
  alpha = 0.0

  IF(myproc == 0) THEN
!-----------------------------------------------------------------------
!
!  Overlaying observations
!
!-----------------------------------------------------------------------

    IF(ovrobs == 1) THEN
      CALL getunit(obunit)
      iob=1
      nobs=0
      DO ifile=1,nsfcobfl
        lsfcobfl=LEN_TRIM(sfcobfl(ifile))
        INQUIRE(FILE=TRIM(sfcobfl(ifile)), EXIST = fexist )
        IF( .NOT.fexist) THEN
          WRITE(6,'(3a)') 'WARNING: Surface obs file ',                &
            TRIM(sfcobfl(ifile)),                                      &
            ' not found. ARPSPLT program will continue.'
          CYCLE
        END IF
        CALL read_surface_obs(sfcobfl(ifile),mxsfcob,atime,n_meso_g,      &
          n_meso_pos,n_sao_g,n_sao_pos_g,n_sao_b,n_sao_pos_b,n_obs_g,     &
          n_obs_pos_g,n_obs_b,n_obs_pos_b,                                &
          stn(iob),obstype(iob),latob(iob),lonob(iob),elevob(iob),wx(iob),&
          tob(iob),tdob(iob),dd(iob),ff(iob),ddg(iob),ffg(iob),pstn(iob), &
          pmsl(iob),alt(iob),kloud(iob),ceil(iob),lowcld(iob),cover(iob), &
          rad(iob),idp3(iob),store_emv(iob,1),store_amt(iob,1),           &
          store_hgt(iob,1),vis(iob),obstime(iob),istatus)

        nobs=nobs+n_obs_b
        WRITE(6,'(2x,a,i8)') ' N sfc obs: ',nobs
        iob=iob+n_obs_b
      END DO

    END IF

!-----------------------------------------------------------------------
!
!  Overlaying airport location
!
!-----------------------------------------------------------------------
!
    IF(ovrstaopt /= 0) THEN
      lstalofl=LEN_TRIM(stalofl)

      INQUIRE(FILE=TRIM(stalofl), EXIST = fexist )
      IF( .NOT.fexist) THEN
        WRITE(6,'(3a)') 'WARNING: Station location file ',             &
          TRIM(stalofl),' not found. ARPSPLT program will continue.'
      ELSE
        CALL getunit(stunit)
        CALL read_station(stalofl,mxstalo,latsta,lonsta,nstatyp,       &
                    nstapro,nsta,s_name,s_state,s_site,s_elev)
        IF(nsta > 0) staset=1
      END IF
    END IF

  END IF   ! myproc == 0
  CALL mpupdatei(n_obs_b,1)
  CALL mpupdatei(staset,1)

!-----------------------------------------------------------------------
!
! Plot beginning below
!
!-----------------------------------------------------------------------

  zxplot_called=0

  DO nf=1, nhisfile

    sigplt(:) = 0

    lenfil = len_trim(hisfile(nf))

    lengbf = len_trim(grdbasfn)

    second1= f_walltime()
    cpu1   = f_cputime()

    IF (rdallhist) THEN
      CALL readarpsmp(ncompressx,ncompressy,nproc_node,                 &
              nprocx_lw,nprocy_lw,hinfmt,                               &
              grdbasfn,lengbf,hisfile(nf),lenfil,nx,ny,nz,nzsoil,nstyps,&
              time,x,y,z,zp,zpsoil,uprt,vprt,wprt,ptprt,pprt,qvprt,     &
              qscalar,tke,kmh,kmv,                                      &
              ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,                   &
              soiltyp,stypfrct,vegtyp,lai,roufns,veg,                   &
              tsoil,qsoil,wetcanp,snowdpth,                             &
              raing,rainc,prcrate,radfrc,radsw,rnflx,radswnet,radlwin,  &
              usflx,vsflx,ptsflx,qvsflx,istatus,tem1,tem2,tem3)
    ELSE
      IF (nf == 1) CALL get_gridxyzzp(nx,ny,nz,grdbasfn,hinfmt,         &
              ncompressx*nproc_x, ncompressy*nproc_y,x,y,z,zp,istatus)
      !
      ! Get file time
      !
      idxfmt = INDEX(hisfile(nf),fmtstr(hinfmt),.TRUE.)  ! backward

      IF(idxfmt > 0) THEN
        READ(hisfile(nf)((idxfmt+3):(idxfmt+8)),'(i6)') filetm
        time=FLOAT(filetm)
        curtim=time
      END IF

      filetm = INDEX(hisfile(nf),'/',.TRUE.)

      idxfmt = INDEX(hisfile(nf),fmtstr(hinfmt))    ! forward
      IF (filetm > idxfmt-2) THEN
        WRITE(6,'(/,1x,3a,/)') 'ERROR: File name ',hisfile(nf),         &
           '. Please check namelist block &history_data and try again.'
        CALL arpsstop('ERROR: file name.',1)
      END IF
      runname = hisfile(nf)(filetm+1:idxfmt-2)    ! decode runname from file name
      lfnkey  = idxfmt-2 - filetm
    END IF

    IF (nf == 1) THEN
      time1c = time
      time1g = time
      time1t = time
    END IF
    time2c = time
    time2g = time
    time2t = time

    cpu2    = f_cputime()
    second2 = f_walltime()

    !print*,'!!!!  total clock time for read data:',second2-second1
    !print*,'!!!!  total cpu time for read data  :', cpu2-cpu1

    IF(zxplot_called == 0 .AND. myproc == 0) THEN
!
!-----------------------------------------------------------------------
!
!  Initialize ZXPLOT plotting package
!
!-----------------------------------------------------------------------
!
! The following two ZXPLOT routines, if called, should be called before xdevic
!
      zxout_unit=2

      IF (LEN_TRIM(outfilename) < 1) THEN
        WRITE(outfilename,'(a)') TRIM(dirname)//runname(1:lfnkey)
      ELSE
        WRITE(outfilename,'(a)') TRIM(dirname)//TRIM(outfilename)
      END IF
!       CALL xpsfn(TRIM(outfilename)//'.ps', zxout_unit)
      CALL xpaprlnth( paprlnth )

      exchange = iwtype/100
      iwtype = MOD(iwtype,100)

      CALL xdevic_new(iwtype,outfilename,LEN_TRIM(outfilename),zxout_unit)

      CALL xstctfn(color_map)

      CALL xsetclrs_new(col_table,exchange)

!       CALL xafstyl(1)
      CALL xcolor(lbcolor)

      CALL xartyp(2)

      CALL xlnmag(lnmag)
      CALL xcfont(fontopt)
      CALL xlabmask(lbmaskopt)

      IF(ctrlbopt == 0) CALL xcltyp(0)
      IF(ctrlbopt == 1) CALL xclfmt('(f10.1)')
      IF(ctrlbopt == 2) CALL xclfmt('(I10)')

      CALL xcmixl
      IF( ctrstyle == 2) THEN
        CALL xcfull
      ELSE IF( ctrstyle == 3) THEN
        CALL xcdash
      END IF

      CALL xclfrq( ctrlbfrq )

      CALL xctrbadv( 1 )  ! Turn on missing value checking for contouring
      CALL xvtrbadv( 1 )  ! Turn on missing value checking for vector plotting
!
!  Default value of -9999.0 for the flag is used.
!
      CALL xbadval ( -9999.0 ) ! Set the missing value flag to -9999.0
!
!  above three routines were not available with older version of ZXPLOT.
!
      CALL xdspac(0.9*winsiz)

      IF( nxpic == 1 .AND. nypic == 1) CALL xdspac( 0.85*winsiz)
      IF( layout == 1 ) THEN
        angl = 00.0
      ELSE
        angl = 90.0
      END IF
      CALL xspace(nxpic, nypic, angl , xlimit,ylimit)
      IF(axlbfmt == (-1)) THEN
        CALL xaxfmt( '*' )
      ELSE IF(axlbfmt == 0) THEN
        CALL xaxfmt( '(i5)' )
      ELSE
        WRITE(stem1,'(i1)')     axlbfmt
        WRITE(stem2,'(a,a1,a)') '(f8.',stem1,')'
        CALL xaxfmt( stem2 )
      END IF

      layover = 0
      CALL xpmagn(margnx*xlimit, margny*ylimit)

    END IF  ! Initalization of plotting package and parameters
!
!   Somewhere in one of the above calls, NCH gets set.  Pass it around.
!
    CALL mpupdatei(nch,1)

    zxplot_called=1

    IF(myproc == 0) WRITE(6,'(/a,F8.0,a)')                              &
         ' ------ Processing file time - ',curtim,' ------'
!
!-----------------------------------------------------------------------
!
!  If islice=-2, plot the y-z cross-sections through w maximum.
!  If jslice=-2, plot the x-z cross-sections through w maximum.
!  If kslice=-2, plot the x-y cross-sections through w maximum.
!
!-----------------------------------------------------------------------
!
    CALL a3dmax(wprt,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,2,nz-1,             &
                wmax,wmin, imax,jmax,kmax, imin,jmin,kmin)

    CALL mpupdatei(imax,1)   ! only processor 0 contains the right indices
    CALL mpupdatei(jmax,1)
    CALL mpupdatei(kmax,1)

    CALL mpupdatei(imin,1)
    CALL mpupdatei(jmin,1)
    CALL mpupdatei(kmin,1)

!  write(6,'(/2(1x,a,f9.3,3(a,i3))/)')
!    :     'wmin =',wmin,' at i=',imin,', j=',jmin,', k=',kmin,
!    :     'wmax =',wmax,' at i=',imax,', j=',jmax,', k=',kmax

    CALL set_coord(iorig,xgrdorg,ygrdorg,xorig,yorig,nx,ny,nz,nzsoil,nzmax, &
                   x,y,z,zp,zpsoil,xc,yc,zc,zpc,zpsoilc,                &
                   dxinv,dyinv,dzinv,dxkm,dykm,dzkm,dzsoilcm,           &
                   xbgn,xend,ybgn,yend,zbgn,zend,zsoilbgn,zsoilend,     &
                   iskip,jskip,                                         &
                   ibgn,iend,jbgn,jend,kbgn,kend,ksoilbgn,ksoilend,     &
                   xpbgn,xpend,ypbgn,ypend,ibgnl,iendl,jbgnl,jendl,     &
                   istride,jstride,kstride,ist,jst,kst,ibgnl2,jbgnl2,   &
                   xr, yr, x1,y1,x2,y2, zr,z1,z2, zsoilr,zsoil1,zsoil2, &
                   zmin,zmax,zsoilmin,zsoilmax,lvldbg,istatus)

!
! Amount of coordinate shift for trajectories
!
    xshift = 0.0
    yshift = 0.0
    IF( iorig == 1) THEN
      xshift = xgrdorg
      yshift = ygrdorg
    ELSE IF( iorig == 2) THEN
      xshift = -(x(2)+xgrdorg) + xorig*1000.0
      yshift = -(y(2)+ygrdorg) + yorig*1000.0
    ELSE IF( iorig == 3) THEN
      xshift = + xorig*1000.0
      yshift = + yorig*1000.0
    ENDIF

    xpbgn1 = xpbgn
    xpend1 = xpend
    ypbgn1 = ypbgn
    ypend1 = ypend

    ibgnl1 = ibgnl
    iendl1 = iendl

    ibgnl21 = ibgnl2

    ibgn1 = ibgn
    iend1 = iend
    xr2 = xr                  ! save the window to be reseted later.
    yr2 = yr
    x11 = x1
    y11 = y1

!-----------------------------------------------------------------------
!
! Calculate ARPS derived variables
!
!-----------------------------------------------------------------------

    IF (rdallhist) THEN
    CALL cal_arpsderived(imove,umove,vmove,umove_readin,vmove_readin,   &
              presaxis_no,nx,ny,nz,zp,zpc,x,y,mapfct,                   &
              ubar,vbar,wbar,qvbar,ptbar,pbar,                          &
              uprt,vprt,wprt,qvprt,ptprt,pprt,                          &
              u,v,w,qv,pt,algpzc,tz,t700,zlevel,psl,td,thet,            &
              capeplt,cinplt,capsplt,liplt,brnplt,ctcplt,               &
              li,cape,cin,capst,ct,                                     &
              heliplt,uhplt,brnuplt,srlfplt,srmfplt,blcoplt,strmplt,    &
              uhmnhgt,uhmxhgt,                                          &
              ustrm,vstrm,srlfl,srmfl,heli,uh,brn,brnu,blcon,           &
              ibgnl,iendl,jbgnl,jendl,kbgn,kend,                        &
              hterain,ztmax,ztmin,u1,v1,                                &
              tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,tem9,             &
              istatus)
    END IF

!
!-----------------------------------------------------------------------
!
!  Re-Initialize the plotting space parameters, which might have
!  been set differently for the vertical profile plotting.
!
!-----------------------------------------------------------------------
!
    IF (profopt == 1)  THEN

      CALL xdspac(0.9*winsiz)

      IF( nxpic == 1 .AND. nypic == 1) CALL xdspac(0.85*winsiz)

      CALL xspace(nxpic, nypic, angl , xlimit,ylimit)
      IF(axlbfmt == -1) THEN
        CALL xaxfmt( '*' )
      ELSE IF(axlbfmt == 0) THEN
        CALL xaxfmt( '(i5)' )
      ELSE
        WRITE(stem1,'(i1)')axlbfmt
        WRITE(stem2,'(a,a1,a)') '(f8.',stem1,')'
        CALL xaxfmt( stem2 )
      END IF

      layover = 0
      CALL xpmagn(margnx*xlimit, margny*ylimit)

    END IF
!
!-----------------------------------------------------------------------
!
!  Set map
!
!-----------------------------------------------------------------------
!
!  IF( ovrmap .eq. 1 .or. ovrobs .eq.1 .or. blcoplt.ne.0 .or.
!    :    ovrstaopt.ne.0 ) THEN

    xl = (nxlg-3)*dxkm
    yl = (nylg-3)*dykm

    IF (myproc == 0) THEN
      xorig_0 = (xc(1,2,2)+xc(2,2,2))*0.5
      yorig_0 = (yc(2,1,2)+yc(2,2,2))*0.5
!    END IF
!    CALL mpupdater(xorig_0,1)
!    CALL mpupdater(yorig_0,1)

      CALL xstpjgrd(mapproj,trulat1,trulat2,trulon,                     &
                    ctrlat,ctrlon,xl,yl,xorig_0,yorig_0)
    END IF
!  ENDIF
!
!-----------------------------------------------------------------------
!
! Calculate Coriolis parameter and map factor
!
!-----------------------------------------------------------------------
!
    DO i=1,nx-1
      xs(i) = (x(i)+x(i+1))*0.5
    END DO
    xs(nx) = -xs(nx-2)+2.0*xs(nx-1)

    DO j=1,ny-1
      ys(j) = (y(j)+y(j+1))*0.5
    END DO
    ys(ny) = -ys(ny-2)+2.0*ys(ny-1)

    CALL xxytoll(nx,ny,xs,ys,tem8(1,1,1),tem8(1,1,2))

    r2d = ATAN(1.0)*4.0/180.0
    omega2 = 2.0* omega

    DO j=1,ny-1
      DO i=1,nx-1
        fcorio(i,j) = omega2*SIN( r2d * tem8(i,j,1) )
      END DO
    END DO

    DO i=1,nx-1
      xp(i) = 0.5*(x(i)+x(i+1))
    END DO
    DO j=1,ny-1
      yp(j) = 0.5*(y(j)+y(j+1))
    END DO
    CALL xytomf(nx-1,ny-1,xp,yp,mapfct)

!-----------------------------------------------------------------------
!  Initialize microphysical varialbes.
!  Added by Youngsun Jung
!-----------------------------------------------------------------------

  mphyopt = rfopt

  IF (mphyopt >= 205) THEN  ! Map COAMPS mphysics to ARPS option number
    mphyopt = mphyopt - 200 + 3
  END IF

  pi = 4.0 * atan(1.0)

  cx(1) = (pi/6.)*rhow
  cx(2) = (pi/6.)*rhow
  cx(3) = 440.0      ! Constant value used in MY scheme for ice
  cx(4) = (pi/6.)*rhosnow
  cx(5) = (pi/6.)*rhogrpl
  cx(6) = (pi/6.)*rhohail

! JYS added
  IF(dualpol == 2) THEN
    SELECT CASE (mphyopt)
    CASE (2:4)
      alpha(:,:,:,2) = alpharain
      alpha(:,:,:,3) = alphaice
      alpha(:,:,:,4) = alphasnow
      alpha(:,:,:,5) = alphahail
      alpha(:,:,:,6) = alphahail
    CASE (5:10,12,205:207,209)
      alpha(:,:,:,2) = alpharain
      alpha(:,:,:,3) = alphaice
      alpha(:,:,:,4) = alphasnow
      alpha(:,:,:,5) = alphagrpl
      alpha(:,:,:,6) = alphahail
    CASE (11,208)
      DO nq=2,6
!        CALL solve_alpha(nx,ny,nz,rhobar,cx(nq),              &
!                 qscalar(:,:,:,nq),qscalar(:,:,:,nq+6),       &
!                 qscalar(:,:,:,nq+11),alpha(:,:,:,nq))
      ENDDO
    END SELECT
  ENDIF
  IF (myproc == 0) THEN
    WRITE(6,'(a)')'Model read-in parameters'
    WRITE(6,'(5x,a,e15.5)')'n0rain = ', n0rain
    WRITE(6,'(5x,a,e15.5)')'n0snow = ', n0snow
    WRITE(6,'(5x,a,e15.5)')'n0grpl = ', n0grpl
    WRITE(6,'(5x,a,e15.5)')'n0hail = ', n0hail
    WRITE(6,'(5x,a,f15.5)')'alpharain = ', alpharain
    WRITE(6,'(5x,a,f15.5)')'alphasnow = ', alphasnow
    WRITE(6,'(5x,a,f15.5)')'alphagrpl = ', alphagrpl
    WRITE(6,'(5x,a,f15.5)')'alphahail = ', alphahail
  END IF

!-----------------------------------------------------------------------
!  Read in the trajectory data and determining the starting and ending
!  points of each trajectory that are to be plotted.
!-----------------------------------------------------------------------

    h_follow_trajc=0

    IF( trajc_plt_opt /= 0 ) THEN

      IF (myproc == 0) THEN
        CALL read_trajc(ireturn)

        IF( ireturn /= 0 ) THEN
          trajc_plt_opt = 0
          PRINT*,'Error reading trajectory data files. Trajectory plotting option turned off.'
        ENDIF

        DO k=1,ntimes
          DO i=1,ntrajcs(k)
            DO j=1,npoints_in(k)
              IF( ABS(xtrajc(j,i,k)+99999.0) > 0.001 .AND. (ytrajc(j,i,k)+99999.0) > 0.001 ) THEN
                xtrajc(j,i,k) = xtrajc(j,i,k) + xshift
                ytrajc(j,i,k) = ytrajc(j,i,k) + yshift
              END IF
            ENDDO
          ENDDO
        ENDDO

        CALL select_trajc_points(x1,x2,y1,y2,time)

        k=1  ! support ntimes=1 only
        itrajc_start = min(ntrajc_start, ntrajcs(k))
        itrajc_end   = min(ntrajc_end  , ntrajcs(k))
        IF( ntrajc_end < 0 ) itrajc_end = ntrajcs(k)

        IF( nslice_h < 0 ) THEN ! plotting constant-height planes following trajectory parcel.
          nslice(4)=0
          DO i=itrajc_start, itrajc_end, ntrajc_stride
            nslice(4)=nslice(4)+1
          ENDDO
          iplot(4) = 1
          h_follow_trajc=1
        ENDIF

        IF( trajc_plt_opt == 2 ) THEN
          nslice(12)=0
          DO i=itrajc_start, itrajc_end, ntrajc_stride
            nslice(12)=nslice(12)+1
          ENDDO
          iplot (12) = 1
        ELSE
          nslice(12)= 0
          iplot(12) = 0
        END IF
      END IF
      CALL mpupdatei(itrajc_start,1)
      CALL mpupdatei(itrajc_end,1)
      CALL mpupdatei(npoints_bgn,ntrajcs_max*nmax_times)
      CALL mpupdatei(npoints_end,ntrajcs_max*nmax_times)
      CALL mpupdater(ttrajc,npoints_max)
      CALL mpupdater(xtrajc,npoints_max*ntrajcs_max*nmax_times)
      CALL mpupdater(ytrajc,npoints_max*ntrajcs_max*nmax_times)
      CALL mpupdater(ztrajc,npoints_max*ntrajcs_max*nmax_times)

      CALL mpupdatei(h_follow_trajc,1)
      CALL mpupdatei(nslice,maxslicopt)
      CALL mpupdatei(iplot,maxslicopt)
    END IF
!
!-----------------------------------------------------------------------
!
!  Loop through slicopts, for different types of
!  2-d slices each time. If iplot(slicopt)=0, plotting for the
!  given slicopt is switched off.
!
!-----------------------------------------------------------------------
!
    DO slicopt = 1,maxslicopt  ! Loop over different types of slices

      first_frame = 1
!
! slicopt = 1: xy slice
! slicopt = 2: xz slice
! slicopt = 3: yz slice
! slicopt = 4: constant height
! slicopt = 5: vertical slice through two points
! slicopt = 6: constant pressure
! slicopt = 7: isentropic level
! slicopt = 8: surface xy slice
! slicopt = 9: xy-soil-slice for the soil model
! slicopt = 10: xz-soil-slice for the soil model
! slicopt = 11: yz-soil-slice for the soil model
! slicopt = 12: vertical slices following trajectories

! reset the windows
      iend = iend1
      ibgn = ibgn1
      xr = xr2
      yr = yr2

      x1 = x11
      x2 = x1+xr
      y1 = y11
      y2 = y1+yr
      z1 = zmin
      z2 = zmax
      zsoil1 = zsoilmin
      zsoil2 = zsoilmax

      xpbgn = xpbgn1
      xpend = xpend1
      ypbgn = ypbgn1
      ypend = ypend1

      ibgnl = ibgnl1
      iendl = iendl1
      ibgnl2 = ibgnl21

      SELECT CASE (slicopt)
      CASE (1, 4, 6, 7, 8, 9)
        aspect=xr/yr
        aspratio = aspect
        IF(yxstrch /= 0.0) aspratio=yxstrch
      CASE (2)
        aspect=xr/zr
        aspratio = aspect
        IF(zxstrch /= 0.0) aspratio=zxstrch
      CASE (3)
        aspect=yr/zr
        aspratio = aspect
        IF(zystrch /= 0.0) aspratio=zystrch
      CASE (5, 12)
!       aspect=
        aspratio = aspect
        IF(zhstrch /= 0.0) aspratio=zhstrch
      CASE (10)
        aspect=xr/zsoilr
        aspratio = aspect
      CASE (11)
        aspect=yr/zsoilr
        aspratio = aspect
      CASE DEFAULT
        WRITE(6,'(1x,a,I2,a)') 'ERROR: Unknown slicopt = ',slicopt,'.'
        CALL arpsstop('ERROR: Wrong option of slicopt',1)
      END SELECT

!
!-----------------------------------------------------------------------
!
!  Loop over all slices (nslice(slicopt)) for each slice type
!
!-----------------------------------------------------------------------
!
      DO indxslic = 1, nslice(slicopt)

! Define which trajectory to plot the vertical scross section for
        IF(trajc_plt_opt/=0) itrajc_index=itrajc_start+(indxslic-1)*ntrajc_stride

        agl=0
        islice=-1
        jslice=-1
        kslice=-1
        isoilslice=-1
        jsoilslice=-1
        ksoilslice=-1

        IF (iplot(slicopt) /= 0) THEN

          ! Regular vertical cross-sections, i.e., i or j is fixed
          IF (slicopt == 2) jslice=slice_xz(indxslic)
          IF (slicopt == 3) islice=slice_yz(indxslic)

          IF (jslice == -1) jslice = (nylg-2)/2+1
          IF (jslice == -2) jslice = jmax
          IF (jslice == -3) jslice = jmin
          IF (islice == -1) islice = (nxlg-2)/2+1
          IF (islice == -2) islice = imax
          IF (islice == -3) islice = imin

!
! Horizontal plots, i.e., k is fixed
!
          IF( slicopt == 1) kslice=slice_xy(indxslic)
          IF( kslice == -1) kslice = (nz-2)/2+1
          IF( kslice == -2) kslice = kmax
!
! Other plots
!
          IF(slicopt == 4) THEN
            IF( h_follow_trajc == 1 ) then ! following the trajectory parcel

              jpoint_final = 0
              DO jpoint=npoints_bgn(itrajc_index,1),npoints_end(itrajc_index,1)
                IF( abs(time-ttrajc(jpoint)) < 0.1 ) THEN
                  jpoint_final = jpoint
                  EXIT
                ENDIF
              ENDDO

              IF( jpoint_final == 0 ) then
                zlevel= ztrajc(npoints_bgn(itrajc_index,1),itrajc_index,1)*0.001  ! Set of the z of 1st point
              ELSE
                zlevel=ztrajc(jpoint_final,itrajc_index,1)*0.001
              ENDIF

            ELSE
              zlevel=slice_h (indxslic)
              IF( zlevel < 0.0 ) THEN
                zlevel = abs(zlevel)
                agl = 1
              END IF
            ENDIF
          END IF

          IF( slicopt <= 5 .OR. slicopt==12) THEN
            DO k=1,nz
              DO j=1,ny
                DO i=1,nx
                  zps3d(i,j,k) = zpc(i,j,k)
                END DO
              END DO
            END DO
          END IF

          IF( slicopt == 6) THEN
            zlevel=-ALOG(100.0*slice_p(indxslic))
            DO k=1,nz
              DO j=1,ny
                DO i=1,nx
                  zps3d(i,j,k) = algpzc(i,j,k)
                END DO
              END DO
            END DO
          END IF

          IF( slicopt == 7) THEN
            zlevel=slice_pt(indxslic)
            DO k=1,nz
              DO j=1,ny
                DO i=1,nx
                  zps3d(i,j,k) = ptbar(i,j,k)+ptprt(i,j,k)
                END DO
              END DO
            END DO
          END IF

          IF( slicopt == 5) THEN
            x01=xi1(indxslic)
            y01=yi1(indxslic)
            x02=xi2(indxslic)
            y02=yi2(indxslic)

            CALL clipwd(x01,y01,x02,y02, idisplay )

            IF(idisplay == 0) THEN
              WRITE(6,'(2(/1x,a))')                                     &
                  'The specified vertical slice was outside the plotting ', &
                  'window, no cross-section is plotted.'
              CYCLE
            END IF

            sqrtdxy = SQRT(dxkm*dykm)
            dist    = SQRT((x01-x02)**2+(y01-y02)**2)
            n       = INT(dist/sqrtdxy+0.5)+1

            IF (n > nx+ny) THEN
              WRITE(6,*) 'ERROR: Working arrays are too small for ',    &
                         ' arbitrary vertical slice. '
              WRITE(6,*) '       See restriction for MPI mode in arpsplt.input'
              CALL arpsstop('Local arrays too small',1)
            END IF
            sinaf=(y02-y01)/dist
            cosaf=(x02-x01)/dist

            DO i=1,n
              xp(i)=x01+(FLOAT(i-1))*sqrtdxy*cosaf
              yp(i)=y01+(FLOAT(i-1))*sqrtdxy*sinaf
            END DO

            ibgn = 1
            iend = n

            ibgnl  = 1
            iendl  = n
            ibgnl2 = 1

            xr = (n-1)*sqrtdxy

            z1 = zmin
            z2 = zmax
            IF (zbgn /= zend) THEN
              z1 = MAX(zbgn,zmin)
              z2 = MIN(zend,zmax)
            END IF

            zr = z2-z1

            x1 = 0.0
            x2 = x1+xr

            IF (zhstrch /= 0.0) THEN
              aspratio = zhstrch
            ELSE
              aspratio = xr/zr
            END IF

          END IF ! slicopt=5

          ! slices for soil variables, tsoil and qsoil
          IF (slicopt == 9) THEN

            ksoilslice=slice_xy_soil(indxslic)
            IF (ksoilslice == -1) ksoilslice = (nzsoil-2)/2+1

            DO k=1,nzsoil
              DO j=1,ny
                DO i=1,nx
                  zpsoils3d(i,j,k) = zpsoilc(i,j,k)
                END DO
              END DO
            END DO

          END IF

          IF (slicopt == 10) jsoilslice=slice_xz_soil(indxslic)
          IF (slicopt == 11) isoilslice=slice_yz_soil(indxslic)

          IF (jsoilslice == -1) jsoilslice = (nylg-2)/2+1
          IF (isoilslice == -1) isoilslice = (nxlg-2)/2+1

          SELECT CASE (slicopt)   ! Adjust active PEs and slice index
          !CASE (1,4,6,7,8,9)      ! xpbgn,xpend,ypbgn,ypend - no change
          CASE (2)
            ypbgn = (jslice-2)/(ny-3) + 1
            ypend = ypbgn
            jslice = MOD((jslice-2),(ny-3)) + 2
          CASE (3)
            xpbgn = (islice-2)/(nx-3) + 1
            xpend = xpbgn
            islice = MOD((islice-2),(nx-3)) + 2
          CASE (5, 12)     ! Only root PE do the plotting
            xpbgn = 1
            xpend = 1
            xpbgn = 1
            ypend = 1
          CASE (10)
            ypbgn = (jsoilslice-2)/(ny-3) + 1
            ypend = (jsoilslice-2)/(ny-3) + 1
            jsoilslice = MOD((jsoilslice-2),(ny-3)) + 2
          CASE (11)
            xpbgn = (isoilslice-2)/(nx-3) + 1
            xpend = (isoilslice-2)/(nx-3) + 1
            isoilslice = MOD((isoilslice-2),(nx-3)) + 2
          END SELECT

        END IF

        CALL xlbint(3)
        CALL xczero(3)

!
!-----------------------------------------------------------------------
!
!  Loop for every priority. small number is high priority, doing high
!  priority first , then low priority, finally doing the zero priority.
!  ipriority(): 1-80 3D plot, 81-120 2D plot, 121-150 surface character.
!  151-170 arbitrary 3D plot, 171-189 arbitrary 2D plot
!
!-----------------------------------------------------------------------
!
        DO ipp = 1,nprio   ! do loop priority

          sovrlay = 0
          plotovr = .false.
          ip=iptemp(ipp)

!
!-----------------------------------------------------------------------
!
!   Start plotting for all slicopt
!
!-----------------------------------------------------------------------
!
          IF( iplot(slicopt) /= 0 .AND. (slicopt <= 7 .OR. slicopt==12) ) THEN
!
!-----------------------------------------------------------------------
!
!   Plot height(10m) on constant pressure level
!        z-coor of sacalar point in physical space (10m)
!
!-----------------------------------------------------------------------
!

            IF(ipriority(1) == ip ) THEN
              IF(hplot == 1 .OR. hplot == 2 .OR. hplot == 4 .OR.        &
                 hplot == 5) THEN
                CALL get_mulovrlay('hplot',5,ovrmul_num,ovrmulname,sovrlay)
                IF(slicopt == 6 .OR. slicopt == 7 ) THEN
                  tem9(:,:,:) = 0.0
                  DO k=1,nz-1
                    DO j=1,ny
                      DO i=1,nx
                        tem9(i,j,k)=100.0*zpc(i,j,k)
                      END DO
                    END DO
                  END DO

                  CALL ctrsetup(hinc,hminc,hmaxc,                       &
                       hovr,hhlf,hzro,hcol1,hcol2,'hplot       ')

                  CALL ctr3d(tem9,xc,yc,zps3d,                          &
                      x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                 &
                      nx,ibgnl,iendl,ny,jbgnl,jendl,nz,kbgn,kend,       &
                      'h (10m)',time,slicopt, kslice,jslice,islice,     &
                      n,xp,yp,b1,b2,zs2,                                &
                      runname,1.0,tem1,tem2,tem3,                       &
                      tem4,tem5,tem6,hterain,hplot)
                END IF
              END IF
            END IF
!
!-----------------------------------------------------------------------
!
!   Plot thickness in 10meters between a given pressure level and
!   850mb level
!
!-----------------------------------------------------------------------
!
            IF(ipriority(130) == ip ) THEN
              IF(thkplt == 1 .OR. thkplt == 2 .OR. thkplt == 4 .OR.     &
                    thkplt == 5) THEN
                CALL get_mulovrlay('thkplt',5,ovrmul_num,ovrmulname,sovrlay)
                IF(slicopt == 6 .OR. slicopt == 7 ) THEN

                  CALL hintrp1                                          &
                       (nx,ny,nz,2,nz-2,zpc ,zps3d,zlevel,tem9(1,1,1))

                  pref = 800.0 ! mb
                  IF( slicopt == 6) THEN
                    tem=-ALOG(100.0*pref)
                  END IF
                  IF( slicopt == 7) THEN
                    tem = 300.0 ! K
                  END IF

                  CALL hintrp1                                          &
                       (nx,ny,nz,2,nz-2,zpc ,zps3d,tem,tem9(1,1,2))

                  DO j=1,ny-1
                    DO i=1,nx-1
                      IF( ABS(tem9(i,j,1)+9999.0) < 0.1 .OR.            &
                            ABS(tem9(i,j,2)+9999.0) < 0.1 ) THEN
                        tem9(i,j,1)=-9999.0
                      ELSE
                        tem9(i,j,1)=100.0*(tem9(i,j,1)-tem9(i,j,2))
                      END IF
                    END DO
                  END DO

                  CALL ctrsetup(thkinc,thkminc,thkmaxc,thkovr,          &
                       thkhlf,thkzro,thkcol1,thkcol2,'thkplt      ')

                  WRITE(label,'(i4,''-'',I4,''MB THICKNESS (10M)'')')   &
                        nint(slice_p(indxslic)),nint(pref)

                  CALL ctrsfc(tem9(1,1,1),xc(1,1,2),yc(1,1,2),          &
                       x1,x2,dxkm,y1,y2,dykm,nx,ibgnl,iendl,ny,         &
                       jbgnl,jendl, label(1:27),time, runname,1.0,      &
                       tem1,tem2,tem3,tem4,tem5,hterain,slicopt,thkplt)
                       ! size of tem5 must be >= 6*nx*ny

                END IF
              END IF
            END IF
!
!-----------------------------------------------------------------------
!
!   Plot ageostrophic winds on pressure surfaces
!
!-----------------------------------------------------------------------
!

            IF(ipriority(129) == ip ) THEN
              IF(vagplt == 1 .OR. vagplt == 2 .OR. vagplt == 4 .OR.     &
                    vagplt == 5) THEN
                CALL get_mulovrlay('vagplt',5,ovrmul_num,ovrmulname,sovrlay)
                IF(slicopt == 6 ) THEN ! Pressue surface plot

                  DO k=1,nz-1
                    DO j=1,ny-1
                      DO i=1,nx-1
                        tem1(i,j,k)=(u(i+1,j,k)+u(i,j,k))*0.5
                        tem2(i,j,k)=(v(i,j+1,k)+v(i,j,k))*0.5
                      END DO
                    END DO
                  END DO

                  CALL hintrp1                                          &
                       (nx,ny,nz,2,nz-2,tem1,zps3d,zlevel,tem9(1,1,1))
                  CALL hintrp1                                          &
                       (nx,ny,nz,2,nz-2,tem2,zps3d,zlevel,tem9(1,1,2))
                  CALL hintrp1                                          &
                       (nx,ny,nz,2,nz-2,zpc ,zps3d,zlevel,tem9(1,1,3))

!          do i=1,20
!            CALL smooth9pmv(tem9(1,1,3),nx,ny,2,nx-2,2,ny-2,tem5)
!          enddo

                  DO j=2,ny-2
                    DO i=1,nx-1
                      IF(ABS(tem9(i,j  ,1)+9999.0) < 0.1.OR.            &
                            ABS(tem9(i,j+1,3)+9999.0) < 0.1.OR.         &
                            ABS(tem9(i,j-1,3)+9999.0) < 0.1 )THEN
                        tem8(i,j,1)=-9999.0
                      ELSE
                        tem8(i,j,1)=                                    &
                             tem9(i,j,1)+ mapfct(i,j)*                  &
                            g*(tem9(i,j+1,3)-tem9(i,j-1,3))/(fcorio(i,j)*2*dykm)
                      END IF
                    END DO
                  END DO
                  DO j=1,ny-1
                    DO i=2,nx-2
                      IF(ABS(tem9(i,j  ,1)+9999.0) < 0.1.OR.            &
                            ABS(tem9(i+1,j,3)+9999.0) < 0.1.OR.         &
                            ABS(tem9(i-1,j,3)+9999.0) < 0.1)THEN
                        tem8(i,j,2)=-9999.0
                      ELSE
                        tem8(i,j,2)=                                    &
                                  tem9(i,j,2)- mapfct(i,j)*             &
                             g*(tem9(i+1,j,3)-tem9(i-1,j,3))/(fcorio(i,j)*2*dxkm)
                      END IF
                    END DO
                  END DO

                  CALL vtrunt ( vagunit )
                  IF(plotovr) CALL overlay ( 1 )
                  IF(.NOT.plotovr)  CALL overlay( vagovr )
                  CALL ctrcol(vagcol1,vagcol2)
                  CALL varplt( 'vagplt      ' )
                  CALL ctrvtr( vagunits, vagtype )

                  WRITE(label,'(''AG AT '',I4,''(HPA)'')')              &
                                nint(slice_p(indxslic))

                  CALL vtrsfc( tem8(1,1,1),tem8(1,1,2),xc(1,1,2),yc(1,1,2), &
                      x1,x2,dxkm,y1,y2,dykm, nx,ibgnl2,iendl,ist,           &
                      ny,jbgnl2,jendl,jst,label(1:15),time, runname,1.0,    &
                      slicopt,tem1,tem2,tem3,tem4,tem5,tem6,hterain)
                      ! size of tem6 must be >= 5*nx*ny
                  sigplt(129) = 1

                END IF
              END IF
            END IF

!
!-----------------------------------------------------------------------
!
!   Plot temperature (C)
!
!-----------------------------------------------------------------------
!
            IF(ipriority(2) == ip ) THEN

              IF( tplot > 0) THEN
                CALL cal_t(tem9,tz,nx, ny, nz,tob, label,length,tunits )
              END IF

              IF(tplot == 1 .OR. tplot == 2 .OR. tplot == 4 .OR.        &
                    tplot == 5 ) THEN
                CALL get_mulovrlay('tplot',5,ovrmul_num,ovrmulname,sovrlay)

!         pltvar = 'tplot'
!         CALL spltpara(tinc, tminc, tmaxc, tovr, thlf, tzro,
!  :                    tcol1,tcol2,pltvar)

                CALL ctrsetup(tinc,tminc,tmaxc,                         &
                     tovr,thlf,tzro,tcol1,tcol2,'tplot       ')

                CALL xcfull

                CALL ctr3d(tem9,xc,yc,zps3d,                            &
                    x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                   &
                    nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,       &
                    label(1:length),time,slicopt, kslice,jslice,islice, &
                    n,xp,yp,b1,b2,zs2,                                  &
                    runname,1.0,tem1,tem2,tem3,                         &
                    tem4,tem5,tem6,hterain,tplot)

                CALL xcmixl

              END IF

            END IF
!
!-----------------------------------------------------------------------
!
!   Plot wind components
!
!-----------------------------------------------------------------------
!
            IF(ipriority(3) == ip ) THEN
              IF(uplot == 1 .OR. uplot == 2 .OR. uplot == 4 .OR.        &
                    uplot == 5) THEN
                CALL get_mulovrlay('uplot',5,ovrmul_num,ovrmulname,sovrlay)

                IF(nobs > 0 .AND. slicopt == 1 .AND. kslice == 2) THEN
                  DO iob=1,nobs
                    IF(dd(iob) >= 0. .AND. dd(iob) < 360. .AND.         &
                          ff(iob) >= 0. .AND. ff(iob) < 60.) THEN
                      CALL ddrotuv(1,lonob(iob),dd(iob),ff(iob),        &
                                   drot,obs1(iob),obs2(iob))
                      obs1(iob)=0.51444*obs1(iob)
                    ELSE
                      obs1(iob)=-999.
                    END IF
                  END DO
                  obsset=1
                END IF

                onvf = 0
                CALL avgx(u , onvf,                                     &
                    nx,ny,nz, 1,nx-1, 1,ny-1, 1,nz-1, tem9)

                CALL ctrsetup(uinc,uminc,umaxc,                         &
                     uovr,uhlf,uzro,ucol1,ucol2,'uplot       ')

                CALL ctr3d( tem9, xc,yc,zps3d,                          &
                    x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                   &
                    nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,       &
                    'u (m/s)',time,slicopt, kslice,jslice,islice,       &
                    n,xp,yp,b1,b2,zs2,                                  &
                    runname,1.0,tem1,tem2,tem3,                         &
                    tem4,tem5,tem6,hterain,uplot)
              END IF

            END IF

            IF(ipriority(4) == ip ) THEN
              IF(vplot == 1 .OR. vplot == 2 .OR. vplot == 4 .OR.        &
                    vplot == 5 ) THEN
                CALL get_mulovrlay('vplot',5,ovrmul_num,ovrmulname,sovrlay)

                IF(nobs > 0 .AND. slicopt == 1 .AND. kslice == 2) THEN
                  DO iob=1,nobs
                    IF(dd(iob) >= 0. .AND. dd(iob) < 360. .AND.         &
                          ff(iob) >= 0. .AND. ff(iob) < 60.) THEN
                      CALL ddrotuv(1,lonob(iob),dd(iob),ff(iob),        &
                                   drot,obs1(iob),obs2(iob))
                      obs1(iob)=0.51444*obs2(iob)
                    ELSE
                      obs1(iob)=-999.
                    END IF
                  END DO
                  obsset=1
                END IF

                onvf = 0
                CALL avgy(v,onvf,nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,tem9)

                CALL ctrsetup(vinc,vminc,vmaxc,                         &
                     vovr,vhlf,vzro,vcol1,vcol2,'vplot       ')

                CALL ctr3d(tem9, xc,yc,zps3d,                           &
                    x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                   &
                    nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,       &
                    'v (m/s)',time,slicopt, kslice,jslice,islice,       &
                    n,xp,yp,b1,b2,zs2,                                  &
                    runname,1.0,tem1,tem2,tem3,                         &
                    tem4,tem5,tem6,hterain,vplot)
              END IF

            END IF

            IF(ipriority(5) == ip ) THEN

              IF(vhplot > 0) THEN
                CALL cal_vh(tem9,u,v,nx,ny,nz,vhunits,label,length,tem8)
              END IF

              IF(vhplot == 1 .OR. vhplot == 2 .OR. vhplot == 4 .OR.     &
                    vhplot == 5 .OR.vhplot == 6 ) THEN

                CALL get_mulovrlay('vhplot',6,ovrmul_num,ovrmulname,sovrlay)

                CALL ctrsetup(vhinc,vhminc,vhmaxc,                      &
                     vhovr,vhhlf,vhzro,vhcol1,vhcol2,'vhplot      ')

                CALL ctr3d(tem9, xc,yc,zps3d,                           &
                    x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                   &
                    nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,       &
                    label(1:length),time,slicopt, kslice,jslice,islice, &
                    n,xp,yp,b1,b2,zs2,                                  &
                    runname,1.0,tem1,tem2,tem3,                         &
                    tem4,tem5,tem6,hterain,vhplot)
              END IF

            END IF

            IF(ipriority(6) == ip ) THEN
              IF(wplot == 1 .OR. wplot == 2 .OR. wplot == 4 .OR.        &
                    wplot == 5   ) THEN
                CALL get_mulovrlay('wplot',5,ovrmul_num,ovrmulname,sovrlay)

                onvf = 0
                CALL avgz(w,onvf,nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,tem9)

                CALL ctrsetup(winc,wminc,wmaxc,                         &
                     wovr,whlf,wzro,wcol1,wcol2,'wplot       ')

                CALL ctr3d(tem9,xc,yc,zps3d,                            &
                    x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                   &
                    nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,       &
                    'w (m/s)',time,slicopt, kslice,jslice,islice,       &
                    n,xp,yp,b1,b2,zs2,                                  &
                    runname,1.0,tem1,tem2,tem3,                         &
                    tem4,tem5,tem6,hterain,wplot)
              END IF

            END IF

!
!-----------------------------------------------------------------------
!
!   Plot scalars
!
!-----------------------------------------------------------------------
!
            IF(ipriority(7) == ip ) THEN
              IF(ptplot == 1 .OR. ptplot == 2 .OR. ptplot == 4 .OR.     &
                    ptplot == 5 ) THEN
                CALL get_mulovrlay('ptplot',6,ovrmul_num,ovrmulname,sovrlay)

                IF(ovrobs == 1 .AND. nobs > 0) THEN
                  DO iob=1,nobs
!
!  Station pressure (mb)
!
                    IF(pstn(iob) > 0.) THEN
                      psta=mbtopa*pstn(iob)
                    ELSE IF(alt(iob) > 0.) THEN
                      psta=mbtopa*alttostpr(alt(iob),elevob(iob))
                    ELSE
                      psta=mbtopa*alttostpr(1013.,elevob(iob))
                    END IF
!
!  Potential temperature (K)
!
                    IF(tob(iob) > -98.) THEN
                      tk=(5.*(tob(iob)-32.)/9.) + 273.15
                      obs1(iob)=tk*((p0/psta)**rddcp)
                    ELSE
                      obs1(iob)=-999.
                    END IF
                  END DO
                  obsset=1
                END IF

                CALL ctrsetup(ptinc,ptminc,ptmaxc,                      &
                     ptovr,pthlf,ptzro,ptcol1,ptcol2,'ptplot      ')

                CALL ctr3d( pt ,xc,yc,zps3d,                            &
                    x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                   &
                    nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,       &
                    'pt (K)',time,slicopt, kslice,jslice,islice,        &
                    n,xp,yp,b1,b2,zs2,                                  &
                    runname,1.0,tem1,tem2,tem3,                         &
                    tem4,tem5,tem6,hterain,ptplot)
              END IF

            END IF

            IF(ipriority(8) == ip ) THEN
              IF(pplot == 1 .OR. pplot == 2 .OR. pplot == 4 .OR.        &
                    pplot == 5 ) THEN
                CALL get_mulovrlay('pplot',5,ovrmul_num,ovrmulname,sovrlay)
                IF( ovrobs == 1 .AND. nobs > 0) THEN
!
                  DO iob=1,nobs
!
!  Station pressure (mb)
!
                    IF(pstn(iob) > 0.) THEN
                      obs1(iob)=pstn(iob)-100.*(MOD(INT(pstn(iob)),100))
                    ELSE IF(alt(iob) > 0.) THEN
                      pmb=alttostpr(alt(iob),elevob(iob))
                      obs1(iob)=pmb -100.*(MOD(INT(pmb),100))
                    ELSE
                      obs1(iob)=-999.
                    END IF
                  END DO
                  obsset=1
                END IF
!
                DO k=1,nz-1
                  DO j=1,ny-1
                    DO i=1,nx-1
                      tem9(i,j,k)=pbar(i,j,k)+pprt (i,j,k)
                    END DO
                  END DO
                END DO

                CALL ctrsetup(pinc,pminc,pmaxc,                         &
                     povr,phlf,pzro,pcol1,pcol2,'pplot       ')

                CALL ctr3d( tem9 ,xc,yc,zps3d,                          &
                    x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                   &
                    nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,       &
                    'p (Pa)',time,slicopt, kslice,jslice,islice,        &
                    n,xp,yp,b1,b2,zs2,                                  &
                    runname,1.0,tem1,tem2,tem3,                         &
                    tem4,tem5,tem6,hterain,pplot)

              END IF
            END IF
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!

            IF(ipriority(9) == ip) THEN
              IF(qvplot == 1 .OR. qvplot == 2 .OR. qvplot == 4 .OR.     &
                    qvplot == 5 ) THEN
                CALL get_mulovrlay('qvplot',6,ovrmul_num,ovrmulname,sovrlay)
                IF( ovrobs == 1 .AND. nobs > 0 ) THEN
                  DO iob=1,nobs
!
!  Station pressure
!
                    IF(pstn(iob) > 0.) THEN
                      psta=mbtopa*pstn(iob)
                    ELSE IF(alt(iob) > 0.) THEN
                      psta=mbtopa*alttostpr(alt(iob),elevob(iob))
                    ELSE
                      psta=mbtopa*alttostpr(1013.,elevob(iob))
                    END IF
!
!  Specific humidity
!
                    IF(tdob(iob) > -90.) THEN
                      tdk=(5.*(tdob(iob)-32.)/9.) + 273.15
                      obs1(iob) = f_qvsat( psta,tdk ) *1000.
                    ELSE
                      obs1(iob)=-999.
                    END IF
                  END DO
                  obsset=1
                END IF

                CALL ctrsetup(qvinc,qvminc,qvmaxc,                      &
                     qvovr,qvhlf,qvzro,qvcol1,qvcol2,'qvplot      ')

                tem7 = qv*1000.0
                CALL ctr3d( tem7,xc,yc,zps3d,                             &
                    x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                   &
                    nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,       &
                    'qv (g/kg)',time,slicopt,kslice,jslice,islice,      &
                    n,xp,yp,b1,b2,zs2,                                  &
                    runname,1.0,tem1,tem2,tem3,                      &
                    tem4,tem5,tem6,hterain,qvplot)
              END IF
            END IF
!
!  Plotting of scalar array (DTD 09/18/06)
!
            factor = 1.0

            DO nq=1,nscalarmax
              IF(ipriority(nq+9) == ip) THEN
                !print*,'ipriority(nq+9),ip',ipriority(nq+9),ip
                IF(qscalarplot(nq) == 1 .OR. qscalarplot(nq) == 2 .OR.  &
                   qscalarplot(nq) == 4 .OR. qscalarplot(nq) == 5 ) THEN

                    CALL get_mulovrlay(TRIM(qnames(nq))//'plot',6,      &
                                       ovrmul_num,ovrmulname,sovrlay)

                    CALL ctrsetup(qscalarinc(nq),qscalarminc(nq),       &
                                  qscalarmaxc(nq), qscalarovr(nq),      &
                                  qscalarhlf(nq),qscalarzro(nq),        &
                                  qscalarcol1(nq),qscalarcol2(nq),      &
                         TRIM(qnames(nq))//'plot      ')

                    !print*,'ipriority(nq+9),ip',ipriority(nq+9),ip
                    !print*,'nq,qscalarovr(nq)',nq,qscalarovr(nq)

                    ! Determine scaling factor based on whether variable to
                    ! be plotted is mixing ratio (Q), number concentration (N)
                    ! or reflectivity (Z)

                    IF(qnames(nq)(1:1)=='q') THEN
                      tem7(:,:,:) = qscalar(:,:,:,nq)*1000.0
                      factor=1.0
                    ELSE IF(qnames(nq)(1:1)=='n') THEN
                      tem7(:,:,:) = qscalar(:,:,:,nq)
                      factor=1.0
                    ELSE IF(qnames(nq)(1:1)=='z') THEN
                      tem7(:,:,:) = qscalar(:,:,:,nq)
                      factor=1.0e18
                    END IF

                    ! Reset zero values to a small negative number to allow meaningful
                    ! zero contours
                    WHERE (tem7 <= 0.0) tem7 = -1.0e-20


                    CALL ctr3d( tem7,xc,yc,zps3d,                       &
                        x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,               &
                        nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,   &
                        qdescp(nq),time,slicopt,kslice,jslice,islice,   &
                        n,xp,yp,b1,b2,zs2,                              &
                        runname,factor,tem1,tem2,tem3,                  &
                        tem4,tem5,tem6,hterain,qscalarplot(nq))

                END IF
              END IF
            END DO

            IF(ipriority(41) == ip ) THEN

              IF(qwplot > 0) THEN
                CALL cal_qw(tem9,qscalar, nx,ny,nz)
              END IF

              IF(qwplot == 1 .OR. qwplot == 2 .OR. qwplot == 4 .OR.     &
                    qwplot == 5 ) THEN

                CALL get_mulovrlay('qwplot',6,ovrmul_num,ovrmulname,sovrlay)

                CALL ctrsetup(qwinc,qwminc,qwmaxc,                      &
                     qwovr,qwhlf,qwzro,qwcol1,qwcol2,'qwplot      ')

                WHERE (tem9 <= 0.0) tem9 = -1.0e-20

                tem9 = tem9 * 1000.0
                CALL ctr3d( tem9,xc,yc,zps3d,                           &
                     x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                  &
                     nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,      &
                     'Total water (g/kg)',time,slicopt,kslice,jslice,islice, &
                     n,xp,yp,b1,b2,zs2,                                 &
                     runname,1.0,tem1,tem2,tem3,                     &
                     tem4,tem5,tem6,hterain,qwplot)
              END IF

            END IF
!
!-----------------------------------------------------------------------
!
!    Calculate relative humidity, where tem1 = temperature,
!    tem2 = saturation qv.
!
!-----------------------------------------------------------------------
!
            IF(ipriority(42) == ip ) THEN
              IF(kmhplt == 1 .OR. kmhplt == 2 .OR. kmhplt == 4 .OR.     &
                    kmhplt == 5 ) THEN
                CALL get_mulovrlay('kmhplt',6,ovrmul_num,ovrmulname,sovrlay)
                CALL ctrsetup(kmhinc,kmhminc,kmhmaxc,                   &
                     kmhovr,kmhhlf,kmhzro,kmhcol1,kmhcol2,'kmhplt       ')
                CALL ctr3d( kmh,xc,yc,zps3d,                            &
                     x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                  &
                     nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,      &
                     'kmh (m**2/s)',time,slicopt,kslice,jslice,islice,  &
                     n,xp,yp,b1,b2,zs2,                                 &
                     runname,1.0,tem1,tem2,tem3,                        &
                     tem4,tem5,tem6,hterain,kmhplt)
              END IF
            END IF

            IF(ipriority(43) == ip ) THEN
              IF(kmvplt == 1 .OR. kmvplt == 2 .OR. kmvplt == 4 .OR.     &
                    kmvplt == 5 ) THEN
                CALL get_mulovrlay('kmvplt',6,ovrmul_num,ovrmulname,sovrlay)
                CALL ctrsetup(kmvinc,kmvminc,kmvmaxc,                   &
                     kmvovr,kmvhlf,kmvzro,kmvcol1,kmvcol2,'kmvplt       ')
                CALL ctr3d( kmv,xc,yc,zps3d,                            &
                     x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                  &
                     nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,      &
                     'kmv (m**2/s)',time,slicopt,kslice,jslice,islice,  &
                     n,xp,yp,b1,b2,zs2,                                 &
                     runname,1.0,tem1,tem2,tem3,                        &
                     tem4,tem5,tem6,hterain,kmvplt)
              END IF
            END IF

            IF(ipriority(44) == ip ) THEN
              IF(tkeplt == 1 .OR. tkeplt == 2 .OR. tkeplt == 4 .OR.     &
                    tkeplt == 5 ) THEN
                CALL get_mulovrlay('tkeplt',6,ovrmul_num,ovrmulname,sovrlay)
                CALL ctrsetup(tkeinc,tkeminc,tkemaxc,                   &
                     tkeovr,tkehlf,tkezro,tkecol1,tkecol2,'tkeplt       ')
                CALL ctr3d( tke,xc,yc,zps3d,                            &
                     x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                  &
                     nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,      &
                     'tke ((m/s)**2)',time,slicopt,kslice,jslice,islice, &
                     n,xp,yp,b1,b2,zs2,                                 &
                     runname,1.0,tem1,tem2,tem3,                        &
                     tem4,tem5,tem6,hterain,tkeplt)
              END IF
            END IF

            IF(ipriority(45) == ip) THEN

              IF(rhplot > 0) THEN
                CALL cal_rh(tem9,pt, pprt ,pbar,qv,tem1,tem2,nx,ny,nz)
              END IF

              IF(rhplot == 1 .OR. rhplot == 2 .OR. rhplot == 4 .OR.     &
                    rhplot == 5 ) THEN

                CALL get_mulovrlay('rhplot',6,ovrmul_num,ovrmulname,sovrlay)
                CALL ctrsetup(rhinc,rhminc,rhmaxc,                      &
                     rhovr,rhhlf,rhzro,rhcol1,rhcol2,'rhplot      ')
                CALL ctr3d( tem9,xc,yc,zps3d,                           &
                     x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                  &
                     nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,      &
                     'RH (%)',time,slicopt,kslice,jslice,islice,        &
                     n,xp,yp,b1,b2,zs2,                                 &
                     runname,1.0,tem1,tem2,tem3,                        &
                     tem4,tem5,tem6,hterain,rhplot)
              END IF

            END IF

!
!-----------------------------------------------------------------------
!
!  Plot dew-point temperature td tdplot
!
!-----------------------------------------------------------------------

            IF(ipriority(46) == ip ) THEN

              IF(tdplot > 0) THEN
                CALL cal_td(tem9,td,nx,ny,nz,tdunits,label, length)
              END IF

              IF(tdplot == 1 .OR. tdplot == 2 .OR. tdplot == 4 .OR.     &
                    tdplot == 5 ) THEN
                CALL get_mulovrlay('tdplot',6,ovrmul_num,ovrmulname,sovrlay)
                CALL cal_tdobs(tdob,tdunits)

                CALL ctrsetup(tdinc,tdminc,tdmaxc,                      &
                     tdovr,tdhlf,tdzro,tdcol1,tdcol2,'tdplot      ')

                CALL ctr3d( tem9,xc,yc,zps3d,                           &
                     x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                  &
                     nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,      &
                     label(1:length),time,slicopt,kslice,jslice,islice, &
                     n,xp,yp,b1,b2,zs2,                                 &
                     runname,1.0,tem1,tem2,tem3,                        &
                     tem4,tem5,tem6,hterain,tdplot)
              END IF
            END IF

!-----------------------------------------------------------------------
!
!  Calculate alpha for DSD to calculte polarimetric variables
!  Added by Youngsun Jung
!-----------------------------------------------------------------------

            IF(ipriority(47) == ip ) THEN
              IF(rfplot == 1 .OR. rfplot == 2 .OR.                      &
                 rfplot == 4 .OR. rfplot == 5 ) THEN
                CALL get_mulovrlay('rfplot',6,ovrmul_num,ovrmulname,sovrlay)

                IF(dualpol >= 1) THEN
!                    CALL ZhhFromDualPolCG(nx,ny,nz,rhobar,tz,qscalar,      &
!                                          tem9,tem8,dualpol,wavelen,alpha)
                ELSE

                  IF (mphyopt < 100) THEN  ! ARPS microphysical schemes

                    SELECT CASE (mphyopt)
                    CASE (8:12)  ! Milbrandt and Yau 3-moment scheme
                      CALL reflec_MM(nx,ny,nz,rhobar,qscalar,tz,tem9)

                    CASE (2:7)   ! Lin,Schultz,Straka Lin,or WSM6 schemes
                      CALL reflec_ferrier(nx,ny,nz, rhobar, qscalar, tz, tem9)

                    CASE (1)     ! Warm rain schemes
                      IF(P_QR > 0) THEN
                        CALL reflec_wr(nx,ny,nz, rhobar, qscalar(:,:,:,P_QR),tem9)
                      ELSE
                        tem9 = 0.0
                      END IF
                    END SELECT

                  ELSE IF (mphyopt > 100 .AND. mphyopt < 200) THEN
                    wrfopt = MOD(mphyopt,100)

                    DO k=1,nz-1
                      DO j=1,ny-1
                        DO i=1,nx-1
                          tem1(i,j,k) = pbar(i,j,k )+pprt (i,j,k)
                        END DO
                      END DO
                    END DO
                    CALL edgfill(tem1,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)

                    SELECT CASE (wrfopt)

                    CASE (1,3,4)     ! Kessler, WSM 3-class
                      print *,' Reflectivity formula to be implemented for WRF mpphysics = ',wrfopt,'.'
                      CALL arpsstop('ERROR: no reflectivity formular implemented.',1)
                    CASE (2,6,10,16) ! Lin, WSM 5- and 6- classes)
                      print *,' Calling reflec_wrf ...'
                      CALL reflec_wrf(nx,ny,nz,qv,qscalar(:,:,:,P_QR),  &
                                                  qscalar(:,:,:,P_QS),  &
                                                  qscalar(:,:,:,P_QG),  &
                                      1,tem1,tz,tem9)
                                                           ! mixed-phase
                    CASE (5)     ! Ferrier microphysics scheme (as in WRF_POST)
                      print *,' Calling reflec_ferrier_wrf ...'
                      tem2(:,:,:) = 2.0
                      CALL reflec_ferrier_wrf(nx,ny,nz,qv,qscalar(:,:,:,P_QC), &
                                                          qscalar(:,:,:,P_QR), &
                                                          qscalar(:,:,:,P_QS), &
                                              tem1,tz,tem9,tem8,tem2)
                    CASE (8) ! Thompson microphysics scheme (old version)
                      print *,' Calling CALREF9s ...'
                      CALL CALREF9s(nx,ny,nz,qscalar(:,:,:,P_QR),  &
                                             qscalar(:,:,:,P_QS),  &
                                             qscalar(:,:,:,P_QG),  &
                                    tem1,tz,tem9)

                    CASE (9,17)
                      print *,' Calling reflec_wrf ...'
                      tem2(:,:,:) = qscalar(:,:,:,P_QG)+qscalar(:,:,:,P_QH)
                      CALL reflec_wrf(nx,ny,nz,qv,qscalar(:,:,:,P_QR),  &
                                                  qscalar(:,:,:,P_QS),  &
                                                  tem2,                 &
                                      1,tem1,tz,tem9) ! mixed-phase
                    CASE DEFAULT
                      WRITE(6,'(1x,a,I2,a)') 'ERROR: Unknown mphyopt = ',mphyopt,'.'
                      CALL arpsstop('ERROR: Unsupported mphopt option',1)
                    END SELECT

                  ELSE IF (mphyopt > 200) THEN ! COAMPS original scheme
                    DO k=1,nz-1
                      DO j=1,ny-1
                        DO i=1,nx-1
                          tem1(i,j,k) = pbar(i,j,k )+pprt (i,j,k)
                        END DO
                      END DO
                    END DO
                    CALL edgfill(tem1,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)

                    SELECT CASE (mphyopt)
                    CASE (201)
                      CALL coamps_rrf(tz,tem1,qscalar(:,:,:,P_QR),      &
                                              qscalar(:,:,:,p_QS),      &
                                              qscalar(:,:,:,P_QI),      &
                                              qscalar(:,:,:,P_QC),      &
                                              qscalar(:,:,:,P_QG),      &
                                      qv,nx,ny,nz,tem9)
                    CASE (204)  ! Thompson scheme
                      CALL CALREF9s(nx,ny,nz,qscalar(:,:,:,P_QR),       &
                                             qscalar(:,:,:,P_QS),       &
                                             qscalar(:,:,:,P_QG),       &
                                    tem1,tz,tem9)
                    CASE (205:209) ! MY scheme

                      CALL reflec_MM(nx,ny,nz,rhobar,qscalar,tz,tem9)

                    END SELECT

                  ELSE
                    print*,'Invalid microphysics option, reflectivity set to zero'
                    tem9 = 0.0
                  END IF

                END IF

                CALL ctrsetup(rfinc,rfminc,rfmaxc,                      &
                     rfovr,rfhlf,rfzro,rfcol1,rfcol2,'rfplot      ')

                CALL ctr3d( tem9,xc,yc,zps3d,                           &
                     x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                  &
                     nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,      &
                     'Ref (dBZ)',time,slicopt,kslice,jslice,islice,     &
                     n,xp,yp,b1,b2,zs2,                                 &
                     runname,1.0,tem1,tem2,tem3,                        &
                     tem4,tem5,tem6,hterain,rfplot)
              END IF

            END IF

            IF(ipriority(71) == ip) THEN

              IF(zdrplt == 1 .OR. zdrplt == 2 .OR. zdrplt == 4 .OR.     &
                    zdrplt == 5 ) THEN
                CALL get_mulovrlay('zdrplt',6,ovrmul_num,ovrmulname,sovrlay)
!                CALL ZdrFromDualPolCG(nx,ny,nz,rhobar,tz,qscalar,tem9,     &
!                                      dualpol,wavelen,alpha)

                WHERE (tem9 <= 0.0) tem9 = -1.0e-20

                CALL ctrsetup(zdrinc,zdrminc,zdrmaxc,                      &
                   zdrovr,zdrhlf,zdrzro,zdrcol1,zdrcol2,'zdrplt      ')
                CALL ctr3d( tem9,xc,yc,zps3d,                           &
                   x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                  &
                   nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,      &
                   'Zdr (dB)',time,slicopt,kslice,jslice,islice,     &
                   n,xp,yp,b1,b2,zs2,                                 &
                   runname,1.0,tem1,tem2,tem3,                        &
                   tem4,tem5,tem6,hterain,zdrplt)

              END IF
            END IF

            IF(ipriority(72) == ip) THEN

              IF(kdpplt == 1 .OR. kdpplt == 2 .OR. kdpplt == 4 .OR.     &
                    kdpplt == 5 ) THEN
                CALL get_mulovrlay('kdpplt',6,ovrmul_num,ovrmulname,sovrlay)

!                CALL Kdp(nx,ny,nz,rhobar,tz,qscalar,tem9,dualpol,          &
!                         wavelen,alpha)

                WHERE (tem9 <= 0.00001) tem9 = -1.0e-20

                CALL ctrsetup(kdpinc,kdpminc,kdpmaxc,                      &
                   kdpovr,kdphlf,kdpzro,kdpcol1,kdpcol2,'kdpplt      ')
                CALL ctr3d( tem9,xc,yc,zps3d,                           &
                   x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                  &
                   nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,      &
                   'Kdp (deg/km)',time,slicopt,kslice,jslice,islice,     &
                   n,xp,yp,b1,b2,zs2,                                  &
                   runname,1.0,tem1,tem2,tem3,                        &
                   tem4,tem5,tem6,hterain,kdpplt)

              END IF
            END IF

            IF(ipriority(73) == ip) THEN

              IF(zdpplt == 1 .OR. zdpplt == 2 .OR. zdpplt == 4 .OR.     &
                    zdpplt == 5 ) THEN
                CALL get_mulovrlay('zdpplt',6,ovrmul_num,ovrmulname,sovrlay)
!                CALL ZdpFromDualPolCG(nx,ny,nz,rhobar,tz,qscalar,tem9,     &
!                                      dualpol,wavelen,alpha)

                WHERE (tem9 <= 0.0) tem9 = -1.0e-20  ! More elegant with where statement

                CALL ctrsetup(zdpinc,zdpminc,zdpmaxc,                      &
                   zdpovr,zdphlf,zdpzro,zdpcol1,zdpcol2,'zdpplt      ')
                CALL ctr3d( tem9,xc,yc,zps3d,                           &
                   x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                  &
                   nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,      &
                   'Zdp (mm**6/m**3)',time,slicopt,kslice,jslice,islice,     &
                   n,xp,yp,b1,b2,zs2,                                 &
                   runname,1.0,tem1,tem2,tem3,                        &
                   tem4,tem5,tem6,hterain,zdpplt)

              END IF
            END IF

            IF(ipriority(74) == ip) THEN

              IF(rhvplt == 1 .OR. rhvplt == 2 .OR. rhvplt == 4 .OR.     &
                    rhvplt == 5 ) THEN
                CALL get_mulovrlay('rhvplt',6,ovrmul_num,ovrmulname,sovrlay)
!                CALL rhvFromDualPolCG(nx,ny,nz,rhobar,tz,qscalar,tem9,     &
!                                      dualpol,wavelen,alpha)

                WHERE (tem9 <= 0.0) tem9 = 1.1

                CALL ctrsetup(rhvinc,rhvminc,rhvmaxc,                      &
                   rhvovr,rhvhlf,rhvzro,rhvcol1,rhvcol2,'rhvplt      ')
                CALL ctr3d( tem9,xc,yc,zps3d,                           &
                   x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                  &
                   nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,      &
                   'rho_hv ()',time,slicopt,kslice,jslice,islice,     &
                   n,xp,yp,b1,b2,zs2,                                 &
                   runname,1.0,tem1,tem2,tem3,                        &
                   tem4,tem5,tem6,hterain,rhvplt)

              END IF
            END IF

!
!-----------------------------------------------------------------------
!
!  Calculate composite reflectivity
!
!-----------------------------------------------------------------------
!

            !IF(ipriority(62) == ip .AND. P_QR > 0 .AND. P_QS > 0 .AND. P_QH > 0) THEN
            IF( ipriority(62) == ip .AND. P_QR > 0 .AND. P_QS > 0 ) THEN

              IF(rfcplt > 0) THEN

                IF(dualpol >= 1) THEN
!                    CALL ZhhFromDualPolCG(nx,ny,nz,rhobar,tz,qscalar,      &
!                                          tem9,tem8,dualpol,wavelen,alpha)
                ELSE

                  ! Milbrandt and Yau 3-moment scheme
                  IF (mphyopt >= 8 .and. mphyopt <= 12) THEN
                    CALL reflec_MM(nx,ny,nz,rhobar,qscalar,tz,tem8)

                  ! Lin,Schultz,Straka Lin,or WSM6 schemes
                  ELSE IF (mphyopt >= 2 .and. mphyopt <= 7) THEN
                    CALL reflec_ferrier(nx,ny,nz, rhobar, qscalar, tz, tem8)

                  ! Warm rain schemes
                  !ELSE IF (mphyopt == 1 .or. mphyopt == 12 .or. mphyopt == 13) THEN
                  ELSE IF (mphyopt == 1) THEN

                    IF(P_QR > 0) THEN
                      CALL reflec_wr(nx,ny,nz, rhobar, qscalar(:,:,:,P_QR),tem9)
                    ELSE
                      tem9 = 0.0
                    END IF

                    !IF (P_QH > 0) THEN
                    !  P_GH = P_QH
                    !ELSE IF (P_QG > 0) THEN
                    !  P_GH = P_QG
                    !ELSE
                    !  CALL arpsstop('ERROR: Neither P_QH nor P_QG is available.',1)
                    !END IF
                    !CALL reflec(nx,ny,nz, rhobar,qscalar(:,:,:,P_QR),     &
                    !                             qscalar(:,:,:,P_QS),     &
                    !                             qscalar(:,:,:,P_GH),     &
                    !            tem8)
                  ELSE IF (mphyopt > 100 .AND. mphyopt < 200) THEN

                    wrfopt = MOD(mphyopt,100)

                    DO k=1,nz-1
                      DO j=1,ny-1
                        DO i=1,nx-1
                          tem1(i,j,k) = pbar(i,j,k )+pprt (i,j,k)
                        END DO
                      END DO
                    END DO
                    CALL edgfill(tem1,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)

                    SELECT CASE (wrfopt)

                    CASE (1,3,4) ! Kessler, WSM 3-class
                      print *,' Reflectivity formula to be implemented for WRF mpphysics = ',wrfopt,'.'
                      CALL arpsstop('ERROR: no reflectivity formular implemented.',1)
                    CASE (2,6,10,16) ! (Lin, WSM 5- and 6- classes)
                      print *,' CALL reflec_wrf ...'
                      CALL reflec_wrf(nx,ny,nz,qv,qscalar(:,:,:,P_QR),  &
                                                  qscalar(:,:,:,P_QS),  &
                                                  qscalar(:,:,:,P_QG),  &
                                      1,tem1,tz,tem8)
                    CASE (5) ! Ferrier microphysics scheme (as in WRF_POST)
                      print *,' Calling reflec_ferrier_wrf ...'
                      tem2(:,:,:) = 2.0
                      CALL reflec_ferrier_wrf(nx,ny,nz,qv,qscalar(:,:,:,P_QC), &
                                                          qscalar(:,:,:,P_QR), &
                                                          qscalar(:,:,:,P_QS), &
                                              tem1,tz,tem8,tem9,tem2)
                    CASE (8) ! Thompson microphysics scheme (old version)
                      print *,' Calling CALREF9s ...'
                      CALL CALREF9s(nx,ny,nz,qscalar(:,:,:,P_QR),         &
                                             qscalar(:,:,:,P_QS),         &
                                             qscalar(:,:,:,P_QG),         &
                                    tem1,tz,tem8)
                    CASE (9,17)
                      print *,' Calling reflec_wrf ...'
                      tem2(:,:,:) = qscalar(:,:,:,P_QG)+qscalar(:,:,:,P_QH)
                      CALL reflec_wrf(nx,ny,nz,qv,qscalar(:,:,:,P_QR),    &
                                                  qscalar(:,:,:,P_QS),    &
                                                  tem2,                   &
                                      1,tem1,tz,tem8) ! mixed-phase
                    CASE DEFAULT
                      WRITE(6,'(1x,a,I2,a)') 'ERROR: Unknown mphyopt = ',mphyopt,'.'
                      CALL arpsstop('ERROR: Unsupported mphopt option',1)
                    END SELECT

                  ELSE IF (mphyopt > 200) THEN ! COAMPS original scheme

                    DO k=1,nz-1
                      DO j=1,ny-1
                        DO i=1,nx-1
                          tem1(i,j,k) = pbar(i,j,k )+pprt (i,j,k)
                        END DO
                      END DO
                    END DO
                    CALL edgfill(tem1,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)

                    SELECT CASE (mphyopt)
                    CASE (201)
                      CALL coamps_rrf(tz,tem1,qscalar(:,:,:,P_QR),        &
                                              qscalar(:,:,:,p_QS),        &
                                              qscalar(:,:,:,P_QI),        &
                                              qscalar(:,:,:,P_QC),        &
                                              qscalar(:,:,:,P_QG),        &
                                      qv,nx,ny,nz,tem8)
                    CASE (204)
                      CALL CALREF9s(nx,ny,nz,qscalar(:,:,:,P_QR),         &
                                             qscalar(:,:,:,P_QS),         &
                                             qscalar(:,:,:,P_QG),         &
                                    tem1,tz,tem8)
                    CASE (205:209)
                      CALL reflec_MM(nx,ny,nz,rhobar,qscalar,tz,tem8)
                    END SELECT

                  ELSE
                    print*,'Invalid microphysics option, reflectivity set to zero'
                    tem8 = 0.0
                  END IF

                  CALL cal_rfc(nx, ny, nz, tem8, tem9)

                END IF
              END IF

              IF(rfcplt == 1 .OR. rfcplt == 2 .OR.                      &
                 rfcplt == 4 .OR. rfcplt == 5 ) THEN
                CALL get_mulovrlay('rfcplt',6,ovrmul_num,ovrmulname,sovrlay)

                label = 'Composite Ref (dBZ)'
                CALL ctrsetup(rfcinc,rfcminc,rfcmaxc,                   &
                     rfcovr,rfchlf,rfczro,rfccol1,rfccol2,'rfcplot     ')

                CALL ctrsfc(tem9,xc(1,1,2),yc(1,1,2),                   &
                     x1,x2,dxkm,y1,y2,dykm,                             &
                    nx,ibgnl,iendl, ny,jbgnl,jendl,                     &
                    label(1:19),time, runname, 1.0,tem1,tem2,tem3,      &
                    tem4,tem5,hterain,slicopt,rfcplt)
                       ! size of tem5 must be >= 6*nx*ny
              END IF

            END IF
!
!-----------------------------------------------------------------------
!
!  Calculate equivalent potential temperature.
!
!-----------------------------------------------------------------------
!
            IF(ipriority(48) == ip ) THEN
              IF(pteplt == 1 .OR. pteplt == 2 .OR. pteplt == 4 .OR.     &
                    pteplt == 5 ) THEN
                CALL get_mulovrlay('pteplt',6,ovrmul_num,ovrmulname,sovrlay)

                DO k=1,nz-1
                  DO j=1,ny-1
                    DO i=1,nx-1
                      tem1(i,j,k)=pprt(i,j,k)+pbar(i,j,k)
                    END DO
                  END DO
                END DO

                CALL pt2pte(nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,tem1,pt,qv,tem9)

                label = 'pte (K)'
                CALL ctrsetup(pteinc,pteminc,ptemaxc,                   &
                     pteovr,ptehlf,ptezro,ptecol1,ptecol2,'pteplot     ')
                CALL ctr3d(tem9,xc,yc,zps3d,                            &
                     x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                  &
                    nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,       &
                    label(1:7),time,slicopt,kslice,jslice,islice,       &
                    n,xp,yp,b1,b2,zs2,                                  &
                    runname,1.0,tem1,tem2,tem3,                         &
                    tem4,tem5,tem6,hterain,pteplt)
              END IF

            END IF
!
!-----------------------------------------------------------------------
!
!  Plot perturbation of wind components
!
!-----------------------------------------------------------------------
!
            IF(ipriority(49) == ip ) THEN

              IF(upplot > 0) THEN
                onvf = 0
                CALL avgx(uprt , onvf,                                  &
                    nx,ny,nz, 1,nx-1, 1,ny-1, 1,nz-1, tem9)
              END IF

              IF(upplot == 1 .OR. upplot == 2 .OR. upplot == 4 .OR.     &
                    upplot == 5 ) THEN
                CALL get_mulovrlay('upplot',6,ovrmul_num,ovrmulname,sovrlay)
                label = 'uprt (m/s)'
                CALL ctrsetup(upinc,upminc,upmaxc,                      &
                     upovr,uphlf,upzro,upcol1,upcol2,'upplot      ')
                CALL ctr3d( tem9,xc,yc,zps3d,                           &
                     x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                  &
                    nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,       &
                    label(1:10),time,slicopt,kslice,jslice,islice,      &
                    n,xp,yp,b1,b2,zs2,                                  &
                    runname,1.0,tem1,tem2,tem3,                         &
                    tem4,tem5,tem6,hterain,upplot)
              END IF

            END IF

            IF(ipriority(50) == ip ) THEN

              IF(vpplot > 0) THEN
                onvf = 0
                CALL avgy(vprt , onvf,                                  &
                    nx,ny,nz, 1,nx-1, 1,ny-1, 1,nz-1, tem9)
              END IF

              IF(vpplot == 1 .OR.  vpplot == 2 .OR.  vpplot == 4 .OR.   &
                    vpplot == 5 ) THEN
                CALL get_mulovrlay('vpplot',6,ovrmul_num,ovrmulname,sovrlay)
                label = 'vprt (m/s)'
                CALL ctrsetup(vpinc,vpminc,vpmaxc,                      &
                     vpovr,vphlf,vpzro,vpcol1,vpcol2,'vpplot      ')
                CALL ctr3d( tem9 ,xc,yc,zps3d,                          &
                     x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                  &
                    nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,       &
                    label(1:10),time,slicopt,kslice,jslice,islice,      &
                    n,xp,yp,b1,b2,zs2,                                  &
                    runname,1.0,tem1,tem2,tem3,                         &
                    tem4,tem5,tem6,hterain,vpplot)
              END IF
            END IF

            IF(ipriority(51) == ip ) THEN

              IF(wpplot > 0) THEN
                onvf = 0
                CALL avgz(wprt , onvf,                                  &
                    nx,ny,nz, 1,nx-1, 1,ny-1, 1,nz-1, tem9)
              END IF

              IF(wpplot == 1 .OR. wpplot == 2.OR. wpplot == 4 .OR.      &
                    wpplot == 5) THEN
                CALL get_mulovrlay('wpplot',6,ovrmul_num,ovrmulname,sovrlay)
                label = 'wprt (m/s)'
                CALL ctrsetup(wpinc,wpminc,wpmaxc,                      &
                     wpovr,wphlf,wpzro,wpcol1,wpcol2,'wpplot      ')
                CALL ctr3d( tem9,xc,yc,zps3d,                           &
                     x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                  &
                    nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,       &
                    label(1:10),time,slicopt,kslice,jslice,islice,      &
                    n,xp,yp,b1,b2,zs2,                                  &
                    runname,1.0,tem1,tem2,tem3,                         &
                    tem4,tem5,tem6,hterain,wpplot)
              END IF

            END IF
!
!-----------------------------------------------------------------------
!
!  Plot perturbation scalars, calculating and storing
!  perturbations in tem4, where necessary
!
!-----------------------------------------------------------------------
!


            IF(ipriority(52) == ip ) THEN
              IF(ptpplt == 1 .OR. ptpplt == 2.OR. ptpplt == 4 .OR.      &
                    ptpplt == 5 ) THEN

              if(.true.) then  ! plot ptprt or buoyancy including water loading

                CALL get_mulovrlay('ptpplt',6,ovrmul_num,ovrmulname,sovrlay)
                label = 'ptprt (K)'
                CALL ctrsetup(ptpinc,ptpminc,ptpmaxc,                   &
                     ptpovr,ptphlf,ptpzro,ptpcol1,ptpcol2,'ptpplot     ')

                CALL ctr3d(ptprt,xc,yc,zps3d,                           &
                     x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                  &
                    nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,       &
                    label(1:10),time,slicopt,kslice,jslice,islice,      &
                    n,xp,yp,b1,b2,zs2,                                  &
                    runname,1.0,tem1,tem2,tem3,                         &
                    tem4,tem5,tem6,hterain,ptpplt)
              else
                call BUOYCY_plt(nx,ny,nz,ptprt,pprt,qv,qscalar,         &
                                ptbar,pbar,rhobar,qvbar, tem6, tem1)

                do k=1,nz-1
                  do j=1,ny-1
                    do i=1,nx-1
                      tem6(i,j,k) = tem6(i,j,k)/(rhobar(i,j,k)*g) *ptbar(i,j,k)
                    enddo
                  enddo
                enddo

                CALL get_mulovrlay('Buoyancy',6,ovrmul_num,ovrmulname,sovrlay)
                label = 'T-Equivalent Buoyancy (K)'
                CALL ctrsetup(ptpinc,ptpminc,ptpmaxc, &
                     ptpovr,ptphlf,ptpzro,ptpcol1,ptpcol2,'ptpplot     ')

                CALL ctr3d(tem6,xc,yc,zps3d,x1,x2,dx,y1,y2,dy,z1,z2,dz, &
                     nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,      &
                     label(1:17),time,slicopt,kslice,jslice,islice,     &
                     n,xp,yp,b1,b2,zs2,                                 &
                     runname,1.0,tem1,tem2,tem3,                        &
                     tem4,tem5,tem6,hterain,ptpplt)
              endif

              END IF

            END IF

            IF(ipriority(53) == ip ) THEN
              IF(ppplot == 1 .OR. ppplot == 2 .OR. ppplot == 4 .OR.     &
                    ppplot == 5 ) THEN
                CALL get_mulovrlay('ppplot',6,ovrmul_num,ovrmulname,sovrlay)
                label = 'pprt (Pa)'
                CALL ctrsetup(ppinc,ppminc,ppmaxc,                      &
                     ppovr,pphlf,ppzro,ppcol1,ppcol2,'ppplot      ')
                CALL ctr3d(pprt ,xc,yc,zps3d,                           &
                     x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                  &
                    nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,       &
                    label(1:10),time,slicopt,kslice,jslice,islice,      &
                    n,xp,yp,b1,b2,zs2,                                  &
                    runname,1.0,tem1,tem2,tem3,                         &
                    tem4,tem5,tem6,hterain,ppplot)
              END IF
            END IF

            IF(ipriority(54) == ip ) THEN
              IF(qvpplt == 1 .OR. qvpplt == 2 .OR. qvpplt == 4 .OR.     &
                    qvpplt == 5 ) THEN
                CALL get_mulovrlay('qvpplt',6,ovrmul_num,ovrmulname,sovrlay)
                label = 'qvprt (g/kg)'
                CALL ctrsetup(qvpinc,qvpminc,qvpmaxc,                   &
                     qvpovr,qvphlf,qvpzro,qvpcol1,qvpcol2,'qvpplt      ')
                CALL ctr3d(qvprt,xc,yc,zps3d,                           &
                     x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                  &
                    nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,       &
                    label(1:12),time,slicopt,kslice,jslice,islice,      &
                    n,xp,yp,b1,b2,zs2,                                  &
                    runname,1000.0,tem1,tem2,tem3,                      &
                    tem4,tem5,tem6,hterain,qvpplt)
              END IF
            END IF

            IF(ipriority(55) == ip ) THEN

              IF(vorpplt > 0) THEN
                CALL cal_vorp(tem9,u,v,x,y,nx,ny,nz,tem1)
              END IF

              IF(vorpplt == 1 .OR. vorpplt == 2 .OR. vorpplt == 4 .OR.  &
                    vorpplt == 5 ) THEN
                CALL get_mulovrlay('vorpplt',7,ovrmul_num,ovrmulname,   &
                    sovrlay)
                label = 'Vort*10^5 (1/s)'
                CALL ctrsetup(vorpinc,vorpminc,vorpmaxc,                &
                    vorpovr,vorphlf,vorpzro,vorpcol1,vorpcol2,'vorpplt     ')
                CALL ctr3d(tem9,xc,yc,zps3d,                            &
                     x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                  &
                    nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,       &
                    label(1:15),time,slicopt,kslice,jslice,islice,      &
                    n,xp,yp,b1,b2,zs2,                                  &
                    runname,1.0,tem1,tem2,tem3,tem4,                    &
                    tem5,tem6,hterain,vorpplt)
              END IF

            END IF

            IF(ipriority(56) == ip ) THEN

              IF(divpplt > 0) THEN
                CALL cal_div(tem9,u,v,x,y,nx,ny,nz,tem1)
              END IF

              IF( divpplt == 1 .OR. divpplt == 2  .OR. divpplt == 4 .OR. &
                    divpplt == 5 ) THEN

                label = 'Div.*1000 (1/s)'
                CALL ctrsetup(divpinc,divpminc,divpmaxc,                &
                    divpovr,divphlf,divpzro,divpcol1,divpcol2,'divpplt     ')
                CALL ctr3d(tem9,xc,yc,zps3d,                            &
                     x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                  &
                    nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,       &
                    label(1:15),time,slicopt,kslice,jslice,islice,      &
                    n,xp,yp,b1,b2,zs2,                                  &
                    runname,1.0,tem1,tem2,tem3,tem4,                    &
                    tem5,tem6,hterain,divpplt)
              END IF

            END IF
!
!-----------------------------------------------------------------------
!
!   Plot divergence of moist (qv*u )
!
!-----------------------------------------------------------------------
!
            IF(ipriority(57) == ip ) THEN

              IF(divqplt > 0) THEN
                CALL cal_divq(tem9,u,v,qv,x,y,nx,ny,nz,tem1)
              END IF

              IF (divqplt == 1 .OR. divqplt == 2 .OR. divqplt == 4 .OR. &
                    divqplt == 5 ) THEN
                CALL get_mulovrlay('divqplt',7,ovrmul_num,ovrmulname,   &
                    sovrlay)
                label = 'Moist Conv.*1000. (g/kg/s)'
                CALL ctrsetup(divqinc,divqminc,divqmaxc,divqovr,        &
                       divqhlf,divqzro,divqcol1,divqcol2,'divqplt     ')
                CALL ctr3d(tem9,xc,yc,zps3d,                            &
                    x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                   &
                    nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,       &
                    label(1:26),time,slicopt,kslice,jslice,islice,      &
                    n,xp,yp,b1,b2,zs2,                                  &
                    runname,1.0,tem1,tem2,tem3,tem4,                    &
                    tem5,tem6,hterain,divqplt)
              END IF

            END IF
!
!-----------------------------------------------------------------------
!
!  Calculate and plot perturbation wind components
!
!-----------------------------------------------------------------------
!
            IF(ipriority(59) == ip ) THEN
              IF(vtpplt == 1) THEN
                CALL get_mulovrlay('vtpplt',6,ovrmul_num,ovrmulname,sovrlay)

                CALL vtrunt ( vtpunit )
                CALL overlay(vtpovr )
                CALL ctrcol(vtpcol1,vtpcol2)
                CALL varplt( 'vtpplt      ' )
                CALL ctrvtr( vtpunits, vtptype )

                CALL cal_vtp(tem7,tem8,tem9,uprt,vprt,wprt,nx,ny,nz,    &
                     vtpunits,label,length)

                CALL vtr3d(tem7, tem8, tem9, xc,yc,zps3d,               &
                    x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                   &
                    nx,ibgnl2,iendl,ist,ny,jbgnl2,jendl,jst,            &
                    nz,kbgn,kend,kst,kslice,jslice,islice,              &
                    label(1:length),time, runname,1.0,                  &
                    slicopt,n,xp,yp,zs2,u1,v1,u2,v2,w2,                 &
                    tem1,tem2,tem3,tem4,tem5,tem6,hterain)

              END IF
            END IF
!
!-----------------------------------------------------------------------
!
!   Plot streamline field
!
!-----------------------------------------------------------------------
!
            IF(ipriority(60) == ip ) THEN
              IF(vtrstrm == 1) THEN
                CALL get_mulovrlay('vtrstrm',7,ovrmul_num,ovrmulname,   &
                    sovrlay)

                CALL cal_vtrstrm(tem7,tem8,tem9,u,v,w,nx,ny,nz,aspratio)

                CALL overlay( vtrstmovr )
                CALL ctrcol(vtrstmcol1,vtrstmcol2)
                CALL varplt( 'vtrstrm     ' )

                CALL strm3d( tem7 , tem8 , tem9, xc,yc,zps3d,           &
                    x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm, nx,ibgn,iend,ist, &
                    ny,jbgn,jend,jst, nz,kbgn,kend,kst,                 &
                    kslice,jslice,islice, time, runname,1.0, slicopt,   &
                    n,xp,yp,zs2,u1,v1,u2,v2,w2,                         &
                    tem1,tem2,tem3,tem4,tem5,tem6,hterain)

              END IF
            END IF
!
!-----------------------------------------------------------------------
!
!   Calculate and plot the perturbation of wind streamlines
!
!-----------------------------------------------------------------------
!
            IF(ipriority(61) == ip ) THEN
              IF(vtpstrm == 1) THEN
                CALL get_mulovrlay('vtpstrm',7,ovrmul_num,ovrmulname,   &
                    sovrlay)

                CALL cal_vtpstrm(tem7,tem8,tem9,uprt,vprt,wprt,nx,ny,nz, &
                    aspratio)

                CALL overlay( vtpstmovr )
                CALL varplt( 'vtpstrm     ' )

                CALL strm3d( tem7 , tem8 , tem9, xc,yc,zps3d,           &
                    x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,nx,ibgn,iend,ist,  &
                    ny,jbgn,jend,jst,nz,kbgn,kend,kst,                  &
                    kslice,jslice,islice, time, runname,1.0, slicopt,   &
                    n,xp,yp,zs2,u1,v1,u2,v2,w2,                         &
                    tem1,tem2,tem3,tem4,tem5,tem6,hterain)

              END IF
            END IF
!
!-----------------------------------------------------------------------
!
!   Plot vertical wind shear, SQRT[(du/dz)^2 + (dv/dz)^2].
!
!-----------------------------------------------------------------------
!
!
            IF(ipriority(63) == ip ) THEN

              IF(vsplot > 0) THEN
                CALL cal_vs(tem9,u,v,zp,tem7,tem8,nx,ny,nz)
              END IF

              IF(vsplot == 1 .OR. vsplot == 2 .OR. vsplot == 4 .OR.     &
                    vsplot == 5 ) THEN
                CALL get_mulovrlay('vsplot',6,ovrmul_num,ovrmulname,sovrlay)

                label = 'Vertical wind shear*1000(1/s)'
                CALL ctrsetup(vsinc,vsminc,vsmaxc,                      &
                     vsovr,vshlf,vszro,vscol1,vscol2,'vsplot      ')

                CALL ctr3d( tem9, xc,yc,zps3d,                          &
                     x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                  &
                    nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,       &
                    label(1:30),time,slicopt,kslice,jslice,islice,      &
                    n,xp,yp,b1,b2,zs2,                                  &
                    runname,1.0,tem1,tem2,tem3,                         &
                    tem4,tem5,tem6,hterain,vsplot)
              END IF

            END IF
!
!-----------------------------------------------------------------------
!
!   Plot gradient Richardson Number  g / theta * d(theta)/dz
!                                    -----------------------
!                                          (dV/dz)^2g
!
!-----------------------------------------------------------------------
!
!
            IF(ipriority(64) == ip ) THEN

              IF(gricplt > 0) THEN
                CALL cal_gric(tem9,u,v,zp,pt,tem7,tem8,nx,ny,nz)
              END IF
              IF(gricplt == 1 .OR. gricplt == 2 .OR. gricplt == 4 .OR.  &
                    gricplt == 5 ) THEN
                CALL get_mulovrlay('gricplt',7,ovrmul_num,ovrmulname,   &
                    sovrlay)

                label = 'Richardson Number'

                CALL ctrsetup(gricinc,gricminc,gricmaxc,                &
                    gricovr,grichlf,griczro,griccol1,griccol2,'gricplt     ')

                CALL ctr3d( tem9, xc,yc,zps3d,                          &
                     x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                  &
                    nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,       &
                    label(1:17),time,slicopt,kslice,jslice,islice,      &
                    n,xp,yp,b1,b2,zs2,                                  &
                    runname,1.0,tem1,tem2,tem3,                         &
                    tem4,tem5,tem6,hterain,gricplt)
              END IF

            END IF
!
!-----------------------------------------------------------------------
!
!   Plot absolute vorticity
!
!-----------------------------------------------------------------------
!
            IF(ipriority(65) == ip ) THEN

              IF(avorplt > 0) THEN
                CALL cal_avor(tem9,u,v,x,y,nx,ny,nz,slicopt,flagsin,omega, &
                    sinlat,tem1,tem2,tem3)
              END IF

              IF( avorplt == 1 .OR. avorplt == 2 .OR. avorplt == 3      &
                    .OR. avorplt == 4 .OR. avorplt == 5  ) THEN

                CALL get_mulovrlay('avorplt',7,ovrmul_num,ovrmulname,   &
                    sovrlay)
                label = 'Absolute Vort*10^5 (1/s)'
                CALL ctrsetup(avorinc,avorminc,avormaxc,                &
                    avorovr,avorhlf,avorzro,avorcol1,avorcol2,'avorplt    ')
                CALL ctr3d(tem9,xc,yc,zps3d,                            &
                     x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                  &
                    nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,       &
                    label(1:24),time,slicopt,kslice,jslice,islice,      &
                    n,xp,yp,b1,b2,zs2,                                  &
                    runname,1.0,tem1,tem2,tem3,tem4,                    &
                    tem5,tem6,hterain,avorplt)
              END IF

            END IF
!
!-----------------------------------------------------------------------
!
!   Plot scalars for qtplot ,similar to qwplot, except that the
!   vapor component is included.
!
!-----------------------------------------------------------------------
!
            IF(ipriority(67) == ip ) THEN

              IF(qtplot > 0) THEN
                CALL cal_qt(tem9,qv,qscalar,nx,ny,nz)
              END IF

              IF( qtplot == 1 .OR. qtplot == 2 .OR. qtplot == 4         &
                    .OR. qtplot == 5 ) THEN

                CALL get_mulovrlay('qtplot',6,ovrmul_num,ovrmulname,sovrlay)
                label = 'Total water & vapor (g/kg)'
                CALL ctrsetup(qtinc,qtminc,qtmaxc,                      &
                     qtovr,qthlf,qtzro,qtcol1,qtcol2,'qtplot      ')

                CALL ctr3d( tem9,xc,yc,zps3d,                           &
                     x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                  &
                    nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,       &
                    label(1:26),time,slicopt,kslice,jslice,islice,      &
                    n,xp,yp,b1,b2,zs2,                                  &
                    runname,1000.0,tem1,tem2,tem3,                      &
                    tem4,tem5,tem6,hterain,qtplot)
              END IF

            END IF
!
!-----------------------------------------------------------------------
!
!   Plot  RHI relative humidity with ice phase
!
!-----------------------------------------------------------------------
!

            IF(ipriority(68) == ip) THEN

              IF(rhiplot > 0) THEN
                CALL cal_rhi(tem9,pt, pprt ,pbar,qv,tem1,tem2,nx,ny,nz)
              END IF

              IF(rhiplot == 1 .OR. rhiplot == 2 .OR. rhiplot == 4 .OR.  &
                    rhiplot == 5 ) THEN

                CALL get_mulovrlay('rhiplot',7,ovrmul_num,ovrmulname,   &
                    sovrlay)
                label = 'RHI'
                CALL ctrsetup(rhiinc,rhiminc,rhimaxc,                   &
                     rhiovr,rhihlf,rhizro,rhicol1,rhicol2,'rhiplot     ')
                CALL ctr3d( tem9,xc,yc,zps3d,                           &
                     x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                  &
                    nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,       &
                    label(1:3),time,slicopt,kslice,jslice,islice,       &
                    n,xp,yp,b1,b2,zs2,                                  &
                    runname,1.0,tem1,tem2,tem3,                         &
                    tem4,tem5,tem6,hterain,rhiplot)
              END IF

            END IF

!
!-----------------------------------------------------------------------
!
!   Plot vector field
!
!-----------------------------------------------------------------------
!
            CALL varplt( 'vtrplt      ' )


            IF(varname(1:6) == ovrname(1:6) .AND. sovrlay == 1) plotovr=.true.

            IF( (ipriority(58) == ip .AND. vtrplt == 0)                  &
                .AND. (.NOT.plotovr)) GOTO 1100

            IF( ipriority(58) == ip .OR. plotovr) THEN

              IF( vtrplt == 1 .OR.  (ovrlaymulopt == 0 .OR. plotovr .OR. &
                    (ovrlaymulopt == 1 .AND. .NOT.plotovr)) )THEN

                CALL vtrunt ( vtrunit )
                IF(plotovr) CALL overlay ( 1 )
                IF(.NOT.plotovr)  CALL overlay( vtrovr )
                CALL ctrcol(vtrcol1,vtrcol2)
                CALL varplt( 'vtrplt      ' )
                CALL ctrvtr( vtrunits, vtrtype )
!
                IF(nobs > 0 .AND. slicopt == 1 .AND. kslice == 2) THEN
                  CALL cal_vtrobs(dd,ff,drot, vtrunits)
                END IF

                CALL cal_vtr(tem7,tem8,tem9,u,v,w,nx,ny,nz,vtrunits,    &
                    label,length)

                CALL vtr3d( tem7,tem8,tem9,xc,yc,zps3d,                 &
                    x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                   &
                    nx,ibgnl2,iendl,ist, ny,jbgnl2,jendl,jst,           &
                    nz,kbgn,kend,kst,kslice,jslice,islice,              &
                    label(1:length),time, runname,1.0,                  &
                    slicopt,n,xp,yp,zs2,u1,v1,u2,v2,w2,                 &
                    tem1,tem2,tem3,tem4,tem5,tem6,hterain)

              END IF

              IF(plotovr) THEN
                plotovr=.false.
                sovrlay = 0
              END IF

            END IF

            1100 CONTINUE
!
!-----------------------------------------------------------------------
!
!   Plot (u,v) in cross-section .
!
!-----------------------------------------------------------------------
!
            CALL varplt( 'xuvplt     ' )
            IF(varname(1:6) == ovrname(1:6) .AND. sovrlay == 1) plotovr=.true.
            IF( ipriority(66) == ip .OR. plotovr) THEN
              IF( xuvplt == 1 .AND. (ovrlaymulopt == 0 .OR. plotovr .OR. &
                    (ovrlaymulopt == 1 .AND. .NOT.plotovr) ) )THEN
                CALL vtrunt ( xuvunit )
                IF(plotovr) CALL overlay ( 1 )
                IF(.NOT.plotovr)  CALL overlay( xuvovr )
                CALL ctrcol(xuvcol1,xuvcol2)
                CALL varplt( 'xuvplt     ' )
                CALL ctrvtr( xuvunits, xuvtype )

                CALL cal_xuv(tem7,tem8,tem9,u,v,w,nx,ny,nz,xuvunits,label, &
                    length)

                IF(slicopt == 2 .OR. slicopt == 3 .OR. slicopt == 5) THEN

                  CALL vtr3d( tem7,tem8,tem9,xc,yc,zps3d,                 &
                    x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm, nx,ibgnl2,iendl,ist,&
                    ny,jbgnl2,jendl,jst, nz,kbgn,kend,kst,                &
                    kslice,jslice,islice, label(1:length),time, runname,1.0, &
                    slicopt,n,xp,yp,zs2,u1,v1,u2,v2,w2,                   &
                    tem1,tem2,tem3,tem4,tem5,tem6,hterain)
                END IF
              END IF
              IF(plotovr) THEN
                plotovr=.false.
                sovrlay = 0
              END IF
            END IF

!
!-----------------------------------------------------------------------
!
!   Plot Mongtgomery Streamfunction (10m) on constant theta surfaces
!
!-----------------------------------------------------------------------
!
          IF(slicopt == 7) THEN
            IF(ipriority(127) == ip ) THEN
              IF(msfplt == 1 .OR. msfplt == 2 .OR. msfplt == 4 .OR. msfplt == 5) THEN
                CALL get_mulovrlay('msfplt',6,ovrmul_num,ovrmulname,sovrlay)

                DO k=1,nz-1
                  DO j=1,ny-1
                    DO i=1,nx-1
                      tem9(i,j,k)=(cp*tz(i,j,k)+g*zpc(i,j,k)*1000.0)/g*0.1
                    END DO
                  END DO
                END DO

                CALL ctrsetup(msfinc,msfminc,msfmaxc,msfovr,            &
                     msfhlf,msfzro,msfcol1,msfcol2,'msfplt      ')
                label = 'MSF (10m)'


                CALL ctr3d(tem9,xc,yc,zps3d,                            &
                     x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                  &
                     nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,      &
                     label(1:9),                                        &
                     time,slicopt, kslice,jslice,islice,                &
                     n,xp,yp,b1,b2,zs2,                                 &
                     runname,1.0,tem1,tem2,tem3,                        &
                     tem4,tem5,tem6,hterain,msfplt)

              END IF
            END IF

            IF(ipriority(128) == ip ) THEN
              IF( ipvplt == 1 .OR.  ipvplt == 2 .OR.  ipvplt == 4 .OR.  &
                    ipvplt == 5 ) THEN

                DO k=1,nz-1
                  DO j=1,ny-1
                    DO i=1,nx-1
                      tem1(i,j,k)=(u(i+1,j,k)+u(i,j,k))*0.5
                      tem2(i,j,k)=(v(i,j+1,k)+v(i,j,k))*0.5
                    END DO
                  END DO
                END DO

                DO k=2,nz-2
                  DO j=1,ny-1
                    DO i=1,nx-1
                      tem3(i,j,k)=-g*(zps3d(i,j,k+1)-zps3d(i,j,k-1))/   &
                          (pbar(i,j,k+1)+pprt(i,j,k+1)-pbar(i,j,k-1)-pprt(i,j,k-1))
                    END DO
                  END DO
                END DO

                CALL hintrp1(nx,ny,nz,2,nz-2,tem1,zps3d,zlevel,tem9(1,1,1))
                CALL hintrp1(nx,ny,nz,2,nz-2,tem2,zps3d,zlevel,tem9(1,1,2))
                CALL hintrp1(nx,ny,nz,2,nz-2,tem3,zps3d,zlevel,tem9(1,1,3))

                DO j=2,ny-2
                  DO i=2,nx-2
                    IF(ABS(tem9(i+1,j,2)+9999.0) < 0.1.OR.              &
                          ABS(tem9(i-1,j,2)+9999.0) < 0.1.OR.           &
                          ABS(tem9(i,j+1,1)+9999.0) < 0.1.OR.           &
                          ABS(tem9(i,j-1,1)+9999.0) < 0.1.OR.           &
                          ABS(tem9(i,j,3  )+9999.0) < 0.1) THEN
                      tem9(i,j,4)= -9999.0
                      relvort = -9999.0
                    ELSE
                      relvort =  (tem9(i+1,j,2)-tem9(i-1,j,2))*dxinv    &
                                -(tem9(i,j+1,1)-tem9(i,j-1,1))*dyinv
                      tem9(i,j,4)=(relvort                              &
                                  +fcorio(i,j))*tem9(i,j,3)*1.0E6
                    END IF
                  END DO
                END DO
                DO j=2,ny-2
                  tem9(   1,j,4)=tem9(   2,j,4)
                  tem9(nx-1,j,4)=tem9(nx-2,j,4)
                END DO
                DO j=1,nx-1
                  tem9(i,   1,4)=tem9(i,   2,4)
                  tem9(i,ny-1,4)=tem9(i,ny-2,4)
                END DO

                IF(mp_opt > 0) THEN
                  CALL mpsendrecv2dew(tem9,nx,ny,nz,1,1,0,tem4)
                  CALL mpsendrecv2dns(tem9,nx,ny,nz,1,1,0,tem5)
                END IF

                CALL get_mulovrlay('ipvplt',6,ovrmul_num,ovrmulname,sovrlay)

                WRITE(label,'(''IPV (PVU) AT THETA='',F5.1)') zlevel

                CALL ctrsetup(ipvinc,ipvminc,ipvmaxc,                   &
                     ipvovr,ipvhlf,ipvzro,ipvcol1,ipvcol2,'ipvplt     ')

                CALL ctrsfc(tem9(1,1,4),xc(1,1,2),yc(1,1,2),            &
                     x1,x2,dxkm,y1,y2,dykm,                             &
                     nx,ibgnl,iendl, ny,jbgnl,jendl,                    &
                     label(1:24),time, runname,1.0 ,tem1,tem2,tem3,     &
                     tem4,tem5,hterain,slicopt,ipvplt)
                       ! size of tem5 must be >= 6*nx*ny

              END IF

            END IF

          END IF ! For slicopt=7 only

        END IF ! For slicopt=1 - 7
!
!-----------------------------------------------------------------------
!
!  Start plotting of 2-D surface fields...
!
!-----------------------------------------------------------------------
!
        IF (slicopt == 8 .OR. slicopt == 1) THEN

          IF (slicopt == 8 ) THEN  ! 2-D surface field
            aspect=xr/yr
            IF(yxstrch /= 0.0) THEN
              aspratio=yxstrch
            ELSE
              aspratio= aspect
            END IF
          END IF

          IF(ipriority(81) == ip) THEN
            IF( (trnplt == 1 .OR. trnplt == 2 .OR. trnplt == 4 .OR.     &
                   trnplt == 5 )  .AND. sigplt(81) == 0 ) THEN
              CALL get_mulovrlay('trnplt',6,ovrmul_num,ovrmulname,sovrlay)
              label = 'Terrain height (m)'
              time0 = 0.0
              CALL ctrsetup(trninc,trnminc,trnmaxc,                     &
                   trnovr,trnhlf,trnzro,trncol1,trncol2,'trnplt      ')
              CALL ctrsfc(zp(1,1,2),xc(1,1,2),yc(1,1,2),                &
                   x1,x2,dxkm,y1,y2,dykm,                               &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  label(1:18),time0, runname,1.0 ,tem1,tem2,tem3,       &
                  tem4,tem5,hterain,slicopt,trnplt)
                       ! size of tem5 must be >= 6*nx*ny
              sigplt(81) = 1
            END IF
          END IF

          IF(ipriority(82) == ip ) THEN
            IF( (wetcanplt == 1 .OR. wetcanplt == 2 .OR. wetcanplt == 4 &
                   .OR. wetcanplt == 5 ) .AND. sigplt(82) == 0 ) THEN
              CALL get_mulovrlay('wetcanplt',9,ovrmul_num,ovrmulname,   &
                                 sovrlay)
              label = 'Canopy water amount'
              CALL ctrsetup(wcpinc,wcpminc,wcpmaxc,                     &
                   wcpovr,wcphlf,wcpzro,wcpcol1,wcpcol2,'wetcanplt   ')

              CALL ctrsfc(wetcanp,xc(1,1,2),yc(1,1,2),                  &
                   x1,x2,dxkm,y1,y2,dykm,                               &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  label(1:19),time, runname,1.0 ,tem1,tem2,tem3,        &
                  tem4,tem5,hterain,slicopt,wetcanplt)
                       ! size of tem5 must be >= 6*nx*ny
              sigplt(82) = 1
            END IF
          END IF

          IF(ipriority(83) == ip) THEN
            IF( (raincplt == 1 .OR.raincplt == 2 .OR.raincplt == 4 .OR. &
                   raincplt == 5 ) .AND. sigplt(83) == 0 ) THEN
              CALL get_mulovrlay('raincplt',8,ovrmul_num,ovrmulname,sovrlay)
              IF(racunit == 0) THEN
                label = 'Cumulus Rainfall (mm)'
                DO j = 1,ny
                  DO i=1,nx
                    tem9(i,j,1) = rainc(i,j)
                  END DO
                END DO
              ELSE IF( racunit == 1) THEN
                label = 'Cumulus Rainfall (in)'    ! unit is inch
                DO j = 1,ny
                  DO i=1,nx
                    tem9(i,j,1) = rainc(i,j)*0.039370079   ! (1/25.4)
                  END DO
                END DO
              END IF
              CALL ctrsetup(raincinc,raincminc,raincmaxc,               &
                   racovr,rachlf,raczro,raccol1,raccol2,'raincplt    ')

              CALL ctrsfc(tem9(1,1,1),xc(1,1,2),yc(1,1,2),              &
                  x1,x2,dxkm,y1,y2,dykm, nx,ibgnl,iendl, ny,jbgnl,jendl,&
                  label(1:28),time, runname,1.0 ,tem1,tem2,tem3,        &
                  tem4,tem5,hterain,slicopt,raincplt)
                       ! size of tem5 must be >= 6*nx*ny
              sigplt(83) = 1
            END IF
          END IF


          IF(ipriority(84) == ip) THEN
            IF( (raingplt == 1 .OR. raingplt == 2 .OR. raingplt == 4 .OR. &
                   raingplt == 5 )  .AND. sigplt(84) == 0 ) THEN
              CALL get_mulovrlay('raingplt',8,ovrmul_num,ovrmulname,sovrlay)

              IF(ragunit == 0 ) THEN
                label = 'Gridscale Rainfall (mm) '
                DO j = 1,ny
                  DO i=1,nx
                    tem9(i,j,1) = raing(i,j)
                  END DO
                END DO
              ELSE IF( ragunit == 1) THEN
                label = 'Gridscale Rainfall (in) '
                DO j = 1,ny
                  DO i=1,nx
                    tem9(i,j,1) = raing(i,j)*0.039370079   ! (1/25.4)
                  END DO
                END DO
              END IF

              CALL ctrsetup(rainginc,raingminc,raingmaxc,               &
                   ragovr,raghlf,ragzro,ragcol1,ragcol2,'raingplt    ')

              CALL ctrsfc(tem9(1,1,1),xc(1,1,2),yc(1,1,2),              &
                  x1,x2,dxkm,y1,y2,dykm, nx,ibgnl,iendl, ny,jbgnl,jendl,&
                  label(1:30),time, runname,1.0 ,tem1,tem2,tem3,        &
                  tem4,tem5,hterain,slicopt,raingplt)
                       ! size of tem5 must be >= 6*nx*ny
              sigplt(84) = 1
            END IF
          END IF


          IF(ipriority(85) == ip) THEN
            IF( (raintplt == 1 .OR.raintplt == 2 .OR.raintplt == 4      &
                  .OR. raintplt == 5 )                                  &
                  .AND. sigplt(85) == 0 ) THEN
              CALL get_mulovrlay('raintplt',8,ovrmul_num,ovrmulname,sovrlay)
              IF(ratunit == 0 ) THEN
                label = 'Total Rainfall (mm)'
                DO j = 1,ny
                  DO i=1,nx
                    tem9(i,j,1) = rainc(i,j) + raing(i,j)
                  END DO
                END DO
              ELSE IF (ratunit == 1) THEN
                label = 'Total Rainfall (in)'
                DO j = 1,ny
                  DO i=1,nx
                    tem9(i,j,1) = (rainc(i,j)+raing(i,j))*0.039370079 ! (1/25.4)
                  END DO
                END DO
              END IF
              CALL ctrsetup(raintinc,raintminc,raintmaxc,               &
                   ratovr,rathlf,ratzro,ratcol1,ratcol2,'raintplt    ')

              CALL ctrsfc(tem9,xc(1,1,2),yc(1,1,2),                     &
                   x1,x2,dxkm,y1,y2,dykm,                               &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  label(1:22),time, runname,1.0 ,tem1,tem2,tem3,        &
                  tem4,tem5,hterain,slicopt,raintplt)
                       ! size of tem5 must be >= 6*nx*ny
              sigplt(85) = 1
            END IF
          END IF

          !
          ! Added Accumulated rainfall plots for priority 111, 112, 113
          !
          IF(ipriority(131) == ip) THEN

            IF( (rainicplt == 1 .OR. rainicplt == 2 .OR.                &
                 rainicplt == 4 .OR. rainicplt == 5 )      .AND.        &
                sigplt(131) == 0 ) THEN
              CALL get_mulovrlay('rainicplt',8,ovrmul_num,ovrmulname,   &
                                 sovrlay)

              DO j=1,ny
                DO i=1,nx
                  acc_rainc(i,j,2)=rainc(i,j)
                END DO
              END DO
              timediff=time2c-time1c           ! time in sec

              IF(raicunit == 0) THEN
                WRITE(label,'(f8.1,a)') timediff,                       &
                         ' s  Accumulated Cumulus Rainfall (mm)'
                DO j = 1,ny
                  DO i=1,nx
                    tem9(i,j,1) = acc_rainc(i,j,2)-acc_rainc(i,j,1)
                    IF (tem9(i,j,1) < 0.001) tem9(i,j,1)=0.0
                  END DO
                END DO

              ELSE IF( raicunit == 1) THEN     ! unit is inch
                WRITE(label,'(f8.1,a)') timediff,                       &
                         ' s  Accumulated Cumulus Rainfall (in)'
                DO j = 1,ny
                  DO i=1,nx
                    tem9(i,j,1)=(acc_rainc(i,j,2)-acc_rainc(i,j,1))     &
                                * 0.039370079  ! (1/25.4)
                    IF (tem9(i,j,1) < 0.0001) tem9(i,j,1)=0.0
                  END DO
                END DO
              END IF

              IF (timediff > 0) THEN
                CALL ctrsetup(rainicinc,rainicminc,rainicmaxc,            &
                           raicovr,raichlf,raiczro,raiccol1,raiccol2,     &
                           'rainicplt   ')

                CALL ctrsfc(tem9(1,1,1),xc(1,1,2),yc(1,1,2),x1,x2,dxkm,   &
                           y1,y2,dykm, nx,ibgnl,iendl, ny,jbgnl,jendl,    &
                           TRIM(label),time, runname,1.0 ,tem1,tem2,tem3, &
                           tem4,tem5,hterain,slicopt,rainicplt)
                           ! size of tem5 must be >= 6*nx*ny
              END IF
              sigplt(131) = 1

              DO j=1,ny
                DO i=1,nx
                  acc_rainc(i,j,1)=acc_rainc(i,j,2)
                END DO
              END DO
              time1c=time2c
            END IF
          END IF

          IF(ipriority(132) == ip) THEN

            IF( (rainigplt == 1 .OR. rainigplt == 2 .OR.                &
                 rainigplt == 4 .OR. rainigplt == 5 )      .AND.        &
                sigplt(132) == 0 ) THEN
              CALL get_mulovrlay('rainigplt',8,ovrmul_num,ovrmulname,   &
                                 sovrlay)

              DO j=1,ny
                DO i=1,nx
                  acc_raing(i,j,2)=raing(i,j)
                END DO
              END DO
              timediff=time2g-time1g           ! time in sec

              IF(raigunit == 0) THEN
                WRITE(label,'(f8.1,a)') timediff,                       &
                         ' s  Accumulated Gridscale Rainfall (mm)'
                DO j = 1,ny
                  DO i=1,nx
                    tem9(i,j,1) = acc_raing(i,j,2)-acc_raing(i,j,1)
                    IF (tem9(i,j,1) < 0.001) tem9(i,j,1)=0.0
                  END DO
                END DO

              ELSE IF( raigunit == 1) THEN     ! unit is inch
                WRITE(label,'(f8.1,a)') timediff,                       &
                         ' s  Accumulated Gridscale Rainfall (in)'
                DO j = 1,ny
                  DO i=1,nx
                    tem9(i,j,1)=(acc_raing(i,j,2)-acc_raing(i,j,1))     &
                                * 0.039370079  ! (1/25.4)
                    IF (tem9(i,j,1) < 0.0001) tem9(i,j,1)=0.0
                  END DO
                END DO
              END IF

              IF (timediff > 0) THEN
                CALL ctrsetup(rainiginc,rainigminc,rainigmaxc,          &
                           raigovr,raighlf,raigzro,raigcol1,raigcol2,   &
                           'rainigplt   ')

                CALL ctrsfc(tem9(1,1,1),xc(1,1,2),yc(1,1,2),x1,x2,dxkm,   &
                           y1,y2,dykm, nx,ibgnl,iendl, ny,jbgnl,jendl,    &
                           TRIM(label),time, runname,1.0 ,tem1,tem2,tem3, &
                           tem4,tem5,hterain,slicopt,rainigplt)
                           ! size of tem5 must be >= 6*nx*ny
              END IF
              sigplt(132) = 1

              DO j=1,ny
                DO i=1,nx
                  acc_raing(i,j,1)=acc_raing(i,j,2)
                END DO
              END DO
              time1g=time2g
            END IF
          END IF

          IF(ipriority(133) == ip) THEN

            IF( (rainitplt == 1 .OR. rainitplt == 2 .OR.                &
                 rainitplt == 4 .OR. rainitplt == 5 )      .AND.        &
                sigplt(133) == 0 ) THEN
              CALL get_mulovrlay('rainitplt',8,ovrmul_num,ovrmulname,   &
                                 sovrlay)

              DO j=1,ny
                DO i=1,nx
                  acc_raint(i,j,2)=rainc(i,j) + raing(i,j)
                END DO
              END DO
              timediff=time2t-time1t           ! time in sec

              IF(raitunit == 0) THEN
                WRITE(label,'(f8.1,a)') timediff,                       &
                         ' s  Accumulated Rainfall (mm)'
                DO j = 1,ny
                  DO i=1,nx
                    tem9(i,j,1) = acc_raint(i,j,2)-acc_raint(i,j,1)
                    IF (tem9(i,j,1) < 0.001) tem9(i,j,1)=0.0
                  END DO
                END DO

              ELSE IF( raitunit == 1) THEN     ! unit is inch
                WRITE(label,'(f8.1,a)') timediff,                       &
                         ' s  Accumulated Rainfall (in)'
                DO j = 1,ny
                  DO i=1,nx
                    tem9(i,j,1)=(acc_raint(i,j,2)-acc_raint(i,j,1))     &
                                * 0.039370079  ! (1/25.4)
                    IF (tem9(i,j,1) < 0.0001) tem9(i,j,1)=0.0
                  END DO
                END DO
              END IF

              IF (timediff > 0.0) THEN  ! skip the first time level
                CALL ctrsetup(rainitinc,rainitminc,rainitmaxc,          &
                           raitovr,raithlf,raitzro,raitcol1,raitcol2,   &
                           'rainitplt   ')

                CALL ctrsfc(tem9(1,1,1),xc(1,1,2),yc(1,1,2),x1,x2,dxkm, &
                           y1,y2,dykm, nx,ibgnl,iendl, ny,jbgnl,jendl,  &
                           TRIM(label),time, runname,1.0 ,tem1,tem2,tem3, &
                           tem4,tem5,hterain,slicopt,rainitplt)
                           ! size of tem5 must be >= 6*nx*ny
              END IF
              sigplt(133) = 1

              DO j=1,ny
                DO i=1,nx
                  acc_raint(i,j,1)=acc_raint(i,j,2)
                END DO
              END DO
              time1t=time2t
            END IF
          END IF

          IF(ipriority(86) == ip) THEN
            IF( (pslplt == 1 .OR. pslplt == 2 .OR. pslplt == 4 .OR.     &
                   pslplt == 5 )                                        &
                  .AND. sigplt(86) == 0 ) THEN
              CALL get_mulovrlay('pslplt',6,ovrmul_num,ovrmulname,sovrlay)
              label = 'Sea Level Pressure (mb)'

              IF(ovrobs == 1 .AND. nobs > 0) THEN
                !
                !  Sea-level pressure
                !
                DO iob=1,nobs
                  IF(pmsl(iob) > 0. ) THEN
                    obs1(iob)=mod(nint(pmsl(iob)),100)
                  !  print *, ' pmsl, obs1: ',pmsl(iob),obs1(iob)
                  !ELSE IF (pstn(iob) > 0. .AND. tob(iob) > -98. ) THEN
                  !  tk=(5.*(tob(iob)-32.)/9.) + 273.15
                  !  obs1(iob)=pstn(iob)*((tk+gamma*elevob(iob))/tk)**ex2
                  !  obs1(iob)=mod(nint(obs1(iob)),100)
                  !  print *, ' pstn, obs1: ',pstn(iob),obs1(iob)
                  !ELSE IF( alt(iob) > 0. ) THEN
                  !  IF ( tob(iob) > -98.) THEN
                  !    tk=(5.*(tob(iob)-32.)/9.) + 273.15
                  !    obs1(iob)=alttostpr(alt(iob),elevob(iob))*    &
                  !                        ((tk+gamma*elevob(iob))/tk)**ex2
                  !    obs1(iob)=mod(nint(obs1(iob)),100)
                  !    print *, ' alt1, obs1: ',alt(iob),obs1(iob)
                  !  ELSE
                  !    obs1(iob)=0.01*alt(iob)
                  !    obs1(iob)=mod(nint(obs1(iob)),100)
                  !    print *, ' alt2, obs1: ',alt(iob),obs1(iob)
                  !  END IF
                  ELSE
                    obs1(iob)=-99.
                  END IF
                END DO
                obsset=1
              END IF

              CALL ctrsetup(pslinc,pslminc,pslmaxc,                     &
                   pslovr,pslhlf,pslzro,pslcol1,pslcol2,'pslplt      ')

              CALL ctrsfc(psl,xc(1,1,2),yc(1,1,2),                      &
                  x1,x2,dxkm,y1,y2,dykm,                                &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  label(1:23),time, runname,1.0 ,tem1,tem2,tem3,        &
                  tem4,tem5,hterain,slicopt,pslplt)
                       ! size of tem5 must be >= 6*nx*ny
              sigplt(86) = 1
            END IF
          END IF


          IF(ipriority(87) == ip) THEN
            IF( (liplt == 1 .OR. liplt == 2 .OR. liplt == 4 .OR.  liplt == 5 ) &
                  .AND. sigplt(87) == 0 ) THEN
              CALL get_mulovrlay('liplt',5,ovrmul_num,ovrmulname,sovrlay)
              label = 'Sfc-based LI (K)'

              CALL ctrsetup(liinc,liminc,limaxc,                        &
                   liovr,lihlf,lizro,licol1,licol2,'liplt       ')

              CALL ctrsfc(li,xc(1,1,2),yc(1,1,2),                       &
                   x1,x2,dxkm,y1,y2,dykm,                               &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  label(1:19),time, runname,1.0 ,tem1,tem2,tem3,        &
                  tem4,tem5,hterain,slicopt,liplt)
                       ! size of tem5 must be >= 6*nx*ny
              sigplt(87) = 1
            END IF
          END IF

          IF(ipriority(88) == ip) THEN
            IF( (capeplt == 1 .OR. capeplt == 2 .OR. capeplt == 4 .OR.  &
                 capeplt == 5) .AND. sigplt(88) == 0 ) THEN
              CALL get_mulovrlay('capeplt',7,ovrmul_num,ovrmulname,sovrlay)
              label = 'Surface CAPE (J/kg)'

              CALL ctrsetup(capeinc,capeminc,capemaxc,                  &
                   capovr,caphlf,capzro,capcol1,capcol2,'capeplt      ')

              CALL ctrsfc(cape,xc(1,1,2),yc(1,1,2),                     &
                   x1,x2,dxkm,y1,y2,dykm,                               &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  label(1:19),time, runname,1.0 ,tem1,tem2,tem3,        &
                  tem4,tem5,hterain,slicopt,capeplt)
                       ! size of tem5 must be >= 6*nx*ny
              sigplt(88) = 1
            END IF
          END IF

          IF(ipriority(89) == ip) THEN
            IF( (cinplt == 1 .OR. cinplt == 2 .OR. cinplt == 4 .OR.     &
                 cinplt == 5 ) .AND. sigplt(89) == 0 ) THEN
              CALL get_mulovrlay('cinplt',6,ovrmul_num,ovrmulname,sovrlay)
              label = 'Surface CIN (J/kg)'

              CALL ctrsetup(cininc,cinminc,cinmaxc,                     &
                   cinovr,cinhlf,cinzro,cincol1,cincol2,'cinplt      ')

              CALL ctrsfc(cin,xc(1,1,2),yc(1,1,2),                      &
                  x1,x2,dxkm,y1,y2,dykm,nx,ibgnl,iendl, ny,jbgnl,jendl, &
                  label(1:18),time, runname,1.0 ,tem1,tem2,tem3,        &
                  tem4,tem5,hterain,slicopt,cinplt)
                       ! size of tem5 must be >= 6*nx*ny
              sigplt(89) = 1
            END IF
          END IF

          IF(ipriority(90) == ip ) THEN
            IF( (thetplt == 1 .OR. thetplt == 2 .OR. thetplt == 4 .OR.  &
                 thetplt == 5 ) .AND. sigplt(90) == 0 ) THEN
              CALL get_mulovrlay('thetplt',7,ovrmul_num,ovrmulname,sovrlay)
              label = 'Surface Theta-e (K)'

              CALL ctrsetup(thetinc,thetminc,thetmaxc,                  &
                   theovr,thehlf,thezro,thecol1,thecol2,'thetplt     ')

              CALL ctrsfc(thet,xc(1,1,2),yc(1,1,2),                     &
                  x1,x2,dxkm,y1,y2,dykm,nx,ibgnl,iendl, ny,jbgnl,jendl, &
                  label(1:19),time, runname,1.0 ,tem1,tem2,tem3,        &
                  tem4,tem5,hterain,slicopt,thetplt)
                       ! size of tem5 must be >= 6*nx*ny
              sigplt(90) = 1
            END IF
          END IF

          IF(ipriority(91) == ip ) THEN
            IF( (heliplt == 1 .OR. heliplt == 2 .OR. heliplt == 4 .OR.  &
                 heliplt == 5 ) .AND. sigplt(91) == 0 ) THEN
              CALL get_mulovrlay('heliplt',7,ovrmul_num,ovrmulname,sovrlay)
              label = '0-3km Helicity (m^2/s^2)'

              !IF(mp_opt >0) THEN
              !  CALL mpsend1dew(heli,nx,ny,1,1,0,mptag1,tem5)
              !  CALL mpsend1dns(heli,nx,ny,1,1,0,mptag2,tem5)
              !
              !  CALL mprecv1dew(heli,nx,ny,1,1,0,mptag1,tem5)
              !  CALL mprecv1dns(heli,nx,ny,1,1,0,mptag2,tem5)
              !END IF

              CALL ctrsetup(heliinc,heliminc,helimaxc,                  &
                   helovr,helhlf,helzro,helcol1,helcol2,'heliplt     ')

              CALL ctrsfc(heli,xc(1,1,2),yc(1,1,2),                     &
                  x1,x2,dxkm,y1,y2,dykm,nx,ibgnl,iendl, ny,jbgnl,jendl, &
                  label(1:24),time, runname,1.0 ,tem1,tem2,tem3,        &
                  tem4,tem5,hterain,slicopt,heliplt)
                       ! size of tem5 must be >= 6*nx*ny
              sigplt(91) = 1
            END IF
          END IF

          IF(ipriority(92) == ip ) THEN
            !print *, ' ip,uhprio,ipriority(92): ',ip,uhprio,ipriority(92)
            IF( (uhplt == 1 .OR. uhplt == 2 .OR. uhplt == 4 .OR.        &
                 uhplt == 5 ) .AND. sigplt(92) == 0 ) THEN
              CALL get_mulovrlay('uhplt',5,ovrmul_num,ovrmulname,sovrlay)
              label = 'Updraft Helic (m^2/s^2)'

              !IF(mp_opt >0) THEN
              !  CALL mpsend1dew(heli,nx,ny,1,1,0,mptag1,tem5)
              !  CALL mpsend1dns(heli,nx,ny,1,1,0,mptag2,tem5)
              !
              !  CALL mprecv1dew(heli,nx,ny,1,1,0,mptag1,tem5)
              !  CALL mprecv1dns(heli,nx,ny,1,1,0,mptag2,tem5)
              !END IF

              CALL ctrsetup(uhinc,uhminc,uhmaxc,                        &
                            uhovr,uhhlf,uhzro,uhcol1,uhcol2,'uhplt     ')

              CALL ctrsfc(uh,xc(1,1,2),yc(1,1,2),                       &
                  x1,x2,dxkm,y1,y2,dykm,nx,ibgnl,iendl, ny,jbgnl,jendl, &
                  label(1:24),time, runname,1.0 ,tem1,tem2,tem3,        &
                  tem4,tem5,hterain,slicopt,uhplt)
                       ! size of tem5 must be >= 6*nx*ny
              sigplt(92) = 1
            END IF
          END IF

          IF(ipriority(93) == ip ) THEN
            IF( (ctcplt == 1 .OR. ctcplt == 2 .OR. ctcplt == 4 .OR.     &
                  ctcplt == 5 ) .AND. sigplt(93) == 0 ) THEN
              CALL get_mulovrlay('ctcplt',6,ovrmul_num,ovrmulname,sovrlay)
              label = 'Convective temp (C)'

              CALL ctrsetup(ctcinc,ctcminc,ctcmaxc,                     &
                   ctcovr,ctchlf,ctczro,ctccol1,ctccol2,'ctcplt      ')

              CALL ctrsfc(ct,xc(1,1,2),yc(1,1,2),                       &
                  x1,x2,dxkm,y1,y2,dykm,                                &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  label(1:19),time, runname,1.0 ,tem1,tem2,tem3,        &
                  tem4,tem5,hterain,slicopt,ctcplt)
                       ! size of tem5 must be >= 6*nx*ny
              sigplt(93) = 1
            END IF
          END IF


          IF(ipriority(94) == ip ) THEN
            IF( (brnplt == 1 .OR. brnplt == 2 .OR. brnplt == 4 .OR.     &
                   brnplt == 5 )                                        &
                  .AND. sigplt(94) == 0 ) THEN
              CALL get_mulovrlay('brnplt',6,ovrmul_num,ovrmulname,sovrlay)
              label = 'BRN Shear Denom (m2/s2)'

              !IF(mp_opt >0) THEN
              !  CALL mpsend1dew(brn,nx,ny,1,1,0,mptag1,tem5)
              !  CALL mpsend1dns(brn,nx,ny,1,1,0,mptag2,tem5)
              !
              !  CALL mprecv1dew(brn,nx,ny,1,1,0,mptag1,tem5)
              !  CALL mprecv1dns(brn,nx,ny,1,1,0,mptag2,tem5)
              !END IF

              CALL ctrsetup(brninc,brnminc,brnmaxc,                     &
                   brnovr,brnhlf,brnzro,brncol1,brncol2,'brnplt      ')

              CALL ctrsfc(brn,xc(1,1,2),yc(1,1,2),                      &
                  x1,x2,dxkm,y1,y2,dykm,                                &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  label(1:23),time, runname,1.0 ,tem1,tem2,tem3,        &
                  tem4,tem5,hterain,slicopt,brnplt)
                       ! size of tem5 must be >= 6*nx*ny
              sigplt(94) = 1
            END IF
          END IF

          IF(ipriority(95) == ip ) THEN
            IF( (brnuplt == 1 .OR. brnuplt == 2 .OR. brnuplt == 4 .OR.  &
                   brnuplt == 5)                                        &
                  .AND. sigplt(95) == 0 ) THEN
              CALL get_mulovrlay('brnuplt',7,ovrmul_num,ovrmulname,sovrlay)
              label = 'Bulk Richardson Shear (1/s)'

              !IF(mp_opt >0) THEN
              !  CALL mpsend1dew(brnu,nx,ny,1,1,0,mptag1,tem5)
              !  CALL mpsend1dns(brnu,nx,ny,1,1,0,mptag2,tem5)
              !
              !  CALL mprecv1dew(brnu,nx,ny,1,1,0,mptag1,tem5)
              !  CALL mprecv1dns(brnu,nx,ny,1,1,0,mptag2,tem5)
              !END IF

              CALL ctrsetup(brnuinc,bruminc,brumaxc,                    &
                  brnuovr,brnuhlf,brnuzro,brnucol1,brnucol2,'brnuplt     ')

              CALL ctrsfc(brnu,xc(1,1,2),yc(1,1,2),                     &
                  x1,x2,dxkm,y1,y2,dykm,                                &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  label(1:27),time, runname,1.0 ,tem1,tem2,tem3,        &
                  tem4,tem5,hterain,slicopt,brnuplt)
              sigplt(95) = 1
            END IF
          END IF

          IF(ipriority(96) == ip ) THEN
            IF( (srlfplt == 1 .OR. srlfplt == 2 .OR. srlfplt == 4 .OR.  &
                   srlfplt == 5 ) .AND. sigplt(96) == 0 ) THEN
              CALL get_mulovrlay('srlfplt',7,ovrmul_num,ovrmulname,sovrlay)
              label = '0-2km StRel Flow (m/s)'

              !IF(mp_opt >0) THEN
              !  CALL mpsend1dew(srlfl,nx,ny,1,1,0,mptag1,tem5)
              !  CALL mpsend1dns(srlfl,nx,ny,1,1,0,mptag2,tem5)
              !
              !  CALL mprecv1dew(srlfl,nx,ny,1,1,0,mptag1,tem5)
              !  CALL mprecv1dns(srlfl,nx,ny,1,1,0,mptag2,tem5)
              !END IF

              CALL ctrsetup(srlfinc,srlminc,srlmaxc,                    &
                  srlfovr,srlfhlf,srlfzro,srlfcol1,srlfcol2,'srlfplt     ')

              CALL ctrsfc(srlfl,xc(1,1,2),yc(1,1,2),                    &
                  x1,x2,dxkm,y1,y2,dykm,                                &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  label(1:22),time, runname,1.0 ,tem1,tem2,tem3,        &
                  tem4,tem5,hterain,slicopt,srlfplt)
              sigplt(96) = 1
            END IF
          END IF

          IF(ipriority(97) == ip ) THEN
            IF( (srmfplt == 1 .OR. srmfplt == 2 .OR. srmfplt == 4 .OR.  &
                   srmfplt == 5)  .AND. sigplt(97) == 0 ) THEN
              CALL get_mulovrlay('srmfplt',7,ovrmul_num,ovrmulname,sovrlay)
              label = '2-9km StRel Flow (m/s)'

              !IF(mp_opt >0) THEN
              !  CALL mpsend1dew(srmfl,nx,ny,1,1,0,mptag1,tem5)
              !  CALL mpsend1dns(srmfl,nx,ny,1,1,0,mptag2,tem5)
              !
              !  CALL mprecv1dew(srmfl,nx,ny,1,1,0,mptag1,tem5)
              !  CALL mprecv1dns(srmfl,nx,ny,1,1,0,mptag2,tem5)
              !END IF

              CALL ctrsetup(srmfinc,srmminc,srmmaxc,                    &
                  srmfovr,srmfhlf,srmfzro,srmfcol1,srmfcol2,'srmfplt     ')

              CALL ctrsfc(srmfl,xc(1,1,2),yc(1,1,2),                    &
                  x1,x2,dxkm,y1,y2,dykm,                                &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  label(1:22),time, runname,1.0 ,tem1,tem2,tem3,        &
                  tem4,tem5,hterain,slicopt,srmfplt)
              sigplt(97) = 1
            END IF
          END IF

          IF(ipriority(98) == ip ) THEN
            IF( (capsplt == 1 .OR. capsplt == 2 .OR. capsplt == 4 .OR.  &
                   capsplt == 5 )   .AND. sigplt(98) == 0 ) THEN
              CALL get_mulovrlay('capsplt',7,ovrmul_num,ovrmulname,sovrlay)
              label = 'Cap Strength (K)'

              CALL ctrsetup(capsinc,capsminc,capsmaxc,                  &
                  capsovr,capshlf,capszro,capscol1,capscol2,'capsplt     ')

              CALL ctrsfc(capst,xc(1,1,2),yc(1,1,2),                    &
                  x1,x2,dxkm,y1,y2,dykm,                                &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  label(1:16),time, runname,1.0 ,tem1,tem2,tem3,        &
                  tem4,tem5,hterain,slicopt,capsplt)
              sigplt(98) = 1
            END IF
          END IF

          IF(ipriority(99) == ip ) THEN
            IF( (blcoplt == 1 .OR. blcoplt == 2 .OR. blcoplt == 4 .OR.  &
                   blcoplt == 5 )    .AND. sigplt(99) == 0 ) THEN
              CALL get_mulovrlay('blcoplt',7,ovrmul_num,ovrmulname,sovrlay)
              !label = '0-2km Layer Conv.*1000 (1/s)'
              label = '0-2km Conv.*1000 (1/s)'

              !IF(mp_opt >0) THEN
              !  CALL mpsend1dew(blcon,nx,ny,1,1,0,mptag1,tem5)
              !  CALL mpsend1dns(blcon,nx,ny,1,1,0,mptag2,tem5)
              !
              !  CALL mprecv1dew(blcon,nx,ny,1,1,0,mptag1,tem5)
              !  CALL mprecv1dns(blcon,nx,ny,1,1,0,mptag2,tem5)
              !END IF

              CALL ctrsetup(blcoinc,blcominc,blcomaxc,                  &
                  blcoovr,blcohlf,blcozro,blcocol1,blcocol2,'blcoplt     ')

              CALL ctrsfc(blcon,xc(1,1,2),yc(1,1,2),                    &
                  x1,x2,dxkm,y1,y2,dykm,                                &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  label(1:22),time, runname,1.0 ,tem1,tem2,tem3,        &
                  tem4,tem5,hterain,slicopt,blcoplt)
              sigplt(99) = 1
            END IF
          END IF

          IF(ipriority(100) == ip .AND. P_QC > 0) THEN
            IF( (viqcplt == 1 .OR. viqcplt == 2 .OR. viqcplt == 4       &
                  .OR. viqcplt == 5 )                                   &
                  .AND. sigplt(100) == 0 ) THEN
              CALL get_mulovrlay('viqcplt',7,ovrmul_num,ovrmulname,sovrlay)
              label = 'Vert. Integrated qc (kg/m2)'

              CALL ctrsetup(viqcinc,viqcminc,viqcmaxc,                  &
                  viqcovr,viqchlf,viqczro,viqccol1,viqccol2,'viqcplt     ')

              CALL cal_viqc(tem7, qscalar(:,:,:,P_QC),rhobar, zp, nx,ny,nz)

              CALL ctrsfc(tem7,xc(1,1,2),yc(1,1,2),                     &
                   x1,x2,dxkm,y1,y2,dykm,                               &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  label(1:27),time, runname,1.0 ,tem1,tem2,tem3,        &
                  tem4,tem5,hterain,slicopt,viqcplt)
              sigplt(100) = 1
            END IF
          END IF

          IF(ipriority(101) == ip .AND. P_QR > 0 ) THEN
            IF( (viqrplt == 1 .OR. viqrplt == 2 .OR. viqrplt == 4       &
                   .OR. viqrplt == 5 )                                  &
                  .AND. sigplt(101) == 0 ) THEN
              CALL get_mulovrlay('viqrplt',7,ovrmul_num,ovrmulname,sovrlay)
              label = 'Vert. Integrated qr (kg/m2)'

              CALL ctrsetup(viqrinc,viqrminc,viqrmaxc,                  &
                  viqrovr,viqrhlf,viqrzro,viqrcol1,viqrcol2,'viqrplt     ')

              CALL cal_viqr(tem7, qscalar(:,:,:,P_QR), rhobar, zp, nx,ny,nz)

              CALL ctrsfc(tem7,xc(1,1,2),yc(1,1,2),                     &
                  x1,x2,dxkm,y1,y2,dykm,                                &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  label(1:27),time, runname,1.0 ,tem1,tem2,tem3,        &
                  tem4,tem5,hterain,slicopt,viqrplt)
              sigplt(101) = 1
            END IF
          END IF

          IF(ipriority(102) == ip .AND. P_QI > 0) THEN
            IF( (viqiplt == 1 .OR. viqiplt == 2 .OR. viqiplt == 4       &
                   .OR. viqiplt == 5 ) .AND. sigplt(102) == 0 ) THEN
              CALL get_mulovrlay('viqiplt',7,ovrmul_num,ovrmulname,sovrlay)
              label = 'Vert. Integrated qi (kg/m2)'

              CALL ctrsetup(viqiinc,viqiminc,viqimaxc,                  &
                  viqiovr,viqihlf,viqizro,viqicol1,viqicol2,'viqiplt     ')

              CALL cal_viqi(tem7, qscalar(:,:,:,P_QI), rhobar, zp, nx,ny,nz)

              CALL ctrsfc(tem7,xc(1,1,2),yc(1,1,2),                     &
                  x1,x2,dxkm,y1,y2,dykm,                                &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  label(1:27),time, runname,1.0 ,tem1,tem2,tem3,        &
                  tem4,tem5,hterain,slicopt,viqiplt)
              sigplt(102) = 1
            END IF
          END IF

          IF(ipriority(103) == ip .AND. P_QS > 0 ) THEN
            IF( (viqsplt == 1 .OR. viqsplt == 2 .OR. viqsplt == 4       &
                   .OR. viqsplt == 5 )                                  &
                  .AND. sigplt(103) == 0 ) THEN
              CALL get_mulovrlay('viqsplt',7,ovrmul_num,ovrmulname,sovrlay)
              label = 'Vert. Integrated qs (kg/m2)'

              CALL ctrsetup(viqsinc,viqsminc,viqsmaxc,                  &
                  viqsovr,viqshlf,viqszro,viqscol1,viqscol2,'viqsplt     ')

              CALL cal_viqs(tem7, qscalar(:,:,:,P_QS), rhobar, zp, nx,ny,nz)

              CALL ctrsfc(tem7,xc(1,1,2),yc(1,1,2),                     &
                  x1,x2,dxkm,y1,y2,dykm,                                &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  label(1:27),time, runname,1.0 ,tem1,tem2,tem3,        &
                  tem4,tem5,hterain,slicopt,viqsplt)
              sigplt(103) = 1
            END IF
          END IF

          IF(ipriority(104) == ip .AND. P_QH > 0 ) THEN
            IF( (viqhplt == 1 .OR. viqhplt == 2 .OR. viqhplt == 4       &
                   .OR. viqhplt == 5 )                                  &
                  .AND. sigplt(104) == 0 ) THEN
              CALL get_mulovrlay('viqhplt',7,ovrmul_num,ovrmulname,sovrlay)
              label = 'Vert. Integrated qh (kg/m2)'

              CALL ctrsetup(viqhinc,viqhminc,viqhmaxc,                  &
                  viqhovr,viqhhlf,viqhzro,viqhcol1,viqhcol2,'viqhplt     ')

              CALL cal_viqh(tem7, qscalar(:,:,:,P_QH), rhobar, zp, nx,ny,nz)

              CALL ctrsfc(tem7,xc(1,1,2),yc(1,1,2),                     &
                  x1,x2,dxkm,y1,y2,dykm,                                &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  label(1:27),time, runname,1.0 ,tem1,tem2,tem3,        &
                  tem4,tem5,hterain,slicopt,viqhplt)
              sigplt(104) = 1
            END IF
          END IF

          IF(ipriority(105) == ip .AND. P_QC > 0 .AND. P_QR > 0 ) THEN
            IF( (vilplt == 1 .OR. vilplt == 2 .OR. vilplt == 4          &
                   .OR. vilplt == 5 )                                   &
                  .AND. sigplt(105) == 0 ) THEN
              CALL get_mulovrlay('vilplt',6,ovrmul_num,ovrmulname,sovrlay)
              label = 'Vert. Integ Liquid (kg/m2)'

              CALL ctrsetup(vilinc,vilminc,vilmaxc,                     &
                   vilovr,vilhlf,vilzro,vilcol1,vilcol2,'vilplt      ')

              CALL cal_vil(tem7,qscalar(:,:,:,P_QC),qscalar(:,:,:,P_QR),&
                           rhobar,zp, nx,ny,nz,tem6)

              CALL ctrsfc(tem7,xc(1,1,2),yc(1,1,2),                     &
                  x1,x2,dxkm,y1,y2,dykm,                                &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  label(1:26),time, runname,1.0 ,tem1,tem2,tem3,        &
                  tem4,tem5,hterain,slicopt,vilplt)
              sigplt(105) = 1
            END IF
          END IF

          IF(ipriority(106) == ip .AND. P_QI > 0 .AND. P_QS > 0 .AND. P_QH > 0) THEN
            ! Why no P_QG? WYH
            IF( (viiplt == 1 .OR. viiplt == 2 .OR. viiplt == 4          &
                   .OR. viiplt == 5 )                                   &
                  .AND. sigplt(106) == 0 ) THEN
              CALL get_mulovrlay('viiplt',6,ovrmul_num,ovrmulname,sovrlay)
              label = 'Vert. Integrated ice (kg/m2)'

              CALL ctrsetup(viiinc,viiminc,viimaxc,                     &
                   viiovr,viihlf,viizro,viicol1,viicol2,'viiplt      ')

              CALL cal_vii(tem7,qscalar(:,:,:,P_QI),qscalar(:,:,:,P_QS),&
                           qscalar(:,:,:,P_QH),rhobar,zp, nx,ny,nz,tem6)

              CALL ctrsfc(tem7,xc(1,1,2),yc(1,1,2),                     &
                  x1,x2,dxkm,y1,y2,dykm,                                &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  label(1:28),time, runname,1.0 ,tem1,tem2,tem3,        &
                  tem4,tem5,hterain,slicopt,viiplt)
              sigplt(106) = 1
            END IF
          END IF

          IF(ipriority(107) == ip ) THEN
            IF( (vicplt == 1 .OR. vicplt == 2 .OR. vicplt == 4          &
                  .OR. vicplt == 5 )                                    &
                  .AND. sigplt(107) == 0 ) THEN
              CALL get_mulovrlay('vicplt',6,ovrmul_num,ovrmulname,sovrlay)
              label = 'Vert. Integ Condensate (kg/m2)'

              CALL ctrsetup(vicinc,vicminc,vicmaxc,                     &
                   vicovr,vichlf,viczro,viccol1,viccol2,'vicplt      ')

              CALL cal_vic(tem7,qscalar,rhobar,zp,nx,ny,nz,tem6)

              CALL ctrsfc(tem7,xc(1,1,2),yc(1,1,2),                     &
                  x1,x2,dxkm,y1,y2,dykm,                                &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  label(1:30),time, runname,1.0 ,tem1,tem2,tem3,        &
                  tem4,tem5,hterain,slicopt,vicplt)
              sigplt(107) = 1
            END IF
          END IF
!
!-----------------------------------------------------------------------
!
!   Plot storm motion vector
!
!-----------------------------------------------------------------------
!
          IF(ipriority(108) == ip ) THEN
            IF( strmplt == 1 .AND. sigplt(108) == 0 ) THEN
              CALL get_mulovrlay('strmplt',7,ovrmul_num,ovrmulname,sovrlay)
              CALL vtrunt ( strmunit )
              CALL overlay( strmovr )
              CALL ctrcol(strmcol1,strmcol2)
              CALL varplt( 'strmplt     ' )
              CALL ctrvtr( strmunits, strmtype )

              !IF(mp_opt >0) THEN
              !  CALL mpsend1dew(ustrm,nx,ny,1,1,0,mptag1,tem4)
              !  CALL mpsend1dns(ustrm,nx,ny,1,1,0,mptag2,tem5)
              !  CALL mprecv1dew(ustrm,nx,ny,1,1,0,mptag1,tem4)
              !  CALL mprecv1dns(ustrm,nx,ny,1,1,0,mptag2,tem5)
              !
              !  CALL mpsend1dew(vstrm,nx,ny,1,1,0,mptag3,tem2)
              !  CALL mpsend1dns(vstrm,nx,ny,1,1,0,mptag4,tem3)
              !  CALL mprecv1dew(vstrm,nx,ny,1,1,0,mptag3,tem2)
              !  CALL mprecv1dns(vstrm,nx,ny,1,1,0,mptag4,tem3)
              !END IF

              CALL cal_strm(tem7, tem8,ustrm,vstrm,strmunits,nx,ny,nz,  &
                  label,length)

              CALL vtrsfc( tem7(1,1,1),tem8(1,1,1),xc(1,1,2),yc(1,1,2), &
                  x1,x2,dxkm,y1,y2,dykm, nx,ibgnl,iendl,ist,            &
                  ny,jbgnl,jendl,jst,                                   &
                  label(1:length),time, runname,1.0,slicopt,            &
                  tem1,tem2,tem3,tem4,tem5,tem6,hterain)
                  ! size of tem6 must be >= 5*nx*ny
              sigplt(108) = 1
            END IF
          END IF
!
!-----------------------------------------------------------------------
!
!   Plot vertically integrated total water (kg/m**2)
!
!-----------------------------------------------------------------------
!
          IF(ipriority(109) == ip ) THEN
            IF( (vitplt == 1 .OR. vitplt == 2 .OR. vitplt == 4          &
                  .OR. vitplt == 5 )                                    &
                  .AND. sigplt(109) == 0 ) THEN
              CALL get_mulovrlay('vitplt',6,ovrmul_num,ovrmulname,sovrlay)
              label = 'Vert. Integ Total Water (kg/m2)'

              CALL ctrsetup(vitinc,vitminc,vitmaxc,                     &
                   vitovr,vithlf,vitzro,vitcol1,vitcol2,'vitplt      ')

              CALL cal_vit(tem7,qv,qscalar,rhobar,zp,nx,ny,nz,tem6)

              CALL ctrsfc(tem7,xc(1,1,2),yc(1,1,2),                     &
                  x1,x2,dxkm,y1,y2,dykm,                                &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  label(1:31),time, runname,1.0 ,tem1,tem2,tem3,        &
                  tem4,tem5,hterain,slicopt,vitplt)
              sigplt(109) = 1
            END IF
          END IF
!
!-----------------------------------------------------------------------
!
!   Plot precipitable water vapor, cm
!
!-----------------------------------------------------------------------
!
          IF(ipriority(110) == ip ) THEN
            IF( (pwplt == 1 .OR. pwplt == 2 .OR. pwplt == 4 .OR. pwplt == 5 ) &
                  .AND. sigplt(110) == 0 ) THEN
              CALL get_mulovrlay('pwplt',5,ovrmul_num,ovrmulname,sovrlay)
              label = 'Precipitable Water Vapor(cm)'

              CALL ctrsetup(pwinc,pwminc,pwmaxc,                        &
                   pwovr,pwhlf,pwzro,pwcol1,pwcol2,'pwplt       ')

              CALL cal_pw(tem7,qv,rhobar,zp,nx,ny,nz,tem6)

              CALL ctrsfc(tem7,xc(1,1,2),yc(1,1,2),                     &
                  x1,x2,dxkm,y1,y2,dykm,                                &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  label(1:28),time, runname,1.0 ,tem1,tem2,tem3,        &
                  tem4,tem5,hterain,slicopt,pwplt)
              sigplt(110) = 1
            END IF
          END IF
!
!-----------------------------------------------------------------------
!
!   Plot total precipitation rate
!
!-----------------------------------------------------------------------
!
          IF(ipriority(111) == ip ) THEN
            IF( (tprplt == 1 .OR. tprplt == 2 .OR. tprplt == 4          &
                  .OR. tprplt == 5 )                                    &
                  .AND. sigplt(111) == 0 ) THEN
              CALL get_mulovrlay('tprplt',6,ovrmul_num,ovrmulname,sovrlay)
              ! label = 'total Precipitation rate( mm/h)'

              CALL ctrsetup(tprinc,tprminc,tprmaxc,                     &
                   tprovr,tprhlf,tprzro,tprcol1,tprcol2,'tprplt      ')

              CALL cal_tpr(tem7,prcrate(1,1,1),nx,ny,nz,tprunits,       &
                  label,length)

              WHERE (tem7 <= 0.0) tem7 = -1.0e-20

              ! call xcltyp(0)

              CALL ctrsfc(tem7,xc(1,1,2),yc(1,1,2),                     &
                  x1,x2,dxkm,y1,y2,dykm,                                &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  label(1:length),time, runname,1.0 ,tem1,tem2,tem3,    &
                  tem4,tem5,hterain,slicopt,tprplt)
              sigplt(111) = 1

              ! call xcltyp(2)
            END IF
          END IF
!
!-----------------------------------------------------------------------
!
!   Plot grid scale precip. rate
!
!-----------------------------------------------------------------------
!
          IF(ipriority(112) == ip ) THEN
            IF( (gprplt == 1 .OR. gprplt == 2 .OR. gprplt == 4          &
                  .OR. gprplt == 5 )                                    &
                  .AND. sigplt(112) == 0 ) THEN
              CALL get_mulovrlay('gprplt',6,ovrmul_num,ovrmulname,sovrlay)
              !label = 'grid scale precip. rate( mm/h)'

              CALL ctrsetup(gprinc,gprminc,gprmaxc,                     &
                   gprovr,gprhlf,gprzro,gprcol1,gprcol2,'gprplt      ')

              CALL cal_gpr(tem7,prcrate(1,1,2),prcrate(1,1,4), nx,ny,nz,&
                  gprunits,label,length)

              CALL ctrsfc(tem7,xc(1,1,2),yc(1,1,2),                     &
                  x1,x2,dxkm,y1,y2,dykm,                                &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  label(1:length),time, runname,1.0 ,tem1,tem2,tem3,    &
                  tem4,tem5,hterain,slicopt,gprplt)
              sigplt(112) = 1
            END IF
          END IF
!
!-----------------------------------------------------------------------
!
!   Plot Convective precip. rate = prcrate(1,1,3)
!
!-----------------------------------------------------------------------
!
          IF(ipriority(113) == ip ) THEN
            IF( (cprplt == 1 .OR. cprplt == 2 .OR. cprplt == 4          &
                  .OR. cprplt == 5 )                                    &
                  .AND. sigplt(113) == 0 ) THEN
              CALL get_mulovrlay('cprplt',6,ovrmul_num,ovrmulname,sovrlay)
              !label = 'Convective precip. rate( mm/h)'

              CALL ctrsetup(cprinc,cprminc,cprmaxc,                     &
                   cprovr,cprhlf,cprzro,cprcol1,cprcol2,'cprplt      ')

              CALL cal_cpr(tem7,prcrate(1,1,3),nx,ny,nz,cprunits,       &
                  label,length)

              CALL ctrsfc(tem7,xc(1,1,2),yc(1,1,2),                     &
                  x1,x2,dxkm,y1,y2,dykm,                                &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  label(1:length),time, runname,1.0 ,tem1,tem2,tem3,    &
                  tem4,tem5,hterain,slicopt,cprplt)
              sigplt(113) = 1
            END IF
          END IF
!
!-----------------------------------------------------------------------
!
!  Plot surface characteristics.
!
!-----------------------------------------------------------------------
!
          IF(ipriority(121) == ip ) THEN
            IF( (soiltpplt == 1 .OR. soiltpplt == 2 .OR. soiltpplt == 4 &
                   .OR. soiltpplt == 5 )   .AND. sigplt(121) == 0 ) THEN
              WRITE(label,'(''SOIL TYPE '',I1,'' (index)'')') soiltpn
              time0=0.0

              DO j=1,ny
                DO i=1,nx
                  tem9(i,j,1)=soiltyp(i,j,soiltpn)
                END DO
              END DO
              CALL ctrsetup(soiltpinc,soiltpminc,soiltpmaxc,            &
                   styovr,styhlf,styzro,stycol1,stycol2,'soiltpplt   ')

              CALL ctrsfc(tem9,xc(1,1,2),yc(1,1,2),                     &
                  x1,x2,dxkm,y1,y2,dykm,                                &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  TRIM(label),time0, runname,1.0 ,tem1,tem2,tem3,       &
                  tem4,tem5,hterain,slicopt,soiltpplt)
              sigplt(121) = 1
            END IF
          END IF

          IF(ipriority(122) == ip ) THEN
            IF( (vegtpplt == 1 .OR. vegtpplt == 2 .OR. vegtpplt == 4    &
                   .OR. vegtpplt == 5 ) .AND. sigplt(122) == 0) THEN
              label = 'Vegetation type (index)'
              time0=0.0

              DO j=1,ny
                DO i=1,nx
                  tem9(i,j,1)=vegtyp(i,j)
                END DO
              END DO
              CALL ctrsetup(vegtpinc,vegtpminc,vegtpmaxc,               &
                   vtyovr,vtyhlf,vtyzro,vtycol1,vtycol2,'vegtpplt    ')

              CALL ctrsfc(tem9,xc(1,1,2),yc(1,1,2),                     &
                  x1,x2,dxkm,y1,y2,dykm,                                &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  TRIM(label),time0, runname,1.0 ,tem1,tem2,tem3,       &
                  tem4,tem5,hterain,slicopt,vegtpplt)
              sigplt(122)=1
            END IF
          END IF


          IF(ipriority(123) == ip ) THEN
            IF( (laiplt == 1 .OR. laiplt == 2 .OR. laiplt == 4 .OR.     &
                   laiplt == 5 ) .AND. sigplt(123) == 0) THEN
              label = 'Leaf Area Index (index)'
              time0=0.0
              CALL ctrsetup(laiinc,laiminc,laimaxc,                     &
                   laiovr,laihlf,laizro,laicol1,laicol2,'laiplt      ')

              CALL ctrsfc(lai,xc(1,1,2),yc(1,1,2),                      &
                  x1,x2,dxkm,y1,y2,dykm,                                &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  TRIM(label),time0, runname,1.0 ,tem1,tem2,tem3,       &
                  tem4,tem5,hterain,slicopt,laiplt)
              sigplt(123) =1
            END IF
          END IF

          IF(ipriority(124) == ip) THEN
            IF( (rouplt == 1 .OR. rouplt == 2 .OR. rouplt == 4 .OR.     &
                   rouplt == 5) .AND. sigplt(124) == 0) THEN
              label = 'Surface roughness (nounit)'
              time0=0.0
              CALL ctrsetup(rouinc,rouminc,roumaxc,                     &
                   rouovr,rouhlf,rouzro,roucol1,roucol2,'rouplt      ')

              CALL ctrsfc(roufns,xc(1,1,2),yc(1,1,2),                   &
                  x1,x2,dxkm,y1,y2,dykm,                                &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  TRIM(label),time0, runname,1.0 ,tem1,tem2,tem3,       &
                  tem4,tem5,hterain,slicopt,rouplt)
              sigplt(124) = 1
            END IF
          END IF

          IF(ipriority(125) == ip) THEN
            IF( (vegplt == 1 .OR. vegplt == 2 .OR. vegplt == 4  .OR.    &
                   vegplt == 5)  .AND. sigplt(125) == 0) THEN
              label = 'Vegetation fraction (fraction)'
              time0=0.0
              CALL ctrsetup(veginc,vegminc,vegmaxc,                     &
                   vegovr,veghlf,vegzro,vegcol1,vegcol2,'vegplt      ')

              CALL ctrsfc(veg,xc(1,1,2),yc(1,1,2),                      &
                  x1,x2,dxkm,y1,y2,dykm,                                &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  TRIM(label),time0, runname,1.0 ,tem1,tem2,tem3,       &
                  tem4,tem5,hterain,slicopt,vegplt)
              sigplt(125) = 1
            END IF
          END IF

          IF(ipriority(126) == ip) THEN
            IF( (snowdplt == 1 .OR. snowdplt == 2 .OR. snowdplt == 4  .OR. &
                   snowdplt == 5)  .AND. sigplt(126) == 0) THEN
              label = 'Snow depth (m)'
              CALL ctrsetup(snowdinc,snowdminc,snowdmaxc,               &
                   snowdovr,snowdhlf,snowdzro,snowdcol1,                &
                   snowdcol2,'snowdplt      ')

              DO j=1,ny
                DO i=1,nx
                  tem9(i,j,1)=snowdpth(i,j)
                END DO
              END DO

              CALL ctrsfc(tem9,xc(1,1,2),yc(1,1,2),                     &
                  x1,x2,dxkm,y1,y2,dykm,                                &
                  nx,ibgnl,iendl, ny,jbgnl,jendl,                       &
                  TRIM(label),time, runname, 1.0,tem1,tem2,tem3,        &
                  tem4,tem5,hterain,slicopt,snowdplt)
              sigplt(126) = 1
            END IF
          END IF

!-----------------------------------------------------------------------
!
!  Overlay wind vector/barb if necessary  on 2-D surface fields.
!
!-----------------------------------------------------------------------
!
          CALL varplt( 'vtrplt      ' )
          IF(varname(1:6) == ovrname(1:6) .AND. sovrlay == 1) plotovr=.true.
          IF( vtrplt == 1 .AND. plotovr ) THEN
            CALL vtrunt ( vtrunit )
            IF(plotovr) CALL overlay ( 1 )
            IF(.NOT.plotovr)  CALL overlay( vtrovr )
            CALL ctrcol(vtrcol1,vtrcol2)
            CALL varplt( 'vtrplt      ' )
            CALL ctrvtr( vtrunits, vtrtype )

            CALL cal_vtr(tem7,tem8,tem9,u,v,w,nx,ny,nz,vtrunits,        &
                label,length)

            CALL vtrsfc( tem7(1,1,1),tem8(1,1,1),xc(1,1,2),yc(1,1,2),   &
                x1,x2,dxkm,y1,y2,dykm, nx,ibgnl,iendl,ist,              &
                ny,jbgnl,jendl,jst,                                     &
                label(1:length),time, runname,1.0,slicopt,              &
                tem1,tem2,tem3,tem4,tem5,tem6,hterain)
                ! size of tem6 must be >= 5*nx*ny
            IF(plotovr) THEN
              plotovr=.false.
              sovrlay = 0
            END IF
          END IF

        END IF  ! For slicopt=8 - two plots
!
!-----------------------------------------------------------------------
!
!  plot arbitrary variables.
!
!-----------------------------------------------------------------------
!

        IF(arbvaropt == 1) THEN
          IF(slicopt <= 7 .AND. var3dnum > 0) THEN
            DO i = 1,var3dnum
              IF(ipriority(150+i) == ip ) THEN
                IF(var3dplot(i) /= 0) THEN

                  CALL ctrsetup(var3dinc(i),var3dminc(i),var3dmaxc(i),  &
                       var3dovr(i),var3dhlf(i),var3dzro(i),             &
                       var3dcol1(i),var3dcol2(i),var3d(i))

                  CALL arbvar_readvar2(nx,ny,nz, var3dv,                &
                      var3d(i), var_name, varunits,                     &
                      time, runname, dirname3d(i), filename3d(i),       &
                      finfmt3d(i), 1, lvldbg, istatus)
                                ! hard coded root read and split

                  label=var3d(i)(1:6)//'('//TRIM(varunits)//')'

                  posdefflag = 1

                  DO kd=2,nz-2
                    DO jd=2,ny-2
                      DO id=2,nx-2
                        IF(var3dv(id,jd,kd) < 0.0) THEN
                          posdefflag = 0
                          EXIT
                        END IF
                      END DO
                    END DO
                  END DO

                  IF(posdefflag == 1) WHERE (var3dv <= 0.0) var3dv = -1.0e-20

                  CALL ctr3d(var3dv,xc,yc,zps3d,                        &
                       x1,x2,dxkm,y1,y2,dykm,z1,z2,dzkm,                &
                      nx,ibgnl,iendl, ny,jbgnl,jendl, nz,kbgn,kend,     &
                      label,time,slicopt,kslice,jslice,islice,          &
                      n,xp,yp,b1,b2,zs2,                                &
                      runname,1.0,tem1,tem2,tem3,                       &
                      tem4,tem5,tem6,hterain,var3dplot(i))

                END IF
              END IF
            END DO
          END IF

          IF((slicopt == 8 .OR. slicopt == 1) .AND. var2dnum > 0) THEN
            DO i = 1,var2dnum
              IF(ipriority(170+i) == ip ) THEN
                IF(var2dplot(i) /= 0 .AND. sigplt(170+i) == 0) THEN
                  CALL ctrsetup(var2dinc(i),var2dminc(i),var2dmaxc(i),  &
                       var2dovr(i),var2dhlf(i),var2dzro(i),             &
                       var2dcol1(i),var2dcol2(i),var2d(i))

                  CALL arbvar_readvar2(nx,ny,1, var2dv,                 &
                      var2d(i), var_name, varunits,                     &
                      time, runname, dirname2d(i), filename2d(i),       &
                      finfmt2d(i), 1, lvldbg, istatus)
                                ! hard coded root read and split

                  IF (istatus /= 0)   &
                    CALL arpsstop('ERROR: after readvar2 in arpsplt.f90.',1)

                  label=var2d(i)(1:6)//'('//varunits//')'

                  CALL ctrsfc(var2dv,xc(1,1,2),yc(1,1,2),               &
                      x1,x2,dxkm,y1,y2,dykm,                            &
                      nx,ibgnl,iendl, ny,jbgnl,jendl,                   &
                      label,time, runname,1.0 ,tem1,tem2,tem3,  &
                      tem4,tem5,hterain,1,var2dplot(i))
                  sigplt(170+i) = 1

                END IF
              END IF
            END DO
          END IF

          IF((slicopt == 1 .OR. slicopt == 6 .OR. slicopt == 8) .AND. vtr2dnum > 0) THEN
            DO i = 1,vtr2dnum
              IF(ipriority(190+i) == ip ) THEN
                IF((vtraplt(i) /= 0 .OR. magaplt(i) /= 0)  &
                    .AND. sigplt(190+i) == 0) THEN
                  need_uv=.true.
                  IF(vtraplt(i) /= 0 .AND. sigplt(190+i) == 0) THEN
!                    CALL get_mulovrlay('vtraplt',7,ovrmul_num,ovrmulname,sovrlay)

                    CALL arbvar_readvar2(nx,ny,1, var2du,               &
                      vtru2d(i), var_name, varunits,                    &
                      time, runname, diruv2d(i), filenameu2d(i),        &
                      finfmtuv2d(i), 1, lvldbg, istatus)
                                ! hard coded root read and split

                    CALL arbvar_readvar2(nx,ny,1, var2dv,               &
                      vtrv2d(i), var_name, varunits,                    &
                      time, runname, diruv2d(i), filenamev2d(i),        &
                      finfmtuv2d(i), 1, lvldbg, istatus)
                                ! hard coded root read and split

                    need_uv = .false.

                    CALL vtrunt ( vtraunit(i) )
                    CALL overlay( vtraovr(i) )
                    CALL ctrcol(vtracol1(i),vtracol2(i))
                    CALL varplt( 'vtraplt     ' )
                    CALL ctrvtr( vtraunits(i), vtratype(i) )

                    label='TEST'
                    length=4
                    IF (myproc == 0) print *, ' Calling cal_arbvtr'

                    CALL cal_arbvtr(tem7,tem8,var2du,var2dv,vtraunits(i), &
                         nx,ny,nz,vtrv2d(i),label,length)

                    IF (myproc == 0) THEN
                      PRINT *, ' label   = "',label(1:length),'", length = ',length
                      PRINT *, ' Calling vtrsfc with slicopt = ',slicopt
                    END IF

                    iast = ist    ! use istride by default
                    jast = jst
                    IF (iastride(i) /= 0) iast = iastride(i)
                    IF (jastride(i) /= 0) jast = jastride(i)

                    ibgnla = (loc_x-1)*(nx-3)+ibgnl  ! hold global index for ibgnl temporarily
                    IF(MOD(ibgnla-ibgn,iast) /= 0)   THEN
                      ibgnla = ibgnl + (iast- MOD(ibgnla - ibgn,iast)) ! for vector plot
                    ELSE
                      ibgnla = ibgnl
                    END IF
                    jbgnla = (loc_y-1)*(ny-3)+jbgnl  ! hold global index for jbgnl
                    IF(MOD(jbgnla-jbgn,jast) /= 0)   THEN
                      jbgnla = jbgnl + (jast- MOD(jbgnla - jbgn,jast)) ! for vector plot
                    ELSE
                      jbgnla = jbgnl
                    END IF

                    CALL vtrsfc( tem7(1,1,1),tem8(1,1,1),xc(1,1,2),yc(1,1,2), &
                         x1,x2,dxkm,y1,y2,dykm,                         &
                         nx,ibgnla,iendl,iast, ny,jbgnla,jendl,jast,    &
                         label(1:length),time, runname,1.0,slicopt,     &
                         tem1,tem2,tem3,tem4,tem5,tem6,hterain)
                  END IF

                  IF(magaplt(i) /= 0 .AND. sigplt(190+i) == 0) THEN
                    IF(need_uv) THEN
                      CALL arbvar_readvar2(nx,ny,1, var2du,             &
                        vtru2d(i), var_name, varunits,                  &
                        time, runname, diruv2d(i), filenameu2d(i),      &
                        finfmtuv2d(i), 1, lvldbg, istatus)
                                ! hard coded root read and split

                      CALL arbvar_readvar2(nx,ny,1, var2dv,             &
                        vtrv2d(i), var_name, varunits,                  &
                        time, runname, diruv2d(i), filenamev2d(i),      &
                        finfmtuv2d(i), 1, lvldbg, istatus)
                                ! hard coded root read and split
                    END IF

                    DO jj=1,ny
                      DO ii=1,nx
                        tem7(ii,jj,1)=sqrt(var2du(ii,jj)*var2du(ii,jj)+ &
                                           var2dv(ii,jj)*var2dv(ii,jj))
                        tem7(ii,jj,1)=min(tem7(ii,jj,1),500.)
                      END DO
                    END DO

                    CALL ctrsetup(magainc(i),magaminc(i),magamaxc(i),   &
                                  magaovr(i),magahlf(i),magazro(i),     &
                                  magacol1(i),magacol2(i),vtru2d(i))

                    CALL ctrsfc(tem7(1,1,1),xc(1,1,2),yc(1,1,2),        &
                      x1,x2,dxkm,y1,y2,dykm,                            &
                      nx,ibgnl,iendl, ny,jbgnl,jendl,                   &
                      label,time, runname,1.0 ,tem1,tem2,tem3,  &
                      tem4,tem5,hterain,1,magaplt(i))
                  END IF
                  sigplt(190+i) = 1
                END IF
              END IF
            END DO
          END IF

        END IF
!
! 06/06/2002 Zuwen He
!
! plotting soil variables, tsoil and qsoil
!
        IF (iplot(slicopt) /= 0 .AND.  &
            (slicopt == 9 .OR. slicopt == 10 .OR. slicopt == 11)) THEN

          IF(ipriority(69) == ip ) THEN
            IF(tsoilplt == 1 .OR. tsoilplt == 2 .OR. tsoilplt == 4 .OR. &
               tsoilplt == 5) THEN
              CALL get_mulovrlay('tsoilplt',8,ovrmul_num,ovrmulname,sovrlay)
              label = 'Soil temperature (K)'
              CALL ctrsetup(tsoilinc,tsoilminc,tsoilmaxc,               &
                   tsoilovr,tsoilhlf,tsoilzro,tsoilcol1,tsoilcol2,      &
                   'tsoilplt    ')

              CALL ctr3d(tsoil(1,1,1,0),xc,yc,zpsoils3d,                &
                   x1,x2,dxkm,y1,y2,dykm,zsoil1,zsoil2,dzsoilcm,        &
                   nx,ibgnl,iendl,ny,jbgnl,jendl,nzsoil,ksoilbgn,ksoilend, &
                   label(1:25),                                         &
                   time,slicopt,ksoilslice,jsoilslice,isoilslice,       &
                   n,xp,yp,b1,b2,zs2,                                   &
                   runname,1.0,tem1,tem2,tem3,                          &
                   tem4,tem5,tem6,hterain,tsoilplt)

              sigplt(69) = 1
            END IF
          END IF

          IF(ipriority(70) == ip ) THEN
            IF(qsoilplt == 1 .OR. qsoilplt == 2 .OR. qsoilplt == 4 .OR. &
                 qsoilplt == 5) THEN
              CALL get_mulovrlay('qsoilplt',8,ovrmul_num,ovrmulname,sovrlay)
              label = 'Soil moisture (m**3/m**3)'
              CALL ctrsetup(qsoilinc,qsoilminc,qsoilmaxc,               &
                   qsoilovr,qsoilhlf,qsoilzro,qsoilcol1,qsoilcol2,      &
                   'qsoilplt    ')

              CALL ctr3d(qsoil(1,1,1,0),xc,yc,zpsoils3d,                &
                   x1,x2,dxkm,y1,y2,dykm,zsoil1,zsoil2,dzsoilcm,        &
                   nx,ibgnl,iendl, ny,jbgnl,jendl,nzsoil,ksoilbgn,ksoilend, &
                   label(1:25),time,slicopt,                            &
                   ksoilslice,jsoilslice,isoilslice,n,xp,yp,b1,b2,zs2,  &
                   runname,1.0,tem1,tem2,tem3,                          &
                   tem4,tem5,tem6,hterain,qsoilplt)

              sigplt(70) = 1
            END IF
          END IF

        END IF
!
      END DO ! loop over priority index ipp, plotting all fields

!-----------------------------------------------------------------------
!  Overlaying trajectories on each slice.
!-----------------------------------------------------------------------

      IF( slicopt /= 8 .AND. trajc_plt_opt /= 0 ) THEN

        xor_current = x(2)
        yor_current = y(2)
        IF (myproc == 0) CALL plt_trajc(slicopt,curtim,x01,x02,y01,y02, &
                    x1,lblmag,iorig,                                    &
                    xor_current,yor_current,xgrdorg,ygrdorg,xorig,yorig)
      ENDIF
!
!
      END DO ! loop over indxslic, the number of slices for each type

    END DO  ! loop over slicopt, for different types of slices (x-y, x-z, etc).

!
!-----------------------------------------------------------------------
!
!  3-D wireframe isosurface plots.
!
!-----------------------------------------------------------------------
!
    IF( w3dplt == 1) THEN
!
!    Plot w=wisosf isosurface.
!
      CALL wirfrm(w,nx,ibgn,iend,ny,jbgn,jend,nz,kbgn,kend,             &
                  wisosf, tem2)

      IF(myproc == 0) THEN

        CALL xqpspc(x1,x2,y1,y2)
        CALL xchmag(0.02*(y2-y1)*lblmag)
        CALL xqmap (x1,x2,y1,y2)
        CALL xcharl(x1,y1-0.02*(y2-y1),runname)

        WRITE(label,'(a,f4.1,a,f6.1,a)')                                &
            'w=',wisosf,' m/s isosurface at t=',time/60.0,' min'

        CALL xcharl(x1,y1-0.04*(y2-y1),label )

      END IF
!
!    Plot w=-wisosf isosurface.
!
      CALL wirfrm(w,nx,ibgn,iend,ny,jbgn,jend,nz,kbgn,kend,             &
                    -wisosf, tem2)

      IF(myproc == 0) THEN

        CALL xqpspc(x1,x2,y1,y2)
        CALL xchmag(0.02*(y2-y1)*lblmag)
        CALL xqmap (x1,x2,y1,y2)
        CALL xcharl(x1,y1-0.02*(y2-y1),runname)

        WRITE(label,'(a,f4.1,a,f6.1,a)')                                &
            'w=',-wisosf,' m/s isosurface at t=',time/60.0,' min'

        CALL xcharl(x1,y1-0.04*(y2-y1),label )

      END IF

    END IF

    IF( q3dplt == 1 .AND. P_QC > 0 .AND. P_QR > 0) THEN
      DO k=1,nz
        DO j=1,ny
          DO i=1,nx
            tem1(i,j,k)=(qscalar(i,j,k,P_QC)+qscalar(i,j,k,P_QR))*1000.0
          END DO
        END DO
      END DO

      CALL wirfrm(tem1,                                                 &
           nx,ibgn,iend,ny,jbgn,jend,nz,kbgn,kend,qisosf,tem2)

      IF(myproc == 0) THEN

        CALL xqpspc(x1,x2,y1,y2)
        CALL xchmag(0.02*(y2-y1)*lblmag)
        CALL xqmap (x1,x2,y1,y2)
        CALL xcharl(x1,y1-0.02*(y2-y1),runname)

        WRITE(label,'(a,f4.1,a,f6.1,a)')                                &
            'qr+qc=',qisosf,' kg/kg isosurface at t=',time/60.0,' min'

        CALL xcharl(x1,y1-0.04*(y2-y1),label )

      END IF

    END IF

    IF (profopt == 0) GO TO 900
!
!-----------------------------------------------------------------------
!
!  Reinitialize plotting parameters for profile plotting
!
!-----------------------------------------------------------------------
!
    CALL xdspac(0.9*winsiz)
    IF( nxprpic == 1 .AND. nyprpic == 1) CALL xdspac(0.85*winsiz)

    CALL xspace(nxprpic, nyprpic, angl , xlimit,ylimit)

    CALL xpmagn(margnx*xlimit, margny*ylimit)

    CALL xnwfrm

!
!-----------------------------------------------------------------------
!
!  Draw all the profiles based on the user inputs
!
!-----------------------------------------------------------------------
!
    IF (uprof == 1) THEN

      onvf = 0
      CALL avgx(u, onvf, nx,ny,nz, 1, nx-1, 1, ny-1, 1, nz-1, tem1)
      onvf = 1
      nzprofc = nz-1
      CALL vprofil (nx,ny,nz,nzprofc,tem1,xc,yc,zpc,uprmin,uprmax,      &
           xprof, yprof, nprof, zprofbgn, zprofend,                     &
           'u (m/s)', 'height (km)',npicprof,tem2,tem3)
      IF(myproc == 0) CALL runlab (runname)

    END IF

    IF (vprof == 1) THEN

      onvf = 0
      CALL avgy(v, onvf, nx,ny,nz, 1, nx-1, 1, ny-1, 1, nz-1, tem1)
      onvf = 1
      nzprofc = nz-1
      CALL vprofil (nx,ny,nz,nzprofc,tem1,xc,yc,zpc,vprmin,vprmax,      &
           xprof, yprof, nprof, zprofbgn, zprofend,                     &
           'v (m/s)', 'height (km)',npicprof,tem2,tem3)
      IF(myproc == 0) CALL runlab (runname)

    END IF

    IF (wprof == 1) THEN

      onvf = 0
      CALL avgz(w, onvf, nx,ny,nz, 1, nx-1, 1, ny-1, 1, nz-1, tem1)
      onvf = 1
      nzprofc = nz-1
      CALL vprofil (nx,ny,nz,nzprofc,tem1,xc,yc,zpc,wprmin,wprmax,      &
           xprof, yprof, nprof, zprofbgn, zprofend,                     &
           'w (m/s)', 'height (km)',npicprof,tem2,tem3)
      IF(myproc == 0) CALL runlab (runname)

    END IF

    IF (ptprof == 1) THEN

      DO k = 1, nz
        DO j = 1, ny
          DO i = 1, nx
            tem1(i,j,k) = ptprt(i,j,k)+ptbar(i,j,k)
          END DO
        END DO
      END DO
      nzprofc = nz-1
      CALL vprofil (nx,ny,nz,nzprofc,tem1,xc,yc,zpc,ptprmin,ptprmax,    &
           xprof, yprof, nprof, zprofbgn, zprofend,                     &
           'pt (K)', 'height (km)',npicprof,tem2,tem3)
      IF(myproc == 0) CALL runlab (runname)

    END IF

    IF (pprof == 1) THEN

      DO k = 1, nz
        DO j = 1, ny
          DO i = 1, nx
            tem1(i,j,k) = pprt(i,j,k)+pbar(i,j,k)
          END DO
        END DO
      END DO
      nzprofc = nz-1
      CALL vprofil (nx,ny,nz,nzprofc,tem1,xc,yc,zpc,pprmin,pprmax,      &
           xprof, yprof, nprof, zprofbgn, zprofend,                     &
           'p (Pa)','height (km)',npicprof,tem2,tem3)
      IF(myproc == 0) CALL runlab (runname)

    END IF

    IF (qvprof == 1) THEN

      DO k = 1, nz
        DO j = 1, ny
          DO i = 1, nx
            tem1(i,j,k) = qv(i,j,k) * 1000.0
          END DO
        END DO
      END DO

      nzprofc = nz-1
      CALL vprofil (nx,ny,nz,nzprofc,tem1,xc,yc,zpc,qvprmin,qvprmax,    &
           xprof, yprof, nprof, zprofbgn, zprofend,                     &
           'qv (g/kg)','height (km)',npicprof,tem2,tem3)
      IF(myproc == 0) CALL runlab (runname)

    END IF

    IF (qcprof == 1 .AND. P_QC >0) THEN

      DO k = 1, nz
        DO j = 1, ny
          DO i = 1, nx
            tem1(i,j,k) = qscalar(i,j,k,P_QC) * 1000.0
          END DO
        END DO
      END DO

      nzprofc = nz-1
      CALL vprofil (nx,ny,nz,nzprofc,tem1,xc,yc,zpc,qcpmin,qcpmax,      &
           xprof, yprof, nprof, zprofbgn, zprofend,                     &
           'qc (g/kg)','height (km)',npicprof,tem2,tem3)
      IF(myproc == 0) CALL runlab (runname)

    END IF

    IF (qrprof == 1 .AND. P_QR > 0) THEN

      DO k = 1, nz
        DO j = 1, ny
          DO i = 1, nx
            tem1(i,j,k) = qscalar(i,j,k,P_QR) * 1000.0
          END DO
        END DO
      END DO

      nzprofc = nz-1
      CALL vprofil (nx,ny,nz,nzprofc,tem1,xc,yc,zpc,qrpmin,qrpmax,      &
           xprof, yprof, nprof, zprofbgn, zprofend,                     &
           'qr (g/kg)','height (km)',npicprof,tem2,tem3)
      IF(myproc == 0) CALL runlab (runname)

    END IF

    IF (qiprof == 1 .AND. P_QI > 0) THEN

      DO k = 1, nz
        DO j = 1, ny
          DO i = 1, nx
            tem1(i,j,k) = qscalar(i,j,k,P_QI) *  1000.0
          END DO
        END DO
      END DO

      nzprofc = nz-1
      CALL vprofil (nx,ny,nz,nzprofc,tem1,xc,yc,zpc,qipmin,qipmax,      &
           xprof, yprof, nprof, zprofbgn, zprofend,                     &
           'qi (g/kg)','height (km)',npicprof,tem2,tem3)
      IF(myproc == 0) CALL runlab (runname)

    END IF

    IF (qsprof == 1 .AND. P_QS > 0) THEN

      DO k = 1, nz
        DO j = 1, ny
          DO i = 1, nx
            tem1(i,j,k) = qscalar(i,j,k,P_QS) * 1000.0
          END DO
        END DO
      END DO

      nzprofc = nz-1
      CALL vprofil (nx,ny,nz,nzprofc,tem1,xc,yc,zpc,qspmin,qspmax,      &
           xprof, yprof, nprof, zprofbgn, zprofend,                     &
           'qs (g/kg)','height (km)',npicprof,tem2,tem3)
      IF(myproc == 0) CALL runlab (runname)

    END IF

    IF (qhprof == 1 .AND. P_QH > 0) THEN

      DO k = 1, nz
        DO j = 1, ny
          DO i = 1, nx
            tem1(i,j,k) = qscalar(i,j,k,P_QH) * 1000.0
          END DO
        END DO
      END DO

      nzprofc = nz-1
      CALL vprofil (nx,ny,nz,nzprofc,tem1,xc,yc,zpc,qhpmin,qhpmax,      &
           xprof, yprof, nprof, zprofbgn, zprofend,                     &
           'qh (g/kg)','height (km)',npicprof,tem2,tem3)
      IF(myproc == 0) CALL runlab (runname)

    END IF

    IF (kmhprof == 1) THEN

      DO k = 1, nz
        DO j = 1, ny
          DO i = 1, nx
            tem1(i,j,k) = kmh(i,j,k)
          END DO
        END DO
      END DO

      nzprofc = nz-1
      CALL vprofil (nx,ny,nz,nzprofc,tem1,xc,yc,zpc,kmhpmin,kmhpmax,    &
           xprof, yprof, nprof, zprofbgn, zprofend,                     &
           'kmh (m**2/s)','height (km)',npicprof,tem2,tem3)
      IF(myproc == 0) CALL runlab (runname)

    END IF

    IF (kmvprof == 1) THEN

      DO  k = 1, nz
        DO  j = 1, ny
          DO  i = 1, nx
            tem1(i,j,k) = kmv(i,j,k) * 1.0
          END DO
        END DO
      END DO

      nzprofc = nz-1
      CALL vprofil (nx,ny,nz,nzprofc,tem1,xc,yc,zpc,kmvpmin,kmvpmax,    &
           xprof, yprof, nprof, zprofbgn, zprofend,                     &
           'kmv (m**2/s)','height (km)',npicprof,tem2,tem3)
      IF(myproc == 0) CALL runlab (runname)

    END IF

    IF (tkeprof == 1) THEN

      DO k = 1, nz
        DO j = 1, ny
          DO i = 1, nx
            tem1(i,j,k) = tke(i,j,k) * 1.0
          END DO
        END DO
      END DO

      nzprofc = nz-1
      CALL vprofil (nx,ny,nz,nzprofc,tem1,xc,yc,zpc,tkepmin,tkepmax,    &
             xprof, yprof, nprof, zprofbgn, zprofend,                   &
             'tke ((m/s)**2)','height (km)',npicprof,tem2,tem3)
      IF(myproc == 0) CALL runlab (runname)

    END IF

    IF (rhprof == 1) THEN

      CALL temper(nx,ny,nz,pt, pprt,pbar, tem1)
      CALL getqvs(nx,ny,nz, 1,nx-1,1,ny-1,1,nz-1, pbar,tem1,tem2)

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem1(i,j,k) = qv(i,j,k)/tem2(i,j,k)
          END DO
        END DO
      END DO

      nzprofc = nz-1
      CALL vprofil (nx,ny,nz,nzprofc,tem1,xc,yc,zpc,rhpmin,rhpmax,          &
             xprof, yprof, nprof, zprofbgn, zprofend,                       &
             'RH', 'height (km)',npicprof,tem2,tem3)
      IF(myproc == 0) CALL runlab (runname)

    END IF

    IF (rfprof == 1) THEN

      IF (mphyopt > 100 .AND. mphyopt < 200) THEN

        wrfopt = MOD(mphyopt,100)

        DO k=1,nz-1
          DO j=1,ny-1
            DO i=1,nx-1
              tem1(i,j,k) = pbar(i,j,k )+pprt (i,j,k)
            END DO
          END DO
        END DO
        CALL edgfill(tem1,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)

        SELECT CASE (wrfopt)
        CASE (1,3,4) ! Kessler, WSM 3-class
          print *,' Reflectivity formula to be implemented for WRF mpphysics = ',wrfopt,'.'
          CALL arpsstop('ERROR: no reflectivity formular implemented.',1)
        CASE (2,6,10,16) ! Lin, WSM 5- and 6- classes
          print *,' CALL reflec_wrf ...'
          CALL reflec_wrf(nx,ny,nz,qv,qscalar(:,:,:,P_QR),              &
                                      qscalar(:,:,:,P_QS),              &
                                      qscalar(:,:,:,P_QG),              &
                                      1,tem1,tz,tem5) ! mixed-phase
        CASE (5) ! Ferrier microphysics scheme (as in WRF_POST)
          print *,' CALL reflec_ferrier_wrf ...'
          tem2(:,:,:) = 2.0
          CALL reflec_ferrier_wrf(nx,ny,nz,qv,qscalar(:,:,:,P_QC),      &
                                              qscalar(:,:,:,P_QR),      &
                                              qscalar(:,:,:,P_QS),      &
                                  tem1,tz,tem5,tem8,tem2)
        CASE (8) ! Thompson microphysics scheme (old version)
          print *,' CALL CALREF9s ...'
          CALL CALREF9s(nx,ny,nz,qscalar(:,:,:,P_QR),                   &
                                 qscalar(:,:,:,P_QS),                   &
                                 qscalar(:,:,:,P_QG),                   &
                                 tem1,tz,tem5)
        CASE (9,17)
          print *,' Calling reflec_wrf ...'
          tem2(:,:,:) = qscalar(:,:,:,P_QG)+qscalar(:,:,:,P_QH)
          CALL reflec_wrf(nx,ny,nz,qv,qscalar(:,:,:,P_QR),              &
                                      qscalar(:,:,:,P_QS),              &
                                      tem2,                             &
                          1,tem1,tz,tem5) ! mixed-phase
        CASE DEFAULT
          WRITE(6,'(1x,a,I2,a)') 'ERROR: Unknown mphyopt = ',mphyopt,'.'
          CALL arpsstop('ERROR: Unsupported mphopt option',1)
        END SELECT

      ELSE IF (mphyopt > 200) THEN ! COAMPS original scheme

        DO k=1,nz-1
          DO j=1,ny-1
            DO i=1,nx-1
              tem1(i,j,k) = pbar(i,j,k )+pprt (i,j,k)
            END DO
          END DO
        END DO
        CALL edgfill(tem1,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)

        SELECT CASE (mphyopt)
        CASE (201)
          CALL coamps_rrf(tz,tem1,qscalar(:,:,:,P_QR),                  &
                                  qscalar(:,:,:,p_QS),                  &
                                  qscalar(:,:,:,P_QI),                  &
                                  qscalar(:,:,:,P_QC),                  &
                                  qscalar(:,:,:,P_QG),                  &
                          qv,nx,ny,nz,tem5)
        CASE (204)
          CALL CALREF9s(nx,ny,nz,qscalar(:,:,:,P_QR),                   &
                                 qscalar(:,:,:,P_QS),                   &
                                 qscalar(:,:,:,P_QG),                   &
                                 tem1,tz,tem5)
        CASE (205:209)
          CALL reflec_MM(nx,ny,nz,rhobar,qscalar,tz,tem5)
        END SELECT

      ! Milbrandt and Yau 3-moment scheme
      ELSE IF (mphyopt >= 8 .and. mphyopt <= 12) THEN
        CALL reflec_MM(nx,ny,nz,rhobar,qscalar,tz,tem5)

      ! Lin,Schultz,Straka Lin,or WSM6 schemes
      ELSE IF (mphyopt >= 2 .and. mphyopt <= 7) THEN
        CALL reflec_ferrier(nx,ny,nz, rhobar, qscalar, tz, tem5)

      ! Warm rain schemes
      ELSE IF (mphyopt == 1) THEN

        IF(P_QR > 0) THEN
          CALL reflec_wr(nx,ny,nz, rhobar, qscalar(:,:,:,P_QR),tem5)
        ELSE
          tem5 = 0.0
        END IF

      ELSE
        print*,'ERROR: Invalid microphysics option, reflectivity set to zero.'
        tem5 = 0.0

      END IF

      nzprofc = nz-1
      CALL vprofil (nx,ny,nz,nzprofc,tem5,xc,yc,zpc,rfpmin,rfpmax,      &
             xprof, yprof, nprof, zprofbgn, zprofend,                   &
             'ref (dBZ)','height (km)',npicprof,tem2,tem3)
      IF(myproc == 0) CALL runlab (runname)

    END IF

    IF (pteprf == 1) THEN

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem2(i,j,k)=pprt(i,j,k)+pbar(i,j,k)
          END DO
        END DO
      END DO

      CALL pt2pte(nx,ny,nz,1,nx-1,1,ny-1,1,nz-1,tem2,pt,qv,tem1)

      nzprofc = nz-1
      CALL vprofil (nx,ny,nz,nzprofc,tem1,xc,yc,zpc,ptepmin,ptepmax,    &
           xprof, yprof, nprof, zprofbgn, zprofend,                     &
           'pte (K)', 'height (km)',npicprof,tem2,tem3)
      IF(myproc == 0) CALL runlab (runname)

    END IF

    IF (upprof == 1) THEN

      nzprofc = nz-1
      CALL vprofil (nx,ny,nz,nzprofc,uprt,xc,yc,zpc,uppmin,uppmax,      &
           xprof, yprof, nprof, zprofbgn, zprofend,                     &
           'uprt (m/s)','height (km)',npicprof,tem2,tem3)
      IF(myproc == 0) CALL runlab (runname)

    END IF

    IF (vpprof == 1) THEN

      nzprofc = nz-1
      CALL vprofil (nx,ny,nz,nzprofc,vprt,xc,yc,zpc,vppmin,vppmax,      &
           xprof, yprof, nprof, zprofbgn, zprofend,                     &
           'vprt (m/s)','height (km)',npicprof,tem2,tem3)
      IF(myproc == 0) CALL runlab (runname)

    END IF

    IF (wpprof == 1) THEN

      nzprofc = nz-1
      CALL vprofil (nx,ny,nz,nzprofc,wprt,xc,yc,zpc,wppmin,wppmax,      &
           xprof, yprof, nprof, zprofbgn, zprofend,                     &
           'wprt (m/s)','height (km)',npicprof,tem2,tem3)
      IF(myproc == 0) CALL runlab (runname)

    END IF

    IF (ptpprf == 1) THEN

      nzprofc = nz-1
      CALL vprofil (nx,ny,nz,nzprofc,ptprt,xc,yc,zpc,ptppmin,ptppmax,   &
           xprof, yprof, nprof, zprofbgn, zprofend,                     &
           'ptprt (K)','height (km)',npicprof,tem2,tem3)
      IF(myproc == 0) CALL runlab (runname)

    END IF

    IF (ppprof == 1) THEN

      nzprofc = nz-1
      CALL vprofil (nx,ny,nz,nzprofc,pprt,xc,yc,zpc,pppmin,pppmax,      &
           xprof, yprof, nprof, zprofbgn, zprofend,                     &
           'pprt (Pa)','height (km)',npicprof,tem2,tem3)
      IF(myproc == 0) CALL runlab (runname)

    END IF

    IF (qvpprf == 1) THEN

      DO k = 1, nz
        DO j = 1, ny
          DO i = 1, nx
            tem1(i,j,k) = qvprt(i,j,k) * 1000.0
          END DO
        END DO
      END DO

      nzprofc = nz-1
      CALL vprofil (nx,ny,nz,nzprofc,tem1,xc,yc,zpc,qvppmin,qvppmax,    &
           xprof, yprof, nprof, zprofbgn, zprofend,                     &
           'qvprt (g/kg)','height (km)',npicprof,tem2,tem3)
      IF(myproc == 0) CALL runlab (runname)

    END IF

    IF (vorpprf == 1) THEN

      DO k=2,nz-2
        DO j=2,ny-2
          DO i=2,nx-2
            tem1(i,j,k)= 1.0E5*(                                        &
                (v(i+1,j,k)-v(i-1,j,k)+v(i+1,j+1,k)-v(i-1,j+1,k))/      &
                (4*(x(i+1)-x(i)))-                                      &
                (u(i,j+1,k)-u(i,j-1,k)+u(i+1,j+1,k)-u(i+1,j-1,k))/      &
                (4*(y(j+1)-y(j))) )
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

      nzprofc = nz-1
      CALL vprofil (nx,ny,nz,nzprofc,tem1,xc,yc,zpc,vorppmin,vorppmax,  &
           xprof, yprof, nprof, zprofbgn, zprofend,                     &
           'Vort *10^5','height (km)',npicprof,tem2,tem3)
      IF(myproc == 0) CALL runlab (runname)

    END IF

    IF (divpprf == 1) THEN

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            tem1(i,j,k)= 1000.0 * (                                     &
                (u(i+1,j,k)-u(i,j,k))/(x(i+1)-x(i))+                    &
                (v(i,j+1,k)-v(i,j,k))/(y(j+1)-y(j)) )
          END DO
        END DO
      END DO

      nzprofc = nz-1
      CALL vprofil (nx,ny,nz,nzprofc,tem1,xc,yc,zpc,divppmin,divppmax,  &
           xprof,yprof,nprof,zprofbgn,zprofend,                         &
           'Div*1000', 'height (km)',npicprof,tem2,tem3)
      IF(myproc == 0) CALL runlab (runname)

    END IF

    IF (tsoilprof == 1) THEN

      DO k=1,nzsoil
        DO j=1,ny
          DO i=1,nx
            tem1(i,j,k)= tsoil(i,j,k,0)
          END DO
        END DO
      END DO

      nzprofc = nzsoil
      CALL vprofil (nx,ny,nzsoil,nzprofc,tem1,xc,yc,zpsoilc,            &
           tsoilprofmin,tsoilprofmax,                                   &
           xprof,yprof,nprof,zsoilprofbgn,zsoilprofend,                 &
           'tsoil', 'height (cm)',npicprof,tem2,tem3)
      IF(myproc == 0) CALL runlab (runname)

    END IF

    IF (qsoilprof == 1) THEN

      DO k=1,nzsoil
        DO j=1,ny
          DO i=1,nx
            tem1(i,j,k)= qsoil(i,j,k,0)
          END DO
        END DO
      END DO

      nzprofc = nzsoil
      CALL vprofil (nx,ny,nzsoil,nzprofc,tem1,xc,yc,zpsoilc,            &
           qsoilprofmin,qsoilprofmax,                                   &
           xprof,yprof,nprof,zsoilprofbgn,zsoilprofend,                 &
           'qsoil', 'height (cm)',npicprof,tem2,tem3)
      IF(myproc == 0) CALL runlab (runname)

    END IF

!
!-----------------------------------------------------------------------
!
!  Set up a new page if inwfrm = 1, i.e., flush the plot buffer and
!  move to the next plotting page.
!
!-----------------------------------------------------------------------
!
    900   CONTINUE

    IF( inwfrm == 1 ) CALL xnwfrm

    !IF ( ireturn == 0 .AND. hinfmt == 9 ) THEN
    !  GO TO 15              ! read one more time in the same GrADS file
    !ELSE
    !  CLOSE ( nchin )      ! read another data file
    !END IF


  END DO    ! data file loop

  IF(myproc ==0 ) THEN

    CALL xgrend

    WRITE(6,'(1x,a)') 'Program completed.'

    WRITE(6,'(/1x,a,e14.6,a)')                                          &
      'Total CPU time used =',f_cputime()-cpu0,'s.'

    WRITE (6,'(1x,a,F12.6)') 'Maxumum memory allocation (in words): ',  &
               max_memory_use
  END IF

  CALL mpexit(0)
END PROGRAM ARPSPLT
