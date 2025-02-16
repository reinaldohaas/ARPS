PROGRAM arpsdiff
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  PROGRAM ARPSDIFF                    ######
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
!
!  Computes difference between two ARPS history files and writes
!  the difference as a similar history file.
!
!  Reads in a history file produced by ARPS in any ARPS format.
!
!  Parameters grdout,varout,mstout,iceout and trbout should be input
!  with the same values as in the data dump subroutines in the model.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster  OU School of Meteorology. February 1992
!
!  MODIFICATION HISTORY:
!   11 Aug 1992  (KB) changed from arps2.5 to arps3.0
!   19 May 1992  (KB) changed from arps3.0 to arps3.2
!   24 May 1992  (KB) added code to handle staggered grid
!                     e.g., intrpu,intrpv,intrpw calls
!   12 Oct 1993  (KB) modified storage to save memory space
!                Differences are now output to interpolated arrays.
!
!   08/01/1995  (Ming Xue)
!   Changed the difference to be that between two total fields.
!   The output perturbation fields represent the difference in
!   the total fields.  The mean-state fields are not changed.
!
!   09/07/1995  (KB)
!   Added differencing and writing of surface (soil) fields.
!   This change affects the input file.
!
!   09/12/1995 (KB)
!   Cleared some bugs and restored some diagnostic printing.
!
!   03/13/1996 (Ming Xue)
!   Added array tke, atke, vtke. Change km to kmh and kmv.
!
!   04/02/1996 (Keith Brewster)
!   New runname is now read-in instead of output file names.
!   File names are constructed from runname just as in arps.
!
!   04/30/1996 (Keith Brewster)
!   Upgraded interpolation to n-order (1-to-fourth) polynomial
!   based on Gauss forward method.  Required some reorganization
!   to use the routine efficiently.  Added "iorder" as input
!   variable.
!
!   11/07/1996 (Keith Brewster)
!   Modified code for new interpolation scheme.
!
!   12/16/1996 (Yuhe Liu)
!   Corrected the dimension definitions for all verification of soil
!   and vegetation variables. They should be (vnx,vny) instead of
!   (nx,ny).
!
!   04/23/1998 (Yvette Richardson)
!   Corrected a bug affecting results when fdx.ne.fdy when
!   calculating fy0.
!
!   12/14/1998 (Donghai Wang)
!   Added the snow cover.
!
!   10/15/1999 (KB via Eric Kemp)
!   Corrected dimension statement of vxs and vys.  Affects only
!   calculation of vxs(vnx) and vys(vny).
!
!   11/03/1999 (KB via Eric Kemp)
!   Corrected dimension of snow cover and sfc flux variables
!   that had recently been added.
!
!   7 May 2002 (Eric Kemp)
!   Added allocation for array lon, and set nstyp = nstyps to
!   correctly read in soil data.  Also, added code to skip interpolation
!   when both grids are identical.
!
!   05/26/2002 (J. Brotzge)
!   Added/modified tsoil/qsoil for new soil scheme
!
!   1 June 2002 (Eric Kemp)
!   Soil model updates.
!
!   07/17/2002 (Yunheng Wang)
!   Added parameters nbeg, nend, nvbeg, nvend in subroutine intsclrs.
!
!   06/26/2012 (Youngsun Jung)
!   Fixed a problem that base state is added to the difference field.
!   
!-----------------------------------------------------------------------
!
!  DATA ARRAYS READ IN:
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       z coordinate of grid points in physical space (m)
!    zpsoil   z coordinate of grid points in the soil (m)
!
!    u        x component of velocity (m/s)
!    v        y component of velocity (m/s)
!    w        vertical component of velocity in Cartesian
!             coordinates (m/s).
!
!    ptprt    perturbation potential temperature (K)
!    pprt     perturbation pressure (Pascal)
!
!    qv       water vapor mixing ratio (kg/kg)
!    qscalar  Hydrometeor scalars
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
!    tsoil    Soil temperature (K)
!    qsoil    Soil moisture (m**3/m**3)
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
!    uprt    perturbation x component of velocity (m/s)
!    vprt    perturbation y component of velocity (m/s)
!    wprt    perturbation z component of velocity (m/s)
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx       ! Number of grid points in the x-direction
  INTEGER :: ny       ! Number of grid points in the y-direction
  INTEGER :: nz       ! Number of grid points in the z-direction
  INTEGER :: nzsoil   ! Number of grid points in the soil
!
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
!
!-----------------------------------------------------------------------
!
!  Arrays to be read in:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nstyps            ! Maximum number of soil types in each
                               ! grid box
!  PARAMETER (nstyps=4)

  REAL, ALLOCATABLE :: x    (:)          ! The x-coord. of the physical and
                                         ! computational grid. Defined at u-point.
  REAL, ALLOCATABLE :: y    (:)          ! The y-coord. of the physical and
                                         ! computational grid. Defined at v-point.
  REAL, ALLOCATABLE :: z    (:)          ! The z-coord. of the computational grid.
                                         ! Defined at w-point on the staggered grid.
  REAL, ALLOCATABLE :: zp   (:,:,:)      ! The physical height coordinate defined at
                                         ! w-point of the staggered grid.
  REAL, ALLOCATABLE :: zpsoil(:,:,:)     ! The physical height coordinate defined at
                                         ! w-point of the soil.

  REAL, ALLOCATABLE :: uprt   (:,:,:)    ! Perturbation u-velocity (m/s)
  REAL, ALLOCATABLE :: vprt   (:,:,:)    ! Perturbation v-velocity (m/s)
  REAL, ALLOCATABLE :: wprt   (:,:,:)    ! Perturbation w-velocity (m/s)
  REAL, ALLOCATABLE :: ptprt  (:,:,:)    ! Perturbation potential temperature (K)
  REAL, ALLOCATABLE :: pprt   (:,:,:)    ! Perturbation pressure (Pascal)
  REAL, ALLOCATABLE :: qvprt  (:,:,:)    ! Perturbation water vapor specific humidity
  REAL, ALLOCATABLE :: qscalar(:,:,:,:)

  REAL, ALLOCATABLE :: tke    (:,:,:)    ! Turbulent Kinetic Energy ((m/s)**2)
  REAL, ALLOCATABLE :: kmh    (:,:,:)    ! Horizontal turb. mixing coef. for
                                         ! momentum. ( m**2/s )
  REAL, ALLOCATABLE :: kmv    (:,:,:)    ! Vertical turb. mixing coef. for
                                         ! momentum. ( m**2/s )

  REAL, ALLOCATABLE :: ubar   (:,:,:)    ! Base state u-velocity (m/s)
  REAL, ALLOCATABLE :: vbar   (:,:,:)    ! Base state v-velocity (m/s)
  REAL, ALLOCATABLE :: wbar   (:,:,:)    ! Base state w-velocity (m/s)
  REAL, ALLOCATABLE :: ptbar  (:,:,:)    ! Base state potential temperature (K)
  REAL, ALLOCATABLE :: pbar   (:,:,:)    ! Base state pressure (Pascal)
  REAL, ALLOCATABLE :: rhobar (:,:,:)    ! Base state air density (kg/m**3)
  REAL, ALLOCATABLE :: qvbar  (:,:,:)    ! Base state water vapor specific humidity

  INTEGER, ALLOCATABLE :: soiltyp (:,:,:)! Soil type
  REAL, ALLOCATABLE :: stypfrct(:,:,:)   ! Soil type
  INTEGER, ALLOCATABLE :: vegtyp  (:,:)  ! Vegetation type
  REAL, ALLOCATABLE :: lai     (:,:)     ! Leaf Area Index
  REAL, ALLOCATABLE :: roufns  (:,:)     ! Surface roughness
  REAL, ALLOCATABLE :: veg     (:,:)     ! Vegetation fraction

  REAL, ALLOCATABLE :: tsoil  (:,:,:,:)    ! Soil temperature (K)
  REAL, ALLOCATABLE :: qsoil  (:,:,:,:)    ! Soil moisture (m**3/m**3)
  REAL, ALLOCATABLE :: wetcanp(:,:,:)    ! Canopy water amount
  REAL, ALLOCATABLE :: snowdpth(:,:)     ! Snow depth (m)

  REAL, ALLOCATABLE :: raing  (:,:)      ! Cumulus convective rain
  REAL, ALLOCATABLE :: rainc  (:,:)      ! Cumulus convective rain
  REAL, ALLOCATABLE :: prcrate(:,:,:)    ! precipitation rate (kg/(m**2*s))
                                         ! prcrate(1,1,1) = total precip. rate
                                         ! prcrate(1,1,2) = grid scale precip. rate
                                         ! prcrate(1,1,3) = cumulus precip. rate
                                         ! prcrate(1,1,4) = microphysics precip. rate

  REAL, ALLOCATABLE :: radfrc(:,:,:)     ! Radiation forcing (K/s)
  REAL, ALLOCATABLE :: radsw (:,:)       ! Solar radiation reaching the surface
  REAL, ALLOCATABLE :: rnflx (:,:)       ! Net radiation flux absorbed by surface
  REAL, ALLOCATABLE :: radswnet(:,:)     ! Net shortwave radiation
  REAL, ALLOCATABLE :: radlwin(:,:)      ! Incoming longwave radiation

  REAL, ALLOCATABLE :: usflx (:,:)       ! Surface flux of u-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: vsflx (:,:)       ! Surface flux of v-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: ptsflx(:,:)       ! Surface heat flux (K*kg/(m*s**2))
  REAL, ALLOCATABLE :: qvsflx(:,:)       ! Surface moisture flux (kg/(m**2*s))
!
!-----------------------------------------------------------------------
!
!  Verification Arrays
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: vx     (:)     ! The x-coord. of the physical and
                                      ! computational grid. Defined at u-point.
  REAL, ALLOCATABLE :: vy     (:)     ! The y-coord. of the physical and
                                      ! computational grid. Defined at v-point.
  REAL, ALLOCATABLE :: vz     (:)     ! The z-coord. of the computational grid.
                                      ! Defined at w-point on the staggered grid.
  REAL, ALLOCATABLE :: vzp    (:,:,:) ! The physical height coordinate defined at
                                      ! w-point of the staggered grid.
  REAL, ALLOCATABLE :: vzpsoil(:,:,:) ! The physical height coordinate defined at
                                      ! w-point of the soil grid.
  REAL, ALLOCATABLE :: vuprt  (:,:,:) ! Perturbation u-velocity (m/s)
  REAL, ALLOCATABLE :: vvprt  (:,:,:) ! Perturbation v-velocity (m/s)
  REAL, ALLOCATABLE :: vwprt  (:,:,:) ! Perturbation w-velocity (m/s)
  REAL, ALLOCATABLE :: vptprt (:,:,:) ! Perturbation potential temperature (K)
  REAL, ALLOCATABLE :: vpprt  (:,:,:) ! Perturbation pressure (Pascal)
  REAL, ALLOCATABLE :: vqvprt (:,:,:) ! Perturbation water vapor specific humidity
  REAL, ALLOCATABLE :: vqscalar(:,:,:,:)

  REAL, ALLOCATABLE :: vtke   (:,:,:) ! Turbulent Kinetic Energy ((m/s)**2)
  REAL, ALLOCATABLE :: vkmh   (:,:,:) ! Horizontal turb. mixing coef. for
                                      ! momentum. ( m**2/s )
  REAL, ALLOCATABLE :: vkmv   (:,:,:) ! Vertical turb. mixing coef. for
                                      ! momentum. ( m**2/s )
  REAL, ALLOCATABLE :: vubar  (:,:,:) ! Base state u-velocity (m/s)
  REAL, ALLOCATABLE :: vvbar  (:,:,:) ! Base state v-velocity (m/s)
  REAL, ALLOCATABLE :: vwbar  (:,:,:) ! Base state w-velocity (m/s)
  REAL, ALLOCATABLE :: vptbar (:,:,:) ! Base state potential temperature (K)
  REAL, ALLOCATABLE :: vpbar  (:,:,:) ! Base state pressure (Pascal)
  REAL, ALLOCATABLE :: vrhobar(:,:,:) ! Base state air density (kg/m**3)
  REAL, ALLOCATABLE :: vqvbar (:,:,:) ! Base state water vapor specific humidity

  INTEGER, ALLOCATABLE :: vsoiltyp (:,:,:)   ! Soil type
  REAL, ALLOCATABLE :: vstypfrct(:,:,:)      ! Soil type
  INTEGER, ALLOCATABLE :: vvegtyp  (:,:)     ! Vegetation type
  REAL, ALLOCATABLE :: vlai     (:,:)        ! Leaf Area Index
  REAL, ALLOCATABLE :: vroufns  (:,:)        ! Surface roughness
  REAL, ALLOCATABLE :: vveg     (:,:)        ! Vegetation fraction

  REAL, ALLOCATABLE :: vtsoil   (:,:,:,:)    ! Soil temperature (K)
  REAL, ALLOCATABLE :: vqsoil  (:,:,:,:)     ! Soil moisture (m**3/m**3)
  REAL, ALLOCATABLE :: vwetcanp(:,:,:)     ! Canopy water amount
  REAL, ALLOCATABLE :: vsnowdpth(:,:)      ! Snow depth (m)

  REAL, ALLOCATABLE :: vraing(:,:)         ! Grid supersaturation rain
  REAL, ALLOCATABLE :: vrainc(:,:)         ! Cumulus convective rain
  REAL, ALLOCATABLE :: vprcrate(:,:,:)     ! precipitation rate (kg/(m**2*s))
                                           ! prcrate(1,1,1) = total precip. rate
                                           ! prcrate(1,1,2) = grid scale precip. rate
                                           ! prcrate(1,1,3) = cumulus precip. rate
                                           ! prcrate(1,1,4) = microphysics precip. rate

  REAL, ALLOCATABLE :: vradfrc(:,:,:)      ! Radiation forcing (K/s)
  REAL, ALLOCATABLE :: vradsw (:,:)        ! Solar radiation reaching the surface
  REAL, ALLOCATABLE :: vrnflx (:,:)        ! Net radiation flux absorbed by surface
  REAL, ALLOCATABLE :: vradswnet(:,:)      ! Net shortwave radiation
  REAL, ALLOCATABLE :: vradlwin(:,:)       ! Incoming longwave radiation

  REAL, ALLOCATABLE :: vusflx (:,:)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: vvsflx (:,:)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: vptsflx(:,:)        ! Surface heat flux (K*kg/(m*s**2))
  REAL, ALLOCATABLE :: vqvsflx(:,:)        ! Surface moisture flux (kg/(m**2*s))
!
!-----------------------------------------------------------------------
!
!  Verification data interpolated to model grid
!  These arrays also hold difference fields after call to diffgr.
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: auprt  (:,:,:)    ! Perturbation u-velocity (m/s)
  REAL, ALLOCATABLE :: avprt  (:,:,:)    ! Perturbation v-velocity (m/s)
  REAL, ALLOCATABLE :: awprt  (:,:,:)    ! Perturbation w-velocity (m/s)
  REAL, ALLOCATABLE :: aptprt (:,:,:)    ! Perturbation potential temperature (K)
  REAL, ALLOCATABLE :: apprt  (:,:,:)    ! Perturbation pressure (Pascal)
  REAL, ALLOCATABLE :: aqvprt (:,:,:)    ! Perturbation water vapor specific humidity
  REAL, ALLOCATABLE :: aqscalar (:,:,:,:)

  REAL, ALLOCATABLE :: atke   (:,:,:)    ! Turbulent Kinetic Energy ((m/s)**2)
  REAL, ALLOCATABLE :: akmh   (:,:,:)    ! Horizontal turb. mixing coef. for
                                         ! momentum. ( m**2/s )
  REAL, ALLOCATABLE :: akmv   (:,:,:)    ! Vertical turb. mixing coef. for
                                         ! momentum. ( m**2/s )

  REAL, ALLOCATABLE :: aubar  (:,:,:)    ! Base state u-velocity (m/s)
  REAL, ALLOCATABLE :: avbar  (:,:,:)    ! Base state v-velocity (m/s)
  REAL, ALLOCATABLE :: awbar  (:,:,:)    ! Base state w-velocity (m/s)
  REAL, ALLOCATABLE :: aptbar (:,:,:)    ! Base state potential temperature (K)
  REAL, ALLOCATABLE :: apbar  (:,:,:)    ! Base state pressure (Pascal)
  REAL, ALLOCATABLE :: arhobar(:,:,:)    ! Base state air density (kg/m**3)
  REAL, ALLOCATABLE :: aqvbar (:,:,:)    ! Base state water vapor specific humidity

  INTEGER, ALLOCATABLE :: asoiltyp (:,:,:)   ! Soil type
  REAL, ALLOCATABLE :: astypfrct(:,:,:)      ! Soil type
  INTEGER, ALLOCATABLE :: avegtyp  (:,:)     ! Vegetation type
  REAL, ALLOCATABLE :: alai     (:,:)        ! Leaf Area Index
  REAL, ALLOCATABLE :: aroufns  (:,:)        ! Surface roughness
  REAL, ALLOCATABLE :: aveg     (:,:)        ! Vegetation fraction

  REAL, ALLOCATABLE :: atsoil  (:,:,:,:)     ! Soil temperature (K)
  REAL, ALLOCATABLE :: aqsoil  (:,:,:,:)     ! Soil moisture (m**3/m**3)
  REAL, ALLOCATABLE :: awetcanp(:,:,:)     ! Canopy water amount
  REAL, ALLOCATABLE :: asnowdpth(:,:)      ! Snow depth (m)

  REAL, ALLOCATABLE :: araing (:,:)        ! Grid supersaturation rain
  REAL, ALLOCATABLE :: arainc (:,:)        ! Cumulus convective rain
  REAL, ALLOCATABLE :: aprcrate(:,:,:)     ! precipitation rate (kg/(m**2*s))
                                           ! prcrate(1,1,1) = total precip. rate
                                           ! prcrate(1,1,2) = grid scale precip. rate
                                           ! prcrate(1,1,3) = cumulus precip. rate
                                           ! prcrate(1,1,4) = microphysics precip. rate

  REAL, ALLOCATABLE :: aradfrc(:,:,:)     ! Radiation forcing (K/s)
  REAL, ALLOCATABLE :: aradsw (:,:)       ! Solar radiation reaching the surface
  REAL, ALLOCATABLE :: arnflx (:,:)       ! Net radiation flux absorbed by surface
  REAL, ALLOCATABLE :: aradswnet(:,:)     ! Net shortwave radiation
  REAL, ALLOCATABLE :: aradlwin(:,:)      ! Incoming longwave radiation

  REAL, ALLOCATABLE :: ausflx (:,:)       ! Surface flux of u-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: avsflx (:,:)       ! Surface flux of v-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: aptsflx(:,:)       ! Surface heat flux (K*kg/(m*s**2))
  REAL, ALLOCATABLE :: aqvsflx(:,:)       ! Surface moisture flux (kg/(m**2*s))
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
  REAL, ALLOCATABLE :: tem3dsoil(:,:,:)

  REAL, ALLOCATABLE :: vtem1(:,:,:)
  REAL, ALLOCATABLE :: vtem2(:,:,:)
  REAL, ALLOCATABLE :: vtem3(:,:,:)

  REAL, ALLOCATABLE :: xs(:)
  REAL, ALLOCATABLE :: ys(:)
  REAL, ALLOCATABLE :: zps(:,:,:)
  REAL, ALLOCATABLE :: x2d(:,:)
  REAL, ALLOCATABLE :: y2d(:,:)
  REAL, ALLOCATABLE :: lat(:,:),lon(:,:)

  REAL, ALLOCATABLE :: vxs(:)
  REAL, ALLOCATABLE :: vys(:)
  REAL, ALLOCATABLE :: vzps(:,:,:)
  REAL, ALLOCATABLE :: dxfld(:)
  REAL, ALLOCATABLE :: dyfld(:)
  REAL, ALLOCATABLE :: rdxfld(:)
  REAL, ALLOCATABLE :: rdyfld(:)

  INTEGER, ALLOCATABLE :: iloc(:,:),jloc(:,:)
  REAL, ALLOCATABLE :: zpver(:,:,:)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: filename,vfilename,grdbasfn,vgrdbasfn,fcrnam,runnmin
  INTEGER :: lengbf,lenfil
  INTEGER :: ifproj,ivproj
  REAL :: flatnot(2),vlatnot(2)
  REAL :: fscale,ftrulon,fdx,fdy,fx0,fy0
  REAL :: fctrlat,fctrlon
  REAL :: vscale,vtrulon,vdx,vdy,vx0,vy0
  REAL :: vctrlat,vctrlon
  REAL :: time,xctr,yctr
  INTEGER :: i,j,k
  INTEGER :: grdbas
  INTEGER :: iorder,hinfmt,vhinfmt
  INTEGER :: ireturn
  LOGICAL :: comcoord,comcoord2
  INTEGER :: nch

!-----------------------------------------------------------------------
!  The following defines the model dimension parameters for the
!  verification data.
!-----------------------------------------------------------------------
!
  INTEGER :: vnstyps            ! Maximum number of soil types in each
                                ! grid box
  INTEGER :: vnx       ! Number of grid points in the x-direction
  INTEGER :: vny       ! Number of grid points in the y-direction
  INTEGER :: vnz       ! Number of grid points in the z-direction
  INTEGER :: vnzsoil   ! Number of grid points in the soil.


  INTEGER :: vnxy, vnxz, vnyz, vnxyz

  INTEGER :: bnscalar, bnscalarq
  INTEGER :: bP_QC, bP_QR, bP_QI, bP_QS, bP_QH, bP_QG
  INTEGER :: bP_NC, bP_NR, bP_NI, bP_NS, bP_NH, bP_NG
  INTEGER ::        bP_ZR, bP_ZI, bP_ZS, bP_ZH, bP_ZG

  CHARACTER(LEN=40) :: bqnames(20)
  CHARACTER(LEN=40) :: bqdescp(20)

  INTEGER :: vnscalar, vnscalarq
  INTEGER :: vP_QC, vP_QR, vP_QI, vP_QS, vP_QH, vP_QG
  INTEGER :: vP_NC, vP_NR, vP_NI, vP_NS, vP_NH, vP_NG
  INTEGER ::        vP_ZR, vP_ZI, vP_ZS, vP_ZH, vP_ZG

  CHARACTER(LEN=40) :: vqnames(20)
  CHARACTER(LEN=40) :: vqdescp(20)

  INTEGER :: anscalar, anscalarq
  INTEGER :: aP_QC, aP_QR, aP_QI, aP_QS, aP_QH, aP_QG
  INTEGER :: aP_NC, aP_NR, aP_NI, aP_NS, aP_NH, aP_NG
  INTEGER ::        aP_ZR, aP_ZI, aP_ZS, aP_ZH, aP_ZG

  CHARACTER(LEN=40) :: aqnames(20)
  CHARACTER(LEN=40) :: aqdescp(20)

  INTEGER, ALLOCATABLE :: vQindex(:), bQindex(:)

  INTEGER :: nq

!-----------------------------------------------------------------------
  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  namelist Declarations:
!
!-----------------------------------------------------------------------
!
!  NAMELIST /grid_dims/ nx,ny,nz

!  NAMELIST /vgrid_dims/ vnx,vny,vnz

  NAMELIST /frcst_fn/ iorder, hinfmt, grdbasfn, filename

  NAMELIST /vrftn_fn/ vhinfmt, vgrdbasfn, vfilename

  NAMELIST /output/ runnmin, hdmpfmt,                      &
                    grdout,   basout,    varout,           &
                    mstout,   iceout,    trbout,           &
                    sfcout,   rainout,   prcout,           &
                    radout,   flxout,    filcmprs

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!  nstyp = nstyps ! EMK

  WRITE(6,'(/10(/5x,a))')                                                 &
      '#################################################################',&
      '#################################################################',&
      '####                                                         ####',&
      '####                  Welcome to ARPSDIFF                    ####',&
      '####                                                         ####',&
      '####   A program that reads in history files generated       ####',&
      '####   by ARPS and produces difference grids.                ####',&
      '####                                                         ####',&
      '#################################################################',&
      '#################################################################'
!
!-----------------------------------------------------------------------
!
!  Set up the default values for all the variables to be read in
!  using the namelist method. In case the user does not specify a
!  particular value, this value will be used.
!
!-----------------------------------------------------------------------
!
   nx   = 67
   ny   = 67
   nz   = 35
   nzsoil=10

   vnx   = 163
   vny   = 99
   vnz   = 43

   iorder   = 3
   hinfmt   = 10
   grdbasfn = 'may20.grbgrdbas'
   filename = 'may20.grb003600'

   runnmin = 'diffm20'
   hdmpfmt = 10
   grdout = 0
   basout = 0
   varout = 1
   mstout = 1
   iceout = 0
   trbout = 0
   sfcout = 0
   rainout = 0
   prcout = 0
   radout = 0
   flxout = 0
   filcmprs = 0

!-----------------------------------------------------------------------
!
!  Read in filenames, as well as nx,ny,nz and vnx, vny, vnz
!
!----------------------------------------------------------------------

!   READ(5,grid_dims,END=100)
!   WRITE(6,'(a)')'Namelist block grid_dims successfully read in.'

!   READ(5,vgrid_dims,END=100)
!   WRITE(6,'(a)')'Namelist block vgrid_dims successfully read in.'

  READ(5,frcst_fn,END=100)
  WRITE(6,'(a)')'Namelist block frcst_fn successfully read in.'

  lengbf=LEN_trim(grdbasfn)
  WRITE(6,'(/a,a)')' The grid/base name is ', grdbasfn(1:lengbf)

  lenfil = LEN_trim(filename)
  WRITE(6,'(/a,/1x,a)')' The data set name is ', filename(1:lenfil)

  CALL get_dims_from_data(hinfmt,filename(1:lenfil),                    &
                          nx,ny,nz,nzsoil,nstyps, ireturn)

  bnscalar  = nscalar
  bnscalarq = nscalarq
  bP_QC = P_QC
  bP_QR = P_QR
  bP_QI = P_QI
  bP_QS = P_QS
  bP_QG = P_QG
  bP_QH = P_QH
  bP_NC = P_NC
  bP_NR = P_NR
  bP_NI = P_NI
  bP_NS = P_NS
  bP_NG = P_NG
  bP_NH = P_NH
  bP_ZR = P_ZR
  bP_ZI = P_ZI
  bP_ZS = P_ZS
  bP_ZG = P_ZG
  bP_ZH = P_ZH
  bqnames(:) = qnames(:)
  bqdescp(:) = qdescp(:)

  IF (nstyps <= 0) nstyps = 1
  nstyp = nstyps ! Copy to global

  READ(5,vrftn_fn,END=100)
  WRITE(6,'(a)')'Namelist block vrftn_fn successfully read in.'

  lengbf=LEN_trim(vgrdbasfn)
  WRITE(6,'(/a,a)')' The grid/base name is ', vgrdbasfn(1:lengbf)

  lenfil=LEN_trim(vfilename)
  WRITE(6,'(/a,a)')' The data set name is ', vfilename(1:lenfil)

  CALL get_dims_from_data(vhinfmt,vfilename(1:lenfil),                  &
                          vnx,vny,vnz,vnzsoil,vnstyps, ireturn)

  IF (vnstyps <= 0) vnstyps = 1
  nstyp = vnstyps ! Copy to global

  vnscalar  = nscalar
  vnscalarq = nscalarq
  vP_QC = P_QC
  vP_QR = P_QR
  vP_QI = P_QI
  vP_QS = P_QS
  vP_QG = P_QG
  vP_QH = P_QH
  vP_NC = P_NC
  vP_NR = P_NR
  vP_NI = P_NI
  vP_NS = P_NS
  vP_NG = P_NG
  vP_NH = P_NH
  vP_ZR = P_ZR
  vP_ZI = P_ZI
  vP_ZS = P_ZS
  vP_ZG = P_ZG
  vP_ZH = P_ZH
  vqnames(:) = qnames(:)
  vqdescp(:) = qdescp(:)

  vnxy=vnx*vny
  vnxz=vnx*vnz
  vnyz=vny*vnz
  vnxyz=vnxy*vnz

!-----------------------------------------------------------------------
!
! Determine Hydrometeor parameters based on input
!
!-----------------------------------------------------------------------

  anscalar  = 0
  anscalarq = 0
  aP_QC = -1; aP_QR = -1; aP_QI = -1; aP_QS = -1; aP_QH = -1; aP_QG = -1
  aP_NC = -1; aP_NR = -1; aP_NI = -1; aP_NS = -1; aP_NG = -1; aP_NH = -1
              aP_ZR = -1; aP_ZI = -1; aP_ZS = -1; aP_ZG = -1; aP_ZH = -1;
  aqnames(:)= ' '; aqdescp(:)= ' '

  ALLOCATE(vQindex(bnscalar), STAT = istatus )
  ALLOCATE(bQindex(bnscalar), STAT = istatus )
  vQindex(:) = -1
  bQindex(:) = -1

  CALL min_set_scalars(bnscalar,bnscalarq,bqnames,bqdescp,              &
                       bP_QC,bP_QR,bP_QI,bP_QS,bP_QG,bP_QH,             &
                       bP_NC,bP_NR,bP_NI,bP_NS,bP_NG,bP_NH,             &
                             bP_ZR,bP_ZI,bP_ZS,bP_ZG,bP_ZH,             &
                       vnscalar,vnscalarq,vqnames,vqdescp,              &
                       vP_QC,vP_QR,vP_QI,vP_QS,vP_QG,vP_QH,             &
                       vP_NC,vP_NR,vP_NI,vP_NS,vP_NG,vP_NH,             &
                             vP_ZR,vP_ZI,vP_ZS,vP_ZG,vP_ZH,             &
                       anscalar,anscalarq,aqnames,aqdescp,              &
                       aP_QC,aP_QR,aP_QI,aP_QS,aP_QG,aP_QH,             &
                       aP_NC,aP_NR,aP_NI,aP_NS,aP_NG,aP_NH,             &
                             aP_ZR,aP_ZI,aP_ZS,aP_ZG,aP_ZH,             &
                       bQindex,vQindex, istatus )

!-----------------------------------------------------------------------
!  Allocate and Initialize the variables
!----------------------------------------------------------------------

  istatus=0

  ALLOCATE(x(nx),STAT=istatus)
  x=0
  ALLOCATE(y(ny),STAT=istatus)
  y=0
  ALLOCATE(z(nz),STAT=istatus)
  z=0
  ALLOCATE(zp(nx,ny,nz),STAT=istatus)
  zp=0
  ALLOCATE(zpsoil(nx,ny,nzsoil),STAT=istatus)
  zpsoil=0

  ALLOCATE(uprt(nx,ny,nz),STAT=istatus)
  uprt=0
  ALLOCATE(vprt(nx,ny,nz),STAT=istatus)
  vprt=0
  ALLOCATE(wprt(nx,ny,nz),STAT=istatus)
  wprt=0
  ALLOCATE(ptprt(nx,ny,nz),STAT=istatus)
  ptprt=0
  ALLOCATE(pprt(nx,ny,nz),STAT=istatus)
  pprt=0
  ALLOCATE(qvprt(nx,ny,nz),STAT=istatus)
  qvprt=0
  ALLOCATE(qscalar(nx,ny,nz,bnscalar),STAT=istatus)
  qscalar = 0.0
  ALLOCATE(tke(nx,ny,nz),STAT=istatus)
  tke=0
  ALLOCATE(kmh(nx,ny,nz),STAT=istatus)
  kmh=0
  ALLOCATE(kmv(nx,ny,nz),STAT=istatus)
  kmv=0
  ALLOCATE(ubar(nx,ny,nz),STAT=istatus)
  ubar=0
  ALLOCATE(vbar(nx,ny,nz),STAT=istatus)
  vbar=0
  ALLOCATE(wbar(nx,ny,nz),STAT=istatus)
  wbar=0
  ALLOCATE(ptbar(nx,ny,nz),STAT=istatus)
  ptbar=0
  ALLOCATE(pbar(nx,ny,nz),STAT=istatus)
  pbar=0
  ALLOCATE(rhobar(nx,ny,nz),STAT=istatus)
  rhobar=0
  ALLOCATE(qvbar(nx,ny,nz),STAT=istatus)
  qvbar=0
  ALLOCATE(soiltyp(nx,ny,nstyps),STAT=istatus)
  soiltyp=0
  ALLOCATE(stypfrct(nx,ny,nstyps),STAT=istatus)
  stypfrct=0
  ALLOCATE(vegtyp(nx,ny),STAT=istatus)
  vegtyp=0
  ALLOCATE(lai(nx,ny),STAT=istatus)
  lai=0
  ALLOCATE(roufns(nx,ny),STAT=istatus)
  roufns=0
  ALLOCATE(veg(nx,ny),STAT=istatus)
  veg=0
  ALLOCATE(tsoil(nx,ny,nzsoil,0:nstyps),STAT=istatus)
  tsoil=0
  ALLOCATE(qsoil(nx,ny,nzsoil,0:nstyps),STAT=istatus)
  qsoil=0
  ALLOCATE(wetcanp(nx,ny,0:nstyps),STAT=istatus)
  wetcanp=0
  ALLOCATE(snowdpth(nx,ny),STAT=istatus)
  snowdpth=0
  ALLOCATE(raing(nx,ny),STAT=istatus)
  raing=0
  ALLOCATE(rainc(nx,ny),STAT=istatus)
  rainc=0
  ALLOCATE(prcrate(nx,ny,4),STAT=istatus)
  prcrate=0
  ALLOCATE(radfrc(nx,ny,nz),STAT=istatus)
  radfrc=0
  ALLOCATE(radsw(nx,ny),STAT=istatus)
  radsw=0
  ALLOCATE(rnflx(nx,ny),STAT=istatus)
  rnflx=0
  ALLOCATE(radswnet(nx,ny),STAT=istatus)
  radswnet=0
  ALLOCATE(radlwin(nx,ny),STAT=istatus)
  radlwin=0
  ALLOCATE(usflx(nx,ny),STAT=istatus)
  usflx=0
  ALLOCATE(vsflx(nx,ny),STAT=istatus)
  vsflx=0
  ALLOCATE(ptsflx(nx,ny),STAT=istatus)
  ptsflx=0
  ALLOCATE(qvsflx(nx,ny),STAT=istatus)
  qvsflx=0

  ALLOCATE(vx(vnx),STAT=istatus)
  vx=0
  ALLOCATE(vy(vny),STAT=istatus)
  vy=0
  ALLOCATE(vz(vnz),STAT=istatus)
  vz=0
  ALLOCATE(vzp(vnx,vny,vnz),STAT=istatus)
  vzp=0
  ALLOCATE(vzpsoil(vnx,vny,vnzsoil),STAT=istatus)
  vzpsoil=0
  ALLOCATE(vuprt(vnx,vny,vnz),STAT=istatus)
  vuprt=0
  ALLOCATE(vvprt(vnx,vny,vnz),STAT=istatus)
  vvprt=0
  ALLOCATE(vwprt(vnx,vny,vnz),STAT=istatus)
  vwprt=0
  ALLOCATE(vptprt(vnx,vny,vnz),STAT=istatus)
  vptprt=0
  ALLOCATE(vpprt(vnx,vny,vnz),STAT=istatus)
  vpprt=0
  ALLOCATE(vqvprt(vnx,vny,vnz),STAT=istatus)
  vqvprt=0
  ALLOCATE(vqscalar(vnx,vny,vnz,vnscalar),STAT=istatus)
  vqscalar = 0.0
  ALLOCATE(vtke(vnx,vny,vnz),STAT=istatus)
  vtke=0
  ALLOCATE(vkmh(vnx,vny,vnz),STAT=istatus)
  vkmh=0
  ALLOCATE(vkmv(vnx,vny,vnz),STAT=istatus)
  vkmv=0
  ALLOCATE(vubar(vnx,vny,vnz),STAT=istatus)
  vubar=0
  ALLOCATE(vvbar(vnx,vny,vnz),STAT=istatus)
  vvbar=0
  ALLOCATE(vwbar(vnx,vny,vnz),STAT=istatus)
  vwbar=0
  ALLOCATE(vptbar(vnx,vny,vnz),STAT=istatus)
  vptbar=0
  ALLOCATE(vpbar(vnx,vny,vnz),STAT=istatus)
  vpbar=0
  ALLOCATE(vrhobar(vnx,vny,vnz),STAT=istatus)
  vrhobar=0
  ALLOCATE(vqvbar(vnx,vny,vnz),STAT=istatus)
  vqvbar=0
  ALLOCATE(vsoiltyp(vnx,vny,vnstyps),STAT=istatus)
  vsoiltyp=0
  ALLOCATE(vstypfrct(vnx,vny,vnstyps),STAT=istatus)
  vstypfrct=0
  ALLOCATE(vvegtyp(vnx,vny),STAT=istatus)
  vvegtyp=0
  ALLOCATE(vlai(vnx,vny),STAT=istatus)
  vlai=0
  ALLOCATE(vroufns(vnx,vny),STAT=istatus)
  vroufns=0
  ALLOCATE(vveg(vnx,vny),STAT=istatus)
  vveg=0
  ALLOCATE(vtsoil(vnx,vny,vnzsoil,0:vnstyps),STAT=istatus)
  vtsoil=0
  ALLOCATE(vqsoil(vnx,vny,vnzsoil,0:vnstyps),STAT=istatus)
  vqsoil=0
  ALLOCATE(vwetcanp(vnx,vny,0:vnstyps),STAT=istatus)
  vwetcanp=0
  ALLOCATE(vsnowdpth(vnx,vny),STAT=istatus)
  vsnowdpth=0
  ALLOCATE(vraing(vnx,vny),STAT=istatus)
  vraing=0
  ALLOCATE(vrainc(vnx,vny),STAT=istatus)
  vrainc=0
  ALLOCATE(vprcrate(vnx,vny,4),STAT=istatus)
  vprcrate=0
  ALLOCATE(vradfrc(vnx,vny,vnz),STAT=istatus)
  vradfrc=0
  ALLOCATE(vradsw(vnx,vny),STAT=istatus)
  vradsw=0
  ALLOCATE(vrnflx(vnx,vny),STAT=istatus)
  vrnflx=0
  ALLOCATE(vradswnet(vnx,vny),STAT=istatus)
  vradswnet=0
  ALLOCATE(vradlwin(vnx,vny),STAT=istatus)
  vradlwin=0
  ALLOCATE(vusflx(vnx,vny),STAT=istatus)
  vusflx=0
  ALLOCATE(vvsflx(vnx,vny),STAT=istatus)
  vvsflx=0
  ALLOCATE(vptsflx(vnx,vny),STAT=istatus)
  vptsflx=0
  ALLOCATE(vqvsflx(vnx,vny),STAT=istatus)
  vqvsflx=0

  ALLOCATE(auprt(nx,ny,nz),STAT=istatus)
  auprt=0
  ALLOCATE(avprt(nx,ny,nz),STAT=istatus)
  avprt=0
  ALLOCATE(awprt(nx,ny,nz),STAT=istatus)
  awprt=0
  ALLOCATE(aptprt(nx,ny,nz),STAT=istatus)
  aptprt=0
  ALLOCATE(apprt(nx,ny,nz),STAT=istatus)
  apprt=0
  ALLOCATE(aqvprt(nx,ny,nz),STAT=istatus)
  aqvprt=0
  ALLOCATE(aqscalar(nx,ny,nz,anscalar),STAT=istatus)
  aqscalar = 0.0
  ALLOCATE(atke(nx,ny,nz),STAT=istatus)
  atke=0
  ALLOCATE(akmh(nx,ny,nz),STAT=istatus)
  akmh=0
  ALLOCATE(akmv(nx,ny,nz),STAT=istatus)
  akmv=0
  ALLOCATE(aubar(nx,ny,nz),STAT=istatus)
  aubar=0
  ALLOCATE(avbar(nx,ny,nz),STAT=istatus)
  avbar=0
  ALLOCATE(awbar(nx,ny,nz),STAT=istatus)
  awbar=0
  ALLOCATE(aptbar(nx,ny,nz),STAT=istatus)
  aptbar=0
  ALLOCATE(apbar(nx,ny,nz),STAT=istatus)
  apbar=0
  ALLOCATE(arhobar(nx,ny,nz),STAT=istatus)
  arhobar=0
  ALLOCATE(aqvbar(nx,ny,nz),STAT=istatus)
  aqvbar=0
  ALLOCATE(asoiltyp(nx,ny,nstyps),STAT=istatus)
  asoiltyp=0
  ALLOCATE(astypfrct(nx,ny,nstyps),STAT=istatus)
  astypfrct=0
  ALLOCATE(avegtyp(nx,ny),STAT=istatus)
  avegtyp=0
  ALLOCATE(alai(nx,ny),STAT=istatus)
  alai=0
  ALLOCATE(aroufns(nx,ny),STAT=istatus)
  aroufns=0
  ALLOCATE(aveg(nx,ny),STAT=istatus)
  aveg=0
  ALLOCATE(atsoil(nx,ny,nzsoil,0:nstyps),STAT=istatus)
  atsoil=0
  ALLOCATE(aqsoil(nx,ny,nzsoil,0:nstyps),STAT=istatus)
  aqsoil=0
  ALLOCATE(awetcanp(nx,ny,0:nstyps),STAT=istatus)
  awetcanp=0
  ALLOCATE(asnowdpth(nx,ny),STAT=istatus)
  asnowdpth=0
  ALLOCATE(araing(nx,ny),STAT=istatus)
  araing=0
  ALLOCATE(arainc(nx,ny),STAT=istatus)
  arainc=0
  ALLOCATE(aprcrate(nx,ny,4),STAT=istatus)
  aprcrate=0
  ALLOCATE(aradfrc(nx,ny,nz),STAT=istatus)
  aradfrc=0
  ALLOCATE(aradsw(nx,ny),STAT=istatus)
  aradsw=0
  ALLOCATE(arnflx(nx,ny),STAT=istatus)
  arnflx=0
  ALLOCATE(aradswnet(nx,ny),STAT=istatus)
  aradswnet=0
  ALLOCATE(aradlwin(nx,ny),STAT=istatus)
  aradlwin=0
  ALLOCATE(ausflx(nx,ny),STAT=istatus)
  ausflx=0
  ALLOCATE(avsflx(nx,ny),STAT=istatus)
  avsflx=0
  ALLOCATE(aptsflx(nx,ny),STAT=istatus)
  aptsflx=0
  ALLOCATE(aqvsflx(nx,ny),STAT=istatus)
  aqvsflx=0
  ALLOCATE(tem1(nx,ny,nz),STAT=istatus)
  tem1=0
  ALLOCATE(tem2(nx,ny,nz),STAT=istatus)
  tem2=0
  ALLOCATE(tem3(nx,ny,nz),STAT=istatus)
  tem3=0
  ALLOCATE(tem3dsoil(nx,ny,nzsoil),STAT=istatus)
  tem3dsoil=0
  ALLOCATE(vtem1(vnx,vny,vnz),STAT=istatus)
  vtem1=0
  ALLOCATE(vtem2(vnx,vny,vnz),STAT=istatus)
  vtem2=0
  ALLOCATE(vtem3(vnx,vny,vnz),STAT=istatus)
  vtem3=0
  ALLOCATE(xs(nx),STAT=istatus)
  xs=0
  ALLOCATE(ys(ny),STAT=istatus)
  ys=0
  ALLOCATE(zps(nx,ny,nz),STAT=istatus)
  zps=0
  ALLOCATE(x2d(nx,ny),STAT=istatus)
  x2d=0
  ALLOCATE(y2d(nx,ny),STAT=istatus)
  y2d=0
  ALLOCATE(lat(nx,ny),STAT=istatus)
  lat=0
  ALLOCATE(lon(nx,ny),STAT=istatus) ! EMK
  lon=0                             ! EMK
  ALLOCATE(vxs(vnx),STAT=istatus)
  vxs=0
  ALLOCATE(vys(vny),STAT=istatus)
  vys=0
  ALLOCATE(vzps(vnx,vny,vnz),STAT=istatus)
  vzps=0
  ALLOCATE(dxfld(vnx),STAT=istatus)
  dxfld=0
  ALLOCATE(dyfld(vny),STAT=istatus)
  dyfld=0
  ALLOCATE(rdxfld(vnx),STAT=istatus)
  rdxfld=0
  ALLOCATE(rdyfld(vny),STAT=istatus)
  rdyfld=0
  ALLOCATE(iloc(nx,ny),STAT=istatus)
  iloc=0
  ALLOCATE(jloc(nx,ny),STAT=istatus)
  jloc=0
  ALLOCATE(zpver(nx,ny,vnz),STAT=istatus)
  zpver=0

  mgrid = 1
  nestgrd = 0
  grbpkbit = 16

!
!-----------------------------------------------------------------------
!
!  Read all input data arrays
!
!-----------------------------------------------------------------------
!
  lengbf=LEN_trim(grdbasfn)
  lenfil = LEN_trim(filename)

  IF (nstyps <= 0) nstyps = 1
  nstyp = nstyps ! Copy to global for dtaread

  nscalar  = bnscalar
  nscalarq = bnscalarq
  P_QC     = bP_QC
  P_QR     = bP_QR
  P_QI     = bP_QI
  P_QS     = bP_QS
  P_QG     = bP_QG
  P_QH     = bP_QH
  P_NC     = bP_NC
  P_NR     = bP_NR
  P_NI     = bP_NI
  P_NS     = bP_NS
  P_NG     = bP_NG
  P_NH     = bP_NH
  P_ZR     = bP_ZR
  P_ZI     = bP_ZI
  P_ZS     = bP_ZS
  P_ZG     = bP_ZG
  P_ZH     = bP_ZH
  qnames(:) = bqnames(:)

  CALL dtaread(nx,ny,nz,nzsoil,nstyps,                                  &
               hinfmt,nch,grdbasfn(1:lengbf),lengbf,                    &
               filename(1:lenfil),lenfil,time,                          &
               x,y,z,zp,zpsoil,uprt ,vprt ,wprt ,ptprt, pprt ,          &
               qvprt, qscalar, tke,kmh,kmv,                             &
               ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,            &
               soiltyp,stypfrct,vegtyp,lai,roufns,veg,                  &
               tsoil,qsoil,wetcanp,snowdpth,                            &
               raing,rainc,prcrate,                                     &
               radfrc,radsw,rnflx,radswnet,radlwin,                     &
               usflx,vsflx,ptsflx,qvsflx,                               &
               ireturn, tem1,tem2,tem3)
!
!-----------------------------------------------------------------------
!
!  ireturn = 0 for a successful read
!
!-----------------------------------------------------------------------
!
  IF( ireturn == 0 ) THEN   ! successful read
    curtim=time
    fcrnam=runname
    ifproj=mapproj
    fscale=sclfct
    flatnot(1)=trulat1
    flatnot(2)=trulat2
    ftrulon=trulon
    fdx=x(3)-x(2)
    fdy=y(3)-y(2)
    fctrlat=ctrlat
    fctrlon=ctrlon
    CALL setmapr(ifproj,fscale,flatnot,ftrulon)
    CALL lltoxy(1,1,fctrlat,fctrlon,xctr,yctr)
    fx0=xctr-fdx*((nx-3)/2)
    fy0=yctr-fdy*((ny-3)/2)
    CALL setorig(1,fx0,fy0)
!
!-----------------------------------------------------------------------
!
!  Establish coordinate for scalar forecast fields.
!
!-----------------------------------------------------------------------
!
    DO i=1,nx-1
      xs(i)=0.5*(x(i)+x(i+1))
    END DO
    xs(nx)=2.*xs(nx-1)-xs(nx-2)
    DO j=1,ny-1
      ys(j)=0.5*(y(j)+y(j+1))
    END DO
    ys(ny)=2.*ys(ny-1)-ys(ny-2)
    CALL xytoll(nx,ny,xs,ys,lat,lon)
    DO k=1,nz-1
      DO j=1,ny
        DO i=1,nx
          zps(i,j,k)=0.5*(zp(i,j,k)+zp(i,j,k+1))
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Set the gridread parameter to 0 so that the verification
!  grid/base file will be read.
!
!-----------------------------------------------------------------------
!
    CALL setgbrd (0)
!
!-----------------------------------------------------------------------
!
!  Read in the verification data.
!
!-----------------------------------------------------------------------
!
    lengbf=LEN_trim(vgrdbasfn)
    lenfil=LEN_trim(vfilename)

    IF (vnstyps <= 0) vnstyps = 1
    nstyp = vnstyps ! Copy to global for dtaread

    nscalar   = vnscalar
    nscalarq  = vnscalarq
    P_QC      = vP_QC
    P_QR      = vP_QR
    P_QI      = vP_QI
    P_QS      = vP_QS
    P_QG      = vP_QG
    P_QH      = vP_QH
    P_NC      = vP_NC
    P_NR      = vP_NR
    P_NI      = vP_NI
    P_NS      = vP_NS
    P_NG      = vP_NG
    P_NH      = vP_NH
    P_ZR      = vP_ZR
    P_ZI      = vP_ZI
    P_ZS      = vP_ZS
    P_ZG      = vP_ZG
    P_ZH      = vP_ZH
    qnames(:) = vqnames(:)

    CALL dtaread(vnx,vny,vnz,vnzsoil,vnstyps,                           &
                 vhinfmt,nch,vgrdbasfn(1:lengbf),lengbf,                &
                 vfilename(1:lenfil),lenfil,time,                       &
                 vx,vy,vz,vzp,vzpsoil, vuprt ,vvprt ,vwprt ,vptprt,     &
                 vpprt,vqvprt, vqscalar, vtke,vkmh,vkmv,                &
                 vubar, vvbar, vwbar, vptbar, vpbar, vrhobar, vqvbar,   &
                 vsoiltyp,vstypfrct,vvegtyp,vlai,vroufns,vveg,          &
                 vtsoil,vqsoil,vwetcanp,vsnowdpth,                      &
                 vraing,vrainc,vprcrate,                                &
                 vradfrc,vradsw,vrnflx,vradswnet,vradlwin,              &
                 vusflx,vvsflx,vptsflx,vqvsflx,                         &
                 ireturn, vtem1,vtem2,vtem3)

    ivproj=mapproj
    vscale=sclfct
    vlatnot(1)=trulat1
    vlatnot(2)=trulat2
    vtrulon=trulon
    vdx=vx(3)-vx(2)
    vdy=vy(3)-vy(2)
    vctrlat=ctrlat
    vctrlon=ctrlon
    CALL setmapr(ivproj,vscale,vlatnot,vtrulon)
    CALL lltoxy(1,1,vctrlat,vctrlon,xctr,yctr)
    vx0=xctr-vdx*((vnx-3)/2)
    vy0=yctr-vdy*((vny-3)/2)
    CALL setorig(1,vx0,vy0)
!
!-----------------------------------------------------------------------
!
!  Establish coordinate for scalar verification fields.
!
!-----------------------------------------------------------------------
!
    DO i=1,vnx-1
      vxs(i)=0.5*(vx(i)+vx(i+1))
    END DO
    vxs(vnx)=2.*vxs(vnx-1)-vxs(vnx-2)
    DO j=1,vny-1
      vys(j)=0.5*(vy(j)+vy(j+1))
    END DO
    vys(vny)=2.*vys(vny-1)-vys(vny-2)
    DO k=1,vnz-1
      DO j=1,vny
        DO i=1,vnx
          vzps(i,j,k)=0.5*(vzp(i,j,k)+vzp(i,j,k+1))
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Find location of scalar forecast fields in verification grid.
!
!-----------------------------------------------------------------------
!
    IF(fx0 == vx0 .AND. fy0 == vy0 .AND.                                &
          flatnot(1) == vlatnot(1) .AND. flatnot(2) == vlatnot(2) .AND. &
          ftrulon == vtrulon .AND. ifproj == ivproj .AND.               &
          fscale == vscale ) THEN
      comcoord=.true.
      WRITE(6,'(//a//)') ' Grids share a common coordinate system'
      DO j=1,ny
        DO i=1,nx
          x2d(i,j)=xs(i)
        END DO
      END DO
      DO j=1,ny
        DO i=1,nx
          y2d(i,j)=ys(j)
        END DO
      END DO
    ELSE
      comcoord=.false.
      WRITE(6,'(//a,a//)') ' Grid coordinate systems differ',           &
                         ' Will convert coordinates via lat,lon.'
      CALL lltoxy(nx,ny,lat,lon,x2d,y2d)
    END IF
!
!-----------------------------------------------------------------------
!
!  Interpolate verification scalars to forecast grid.
!
!-----------------------------------------------------------------------
!

!EMK Special for same grids
     IF ((comcoord) .AND. (nx == vnx) .AND. (ny == vny) .AND. &
         (nz == vnz) .AND. (nzsoil == vnzsoil)) THEN
       comcoord2 = .TRUE.
       DO k = 1,nz
         IF (z(k) /= vz(k)) comcoord2 = .FALSE.
       END DO
       IF (comcoord2) THEN
         DO k = 1,nzsoil
           DO j = 1,ny-1
             DO i = 1,nx-1
               IF (zpsoil(i,j,k) /= vzpsoil(i,j,k)) comcoord2 = .FALSE.
             END DO
           END DO
         END DO
       END IF
       IF (comcoord2) THEN
         aptprt = vptprt
         apprt = vpprt
         aqvprt = vqvprt
         DO nq = 1, anscalar
           aqscalar(:,:,:,nq) = vqscalar(:,:,:,vQindex(nq))
         END DO
         atke = vtke
         akmh = vkmh
         akmv = vkmv
         aptbar = vptbar
         apbar = vpbar
         arhobar = vrhobar
         aqvbar = vqvbar
!         atsoil = vtsoil
!         aqsoil = vqsoil
!         awetcanp = vwetcanp
         atsoil(:,:,:,0) = vtsoil(:,:,:,0)
         aqsoil(:,:,:,0) = vqsoil(:,:,:,0)
         awetcanp(:,:,0) = vwetcanp(:,:,0)
         araing = vraing
         arainc = vrainc
         aprcrate = vprcrate
         aradfrc = vradfrc
         aradsw = vradsw
         arnflx = vrnflx
         aradswnet = vradswnet
         aradlwin = vradlwin
         ausflx = vusflx
         avsflx = vvsflx
         aptsflx = vptsflx
         aqvsflx = vqvsflx
         awbar = vwbar
         awprt = vwprt
         avbar = vvbar
         avprt = vvprt
         aubar = vubar
         auprt = vuprt

         GO TO 1000 ! Skip interpolation
      END IF
    END IF
!EMK Special for same grids

    CALL setdxdy(vnx,vny,                                               &
                 1,vnx-1,1,vny-1,                                       &
                 vxs,vys,dxfld,dyfld,rdxfld,rdyfld)

    CALL intsclrs(nx, ny, nz,nzsoil, vnx, vny, vnz,vnzsoil,             &
                  1, nx-1,1, ny-1,1, nz-1, 1, nzsoil-1,                 &
                  1,vnx-1,1,vny-1,1,vnz-1, 1, vnzsoil-1,                &
                  iorder, vnscalar, anscalar, vQindex,                  &
                  x2d, y2d, zps,zpsoil,vxs,vys,vzps,vzpsoil,            &
                  vptprt, vpprt,                                        &
                  vqvprt, vqscalar, vtke,vkmh,vkmv,                     &
                  vptbar, vpbar, vrhobar, vqvbar,                       &
                  vtsoil,vqsoil,vwetcanp,                               &
                  vraing,vrainc,vprcrate,                               &
                  vradfrc,vradsw,vrnflx,vradswnet,vradlwin,             &
                  vusflx,vvsflx,vptsflx,vqvsflx,                        &
                  aptprt, apprt,                                        &
                  aqvprt, aqscalar, atke,akmh,akmv,                     &
                  aptbar, apbar, arhobar, aqvbar,                       &
                  atsoil,aqsoil,awetcanp,                               &
                  araing,arainc,aprcrate,                               &
                  aradfrc,aradsw,arnflx,aradswnet,aradlwin,             &
                  ausflx,avsflx,aptsflx,aqvsflx,                        &
                  iloc,jloc,zpver,dxfld,dyfld,rdxfld,rdyfld,            &
                  vtem1,vtem2,vtem3,                                    &
                  ireturn )
!
!-----------------------------------------------------------------------
!
!  Interpolate verification w to forecast grid.
!
!-----------------------------------------------------------------------
!
    CALL intonef(nx, ny, nz, vnx, vny, vnz,                             &
                 1,nx-1,1,ny-1,1,nz,                                    &
                 1,vnx-1,1,vny-1,1,vnz,                                 &
                 iorder,                                                &
                 x2d, y2d, zp, vxs,vys,vzp,                             &
                 vwprt, vwbar, awprt, awbar,                            &
                 iloc,jloc,zpver,dxfld,dyfld,rdxfld,rdyfld,             &
                 vtem1,vtem2,vtem3,                                     &
                 ireturn )
!
!-----------------------------------------------------------------------
!
!  Find location of u forecast field in verification grid.
!
!-----------------------------------------------------------------------
!
    IF(comcoord) THEN
      DO j=1,ny
        DO i=1,nx
          x2d(i,j)=x(i)
        END DO
      END DO
    ELSE
      CALL setmapr(ifproj,fscale,flatnot,ftrulon)
      CALL setorig(1,fx0,fy0)
      CALL xytoll(nx,ny,x,ys,lat,lon)
      CALL setmapr(ivproj,vscale,vlatnot,vtrulon)
      CALL setorig(1,vx0,vy0)
      CALL lltoxy(nx,ny,lat,lon,x2d,y2d)
    END IF
!
!-----------------------------------------------------------------------
!
!  Interpolate verification u to forecast grid.
!
!-----------------------------------------------------------------------
!
    CALL setdxdy(vnx,vny,                                               &
                 1,vnx,1,vny-1,                                         &
                 vx,vys,dxfld,dyfld,rdxfld,rdyfld)
    CALL intonef(nx, ny, nz, vnx, vny, vnz,                             &
                 1,nx,1,ny-1,1,nz-1,                                    &
                 1,vnx,1,vny-1,1,vnz-1,                                 &
                 iorder,                                                &
                 x2d, y2d, zp, vxs,vys,vzp,                             &
                 vuprt, vubar, auprt, aubar,                            &
                 iloc,jloc,zpver,dxfld,dyfld,rdxfld,rdyfld,             &
                 vtem1,vtem2,vtem3,                                     &
                 ireturn )
!
!-----------------------------------------------------------------------
!
!  Find location of v forecast field in verification grid.
!
!-----------------------------------------------------------------------
!
    IF(comcoord) THEN
      DO j=1,ny
        DO i=1,nx
          x2d(i,j)=xs(i)
        END DO
      END DO
      DO j=1,ny
        DO i=1,nx
          y2d(i,j)=y(j)
        END DO
      END DO
    ELSE
      CALL setmapr(ifproj,fscale,flatnot,ftrulon)
      CALL setorig(1,fx0,fy0)
      CALL xytoll(nx,ny,xs,y,lat,lon)
      CALL setmapr(ivproj,vscale,vlatnot,vtrulon)
      CALL setorig(1,vx0,vy0)
      CALL lltoxy(nx,ny,lat,lon,x2d,y2d)
    END IF
!
!-----------------------------------------------------------------------
!
!  Interpolate verification v to forecast grid.
!
!-----------------------------------------------------------------------
!
    CALL setdxdy(vnx,vny,                                               &
                 1,vnx-1,1,vny,                                         &
                 vxs,vy,dxfld,dyfld,rdxfld,rdyfld)
    CALL intonef(nx, ny, nz, vnx, vny, vnz,                             &
                 1,nx-1,1,ny,1,nz-1,                                    &
                 1,vnx-1,1,vny,1,vnz-1,                                 &
                 iorder,                                                &
                 x2d, y2d, zp, vxs,vys,vzp,                             &
                 vvprt, vvbar, avprt, avbar,                            &
                 iloc,jloc,zpver,dxfld,dyfld,rdxfld,rdyfld,             &
                 vtem1,vtem2,vtem3,                                     &
                 ireturn )
!
!-----------------------------------------------------------------------
!
!  Find   difference = forecast - verification
!
!  To reduce memory requirements, the difference fields are
!  written to the same arrays as the interpolated fields.
!
!-----------------------------------------------------------------------
!
    1000 CONTINUE

    CALL diffield(nx,ny,nz,nzsoil,bnscalar,anscalar,bQindex,            &
                  uprt, vprt, wprt, ptprt, pprt,                        &
                  qvprt, qscalar, tke,kmh,kmv,                          &
                  ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,         &
                  tsoil,qsoil,wetcanp,                                  &
                  raing,rainc,prcrate,                                  &
                  radfrc,radsw,rnflx,radswnet,radlwin,                  &
                  usflx,vsflx,ptsflx,qvsflx,                            &
                  auprt, avprt, awprt, aptprt, apprt,                   &
                  aqvprt, aqscalar, atke,akmh,akmv,                     &
                  aubar, avbar, awbar, aptbar, apbar, arhobar, aqvbar,  &
                  atsoil,aqsoil,awetcanp,                               &
                  araing,arainc,aprcrate,                               &
                  aradfrc,aradsw,arnflx,aradswnet,aradlwin,             &
                  ausflx,avsflx,aptsflx,aqvsflx,                        &
                  auprt, avprt, awprt, aptprt, apprt,                   &
                  aqvprt, aqscalar, atke,akmh,akmv,                     &
                  atsoil,aqsoil,awetcanp,                               &
                  araing,arainc,aprcrate,                               &
                  aradfrc,aradsw,arnflx,aradswnet,aradlwin,             &
                  ausflx,avsflx,aptsflx,aqvsflx,                        &
                  tem1,tem3dsoil,                                       &
                  ireturn )
!
!-----------------------------------------------------------------------
!
!  Set output variables to forecast coordinates
!
!-----------------------------------------------------------------------
!
    curtim=time
    mapproj=ifproj
    sclfct=fscale
    trulat1=flatnot(1)
    trulat2=flatnot(2)
    trulon=ftrulon
    ctrlat=fctrlat
    ctrlon=fctrlon
!
!-----------------------------------------------------------------------
!
!  Get output info
!
!-----------------------------------------------------------------------
!
    WRITE(6,'(/1x,a/)') 'The program is ready for generating OUTPUT file.'
!
!-----------------------------------------------------------------------
!
!  Get runname to use for output data.
!
!-----------------------------------------------------------------------
!
    READ(5,output,END=100)
    WRITE(6,'(1x,2a)') 'The output run name is: ', TRIM(runnmin)

    runname=runnmin
!
!-----------------------------------------------------------------------
!
!  Find out the number of characters to be used to construct file
!  names.
!
!-----------------------------------------------------------------------
!
    CALL gtlfnkey( runname, lfnkey )
!
!-----------------------------------------------------------------------
!
!  Find out the number of characters to be used to construct file
!  names.
!
!-----------------------------------------------------------------------
!
    CALL gtlfnkey( runname, lfnkey )
!
!    WRITE(6,'(a)')' Please input the data format flag value 1/2/3/4 '
!
!-----------------------------------------------------------------------
!
!  Set control parameters for
!  grid, base state, moisture, and ice variable dumping.
!
!-----------------------------------------------------------------------
!
    varout=1
!    WRITE(6,'(a/a/a)')                                                  &
!        ' Will it contain any grid information? (1 or 0)',              &
!        ' If it will not, grid and base information will be dumped',    &
!        ' to a separate file (filename input later as grdbasfn).'
!
!    WRITE(6,'(a)')                                                      &
!        ' Will it contain any base state data? (1 or 0)'
!
!    WRITE(6,'(2(/5x,a)/)')                                              &
!         'Do you want to write u, v, w, ptprt and pprt arrays?',        &
!         '(select 0 or 1)'
!
!    WRITE(6,'(a)')                                                      &
!             ' Write moisture fields to the output file?(1 or 0)'
!
!    WRITE(6,'(a)')' Write ice fields to the output file?(1 or 0)'
!
!    WRITE(6,'(a)')                                                      &
!             ' Write turbulence fields to the output file?(1 or 0)'
!
!    WRITE(6,'(a)')                                                      &
!        ' Write surface (soil) fields to the output file?(1 or 0)'
!
!    WRITE(6,'(a)')' Write rain fields to the output file?(1 or 0)'
!
!    WRITE(6,'(a)')' Write precipitation rate to output file?(1 or 0)'
!
!    WRITE(6,'(a)')' Write radiation arrays to output file?(1 or 0)'
!
!    WRITE(6,'(a)')' Write the surface fluxes to output file?(1 or 0)'
!
!    WRITE(6,'(/5x,a/)')                                                 &
!         'Do you want to compress the output data? (select 0 or 1)'
!

    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          auprt(i,j,k)=aubar(i,j,k)+auprt(i,j,k)
          avprt(i,j,k)=avbar(i,j,k)+avprt(i,j,k)
          awprt(i,j,k)=awbar(i,j,k)+awprt(i,j,k)
          aqvprt(i,j,k)=aqvbar(i,j,k)+aqvprt(i,j,k)
        END DO
      END DO
    END DO

    nscalar   = anscalar
    nscalarq  = anscalarq
    P_QC      = aP_QC
    P_QR      = aP_QR
    P_QI      = aP_QI
    P_QS      = aP_QS
    P_QG      = aP_QG
    P_QH      = aP_QH
    P_NC      = aP_NC
    P_NR      = aP_NR
    P_NI      = aP_NI
    P_NS      = aP_NS
    P_NG      = aP_NG
    P_NH      = aP_NH
    P_ZR      = aP_ZR
    P_ZI      = aP_ZI
    P_ZS      = aP_ZS
    P_ZG      = aP_ZG
    P_ZH      = aP_ZH
    qnames(:) = aqnames(:)
    qdescp(:) = aqdescp(:)

    IF (hdmpfmt == 9) GO TO 700

    CALL gtbasfn(runname(1:lfnkey),'./',2,hdmpfmt,mgrid,nestgrd,        &
                 grdbasfn, lengbf)

    WRITE(6,'(/1x,a,a)')                                                &
        'Output grid/base state file is ', grdbasfn(1:lengbf)

    nchdmp = 80
    grdbas = 1      ! Dump out grd and base state arrays only

    CALL dtadump(nx,ny,nz,nzsoil,nstyps,                                &
                 hdmpfmt,nchdmp,grdbasfn(1:lengbf),grdbas,filcmprs,     &
                 auprt,avprt,awprt,aptprt,apprt,                        &
                 aqvprt,aqscalar,atke,akmh,akmv,                        &
                 ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,                &
                 x,y,z,zp,zpsoil,                                       &
                 asoiltyp,astypfrct,avegtyp,alai,aroufns,aveg,          &
                 atsoil,aqsoil,awetcanp,asnowdpth,                      &
                 araing,arainc,aprcrate,                                &
                 aradfrc,aradsw,arnflx,aradswnet,aradlwin,              &
                 ausflx,avsflx,aptsflx,aqvsflx,                         &
                 tem1,tem2,tem3)

!
!-----------------------------------------------------------------------
!
!  Find a unique name hdmpfn(1:ldmpf) for history dump data set
!  at time 'curtim'.
!
!-----------------------------------------------------------------------
!
    grdbas = 0      ! Not just dump out time dependent arrays

    700 CONTINUE

    CALL gtdmpfn(runname(1:lfnkey),'./',2,                              &
                 curtim,hdmpfmt,                                        &
                 mgrid,nestgrd, hdmpfn, ldmpf)

    WRITE(6,'(/1x,a,f10.0,a,a)')                                        &
        'Output file at time ',curtim,' (s) is ', hdmpfn(1:ldmpf)

    CALL dtadump(nx,ny,nz,nzsoil,nstyps,                                &
                 hdmpfmt,nchdmp,hdmpfn(1:ldmpf),grdbas,filcmprs,        &
                 auprt,avprt,awprt,aptprt,apprt,                        &
                 aqvprt,aqscalar,atke,akmh,akmv,                        &
                 ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,                &
                 x,y,z,zp,zpsoil,                                       &
                 asoiltyp,astypfrct,avegtyp,alai,aroufns,aveg,          &
                 atsoil,aqsoil,awetcanp,asnowdpth,                      &
                 araing,arainc,aprcrate,                                &
                 aradfrc,aradsw,arnflx,aradswnet,aradlwin,              &
                 ausflx,avsflx,aptsflx,aqvsflx,                         &
                 tem1,tem2,tem3)

  END IF

  GOTO 101

  100 CONTINUE

  WRITE(6,'(a)')'Namelist block READ in error. Then program will terminated.'

  101 CONTINUE

  STOP

END PROGRAM arpsdiff
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE DIFFIELD                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE diffield(nx,ny,nz,nzsoil,bnscalar,anscalar,Qindex,           &
           uprt, vprt, wprt, ptprt, pprt,                               &
           qvprt, qscalar, tke,kmh,kmv,                                 &
           ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,                &
           tsoil,qsoil,wetcanp,                                         &
           raing,rainc,prcrate,                                         &
           radfrc,radsw,rnflx,radswnet,radlwin,                         &
           usflx,vsflx,ptsflx,qvsflx,                                   &
           auprt, avprt, awprt, aptprt, apprt,                          &
           aqvprt, aqscalar, atke,akmh,akmv,                            &
           aubar, avbar, awbar, aptbar, apbar, arhobar, aqvbar,         &
           atsoil,aqsoil,awetcanp,                                      &
           araing,arainc,aprcrate,                                      &
           aradfrc,aradsw,arnflx,aradswnet,aradlwin,                    &
           ausflx,avsflx,aptsflx,aqvsflx,                               &
           duprt, dvprt, dwprt, dptprt, dpprt,                          &
           dqvprt, dqscalar, dtke,dkmh,dkmv,                            &
           dtsoil,dqsoil,dwetcanp,                                      &
           draing,drainc,dprcrate,                                      &
           dradfrc,dradsw,drnflx,dradswnet,dradlwin,                    &
           dusflx,dvsflx,dptsflx,dqvsflx,                               &
           tem1,tem3dsoil,                                              &
           ireturn)
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Subtract the forecast fields from the interpolated verification
!  fields (names beginning with "a") and output to the difference
!  fields (names beginning with "d").  The input and difference
!  fields may share the same storage location.  For this subroutine
!  it is assumed the forecast and corresponding verification
!  data are at the same physical location, however, the physical
!  location may differ between variables. That is uprt and auprt
!  are at the same location, but that may differ from pprt and apprt.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster  Ou School of Meteorology. April 1992
!
!  MODIFICATION HISTORY:
!   14 May 1992  (KB) changed from arps2.5 to arps3.0
!   03 Aug 1992  (KB) updated to account for changes in arps3.0
!
!   09/07/1995  (KB)
!   Added differencing of surface (soil) fields.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!    nx,ny,nz Array dimensions for forecast field.
!
!    FORECAST FIELDS:
!
!    uprt     perturbation x component of velocity (m/s)
!    vprt     perturbation y component of velocity (m/s)
!    wprt     perturbation vertical component of velocity in Cartesian
!             coordinates (m/s).
!
!    ptprt    perturbation potential temperature (K)
!    pprt     perturbation pressure (Pascal)
!
!    qvprt    perturbation water vapor mixing ratio (kg/kg)
!    qscalar  Hydrometeor scalars
!
!    ubar     Base state x velocity component (m/s)
!    vbar     Base state y velocity component (m/s)
!    wbar     Base state z velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhobar   Base state density (kg/m**3)
!    qvbar    Base state water vapor mixing ratio (kg/kg)
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!    kmh      Horizontal turbulent mixing coefficient (m**2/s)
!    kmv      Vertical turbulent mixing coefficient (m**2/s)
!
!    tsoil    Soil temperature (K)
!    qsoil    Soil temperature (m**3/m**3)
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
!    INTERPOLATED VERIFICATION FIELDS:
!
!    auprt     perturbation x component of velocity (m/s)
!    avprt     perturbation y component of velocity (m/s)
!    awprt     perturbation vertical component of velocity in Cartesian
!              coordinates (m/s).
!
!    aptprt    perturbation potential temperature (K)
!    apprt     perturbation pressure (Pascal)
!
!    aqvprt    perturbation water vapor mixing ratio (kg/kg)
!    aqscalar
!
!    aubar     Base state x velocity component (m/s)
!    avbar     Base state y velocity component (m/s)
!    awbar     Base state z velocity component (m/s)
!    aptbar    Base state potential temperature (K)
!    apbar     Base state pressure (Pascal)
!    arhobar   Base state density (kg/m**3)
!    aqvbar    Base state water vapor mixing ratio (kg/kg)
!
!    atsoil    Soil temperature (K)
!    aqsoil    Soil temperature (m**3/m**3)
!    awetcanp  Canopy water amount
!
!    araing    Grid supersaturation rain
!    arainc    Cumulus convective rain
!    aprcrate  Precipitation rates
!
!    aradfrc   Radiation forcing (K/s)
!    aradsw    Solar radiation reaching the surface
!    arnflx    Net radiation flux absorbed by surface
!    aradswnet Net shortwave radiation, SWin - SWout
!    aradlwin  Incoming longwave radiation
!
!    ausflx    Surface flux of u-momentum (kg/(m*s**2))
!    avsflx    Surface flux of v-momentum (kg/(m*s**2))
!    aptsflx   Surface heat flux (K*kg/(m**2 * s ))
!    aqvsflx   Surface moisture flux of (kg/(m**2 * s))
!
!  OUTPUT :
!
!    DIFFERENCE FIELDS (may share storage with forecast fields
!                       or interpolated fields in calling program):
!
!    duprt     perturbation x component of velocity (m/s)
!    dvprt     perturbation y component of velocity (m/s)
!    dwprt     perturbation vertical component of velocity in Cartesian
!              coordinates (m/s).
!
!    dptprt    perturbation potential temperature (K)
!    dpprt     perturbation pressure (Pascal)
!
!    dqvprt    perturbation water vapor mixing ratio (kg/kg)
!    dqscalar
!
!    dtsoil    Soil temperature (K)
!    dqsoil    Soil moisture (m**3/m**3)
!    dwetcanp  Canopy water amount
!
!    draing    Grid supersaturation rain
!    drainc    Cumulus convective rain
!    dprcrate  Precipitation rates
!
!    dradfrc   Radiation forcing (K/s)
!    dradsw    Solar radiation reaching the surface
!    drnflx    Net radiation flux absorbed by surface
!    dradswnet Net shortwave radiation, SWin - SWout
!    dradlwin  Incoming longwave radiation
!
!    dusflx    Surface flux of u-momentum (kg/(m*s**2))
!    dvsflx    Surface flux of v-momentum (kg/(m*s**2))
!    dptsflx   Surface heat flux (K*kg/(m**2 * s ))
!    dqvsflx   Surface moisture flux of (kg/(m**2 * s))
!
!    tem1      Work array
!    tem3dsoil Work array
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz         ! 3 dimensions of array
  INTEGER :: nzsoil           ! soil array

  INTEGER :: bnscalar, anscalar, Qindex(anscalar)
!
!-----------------------------------------------------------------------
!
!  Model Arrays
!
!-----------------------------------------------------------------------
!
  REAL :: uprt   (nx,ny,nz)    ! Perturbation u-velocity (m/s)
  REAL :: vprt   (nx,ny,nz)    ! Perturbation v-velocity (m/s)
  REAL :: wprt   (nx,ny,nz)    ! Perturbation w-velocity (m/s)
  REAL :: ptprt  (nx,ny,nz)    ! Perturbation potential temperature (K)
  REAL :: pprt   (nx,ny,nz)    ! Perturbation pressure (Pascal)
  REAL :: qvprt  (nx,ny,nz)    ! Perturbation water vapor specific humidity
  REAL :: qscalar(nx,ny,nz,bnscalar)

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
  REAL :: qvbar  (nx,ny,nz)    ! Base state water vapor specific humidity

  REAL :: tsoil  (nx,ny,nzsoil) ! Soil temperature (K)
  REAL :: qsoil (nx,ny,nzsoil)  ! Soil moisture (m**3/m**3)
  REAL :: wetcanp(nx,ny)      ! Canopy water amount

  REAL :: raing (nx,ny)       ! Grid supersaturation rain
  REAL :: rainc (nx,ny)       ! Cumulus convective rain
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
!
!-----------------------------------------------------------------------
!
!  Verification data interpolated to model grid
!
!-----------------------------------------------------------------------
!
  REAL :: auprt  (nx,ny,nz)    ! Perturbation u-velocity (m/s)
  REAL :: avprt  (nx,ny,nz)    ! Perturbation v-velocity (m/s)
  REAL :: awprt  (nx,ny,nz)    ! Perturbation w-velocity (m/s)
  REAL :: aptprt (nx,ny,nz)    ! Perturbation potential temperature (K)
  REAL :: apprt  (nx,ny,nz)    ! Perturbation pressure (Pascal)
  REAL :: aqvprt (nx,ny,nz)    ! Perturbation water vapor specific humidity
  REAL :: aqscalar(nx,ny,nz,anscalar)

  REAL :: atke   (nx,ny,nz)    ! Turbulent Kinetic Energy ((m/s)**2)
  REAL :: akmh   (nx,ny,nz)    ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: akmv   (nx,ny,nz)    ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )

  REAL :: aubar  (nx,ny,nz)    ! Base state u-velocity (m/s)
  REAL :: avbar  (nx,ny,nz)    ! Base state v-velocity (m/s)
  REAL :: awbar  (nx,ny,nz)    ! Base state w-velocity (m/s)
  REAL :: aptbar (nx,ny,nz)    ! Base state potential temperature (K)
  REAL :: arhobar(nx,ny,nz)    ! Base state density (kg/m**3)
  REAL :: apbar  (nx,ny,nz)    ! Base state pressure (Pascal)
  REAL :: aqvbar (nx,ny,nz)    ! Base state water vapor specific humidity

  REAL :: atsoil (nx,ny,nzsoil) ! Soil temperature (K)
  REAL :: aqsoil (nx,ny,nzsoil) ! Soil moisture (m**3/m**3))
  REAL :: awetcanp(nx,ny)      ! Canopy water amount

  REAL :: araing (nx,ny)       ! Grid supersaturation rain
  REAL :: arainc (nx,ny)       ! Cumulus convective rain
  REAL :: aprcrate(nx,ny,4)    ! precipitation rate (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precip. rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulus precip. rate
                               ! prcrate(1,1,4) = microphysics precip. rate

  REAL :: aradfrc(nx,ny,nz)     ! Radiation forcing (K/s)
  REAL :: aradsw (nx,ny)        ! Solar radiation reaching the surface
  REAL :: arnflx (nx,ny)        ! Net radiation flux absorbed by surface
  REAL :: aradswnet(nx,ny)      ! Net shortwave radiation
  REAL :: aradlwin(nx,ny)       ! Incoming longwave radiation

  REAL :: ausflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL :: avsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL :: aptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL :: aqvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))
!
!-----------------------------------------------------------------------
!
!  Difference arrays
!
!-----------------------------------------------------------------------
!
  REAL :: duprt  (nx,ny,nz)    ! perturbation x component of velocity (m/s)
  REAL :: dvprt  (nx,ny,nz)    ! perturbation y component of velocity (m/s)
  REAL :: dwprt  (nx,ny,nz)    ! perturbation vertical component of
                               ! velocity in Cartesian coordinates (m/s)
  REAL :: dptprt (nx,ny,nz)    ! perturbation potential temperature (K)
  REAL :: dpprt  (nx,ny,nz)    ! perturbation pressure (Pascal)
  REAL :: dqvprt (nx,ny,nz)    ! perturbation water vapor mixing ratio (kg/kg)
  REAL :: dqscalar(nx,ny,nz,anscalar)
  REAL :: dtke   (nx,ny,nz)    ! Turbulent Kinetic Energy ((m/s)**2)
  REAL :: dkmh   (nx,ny,nz)    ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: dkmv   (nx,ny,nz)    ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )

  REAL :: dtsoil (nx,ny)       ! Soil temperature at surface (K)
  REAL :: dqsoil (nx,ny)       ! Soil moisture (m**3/m**3)
  REAL :: dwetcanp(nx,ny)      ! Canopy water amount

  REAL :: draing (nx,ny)       ! Grid supersaturation rain
  REAL :: drainc (nx,ny)       ! Cumulus convective rain
  REAL :: dprcrate(nx,ny,4)    ! precipitation rate (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precip. rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulus precip. rate
                               ! prcrate(1,1,4) = microphysics precip. rate

  REAL :: dradfrc(nx,ny,nz)    ! Radiation forcing (K/s)
  REAL :: dradsw (nx,ny)       ! Solar radiation reaching the surface
  REAL :: drnflx (nx,ny)       ! Net radiation flux absorbed by surface
  REAL :: dradswnet(nx,ny)      ! Net shortwave radiation
  REAL :: dradlwin(nx,ny)       ! Incoming longwave radiation

  REAL :: dusflx (nx,ny)       ! Surface flux of u-momentum (kg/(m*s**2))
  REAL :: dvsflx (nx,ny)       ! Surface flux of v-momentum (kg/(m*s**2))
  REAL :: dptsflx(nx,ny)       ! Surface heat flux (K*kg/(m*s**2))
  REAL :: dqvsflx(nx,ny)       ! Surface moisture flux (kg/(m**2*s))

  REAL :: tem1   (nx,ny,nz)    ! A work array
  REAL :: tem3dsoil(nx,ny,nzsoil) ! A work array
  INTEGER :: ireturn, i,j,k, nq
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: is,js,ks,ls,ie,je,ke,le
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  is=1
  js=1
  ks=1
  ls=1
  ie=nx-1
  je=ny-1
  ke=nz-1
  le=nzsoil-1

!-----------------------------------------------------------------------
!
!  Scalars
!
!-----------------------------------------------------------------------

  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        tem1(i,j,k)=0.0
      END DO
    END DO
  END DO

  tem3dsoil = 0.0

  PRINT *, ' ptprt: '
  CALL subtr(nx,ny,nz, ptprt,ptbar,aptprt,aptbar,dptprt,                &
             is,js,ks,ie,je,ke)
  PRINT *, ' pprt: '
  CALL subtr(nx,ny,nz,  pprt, pbar, apprt, apbar, dpprt,                &
             is,js,ks,ie,je,ke)
  PRINT *, ' qvprt: '
  CALL subtr(nx,ny,nz, qvprt,qvbar,aqvprt,aqvbar,dqvprt,                &
             is,js,ks,ie,je,ke)

  DO nq = 1,anscalar
    PRINT *, ' nq = ',nq
    CALL subtr(nx,ny,nz,    qscalar(:,:,:,Qindex(nq)), tem1,           &
               aqscalar(:,:,:,nq),  tem1,   dqscalar(:,:,:,nq),         &
               is,js,ks,ie,je,ke)
  END DO
  PRINT *, ' tke: '
  CALL subtr(nx,ny,nz,    tke, tem1,   atke,  tem1,   dtke,             &
             is,js,ks,ie,je,ke)
  PRINT *, ' kmh: '
  CALL subtr(nx,ny,nz,    kmh, tem1,   akmh,  tem1,   dkmh,             &
             is,js,ks,ie,je,ke)
  PRINT *, ' kmv: '
  CALL subtr(nx,ny,nz,    kmv, tem1,   akmv,  tem1,   dkmv,             &
             is,js,ks,ie,je,ke)

!-----------------------------------------------------------------------
!
!  u wind components
!
!-----------------------------------------------------------------------

  ie=nx
  PRINT *, ' uprt: '
  CALL subtr(nx,ny,nz,uprt,ubar,auprt,aubar,duprt,                      &
             is,js,ks,ie,je,ke)

!-----------------------------------------------------------------------
!
!  v wind components
!
!-----------------------------------------------------------------------

  ie=nx-1
  je=ny
  PRINT *, ' vprt: '
  CALL subtr(nx,ny,nz,vprt,vbar,avprt,avbar,dvprt,                      &
             is,js,ks,ie,je,ke)

!-----------------------------------------------------------------------
!
!  w wind components
!
!-----------------------------------------------------------------------

  je=ny-1
  ke=nz
  CALL subtr(nx,ny,nz,wprt,tem1,awprt,tem1,dwprt,                       &
             is,js,ks,ie,je,ke)

!-----------------------------------------------------------------------
!
!  2-d surface (soil) variables
!
!-----------------------------------------------------------------------

  ie=nx-1
  je=ny-1
  le=nzsoil
  ks=1
  ke=1

!
  PRINT *, ' tsoil:'
  CALL subtr(nx,ny,nzsoil, tsoil,tem3dsoil,atsoil,tem3dsoil, dtsoil,    &
             is,js,ls,ie,je,le)
  PRINT *, ' qsoil:'
  CALL subtr(nx,ny,nzsoil, qsoil,tem3dsoil, aqsoil,tem3dsoil, dqsoil,   &
             is,js,ls,ie,je,le)
  PRINT *, ' wetcanp:'
  CALL subtr(nx,ny,1,wetcanp,tem1,awetcanp,tem1,dwetcanp,               &
             is,js,ks,ie,je,ke)
  PRINT *, ' raing:'
  CALL subtr(nx,ny,1,  raing,tem1,  araing,tem1,  draing,               &
             is,js,ks,ie,je,ke)
  PRINT *, ' rainc:'
  CALL subtr(nx,ny,1,  rainc,tem1,  arainc,tem1,  drainc,               &
             is,js,ks,ie,je,ke)
  PRINT *, ' prcrate1:'
  CALL subtr(nx,ny,1,  prcrate(1,1,1),tem1, aprcrate(1,1,1),tem1,       &
             dprcrate(1,1,1), is,js,ks,ie,je,ke)
  PRINT *, ' prcrate2:'
  CALL subtr(nx,ny,1,  prcrate(1,1,2),tem1, aprcrate(1,1,2),tem1,       &
             dprcrate(1,1,2), is,js,ks,ie,je,ke)
  PRINT *, ' prcrate3:'
  CALL subtr(nx,ny,1,  prcrate(1,1,3),tem1, aprcrate(1,1,3),tem1,       &
             dprcrate(1,1,3), is,js,ks,ie,je,ke)
  PRINT *, ' prcrate4:'
  CALL subtr(nx,ny,1,  prcrate(1,1,4),tem1, aprcrate(1,1,4),tem1,       &
             dprcrate(1,1,4), is,js,ks,ie,je,ke)

  PRINT *, ' radfrc:'
  CALL subtr(nx,ny,nz, radfrc,tem1, aradfrc,tem1,                       &
             dradfrc, is,js,ks,ie,je,ke)
  PRINT *, ' radsw:'
  CALL subtr(nx,ny,1, radsw,tem1, aradsw,tem1,                          &
             dradsw, is,js,ks,ie,je,ke)
  PRINT *, ' rnflx:'
  CALL subtr(nx,ny,1, rnflx,tem1, arnflx,tem1,                          &
             drnflx, is,js,ks,ie,je,ke)
  PRINT *, ' radswnet:'
  CALL subtr(nx,ny,1, radswnet,tem1,aradswnet,tem1,                     &
             dradswnet, is,js,ks,ie,je,ke)
  PRINT *, ' radlwin:'
  CALL subtr(nx,ny,1, radlwin,tem1, aradlwin,tem1,                      &
             dradlwin, is,js,ks,ie,je,ke)

  RETURN
END SUBROUTINE diffield
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE SUBTR                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE subtr(nx,ny,nz, a,abar,b,bbar,c,                             &
           istr,jstr,kstr,iend,jend,kend)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Subtracts 2 three-dimensional arrays, represented by
!  means plus perturbations.
!
!  AUTHOR: Keith Brewster  OU School of Meteorology.  Feb 1992
!
!  MODIFICATION HISTORY:
!   11 Aug 1992  (KB) changed from arps2.5 to arps3.0
!
!   7 May 2002 (Eric Kemp)
!   Minor change to list all non-zero differences.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    a        perturbation data array
!    abar     mean data array
!    b        perturbation data array to subtract from a
!    bbar     mean data array to subtract from a
!
!  OUTPUT:
!    c        difference array a-b
!             (may share storage in calling program with array a or b)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! 3 dimensions of array

  REAL :: a   (nx,ny,nz)       ! data array
  REAL :: abar(nx,ny,nz)       ! base state of data array a
  REAL :: b   (nx,ny,nz)       ! data array to subtract from a
  REAL :: bbar(nx,ny,nz)       ! base state of data arrya b
  REAL :: c   (nx,ny,nz)       ! difference array a-b
  INTEGER :: istr,jstr,kstr
  INTEGER :: iend,jend,kend
  INTEGER :: i,j,k,imid,jmid,kmid
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  imid=nint(0.5*FLOAT(istr+iend))
  jmid=nint(0.5*FLOAT(jstr+jend))
  kmid=nint(0.5*FLOAT(kstr+kend))
!
!-----------------------------------------------------------------------
!
!  Tell us about a sample input point
!
!-----------------------------------------------------------------------
!
  PRINT *, ' sample, a= ',(a(imid,jmid,kmid)+abar(imid,jmid,kmid)),     &
                   ' b= ',(b(imid,jmid,kmid)+bbar(imid,jmid,kmid))
!
!-----------------------------------------------------------------------
!
!  Subtraction
!
!-----------------------------------------------------------------------
!
  DO k=kstr,kend
    DO j=jstr,jend
      DO i=istr,iend
        c(i,j,k)=a(i,j,k)+abar(i,j,k)-(b(i,j,k)+bbar(i,j,k))

        ! Y.J. removed base field
        abar(i,j,k) = 0.0 
        bbar(i,j,k) = 0.0

!          IF (c(i,j,k) /= 0) THEN
!            WRITE(6,*)'Non-zero difference in field! c = ',c(i,j,k), &
!                      ' at i,j,k: ',i,j,k
!          END IF

      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Tell us about a sample output point
!
!-----------------------------------------------------------------------
!
  PRINT *, '         c= ',c(imid,jmid,kmid)
  RETURN
END SUBROUTINE subtr
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE INTSCLRS                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE intsclrs(nx,ny,nz,nzsoil,vnx,vny,vnz,vnzsoil,                &
           ibeg,iend,jbeg,jend,kbeg,kend,nbeg, nend,                    &
           ivbeg,ivend,jvbeg,jvend,kvbeg,kvend,nvbeg, nvend,            &
           iorder,vnscalar, anscalar,Qindex,                            &
           xs2d,ys2d,zps,zpsoil, vxs,vys,vzps,vzpsoil,                  &
           vptprt, vpprt,                                               &
           vqvprt, vqscalar, vtke,vkmh,vkmv,                            &
           vptbar, vpbar, vrhobar, vqvbar,                              &
           vtsoil,vqsoil,vwetcanp,                                      &
           vraing,vrainc,vprcrate,                                      &
           vradfrc,vradsw,vrnflx,vradswnet,vradlwin,                    &
           vusflx,vvsflx,vptsflx,vqvsflx,                               &
           aptprt, apprt,                                               &
           aqvprt, aqscalar, atke,akmh,akmv,                            &
           aptbar, apbar, arhobar, aqvbar,                              &
           atsoil,aqsoil,awetcanp,                                      &
           araing,arainc,aprcrate,                                      &
           aradfrc,aradsw,arnflx,aradswnet,aradlwin,                    &
           ausflx,avsflx,aptsflx,aqvsflx,                               &
           iloc,jloc,zpver,dxfld,dyfld,rdxfld,rdyfld,                   &
           slopey,alphay,betay,                                         &
           ireturn )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Intfield interpolates scalars from a set of fields (the verification
!    fields, "verif") having Cartesian coordinates described by vx,vy,vzp
!    to a second set of fields described by cartesion coordinates x,y,zp.
!    It is assumed that x,y,zp and vx,vy,vzp are monotonically increasing
!    with increasing index.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster  OU School of Meteorology.  Feb 1992
!
!  MODIFICATION HISTORY:
!    12 Aug 1992  (KB) changed from arps2.5 to arps3.0
!    19 May 1993  (KB) changed from arps3.1 to arps3.2
!    24 May 1993  (KB) changed to special version for scalars only.
!
!     9 Sep 1995  (KB) added processing of sfc (soil) fields
!    26 Apr 1996  (KB) Version 2.0 -- Uses Gauss Forward routines for
!                      interpolation rather than piecewise linear.
!    07 Nov 1996  (KB) Replaced interpolation scheme.
!                      Reordered sequence of variables in call.
!    03 Nov 1999  (KB via Eric Kemp) Corrected dimensions of precip
!                      and flux v* (verif) variables that had been
!                      recently added.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    nx       Number of grid pts in the x-direction (east/west)
!    ny       Number of grid pts in the y-direction (north/south)
!    nz       Number of grid pts in the vertical
!    nzsoil   Number of grid pts in the soil
!
!    vnx      Number of verif grid points in the x-direction (east/west)
!    vny      Number of verif grid points in the y-direction (north/south)
!    vnz      Number of verif grid points in the vertical
!    vnzsoil  Number of verif grid points in the soil.
!
!    ibeg,iend   Range of x index to do interpolation
!    jbeg,jend   Range of y index to do interpolation
!    kbeg,kend   Range of z index to do interpolation
!
!    ivbeg,ivend   Range of x index to use in verification array
!    jvbeg,jvend   Range of y index to use in verification array
!    kvbeg,kvend   Range of z index to use in verification array
!
!    iorder   Interpolation parameter.
!             iorder specifies the order of interpolation
!             1 = bi-linear
!             2 = bi-quadratic
!
!    xs2d     x coordinate of scalar grid points in physical space (m)
!    ys2d     y coordinate of scalar grid points in physical space (m)
!    zps      z coordinate of scalar grid points in physical space (m)
!    zpsoil   z coordinate of soil levels (m)
!
!    vxs      x coordinate of verif scalar grid points in physical space (m)
!    vys      y coordinate of verif scalar grid points in physical space (m)
!    vzps     z coordinate of verif scalar grid points in physical space (m)
!    vzpsoil  z coordinate of verif soil levels (m)
!
!    vpt      Potential temperature
!    vpprt    Perturbation pressure  (Pascal)
!    vqv      Water vapor specific humidity  (kg/kg)
!    vqscalar
!
!    vptbar    Base state potential temperature (K)
!    vpbar     Base state pressure (Pascal)
!    vrhobar   Base state density (kg/m**3)
!    vqvbar    Base state water vapor mixing ratio (kg/kg)
!
!    vtsoil    Soil temperature (K)
!    vqsoil    Soil moisture (m**3/m**3)
!    vwetcanp  Canopy water amount
!
!    vraing    Grid supersaturation rain
!    vrainc    Cumulus convective rain
!    vprcrate  Precipitation rates
!
!    vradfrc   Radiation forcing (K/s)
!    vradsw    Solar radiation reaching the surface
!    vrnflx    Net radiation flux absorbed by surface
!    vradswnet Net shortwave radiation, SWin - SWout
!    vradlwin  Incoming longwave radiation
!
!    vusflx    Surface flux of u-momentum (kg/(m*s**2))
!    vvsflx    Surface flux of v-momentum (kg/(m*s**2))
!    vptsflx   Surface heat flux (K*kg/(m**2 * s ))
!    vqvsflx   Surface moisture flux of (kg/(m**2 * s))
!
!  OUTPUT:
!    apt      Interpolated potential temperature
!    apprt    Interpolated perturbation pressure  (Pascal)
!    aqv      Interpolated water vapor specific humidity  (kg/kg)
!    aqscalar
!
!    aptbar   Interpolated base state potential temperature (K)
!    apbar    Interpolated base state pressure (Pascal)
!    arhobar  Interpolated base state density (kg/m**3)
!    aqvbar   Interpolated base state water vapor mixing ratio (kg/kg)
!
!    atsoil    Interpolated temperature (K)
!    aqsoil    Interpolated soil moisture (m**3/m**3)
!    awetcanp  Interpolated canopy water amount
!
!    araing    Interpolated grid supersaturation rain
!    arainc    Interpolated cumulus convective rain
!    aprcrate  Precipitation rates
!
!    aradfrc   Radiation forcing (K/s)
!    aradsw    Solar radiation reaching the surface
!    arnflx    Net radiation flux absorbed by surface
!    aradswnet Net shortwave radiation, SWin - SWout
!    aradlwin  Incoming longwave radiation
!
!    ausflx    Surface flux of u-momentum (kg/(m*s**2))
!    avsflx    Surface flux of v-momentum (kg/(m*s**2))
!    aptsflx   Surface heat flux (K*kg/(m**2 * s ))
!    aqvsflx   Surface moisture flux of (kg/(m**2 * s))
!
!    ireturn
!
!  WORK ARRAYS:
!
!    iloc     i-index of interpolation points in field to be interpolated
!    jloc     j-index of interpolation points in field to be interpolated
!    dxfld    Vector of delta-x (m) of field to be interpolated
!    dyfld    Vector of delta-y (m) of field to be interpolated
!    rdxfld   Vector of 1./delta-x (1/m) of field to be interpolated
!    rdyfld   Vector of 1./delta-y (1/m) of field to be interpolated
!
!    slopey   Piecewise linear df/dy
!    alphay   Coefficient of y-squared term in y quadratic interpolator
!    betay    Coefficient of y term in y quadratic interpolator
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz,nzsoil
  INTEGER :: vnx,vny,vnz,vnzsoil
  INTEGER :: ibeg,iend,jbeg,jend,kbeg,kend,nbeg,nend
  INTEGER :: ivbeg,ivend,jvbeg,jvend,kvbeg,kvend,nvbeg,nvend

  REAL :: xs2d(nx,ny)
  REAL :: ys2d(nx,ny)
  REAL :: zps(nx,ny,nz)
  REAL :: zpsoil(nx,ny,nzsoil)

  REAL :: vxs(vnx)
  REAL :: vys(vny)
  REAL :: vzps(vnx,vny,vnz)
  REAL :: vzpsoil(vnx,vny,vnzsoil)


  INTEGER :: iorder
  INTEGER :: vnscalar, anscalar, Qindex(anscalar)
!
!-----------------------------------------------------------------------
!
!  Original arrays (verification field)
!
!-----------------------------------------------------------------------
!
  REAL :: vptprt(vnx,vny,vnz)
  REAL :: vpprt (vnx,vny,vnz)
  REAL :: vqvprt(vnx,vny,vnz)
  REAL :: vqscalar   (vnx,vny,vnz,vnscalar)
  REAL :: vtke  (vnx,vny,vnz)
  REAL :: vkmh  (vnx,vny,vnz)
  REAL :: vkmv  (vnx,vny,vnz)

  REAL :: vptbar (vnx,vny,vnz)
  REAL :: vrhobar(vnx,vny,vnz)
  REAL :: vpbar  (vnx,vny,vnz)
  REAL :: vqvbar (vnx,vny,vnz)

  REAL :: vtsoil (vnx,vny,vnzsoil) ! Soil Temperature (K)
  REAL :: vqsoil (vnx,vny,vnzsoil) ! Soil moisture
  REAL :: vwetcanp(vnx,vny)    ! Canopy water amount

  REAL :: vraing (vnx,vny)     ! Grid supersaturation rain
  REAL :: vrainc (vnx,vny)     ! Cumulus convective rain
  REAL :: vprcrate(vnx,vny,4)    ! precipitation rates (kg/(m**2*s))
                                 ! prcrate(1,1,1) = total precip. rate
                                 ! prcrate(1,1,2) = grid scale precip. rate
                                 ! prcrate(1,1,3) = cumulus precip. rate
                                 ! prcrate(1,1,4) = microphysics precip. rate

  REAL :: vradfrc(vnx,vny,vnz)   ! Radiation forcing (K/s)
  REAL :: vradsw (vnx,vny)       ! Solar radiation reaching the surface
  REAL :: vrnflx (vnx,vny)       ! Net radiation flux absorbed by surface
  REAL :: vradswnet(vnx,vny)     ! Net shortwave radiation
  REAL :: vradlwin(vnx,vny)      ! Incoming longwave radiation

  REAL :: vusflx (vnx,vny)       ! Surface flux of u-momentum (kg/(m*s**2))
  REAL :: vvsflx (vnx,vny)       ! Surface flux of v-momentum (kg/(m*s**2))
  REAL :: vptsflx(vnx,vny)       ! Surface heat flux (K*kg/(m*s**2))
  REAL :: vqvsflx(vnx,vny)       ! Surface moisture flux (kg/(m**2*s))
!
!-----------------------------------------------------------------------
!
!  Arrays interpolated to model grid
!
!-----------------------------------------------------------------------
!
  REAL :: aptprt(nx,ny,nz)
  REAL :: apprt (nx,ny,nz)
  REAL :: aqvprt(nx,ny,nz)
  REAL :: aqscalar   (nx,ny,nz,anscalar)
  REAL :: atke  (nx,ny,nz)
  REAL :: akmh  (nx,ny,nz)
  REAL :: akmv  (nx,ny,nz)

  REAL :: aptbar (nx,ny,nz)
  REAL :: arhobar(nx,ny,nz)
  REAL :: apbar  (nx,ny,nz)
  REAL :: aqvbar (nx,ny,nz)

  REAL :: atsoil (nx,ny,nzsoil) ! Soil temperature (K)
  REAL :: aqsoil (nx,ny,nzsoil) ! Soil moisture (m**3/m**3)
  REAL :: awetcanp(nx,ny)      ! Canopy water amount

  REAL :: araing (nx,ny)       ! Grid supersaturation rain
  REAL :: arainc (nx,ny)       ! Cumulus convective rain
  REAL :: aprcrate(nx,ny,4)    ! precipitation rate (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precip. rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulus precip. rate
                               ! prcrate(1,1,4) = microphysics precip. rate

  REAL :: aradfrc(nx,ny,nz)     ! Radiation forcing (K/s)
  REAL :: aradsw (nx,ny)        ! Solar radiation reaching the surface
  REAL :: arnflx (nx,ny)        ! Net radiation flux absorbed by surface
  REAL :: aradswnet(nx,ny)      ! Net shortwave radiation
  REAL :: aradlwin(nx,ny)       ! Incoming longwave radiation

  REAL :: ausflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL :: avsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL :: aptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL :: aqvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))
!
!-----------------------------------------------------------------------
!
!  Work arrays
!
!-----------------------------------------------------------------------
!
  INTEGER :: iloc(nx,ny)
  INTEGER :: jloc(nx,ny)
  REAL :: zpver(nx,ny,vnz)
!
  REAL :: dxfld(vnx)
  REAL :: dyfld(vny)
  REAL :: rdxfld(vnx)
  REAL :: rdyfld(vny)
  REAL :: slopey(vnx,vny,vnz)
  REAL :: alphay(vnx,vny,vnz)
  REAL :: betay(vnx,vny,vnz)
!
  INTEGER :: ireturn
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: k,korder, nq
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Find i,j indices in verfication grid of each forecast point
!
!-----------------------------------------------------------------------
!
  CALL setijloc(nx,ny,vnx,vny,xs2d,ys2d,vxs,vys,iloc,jloc)
  CALL setdxdy(vnx,vny,                                                 &
               1,vnx-1,1,vny-1,                                         &
               vxs,vys,dxfld,dyfld,rdxfld,rdyfld)
!
!-----------------------------------------------------------------------
!
!  Interpolate 3-d, 2-d fields
!
!-----------------------------------------------------------------------
!
  korder=MIN(iorder,3)

  CALL fldint3d(nx,ny,nzsoil,vnx,vny,vnzsoil,                           &
                ibeg,iend,jbeg,jend,nbeg,nend,                          &
                ivbeg,ivend,jvbeg,jvend,nvbeg,nvend,                    &
                korder,xs2d,ys2d,zpsoil,vtsoil,                         &
                vxs,vys,vzpsoil,iloc,jloc,                              &
                dxfld,dyfld,rdxfld,rdyfld,                              &
                slopey,alphay,betay,                                    &
                atsoil)

  CALL fldint3d(nx,ny,nzsoil,vnx,vny,vnzsoil,                           &
                ibeg,iend,jbeg,jend,nbeg,nend,                          &
                ivbeg,ivend,jvbeg,jvend,nvbeg,nvend,                    &
                korder,xs2d,ys2d,zpsoil,vqsoil,                         &
                vxs,vys,vzpsoil,iloc,jloc,                              &
                dxfld,dyfld,rdxfld,rdyfld,                              &
                slopey,alphay,betay,                                    &
                aqsoil)

!
  CALL fldint2d(nx,ny,vnx,vny,                                          &
                ibeg,iend,jbeg,jend,                                    &
                ivbeg,ivend,jvbeg,jvend,                                &
                korder,xs2d,ys2d,vwetcanp,vxs,vys,iloc,jloc,            &
                dxfld,dyfld,rdxfld,rdyfld,                              &
                slopey,alphay,betay,                                    &
                awetcanp)
!
  CALL fldint2d(nx,ny,vnx,vny,                                          &
                ibeg,iend,jbeg,jend,                                    &
                ivbeg,ivend,jvbeg,jvend,                                &
                korder,xs2d,ys2d,vraing,vxs,vys,iloc,jloc,              &
                dxfld,dyfld,rdxfld,rdyfld,                              &
                slopey,alphay,betay,                                    &
                araing)
!
  CALL fldint2d(nx,ny,vnx,vny,                                          &
                ibeg,iend,jbeg,jend,                                    &
                ivbeg,ivend,jvbeg,jvend,                                &
                korder,xs2d,ys2d,vrainc,vxs,vys,iloc,jloc,              &
                dxfld,dyfld,rdxfld,rdyfld,                              &
                slopey,alphay,betay,                                    &
                arainc)
!
  DO k=1,4
    CALL fldint2d(nx,ny,vnx,vny,                                        &
                  ibeg,iend,jbeg,jend,                                  &
                  ivbeg,ivend,jvbeg,jvend,                              &
                  korder,xs2d,ys2d,vprcrate(1,1,k),                     &
                  vxs,vys,iloc,jloc,                                    &
                  dxfld,dyfld,rdxfld,rdyfld,                            &
                  slopey,alphay,betay,                                  &
                  aprcrate(1,1,k))
  END DO

  CALL fldint2d(nx,ny,vnx,vny,                                          &
                ibeg,iend,jbeg,jend,                                    &
                ivbeg,ivend,jvbeg,jvend,                                &
                korder,xs2d,ys2d,vradsw,vxs,vys,iloc,jloc,              &
                dxfld,dyfld,rdxfld,rdyfld,                              &
                slopey,alphay,betay,                                    &
                aradsw)
!
  CALL fldint2d(nx,ny,vnx,vny,                                          &
                ibeg,iend,jbeg,jend,                                    &
                ivbeg,ivend,jvbeg,jvend,                                &
                korder,xs2d,ys2d,vrnflx,vxs,vys,iloc,jloc,              &
                dxfld,dyfld,rdxfld,rdyfld,                              &
                slopey,alphay,betay,                                    &
                arnflx)
  CALL fldint2d(nx,ny,vnx,vny,                                          &
                ibeg,iend,jbeg,jend,                                    &
                ivbeg,ivend,jvbeg,jvend,                                &
                korder,xs2d,ys2d,vradswnet,vxs,vys,iloc,jloc,           &
                dxfld,dyfld,rdxfld,rdyfld,                              &
                slopey,alphay,betay,                                    &
                aradswnet)
  CALL fldint2d(nx,ny,vnx,vny,                                          &
                ibeg,iend,jbeg,jend,                                    &
                ivbeg,ivend,jvbeg,jvend,                                &
                korder,xs2d,ys2d,vradlwin,vxs,vys,iloc,jloc,            &
                dxfld,dyfld,rdxfld,rdyfld,                              &
                slopey,alphay,betay,                                    &
                aradlwin)
!
  CALL fldint2d(nx,ny,vnx,vny,                                          &
                ibeg,iend,jbeg,jend,                                    &
                ivbeg,ivend,jvbeg,jvend,                                &
                korder,xs2d,ys2d,vusflx,vxs,vys,iloc,jloc,              &
                dxfld,dyfld,rdxfld,rdyfld,                              &
                slopey,alphay,betay,                                    &
                ausflx)
!
  CALL fldint2d(nx,ny,vnx,vny,                                          &
                ibeg,iend,jbeg,jend,                                    &
                ivbeg,ivend,jvbeg,jvend,                                &
                korder,xs2d,ys2d,vvsflx,vxs,vys,iloc,jloc,              &
                dxfld,dyfld,rdxfld,rdyfld,                              &
                slopey,alphay,betay,                                    &
                avsflx)
!
  CALL fldint2d(nx,ny,vnx,vny,                                          &
                ibeg,iend,jbeg,jend,                                    &
                ivbeg,ivend,jvbeg,jvend,                                &
                korder,xs2d,ys2d,vptsflx,vxs,vys,iloc,jloc,             &
                dxfld,dyfld,rdxfld,rdyfld,                              &
                slopey,alphay,betay,                                    &
                aptsflx)
!
  CALL fldint2d(nx,ny,vnx,vny,                                          &
                ibeg,iend,jbeg,jend,                                    &
                ivbeg,ivend,jvbeg,jvend,                                &
                korder,xs2d,ys2d,vqvsflx,vxs,vys,iloc,jloc,             &
                dxfld,dyfld,rdxfld,rdyfld,                              &
                slopey,alphay,betay,                                    &
                aqvsflx)
!
!-----------------------------------------------------------------------
!
!  Create array of verification heights at
!  forecast x,y locations
!
!-----------------------------------------------------------------------
!
  DO k=1,vnz-1
    CALL fldint2d(nx,ny,vnx,vny,                                        &
                  ibeg,iend,jbeg,jend,                                  &
                  ivbeg,ivend,jvbeg,jvend,                              &
                  korder,xs2d,ys2d,vzps(1,1,k),vxs,vys,iloc,jloc,       &
                  dxfld,dyfld,rdxfld,rdyfld,                            &
                  slopey,alphay,betay,                                  &
                  zpver(1,1,k))
  END DO
!
!-----------------------------------------------------------------------
!
!  Interpolate 3d scalar fields
!
!-----------------------------------------------------------------------
!
  CALL fldint3d(nx,ny,nz,vnx,vny,vnz,                                   &
                ibeg,iend,jbeg,jend,kbeg,kend,                          &
                ivbeg,ivend,jvbeg,jvend,kvbeg,kvend,                    &
                korder,xs2d,ys2d,zps,vptprt,                            &
                vxs,vys,zpver,iloc,jloc,                                &
                dxfld,dyfld,rdxfld,rdyfld,                              &
                slopey,alphay,betay,                                    &
                aptprt)
!
  CALL fldint3d(nx,ny,nz,vnx,vny,vnz,                                   &
                ibeg,iend,jbeg,jend,kbeg,kend,                          &
                ivbeg,ivend,jvbeg,jvend,kvbeg,kvend,                    &
                korder,xs2d,ys2d,zps,vpprt,                             &
                vxs,vys,zpver,iloc,jloc,                                &
                dxfld,dyfld,rdxfld,rdyfld,                              &
                slopey,alphay,betay,                                    &
                apprt)
!
  CALL fldint3d(nx,ny,nz,vnx,vny,vnz,                                   &
                ibeg,iend,jbeg,jend,kbeg,kend,                          &
                ivbeg,ivend,jvbeg,jvend,kvbeg,kvend,                    &
                korder,xs2d,ys2d,zps,vqvprt,                            &
                vxs,vys,zpver,iloc,jloc,                                &
                dxfld,dyfld,rdxfld,rdyfld,                              &
                slopey,alphay,betay,                                    &
                aqvprt)

  DO nq = 1, anscalar
    CALL fldint3d(nx,ny,nz,vnx,vny,vnz,                                 &
                  ibeg,iend,jbeg,jend,kbeg,kend,                        &
                  ivbeg,ivend,jvbeg,jvend,kvbeg,kvend,                  &
                  korder,xs2d,ys2d,zps,vqscalar(:,:,:,Qindex(nq)),      &
                  vxs,vys,zpver,iloc,jloc,                              &
                  dxfld,dyfld,rdxfld,rdyfld,                            &
                  slopey,alphay,betay,                                  &
                  aqscalar(:,:,:,nq))
  END DO

  CALL fldint3d(nx,ny,nz,vnx,vny,vnz,                                   &
                ibeg,iend,jbeg,jend,kbeg,kend,                          &
                ivbeg,ivend,jvbeg,jvend,kvbeg,kvend,                    &
                korder,xs2d,ys2d,zps,vtke,                              &
                vxs,vys,zpver,iloc,jloc,                                &
                dxfld,dyfld,rdxfld,rdyfld,                              &
                slopey,alphay,betay,                                    &
                atke)
!
  CALL fldint3d(nx,ny,nz,vnx,vny,vnz,                                   &
                ibeg,iend,jbeg,jend,kbeg,kend,                          &
                ivbeg,ivend,jvbeg,jvend,kvbeg,kvend,                    &
                korder,xs2d,ys2d,zps,vkmh,                              &
                vxs,vys,zpver,iloc,jloc,                                &
                dxfld,dyfld,rdxfld,rdyfld,                              &
                slopey,alphay,betay,                                    &
                akmh)
!
  CALL fldint3d(nx,ny,nz,vnx,vny,vnz,                                   &
                ibeg,iend,jbeg,jend,kbeg,kend,                          &
                ivbeg,ivend,jvbeg,jvend,kvbeg,kvend,                    &
                korder,xs2d,ys2d,zps,vkmv,                              &
                vxs,vys,zpver,iloc,jloc,                                &
                dxfld,dyfld,rdxfld,rdyfld,                              &
                slopey,alphay,betay,                                    &
                akmv)
!
  CALL fldint3d(nx,ny,nz,vnx,vny,vnz,                                   &
                ibeg,iend,jbeg,jend,kbeg,kend,                          &
                ivbeg,ivend,jvbeg,jvend,kvbeg,kvend,                    &
                korder,xs2d,ys2d,zps,vptbar,                            &
                vxs,vys,zpver,iloc,jloc,                                &
                dxfld,dyfld,rdxfld,rdyfld,                              &
                slopey,alphay,betay,                                    &
                aptbar)
!
  CALL fldint3d(nx,ny,nz,vnx,vny,vnz,                                   &
                ibeg,iend,jbeg,jend,kbeg,kend,                          &
                ivbeg,ivend,jvbeg,jvend,kvbeg,kvend,                    &
                korder,xs2d,ys2d,zps,vpbar,                             &
                vxs,vys,zpver,iloc,jloc,                                &
                dxfld,dyfld,rdxfld,rdyfld,                              &
                slopey,alphay,betay,                                    &
                apbar)
!
  CALL fldint3d(nx,ny,nz,vnx,vny,vnz,                                   &
                ibeg,iend,jbeg,jend,kbeg,kend,                          &
                ivbeg,ivend,jvbeg,jvend,kvbeg,kvend,                    &
                korder,xs2d,ys2d,zps,vrhobar,                           &
                vxs,vys,zpver,iloc,jloc,                                &
                dxfld,dyfld,rdxfld,rdyfld,                              &
                slopey,alphay,betay,                                    &
                arhobar)
!
  CALL fldint3d(nx,ny,nz,vnx,vny,vnz,                                   &
                ibeg,iend,jbeg,jend,kbeg,kend,                          &
                ivbeg,ivend,jvbeg,jvend,kvbeg,kvend,                    &
                korder,xs2d,ys2d,zps,vqvbar,                            &
                vxs,vys,zpver,iloc,jloc,                                &
                dxfld,dyfld,rdxfld,rdyfld,                              &
                slopey,alphay,betay,                                    &
                aqvbar)
!
  CALL fldint3d(nx,ny,nz,vnx,vny,vnz,                                   &
                ibeg,iend,jbeg,jend,kbeg,kend,                          &
                ivbeg,ivend,jvbeg,jvend,kvbeg,kvend,                    &
                korder,xs2d,ys2d,zps,vradfrc,                           &
                vxs,vys,zpver,iloc,jloc,                                &
                dxfld,dyfld,rdxfld,rdyfld,                              &
                slopey,alphay,betay,                                    &
                aradfrc)
!
  RETURN
END SUBROUTINE intsclrs
