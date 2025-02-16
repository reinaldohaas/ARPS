!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE grafclose                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE grafclose
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  A dummy routine to be called instead of the functional routine
!  svidump when the savi3D library is not available.
!
!-----------------------------------------------------------------------
!
  RETURN

END SUBROUTINE grafclose
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE  MCLOSEDATASET             ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mclosedataset(gridid, ierr)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  A dummy routine to be called instead of the functional routine
!  svidump when the savi3D library is not available.
!
!-----------------------------------------------------------------------
!
  INTEGER :: gridid, ierr

  RETURN

END SUBROUTINE mclosedataset
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE  MCLOSESCHEME                 ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mclosescheme(gridid, ierr)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  A dummy routine to be called instead of the functional routine
!  svidump when the savi3D library is not available.
!
!-----------------------------------------------------------------------
!
  INTEGER :: gridid, ierr

  RETURN

END SUBROUTINE mclosescheme
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE SVIDUMP                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE svidump(nx,ny,nz,nzsoil,nstyps, nchout,graffn, grdbas,       &
           u,v,w,ptprt,pprt,qv,qc,qr,qi,qs,qh,tke,kmh,kmv,              &
           ubar,vbar,wbar,ptbar,pbar,rhobar,qvbar,                      &
           x,y,z,zp,zpsoil,                                             &
           soiltyp,stypfrct,vegtyp,lai,roufns,veg,                      &
           tsoil,qsoil,wetcanp,snowdpth,                                &
           raing,rainc,prcrate,                                         &
           radfrc,radsw,rnflx,                                          &
           usflx,vsflx,ptsflx,qvsflx,                                   &
           tem1, tem2)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Dump a data file for the visualization program Savi3D.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Jason J. Levit
!  07/20/92
!
!  MODIFICATION HISTORY:
!
!  08/02/92 (J. Levit)
!  Added full documentation and performed a clean-up.
!
!  08/04/92 (M. Xue)
!  Subroutine streamlined to conform to data dump format standards.
!
!  8/23/92 (M. Xue)
!  Modify to perform the dumping of both base and t-dependent arrays
!  and added control on grid staggering.
!
!  9/18/92 (J. Levit, M. Xue)
!  Added code to produce a configuration file for Savi3D.
!
!  11/2/92 (M. Xue)
!
!  Major overhaul. grafwritescalarpoint is called rather than
!  grafwritesclararray. This elliminates the need to define the
!  grid work array and make the data writing more flexible.
!
!  The capability to write out part of a data array implemented.
!
!  09/02/94 (J. Levit & Y. Lu)
!  Cleaned up documentation.
!
!  11/10/94 (Liping Sun & Min Zou)
!  Upgraded to version 1.2.2 of Savi3D. The file format was changed
!  from GRAF to MeRAF.
!
!  11/10/94 (Y. Liu)
!  Merged the upgraded version with the documentation cleaned up
!  version.
!
!  12/09/1998 (Donghai Wang)
!  Added the snow cover.
!
!  05/23/2002 (J. Brotzge)
!  Added variables for soil levels and new soil schemes.  
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
!    graffn   Name of the Savi3D MeRAF file.
!    grdbas   Flag indicating if this is a call for the data dump
!             of grid and base state arrays only. If so, grdbas=1
!             (not used in this routine).
!
!    u        x component of velocity at a given time level (m/s)
!    v        y component of velocity at a given time level (m/s)
!    w        Vertical component of Cartesian velocity at a given
!             time level (m/s)
!    ptprt    Perturbation potential temperature at a given time
!             level (K)
!    pprt     Perturbation pressure at  a given time level (Pascal)
!    qv       Water vapor specific humidity at a given time level
!             (kg/kg)
!    qc       Cloud water mixing ratio at a given time level (kg/kg)
!    qr       Rainwater mixing ratio at a given time level (kg/kg)
!    qi       Cloud ice mixing ratio at a given time level (kg/kg)
!    qs       Snow mixing ratio at a given time level (kg/kg)
!    qh       Hail mixing ratio at a given time level (kg/kg)
!
!    tke      Turbulent Kinetic Energy ((m/s)**2)
!    kmh      Horizontal turb. mixing coef. for momentum ( m**2/s )
!    kmv      Vertical turb. mixing coef. for momentum ( m**2/s )
!
!    ubar     Base state zonal velocity component (m/s)
!    vbar     Base state meridional velocity component (m/s)
!    wbar     Base state vertial velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhobar   Base state density (kg/m**3)
!    qvbar    Base state water vapor specific humidity (kg/kg)
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!    zp       Vertical coordinate of grid points in physical space(m)
!    zpsoil   Vertical coordinate of grid points in the soil (m)
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
!    grafgrid Passed from a dummy variable, used to define the grid
!             in Savi3D
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
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of grid points in the soil  
!
  CHARACTER (LEN=*     ) :: graffn ! Name of the Savi3D MeRAF file
  INTEGER :: grdbas            ! If this is a grid/base state dump
!
  REAL :: u     (nx,ny,nz)     ! Total u-velocity (m/s)
  REAL :: v     (nx,ny,nz)     ! Total v-velocity (m/s)
  REAL :: w     (nx,ny,nz)     ! Total w-velocity (m/s)
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)
!
  REAL :: qv    (nx,ny,nz)     ! Water vapor specific humidity (kg/kg)
  REAL :: qc    (nx,ny,nz)     ! Cloud water mixing ratio (kg/kg)
  REAL :: qr    (nx,ny,nz)     ! Rain water mixing ratio (kg/kg)
  REAL :: qi    (nx,ny,nz)     ! Cloud ice mixing ratio (kg/kg)
  REAL :: qs    (nx,ny,nz)     ! Snow mixing ratio (kg/kg)
  REAL :: qh    (nx,ny,nz)     ! Hail mixing ratio (kg/kg)
  REAL :: tke   (nx,ny,nz)     ! Turbulent Kinetic Energy ((m/s)**2)
!
  REAL :: kmh   (nx,ny,nz)     ! Horizontal turb. mixing coef. for
                               ! momentum. ( m**2/s )
  REAL :: kmv   (nx,ny,nz)     ! Vertical turb. mixing coef. for
                               ! momentum. ( m**2/s )
!
  REAL :: ubar  (nx,ny,nz)     ! Base state u-velocity (m/s)
  REAL :: vbar  (nx,ny,nz)     ! Base state v-velocity (m/s)
  REAL :: wbar  (nx,ny,nz)     ! Base state w-velocity (m/s)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL :: rhobar(nx,ny,nz)     ! Base state air density (kg/m**3)
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific
                               ! humidity (kg/kg)

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.
  REAL :: zp    (nx,ny,nz)     ! The physical height coordinate defined
                               ! at w-point of the staggered grid.
  REAL :: zpsoil (nx,ny,nzsoil) ! The physical height coordinate defined
                               ! at w-point of the soil.


  INTEGER :: nstyps
  INTEGER :: soiltyp (nx,ny,nstyps)    ! Soil type
  REAL :: stypfrct(nx,ny,nstyps)       ! Soil type
  INTEGER :: vegtyp(nx,ny)             ! Vegetation type
  REAL :: lai    (nx,ny)            ! Leaf Area Index
  REAL :: roufns (nx,ny)            ! Surface roughness
  REAL :: veg    (nx,ny)            ! Vegetation fraction

  REAL :: tsoil  (nx,ny,nzsoil,0:nstyps) ! Deep soil temperature (K)
  REAL :: qsoil  (nx,ny,nzsoil,0:nstyps) ! Deep soil moisture
  REAL :: wetcanp(nx,ny,0:nstyps)      ! Canopy water amount
  REAL :: snowdpth(nx,ny)              ! Snow depth (m)

  REAL :: raing(nx,ny)         ! Grid supersaturation rain
  REAL :: rainc(nx,ny)         ! Cumulus convective rain
  REAL :: prcrate(nx,ny,4)     ! precipitation rates (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precip. rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulus precip. rate
                               ! prcrate(1,1,4) = microphysics precip. rate

  REAL :: radfrc(nx,ny,nz)     ! Radiation forcing (K/s)
  REAL :: radsw (nx,ny)        ! Solar radiation reaching the surface
  REAL :: rnflx (nx,ny)        ! Net radiation flux absorbed by surface

  REAL :: usflx (nx,ny)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL :: vsflx (nx,ny)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL :: ptsflx(nx,ny)        ! Surface heat flux (K*kg/(m*s**2))
  REAL :: qvsflx(nx,ny)        ! Surface moisture flux (kg/(m**2*s))

  REAL :: tem1  (nx*ny*nz)     ! Temporary work array
  REAL :: tem2  (nx,ny,nz)     ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
!  include 'meraf.inc'
!  include 'globcst.inc'
!  include 'grafibm.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nxout,nyout,nzout ! The size of array to be written out.

  INTEGER :: ist ,ind ,isk ,jst ,jnd ,jsk ,kst ,knd ,ksk
  INTEGER :: ist1,ind1,isk1,jst1,jnd1,jsk1,kst1,knd1,ksk1

  INTEGER :: uid, vid, wid, uprtid, vprtid, wprtid        !*****
  INTEGER :: qvprtid, qcid, qrid, qwid, qiid, qsid, qhid  !*****
  INTEGER :: vortid, divid, ubarid, vbarid, wbarid        !*****
  INTEGER :: pbarid, rhobarid, qvbarid   !*****
  INTEGER :: windid, totalwindid         !*****
  INTEGER :: ptprtid, ptbarid, pprtid    !*****
  INTEGER :: value             ! *****
  INTEGER :: frame             ! *****

  INTEGER :: ierr              ! Used as an int'l error code by Savi3D.
  INTEGER :: i,j,k             ! Used by do loops.
  INTEGER :: ii,jj,kk
  INTEGER :: nchout            ! Unused.
  INTEGER :: ishf, jshf, kshf
  CHARACTER (LEN=50) :: configname   ! Used to create Savi3D config file.

  CHARACTER (LEN=20) :: schemename 
  CHARACTER (LEN=40) :: errorstring

  REAL*8 xbase,ybase,zbase  ! ****

  REAL :: conx, cony       ! Used to create Savi3D config file.
  INTEGER :: gbwrtn            ! See if grid and base state
                               ! parameter/arrays have been written
                               ! into the data file
  INTEGER :: ncalls
  DATA gbwrtn,ncalls /0,0/
  SAVE gbwrtn,ncalls
  CHARACTER (LEN=7) :: chtem2
  CHARACTER (LEN=7) :: chtem1
                                ! Used to create Savi3D config file.
  CHARACTER (LEN=6) :: timhms

!  integer year,month,day,hour,minute,second,node ! def. in globcst.inc
  INTEGER :: node
  REAL :: tem
  REAL :: xcord,ycord,zcord
  REAL*8 second1            ! Used for Savi3D
  INTEGER :: nchout0         ! Used to open Savi3D config file.

  INTEGER :: setdomn,setskip
  SAVE setdomn, setskip
  SAVE ist,ind,isk,jst,jnd,jsk,kst,knd,ksk
  DATA setdomn/0/, setskip /0/

  WRITE (6,'(/a/a,a/)')                                                 &
      'Option -io svi was not used when doing makearps.',               &
      'No Vis5d output was produced. ',                                 &
      'Re-do makearps with -io v5d option.',                            &
      'You also need to have Savi3D data dump library properly set up.'

  CALL arpsstop('Dummy version of SVIDUMP', 1)

  RETURN
END SUBROUTINE svidump
