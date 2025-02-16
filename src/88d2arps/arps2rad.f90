PROGRAM arps2rad
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   PROGRAM ARPS2RAD                   ######
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
!  Generate fake radar data from an ARPS history file.
!
!  It shares with the model the include file 'globcst.inc'
!  for storage parameters.
!
!  It reads in a history file produced by ARPS 4.0 in a user
!  specified format.
!
!  Parameters grdin,basin,mstin,icein,trbin are read in from the
!  data file itself, therefore are determined internally.
!  Arrays that are not read in retain their initial zero values.
!  These parameters are passed among subroutines through
!  a common block defined in 'indtflg.inc'.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS, December, 1995
!
!  MODIFICATION HISTORY:
!
!  12/06/2011 (Y. Wang)
!  Upgraded to arps5.3. Updated the reflectivity formula and output
!  interface.
!
!-----------------------------------------------------------------------
!
!  DATA ARRAYS READ IN:
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (km)
!    zp       z coordinate of grid points in computational space (m)
!    zpsoil   z coordinate of soil levels (m)
!    uprt     x component of perturbation velocity (m/s)
!    vprt     y component of perturbation velocity (m/s)
!    wprt     Vertical component of perturbation velocity in Cartesian
!             coordinates (m/s).
!
!    ptprt    Perturbation potential temperature (K)
!    pprt     Perturbation pressure (Pascal)
!
!    qvprt    Perturbation water vapor mixing ratio (kg/kg)
!    qc       Cloud water mixing ratio (kg/kg)
!    qr       Rainwater mixing ratio (kg/kg)
!    qi       Cloud ice mixing ratio (kg/kg)
!    qs       Snow mixing ratio (kg/kg)
!    qh       Hail mixing ratio (kg/kg)
!
!    ubar     Base state x velocity component (m/s)
!    vbar     Base state y velocity component (m/s)
!    wbar     Base state z velocity component (m/s)
!    ptbar    Base state potential temperature (K)
!    pbar     Base state pressure (Pascal)
!    rhobar   Base state air density (kg/m**3)
!    qvbar    Base state water vapor mixing ratio (kg/kg)
!
!  CALCULATED DATA ARRAYS:
!
!    u        x component of velocity (m/s)
!    v        y component of velocity (m/s)
!    w        z component of velocity (m/s)
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!
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
  INCLUDE 'grid.inc'
  INCLUDE 'indtflg.inc'
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  Arrays to be read in:
!
!-----------------------------------------------------------------------
!

  REAL, ALLOCATABLE  :: x     (:)         ! The x-coord. of the physical and
                                          ! computational grid. Defined at u-point.
  REAL, ALLOCATABLE  :: y     (:)         ! The y-coord. of the physical and
                                          ! computational grid. Defined at v-point.
  REAL, ALLOCATABLE  :: z     (:)         ! The z-coord. of the computational grid.
                                          ! Defined at w-point on the staggered grid.
  REAL, ALLOCATABLE  :: zp    (:,:,:)     ! The physical height coordinate defined at
                                          ! w-point of the staggered grid.
  REAL, ALLOCATABLE  :: zpsoil(:,:,:)     ! Soil level depth.
  REAL, ALLOCATABLE  :: j1    (:,:,:)     ! Coordinate transformation Jacobian defined
                                          ! as - d( zp )/d( x )
  REAL, ALLOCATABLE  :: j2    (:,:,:)     ! Coordinate transformation Jacobian defined
                                          ! as - d( zp )/d( y )
  REAL, ALLOCATABLE  :: j3    (:,:,:)     ! Coordinate transformation Jacobian defined
                                          ! as d( zp )/d( z )
  REAL, ALLOCATABLE  :: j3soil(:,:,:)     ! Coordinate transformation Jacobian defined
                                          ! as d( zp )/d( z )

  REAL, ALLOCATABLE  :: u     (:,:,:)     ! Total u-velocity (m/s)
  REAL, ALLOCATABLE  :: v     (:,:,:)     ! Total v-velocity (m/s)
  REAL, ALLOCATABLE  :: w     (:,:,:)     ! Total w-velocity (m/s)

  REAL, ALLOCATABLE  :: qscalar(:,:,:,:)

  REAL, ALLOCATABLE  :: tke   (:,:,:)     ! Turbulent kinetic energy
  REAL, ALLOCATABLE  :: kmh   (:,:,:)     ! The turbulent mixing coefficient for
                                          ! momentum. ( m**2/s ) horizontal
  REAL, ALLOCATABLE  :: kmv   (:,:,:)     ! The turbulent mixing coefficient for
                                          ! momentum. ( m**2/s ) vertical

  REAL, ALLOCATABLE  :: ubar  (:,:,:)     ! Base state u-velocity (m/s)
  REAL, ALLOCATABLE  :: vbar  (:,:,:)     ! Base state v-velocity (m/s)
  REAL, ALLOCATABLE  :: wbar  (:,:,:)     ! Base state w-velocity (m/s)
  REAL, ALLOCATABLE  :: ptbar (:,:,:)     ! Base state potential temperature (K)
  REAL, ALLOCATABLE  :: rhobar(:,:,:)     ! Base state air density (kg/m**3)
  REAL, ALLOCATABLE  :: pbar  (:,:,:)     ! Base state pressure (Pascal)
  REAL, ALLOCATABLE  :: qvbar (:,:,:)     ! Base state water vapor specific humidity
                                          ! (kg/kg)

  INTEGER, ALLOCATABLE  :: soiltyp(:,:,:)   ! Soil type
  INTEGER, ALLOCATABLE  :: vegtyp(:,:)    ! Vegetation type
  REAL, ALLOCATABLE  :: stypfrct(:,:,:)     ! Soil type fraction
  REAL, ALLOCATABLE  :: lai    (:,:)      ! Leaf Area Index
  REAL, ALLOCATABLE  :: roufns (:,:)      ! Surface roughness
  REAL, ALLOCATABLE  :: veg    (:,:)      ! Vegetation fraction

  REAL, ALLOCATABLE :: tsoil  (:,:,:,:)   ! soil temperature (K)
  REAL, ALLOCATABLE :: qsoil  (:,:,:,:)   ! soil moisture (g/kg)

  REAL, ALLOCATABLE  :: wetcanp(:,:)       ! Canopy water amount
  REAL, ALLOCATABLE  :: snowdpth(:,:)      ! Snow depth (m)

  REAL, ALLOCATABLE  :: raing  (:,:)       ! Grid supersaturation rain
  REAL, ALLOCATABLE  :: rainc  (:,:)       ! Cumulus convective rain
  REAL, ALLOCATABLE  :: prcrate(:,:,:)     ! precipitation rate (kg/(m**2*s))
                                           ! prcrate(1,1,1) = total precip. rate
                                           ! prcrate(1,1,2) = grid scale precip. rate
                                           ! prcrate(1,1,3) = cumulus precip. rate
                                           ! prcrate(1,1,4) = microphysics precip. rate

  REAL, ALLOCATABLE  :: radfrc(:,:,:)      ! Radiation forcing (K/s)
  REAL, ALLOCATABLE  :: radsw (:,:)        ! Solar radiation reaching the surface
  REAL, ALLOCATABLE  :: rnflx (:,:)        ! Net radiation flux absorbed by surface
  REAL, ALLOCATABLE :: radswnet (:,:) ! Net solar radiation, SWin - SWout
  REAL, ALLOCATABLE :: radlwin  (:,:) ! Incoming longwave radiation

  REAL, ALLOCATABLE  :: usflx (:,:)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE  :: vsflx (:,:)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE  :: ptsflx(:,:)        ! Surface heat flux (K*kg/(m*s**2))
  REAL, ALLOCATABLE  :: qvsflx(:,:)        ! Surface moisture flux (kg/(m**2*s))

  REAL, ALLOCATABLE  :: uprt  (:,:,:)     ! Perturbation u-velocity (m/s)
  REAL, ALLOCATABLE  :: vprt  (:,:,:)     ! Perturbation v-velocity (m/s)
  REAL, ALLOCATABLE  :: wprt  (:,:,:)     ! Perturbation w-velocity (m/s)
  REAL, ALLOCATABLE  :: ptprt (:,:,:)     ! Perturbation potential temperature (K)
  REAL, ALLOCATABLE  :: pprt  (:,:,:)     ! Perturbation pressure (Pascal)
  REAL, ALLOCATABLE  :: rhoprt(:,:,:)     ! Perturbation air density (kg/m**3)
  REAL, ALLOCATABLE  :: qvprt (:,:,:)     ! Perturbation water vapor specific
                                          ! humidity (kg/kg)

!
!-----------------------------------------------------------------------
!
!  Other data variables
!
!-----------------------------------------------------------------------
!
  REAL :: time
  REAL, ALLOCATABLE  :: tk (:,:,:)     ! Temperature (K)
!
!-----------------------------------------------------------------------
!
!  Temporary working arrays:
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE  :: tem1(:,:,:)
  REAL, ALLOCATABLE  :: tem2(:,:,:)
  REAL, ALLOCATABLE  :: tem3(:,:,:)
!
!-----------------------------------------------------------------------
!
!  ARPS dimensions:
!
!  nx, ny, nz: Dimensions of computatioal grid. When run on MPP
!              with PVM or MPI, they represent of the size of the
!              subdomains. See below.
!
!              Given nx, ny and nz, the physical domain size will be
!              xl=(nx-3)*dx by yl=(ny-3)*dy by zh=(nz-3)*dz. The
!              variables nx, ny, nz, dx, dy and dz are read in from
!              the input file by subroutine INITPARA.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx       ! Number of grid points in the x-direction
  INTEGER :: ny       ! Number of grid points in the y-direction
  INTEGER :: nz       ! Number of grid points in the z-direction
  INTEGER :: nzsoil   ! Number of soil levels
  INTEGER :: nstyps   ! Number of soil types
!
!-----------------------------------------------------------------------
!
!  User request stuff
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=4) :: radid
  CHARACTER (LEN=256) :: grdbasfn,filename
  CHARACTER (LEN=256) :: rfname
  REAL :: latrad,lonrad,elvrad
  REAL :: thref,elvmin,elvmax,rngmin,rngmax
!
  INTEGER :: hinfmt
  INTEGER :: dmpfmt,hdf4cmpr,rfopt
  INTEGER :: isource
  REAL :: refelvmin,refelvmax
  REAL :: refrngmin,refrngmax
!
!-----------------------------------------------------------------------
!
!  Scalar grid location variables
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE  :: lat(:,:),lon(:,:)

!-----------------------------------------------------------------------
!
!  Radar variables
!
!-----------------------------------------------------------------------
!
  INCLUDE 'remap.inc'

  REAL, ALLOCATABLE :: gridvel(:,:,:)
  REAL, ALLOCATABLE :: gridref(:,:,:)
  !REAL, ALLOCATABLE :: gridnyq(:,:,:)
  !REAL, ALLOCATABLE :: gridtim(:,:,:)
  REAL, ALLOCATABLE :: xs(:),ys(:)
  REAL, ALLOCATABLE :: zps(:,:,:)
!
!-----------------------------------------------------------------------
!
!  Radar angles
!
!-----------------------------------------------------------------------
!
  INTEGER :: maxtilt
  PARAMETER (maxtilt=14)
  REAL :: hlfwid
  PARAMETER (hlfwid=0.45)
  INTEGER :: ntilt
  REAL :: ang(maxtilt)
  REAL :: ang11(maxtilt),ang12(maxtilt)
  REAL :: ang31(maxtilt),ang32(maxtilt)
  DATA ang11 / 0.5,1.5,2.4,3.4,4.3,5.3,6.2,7.5,8.7,                     &
               10.0,12.0,14.0,16.7,19.5/
  DATA ang12 / 0.5,1.5,2.4,3.4,4.3,6.0,9.9,14.6,19.5,                   &
               19.5,19.5,19.5,19.5,19.5/
  DATA ang31 / 0.5,1.5,2.5,3.5,4.5,4.5,4.5,4.5,4.5,                     &
                4.5, 4.5, 4.5, 4.5, 4.5/
  DATA ang32 / 0.5,1.5,2.5,3.5,4.5,4.5,4.5,4.5,4.5,                     &
                4.5, 4.5, 4.5, 4.5, 4.5/
!
!-----------------------------------------------------------------------
!
! Variables for terminal velocity estimate
!
!-----------------------------------------------------------------------
!
  REAL :: denom,refz,rhofact,s1,s2,vt
  REAL, PARAMETER :: zfrez=3000.      ! freezing level (m)
  REAL, PARAMETER :: zice=8000.       ! level above which entirely ice
  REAL, PARAMETER :: h0=7000.
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: dmpzero = 0
  INTEGER, PARAMETER :: dualpout= 0
  CHARACTER (LEN=60) :: vcptxt
  INTEGER :: i,j,k,itilt,ireturn,mode
  INTEGER :: ngchan,nchanl,lenfil,lengbf,nchin
  INTEGER :: iradfmt,vcpnum,myr
  INTEGER :: knt,knttry,kntref,kntrng,kntang
  INTEGER :: abstsec,iyr,imon,iday,ihr,imin,isec
  REAL :: latnot(2)
  REAL :: gdx,gdy,xctr,yctr,x0sc,y0sc,radarx,radary,delx,dely
  REAL :: azim,dist,delz,eleva,range,dhdr,dsdr
  REAL :: ddrot,uazmrad,vazmrad
  REAL :: perc,prctry,prcref,prcdst,prcang
  INTEGER :: istatus
  LOGICAL :: angmatch

  INTEGER :: nxny, idat, ioffset, joffset, ncoltot
  INTEGER :: wrfopt

  LOGICAL, ALLOCATABLE :: havdat(:,:)

  INTEGER, ALLOCATABLE :: icolp(:)
  INTEGER, ALLOCATABLE :: jcolp(:)
  REAL,    ALLOCATABLE :: xcolp(:)
  REAL,    ALLOCATABLE :: ycolp(:)
  REAL,    ALLOCATABLE :: zcolp(:,:)

  REAL,    ALLOCATABLE :: colvel(:,:)
  REAL,    ALLOCATABLE :: colref(:,:)
  REAL,    ALLOCATABLE :: colnyq(:,:)
  REAL,    ALLOCATABLE :: coltim(:,:)

!-----------------------------------------------------------------------
!
! NAMELIST DECLARATION
!
!-----------------------------------------------------------------------

  NAMELIST /file_name/ hinfmt,grdbasfn,filename
  NAMELIST /radar_info/ radid, latrad, lonrad, elvrad, vcpnum, rngmin,   &
                        rngmax, rfopt, thref, dmpfmt, hdf4cmpr

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  WRITE(6,'(/11(/5x,a)/)')                                              &
     '###############################################################', &
     '###############################################################', &
     '###                                                         ###', &
     '###                  Welcome to ARPS2RAD                    ###', &
     '###      This program reads in the history dump data        ###', &
     '###      sets generated by ARPS, and generates a fake       ###', &
     '###      radar data file.                                   ###', &
     '###                                                         ###', &
     '###############################################################', &
     '###############################################################'

!-----------------------------------------------------------------------
!
!  Get the name of the input data set.
!
!-----------------------------------------------------------------------
!
  READ(5,file_name,END=100)
  WRITE(6,'(/a)')' Namelist file_name sucessfully read in.'

  lengbf=LEN_trim(grdbasfn)
  lenfil=LEN_trim(filename)

  WRITE(6,'(/a,a)')' Grid-base file name is ', grdbasfn(1:lengbf)
  WRITE(6,'(/a,a)')' The  data  set name is ', filename(1:lenfil)

  CALL get_dims_from_data(hinfmt,filename(1:lenfil),                    &
                          nx,ny,nz,nzsoil,nstyps, ireturn)

!-----------------------------------------------------------------------
!
!  Initialize variables.
!
!-----------------------------------------------------------------------
!
  istatus=0
  dmpfmt=3
  hdf4cmpr=1
  nstyp=nstyps
  denom=1./(zice-zfrez)

  ALLOCATE(x(nx))
  ALLOCATE(y(ny))
  ALLOCATE(z(nz))
  ALLOCATE(zp(nx,ny,nz))
  ALLOCATE(zpsoil(nx,ny,nzsoil))
  ALLOCATE(j1(nx,ny,nz))
  ALLOCATE(j2(nx,ny,nz))
  ALLOCATE(j3(nx,ny,nz))
  ALLOCATE(j3soil(nx,ny,nzsoil))
  x=0
  y=0
  z=0
  zp=0
  zpsoil=0
  j1=0
  j2=0
  j3=0
  j3soil=0

  ALLOCATE(u(nx,ny,nz))
  u=0
  ALLOCATE(v(nx,ny,nz))
  v=0
  ALLOCATE(w(nx,ny,nz))
  w=0
  ALLOCATE(qscalar(nx,ny,nz,nscalar), STAT = istatus)
  qscalar = 0.0
  ALLOCATE(tke(nx,ny,nz))
  tke=0
  ALLOCATE(kmh(nx,ny,nz))
  kmh=0
  ALLOCATE(kmv(nx,ny,nz))
  kmv=0
  ALLOCATE(ubar(nx,ny,nz))
  ubar=0
  ALLOCATE(vbar(nx,ny,nz))
  vbar=0
  ALLOCATE(wbar(nx,ny,nz))
  wbar=0
  ALLOCATE(ptbar(nx,ny,nz))
  ptbar=0
  ALLOCATE(rhobar(nx,ny,nz))
  rhobar=0
  ALLOCATE(pbar(nx,ny,nz))
  pbar=0
  ALLOCATE(qvbar(nx,ny,nz))
  qvbar=0
  ALLOCATE(soiltyp(nx,ny,nstyps))
  soiltyp=0
  ALLOCATE(stypfrct(nx,ny,nstyps))
  stypfrct=0
  ALLOCATE(vegtyp(nx,ny))
  vegtyp=0
  ALLOCATE(lai(nx,ny))
  lai=0
  ALLOCATE(roufns(nx,ny))
  roufns=0
  ALLOCATE(veg(nx,ny))
  veg=0
  ALLOCATE(tsoil   (nx,ny,nzsoil,0:nstyps))
  tsoil = 0
  ALLOCATE(qsoil   (nx,ny,nzsoil,0:nstyps))
  qsoil = 0

  ALLOCATE(wetcanp(nx,ny))
  wetcanp=0
  ALLOCATE(snowdpth(nx,ny))
  snowdpth=0
  ALLOCATE(raing(nx,ny))
  raing=0
  ALLOCATE(rainc(nx,ny))
  rainc=0
  ALLOCATE(prcrate(nx,ny,4))
  prcrate=0
  ALLOCATE(radfrc(nx,ny,nz))
  radfrc=0
  ALLOCATE(radsw(nx,ny))
  radsw=0
  ALLOCATE(rnflx(nx,ny))
  rnflx=0
  ALLOCATE(radswnet (nx,ny))
  radswnet = 0
  ALLOCATE(radlwin (nx,ny))
  radlwin = 0

  ALLOCATE(usflx(nx,ny))
  usflx=0
  ALLOCATE(vsflx(nx,ny))
  vsflx=0
  ALLOCATE(ptsflx(nx,ny))
  ptsflx=0
  ALLOCATE(qvsflx(nx,ny))
  qvsflx=0
  ALLOCATE(uprt(nx,ny,nz))
  uprt=0
  ALLOCATE(vprt(nx,ny,nz))
  vprt=0
  ALLOCATE(wprt(nx,ny,nz))
  wprt=0
  ALLOCATE(ptprt(nx,ny,nz))
  ptprt=0
  ALLOCATE(pprt(nx,ny,nz))
  pprt=0
  ALLOCATE(rhoprt(nx,ny,nz))
  rhoprt=0
  ALLOCATE(qvprt(nx,ny,nz))
  qvprt=0
  ALLOCATE(tk(nx,ny,nz))
  tk=0

  ALLOCATE(tem1(nx,ny,nz))
  tem1=0
  ALLOCATE(tem2(nx,ny,nz))
  tem2=0
  ALLOCATE(tem3(nx,ny,nz))
  tem3=0

  ALLOCATE(lat(nx,ny))
  lat=0
  ALLOCATE(lon(nx,ny))
  lon=0

!-------------------------------------------------------------
!
! Allocate for variables from remap_d.inc
!
!-------------------------------------------------------------

  ALLOCATE (gridvel(nx,ny,nz))
  gridvel=velmis
  ALLOCATE (gridref(nx,ny,nz))
  gridref=refmis
  !ALLOCATE (gridnyq(nx,ny,nz))
  !gridnyq=100.
  !ALLOCATE (gridtim(nx,ny,nz))
  !gridtim=0.
  ALLOCATE (xs(nx),ys(ny))
  ALLOCATE (zps(nx,ny,nz))

  ALLOCATE(havdat(nx,ny), STAT = istatus)
!
!------------------------------------------------------------------------
!
! Get the radar information
!
!------------------------------------------------------------------------

  READ(5,radar_info,END=100)
  WRITE(6,'(/a)')' Namelist radar_info sucessfully read in.'

!-----------------------------------------------------------------------
!
!  Set elevation angles
!
!-----------------------------------------------------------------------

  IF(vcpnum == 11) THEN
    ntilt=14
    vcptxt='Storm mode  14 tilts 0.5-19.5 deg'
    DO itilt=1,ntilt
      ang(itilt)=ang11(itilt)
    END DO
  ELSE IF (vcpnum == 12) THEN
    ntilt=9
    vcptxt='Storm mode   9 tilts 0.5-19.5 deg'
    DO itilt=1,ntilt
      ang(itilt)=ang11(itilt)
    END DO
  ELSE IF (vcpnum == 31) THEN
    ntilt=5
    vcptxt='Clear-air    5 tilts 0.5- 4.5 deg'
    DO itilt=1,ntilt
      ang(itilt)=ang11(itilt)
    END DO
  ELSE IF (vcpnum == 32) THEN
    ntilt=5
    vcptxt='Clear-air    5 tilts 0.5- 4.5 deg'
    DO itilt=1,ntilt
      ang(itilt)=ang11(itilt)
    END DO
  ELSE
    WRITE(6,*) vcpnum,' is a bogus vcpnum'
    STOP
  END IF
  elvmin=ang(1)-hlfwid
  elvmax=ang(ntilt)+hlfwid
  refelvmin=ang(1)
  refelvmax=ang(ntilt)
!
!-----------------------------------------------------------------------
!
!  Read all input data arrays
!
!-----------------------------------------------------------------------
!
  CALL dtaread(nx,ny,nz,nzsoil,nstyps,                                  &
               hinfmt,nchin,grdbasfn(1:lengbf),lengbf,                  &
               filename(1:lenfil),lenfil,time,                          &
               x,y,z,zp,zpsoil,uprt,vprt,wprt,ptprt,pprt,               &
               qvprt, qscalar, tke, kmh, kmv,                           &
               ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,            &
               soiltyp,stypfrct,vegtyp,lai,roufns,veg,                  &
               tsoil,qsoil,wetcanp,snowdpth,                            &
               raing,rainc,prcrate,                                     &
               radfrc,radsw,rnflx,radswnet,radlwin,                     &
               usflx,vsflx,ptsflx,qvsflx,                               &
               ireturn, tem1,tem2,tem3)

  curtim = time
!
!-----------------------------------------------------------------------
!
!  ireturn = 0 for a successful read
!
!-----------------------------------------------------------------------
!
  IF( ireturn == 0 ) THEN   ! successful read

!
!-----------------------------------------------------------------------
!
!  Set map projection parameters
!
!-----------------------------------------------------------------------
!
    latnot(1)=trulat1
    latnot(2)=trulat2
    CALL setmapr(mapproj,sclfct,latnot,trulon)
    CALL lltoxy(1,1,ctrlat,ctrlon,xctr,yctr)
!
    gdx=x(2)-x(1)
    x0sc=xctr - 0.5*(nx-3)*gdx
    gdy=y(2)-y(1)
    y0sc=yctr - 0.5*(ny-3)*gdy
    CALL setorig(1,x0sc,y0sc)
!
!-----------------------------------------------------------------------
!
!  Calculate lat,lon locations of the scalar grid points
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
    CALL lltoxy(1,1,latrad,lonrad,radarx,radary)
!
!-----------------------------------------------------------------------
!
!  Move z field onto the scalar grid.
!
!-----------------------------------------------------------------------
!
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          zps(i,j,k)=0.5*(zp(i,j,k)+zp(i,j,k+1))
        END DO
      END DO
    END DO
    DO j=1,ny-1
      DO i=1,nx-1
        zps(i,j,nz)=(2.*zps(i,j,nz-1))                                &
                       -zps(i,j,nz-2)
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Combine wind perturbations and mean fields and put them
!  on scalar grid.
!
!-----------------------------------------------------------------------
!
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          u(i,j,k)=0.5*(ubar(i,j,k)  +uprt(i,j,k) +                     &
                        ubar(i+1,j,k)+uprt(i+1,j,k))
        END DO
      END DO
    END DO
!
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          v(i,j,k)=0.5*(vbar(i,j,k)  +vprt(i,j,k) +                     &
                        vbar(i,j+1,k)+vprt(i,j+1,k))
        END DO
      END DO
    END DO
!
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          w(i,j,k)=0.5*(wbar(i,j,k)  +wprt(i,j,k) +                     &
                        wbar(i,j,k+1)+wprt(i,j,k+1))
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Get reflectivity from qr, qs and qh.
!
!-----------------------------------------------------------------------
!
    DO k=1, nz-1
      DO j=1, ny-1
        DO i=1, nx-1
          tem2(i,j,k) = ptbar(i,j,k)+ptprt(i,j,k)   ! tem2 is pressure
          tk(i,j,k)   = tem2(i,j,k) * ((tem2(i,j,k)/p0)**rddcp)
        END DO
      END DO
    END DO

    !IF (rfopt == 2) THEN
    !  CALL reflec_ferrier(nx,ny,nz, rhobar, qr, qs, qh, tk, tem1)
    !ELSE
    !  CALL reflec(nx,ny,nz, rhobar, qr, qs, qh, tem1)
    !ENDIF

    IF (rfopt >= 8 .AND. rfopt <= 12) THEN
      CALL reflec_MM(nx,ny,nz,rhobar,qscalar,tk,tem1)

    ! Lin,Schultz,Straka Lin,or WSM6 schemes
    ELSE IF (rfopt >= 2 .AND. rfopt <= 7) THEN
      CALL reflec_ferrier(nx,ny,nz, rhobar, qscalar, tk, tem1)

    ! Warm rain schemes
    ELSE IF (rfopt == 1) THEN
      IF(P_QR > 0) THEN
        CALL reflec_wr(nx,ny,nz, rhobar, qscalar(:,:,:,P_QR),tem1)
      ELSE
        tem1 = 0.0
      END IF

    !ELSE IF (rfopt > 100 .AND. rfopt < 200) THEN
    !  wrfopt = MOD(rfopt,100)
    !
    !  SELECT CASE (wrfopt)
    !
    !  CASE (1,3)     ! Kessler, WSM 3-class
    !    print *,' Calling reflec_wrf ...'
    !    CALL reflec_wrf(nx,ny,nz,qvprt,                                   &
    !         qscalar(:,:,:,P_QR),qscalar(:,:,:,P_QS),qscalar(:,:,:,P_QH), &
    !         0,tem2,tk,tem1)
    !  CASE (2,4,6) ! Lin, WSM 5- and 6- classes)
    !    print *,' Calling reflec_wrf ...'
    !    CALL reflec_wrf(nx,ny,nz,qvprt,                                   &
    !         qscalar(:,:,:,P_QR),qscalar(:,:,:,P_QS),qscalar(:,:,:,P_QH), &
    !         1,tem2,tk,tem1)
    !                                         ! mixed-phase
    !  CASE (5)     ! Ferrier microphysics scheme (as in WRF_POST)
    !    print *,' Calling reflec_ferrier_wrf ...'
    !    tem3(:,:,:) = 2.0
    !    CALL reflec_ferrier_wrf(nx,ny,nz,qvprt,                           &
    !         qscalar(:,:,:,P_QC),qscalar(:,:,:,P_QR),qscalar(:,:,:,P_QS), &
    !         tem2,tk,tem1,pprt,tem3)  ! pprt is temperory array here
    !  CASE (8) ! Thompson microphysics scheme (old version)
    !    print *,' Calling CALREF9s ...'
    !    CALL CALREF9s(nx,ny,nz,                                           &
    !         qscalar(:,:,:,P_QR),qscalar(:,:,:,P_QS),qscalar(:,:,:,P_QH), &
    !         tem2,tk,tem1)
    !  END SELECT

    ELSE
      print*,'Invalid microphysics option, reflectivity set to zero'
      tem1 = 0.0
    END IF
!
!-----------------------------------------------------------------------
!
!  Get radial velocity from 3-d velocity components
!
!-----------------------------------------------------------------------
!
    havdat(:,:) = .FALSE.

    knt=0
    kntref=0
    kntang=0
    kntrng=0
    DO j=2,ny-2
      DO i=2,nx-2
        delx=xs(i)-radarx
        dely=ys(j)-radary
        dist=sqrt(delx*delx+dely*dely)
        DO k=2,nz-2
          delz=zps(i,j,k)-elvrad
          CALL beamelv(delz,dist,eleva,range)
          IF(range > rngmin .AND. range < rngmax ) THEN
            IF(eleva > elvmin .AND. eleva < elvmax ) THEN

              angmatch=.false.
              DO itilt=1,ntilt
                IF(ABS(eleva-ang(itilt)) < hlfwid) THEN
                  angmatch=.true.
                  EXIT
                END IF
              END DO

              IF( angmatch ) THEN
                IF(tem1(i,j,k) > thref) THEN
                  havdat(i,j) = .TRUE.
                  knt=knt+1
                  gridref(i,j,k)=tem1(i,j,k)
!
! Find precipitation terminal velocity
! For now a fixed freezing level is used, reversible with rmvterm
!
                  refz=10.**(0.1*gridref(i,j,k))
                  rhofact=EXP(0.4*zps(i,j,k)/h0)
                  IF(zps(i,j,k) < zfrez) THEN
                    vt=2.6*(refz**0.107)*rhofact
                  ELSE IF(zps(i,j,k) < zice) THEN
                    s1=(zice-zps(i,j,k))*denom
                    s2=2.*(zps(i,j,k)-zfrez)*denom
                    vt=s1*2.6*(refz**0.107)*rhofact + s2
                  ELSE
                    vt=2.0
                  END IF
!
! Find local viewing angles
!
                  CALL dhdrange(eleva,range,dhdr)
                  dsdr=SQRT(AMAX1(0.,(1.-dhdr*dhdr)))
                  uazmrad=delx/dist
                  vazmrad=dely/dist
                  gridvel(i,j,k)=u(i,j,k)*uazmrad*dsdr                &
                                +v(i,j,k)*vazmrad*dsdr                &
                                +(w(i,j,k)-vt)*dhdr
                ELSE
                  kntref=kntref+1
                END IF
              ELSE
                kntang=kntang+1
              END IF

            ELSE
              kntang=kntang+1
            END IF
          ELSE
            kntrng=kntrng+1
          END IF
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Write the data counters.
!
!-----------------------------------------------------------------------
!
    knttry=(nx-3)*(ny-3)*(nz-3)
    perc  =100.*FLOAT(knt)/FLOAT(knttry)
    prcdst=100.*FLOAT(kntrng)/FLOAT(knttry)
    prcang=100.*FLOAT(kntang)/FLOAT(knttry)
    prcref=100.*FLOAT(kntref)/FLOAT(knttry)
!
    WRITE(6,'(a,i4/,a/)')       ' Used VCP Number',  vcpnum,vcptxt
    WRITE(6,'(a,i12)')          ' Total points:   ', knttry
    WRITE(6,'(a,i12,f6.1,a)')   ' Generated data: ', knt,perc,' %'
    WRITE(6,'(a,i12,f6.1,a)')   ' Range rejected: ', kntrng,prcdst,' %'
    WRITE(6,'(a,i12,f6.1,a)')   ' Angle samp rej: ', kntang,prcang,' %'
    WRITE(6,'(a,i12,f6.1,a//)') ' Ref thresh rej: ', kntref,prcref,' %'
!
!-----------------------------------------------------------------------
!
!  Write the radar data out.
!
!-----------------------------------------------------------------------
!
    iradfmt=1
    CALL ctim2abss(year,month,day,hour,minute,second, abstsec)
    abstsec=abstsec+curtim
    CALL abss2ctim( abstsec, iyr, imon, iday, ihr, imin, isec )

    myr=MOD(iyr,100)
    WRITE(rfname,'(a,a,i2.2,i2.2,i2.2,a,i2.2,i2.2)')                    &
        radid,'.',myr,imon,iday,'.',ihr,imin
    PRINT *, ' Writing data into ',TRIM(rfname)

    !CALL wtradcol(nx,ny,nz,dmpfmt,iradfmt,hdf4cmpr,rfname,              &
    !              radid,latrad,lonrad,elvrad,                           &
    !              iyr,imon,iday,ihr,imin,isec,vcpnum,isource,           &
    !              refelvmin,refelvmax,rngmin,rngmax,                    &
    !              xs,ys,zps,gridvel,gridref,gridnyq,gridtim,tem2d);
    !

    nxny=nx*ny
    ALLOCATE(icolp(nxny),stat=istatus)
    CALL check_alloc_status(istatus,'f88d2arps:icolp')
    ALLOCATE(jcolp(nxny),stat=istatus)
    CALL check_alloc_status(istatus,'f88d2arps:jcolp')
    ALLOCATE(xcolp(nxny),stat=istatus)
    CALL check_alloc_status(istatus,'f88d2arps:xcolp')
    ALLOCATE(ycolp(nxny),stat=istatus)
    CALL check_alloc_status(istatus,'f88d2arps:ycolp')
    ALLOCATE(zcolp(nz,nxny),stat=istatus)
    CALL check_alloc_status(istatus,'f88d2arps:zcolp')

    ALLOCATE(colref(nz,nxny),stat=istatus)
    CALL check_alloc_status(istatus,'f88d2arps:colref')
    ALLOCATE(colvel(nz,nxny),stat=istatus)
    CALL check_alloc_status(istatus,'f88d2arps:colvel')
    ALLOCATE(colnyq(nz,nxny),stat=istatus)
    CALL check_alloc_status(istatus,'f88d2arps:colnyq')
    ALLOCATE(coltim(nz,nxny),stat=istatus)
    CALL check_alloc_status(istatus,'f88d2arps:coltim')

    colref(:,:) = refmis
    colvel(:,:) = velmis
    colnyq(:,:) = -999.0
    coltim(:,:) = -999.0

    idat=0
    ioffset=0
    joffset=0

    DO j=2,ny-2
      DO i=2,nx-2
        IF(havdat(i,j)) THEN
          idat=idat+1
          icolp(idat)=i+ioffset
          jcolp(idat)=j+joffset
          xcolp(idat)=xs(i)
          ycolp(idat)=ys(j)
          DO k=1,nz
            zcolp(k,idat)=zps(i,j,k)
            colref(k,idat) = gridref(i,j,k)
            colvel(k,idat) = gridvel(i,j,k)
            colnyq(k,idat) = 100.0
            coltim(k,idat) = 0.0
          END DO
        END IF
      END DO
    END DO
    ncoltot = idat

    CALL wrtradcol(nx,ny,nxny,nz,ncoltot,ncoltot,                       &
                  dmpfmt,iradfmt,hdf4cmpr,dmpzero,dualpout,             &
                  rfname,radid,latrad,lonrad,elvrad,                    &
                  iyr,imon,iday,ihr,imin,isec,vcpnum,isource,           &
                  refelvmin,refelvmax,rngmin,rngmax,                    &
                  icolp,jcolp,xcolp,ycolp,zcolp,                        &
                  colref,colvel,colnyq,coltim,istatus)

  END IF                                       ! successful read

  GOTO 110

  100 CONTINUE

  Write(6,'(1x,a)')'Namelist read end encountered. Program stopped.'

  110 CONTINUE

  STOP
END PROGRAM arps2rad
