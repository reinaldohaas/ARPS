PROGRAM arps2gem
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  PROGRAM ARSP2GEM                    ######
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
!  Converts ARPS history files to a GEMPAK format file.
!
!  Reads in a history file produced by ARPS in any ARPS format.
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  February, 1995
!
!  MODIFICATION HISTORY:
!  Added option for pressure-level output.
!  March, 1995   KB
!
!  Added sea-level pressure and w, combined rain and liquid water
!  quantities.   Added stability indices.  April, 1995   KB
!
!  Upgraded for new I/O routines including tke,kmh,kmv.
!  May, 1996, KB
!
!  3/19/98:
!  -Added option for specifying GEMPAK gdfile in namelist history_data
!  (Jonathan Case)
!  -Added igempr,igemz, and kintvl to namelist outopts, JC
!  -Added arps2gem.inc file for specfying pressure level output
!  -Fixed bugs in extrph subroutine, JC
!  -modified extrapolation routines such that all scalars are set
!   to zero for below ground values, except height and temp., JC
!  -Added data packing in GRIB/GEMPAK format (reduces file size)
!  -Added zps_km array for stability calculations (height must be in
!   km for use in arps_be subroutine), fixed bug in arps_be call, JC
!  -created separate arrays for grid and convective rains, JC
!  -Added relative humidity, surface pressure gridded output, JC
!  -Replaced UTM projection with MER projection, fixed bug in
!   setting up polar stereographic projection on GEMPAK grid, JC
!  -Kept qr, qc grids separate, Lumped QI,QS, and QH together
!   into QICE grid, JC.
!
!  1 June 2002 Eric Kemp
!  Soil variable updates.
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
!    w        vertical component of velocity in Cartesian
!             coordinates (m/s).
!
!    ptprt    perturbation potential temperature (K)
!    pprt     perturbation pressure (Pascal)
!    uprt     perturbation x velocity component (m/s)
!    vprt     perturbation y velocity component (m/s)
!    wprt     perturbation z velocity component (m/s)
!
!    qv       water vapor mixing ratio (kg/kg)
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
!    qsoil    soil moisture
!    wetcanp  Canopy water amount
!
!    rain     Total rain (raing + rainc)
!    raing    Grid supersaturation rain
!    rainc    Cumulus convective rain
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!  implicit none is commented out except for test compilations
!  due to the GEMPAK include file.   UNIDATA has been harrased
!  on that issue.
!
!-----------------------------------------------------------------------
!
! implicit none

  INTEGER :: nx,ny,nz          ! Grid dimensions.
  INTEGER :: nzsoil            ! soil levels
  INTEGER :: nstyps            ! Maximum number of soil types.

  INTEGER :: hinfmt,nhisfile_max,nhisfile,lengbf,nf,lenfil
  PARAMETER (nhisfile_max=200)
  CHARACTER (LEN=256) :: grdbasfn,hisfile(nhisfile_max)
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'arps2gem.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'GEMINC:GEMPRM.PRM'
!
!-----------------------------------------------------------------------
!
!  Arrays to be read in:
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: x(:)       ! The x-coord. of the physical and
                                  ! computational grid.
                                  ! Defined at u-point.
  REAL, ALLOCATABLE :: y(:)       ! The y-coord. of the physical and
                                  ! computational grid.
                                  ! Defined at v-point.
  REAL, ALLOCATABLE :: z(:)       ! The z-coord. of the computational
                                  ! grid. Defined at w-point.

  REAL, ALLOCATABLE :: zp(:,:,:)  ! The height of the terrain.
  REAL, ALLOCATABLE :: zpsoil(:,:,:)  ! The height of the soil level.

  REAL, ALLOCATABLE :: uprt   (:,:,:) ! Perturbation u-velocity (m/s)
  REAL, ALLOCATABLE :: vprt   (:,:,:) ! Perturbation v-velocity (m/s)
  REAL, ALLOCATABLE :: wprt   (:,:,:) ! Perturbation w-velocity (m/s)
  REAL, ALLOCATABLE :: ptprt  (:,:,:) ! Perturbation potential temperature (K)
  REAL, ALLOCATABLE :: pprt   (:,:,:) ! Perturbation pressure (Pascal)
  REAL, ALLOCATABLE :: qvprt  (:,:,:) ! Perturbation water vapor specific
                                      ! humidity (kg/kg)
  REAL, ALLOCATABLE :: qc     (:,:,:) ! Cloud water mixing ratio (kg/kg)
  REAL, ALLOCATABLE :: qr     (:,:,:) ! Rain water mixing ratio (kg/kg)
  REAL, ALLOCATABLE :: qi     (:,:,:) ! Cloud ice mixing ratio (kg/kg)

  REAL, ALLOCATABLE :: qscalar (:,:,:,:) ! Microphysics scalar array

  REAL, ALLOCATABLE :: tke    (:,:,:) ! Turbulent Kinetic Energy ((m/s)**2)
  REAL, ALLOCATABLE :: kmh    (:,:,:) ! Horizontal turb. mixing coef. for
                                      ! momentum. ( m**2/s )
  REAL, ALLOCATABLE :: kmv    (:,:,:) ! Vertical turb. mixing coef. for
                                      ! momentum. ( m**2/s )
  REAL, ALLOCATABLE :: ubar   (:,:,:) ! Base state u-velocity (m/s)
  REAL, ALLOCATABLE :: vbar   (:,:,:) ! Base state v-velocity (m/s)
  REAL, ALLOCATABLE :: wbar   (:,:,:) ! Base state w-velocity (m/s)
  REAL, ALLOCATABLE :: ptbar  (:,:,:) ! Base state potential temperature (K)
  REAL, ALLOCATABLE :: pbar   (:,:,:) ! Base state pressure (Pascal)
  REAL, ALLOCATABLE :: rhobar (:,:,:) ! Base state air density (kg/m**3)
  REAL, ALLOCATABLE :: qvbar  (:,:,:) ! Base state water vapor specific
                                      ! humidity (kg/kg)

  REAL, ALLOCATABLE :: u      (:,:,:) ! Total u-velocity (m/s)
  REAL, ALLOCATABLE :: v      (:,:,:) ! Total v-velocity (m/s)
  REAL, ALLOCATABLE :: w      (:,:,:) ! Total w-velocity (m/s)
  REAL, ALLOCATABLE :: qv     (:,:,:) ! Water vapor specific humidity (kg/kg)

  INTEGER, ALLOCATABLE :: soiltyp (:,:,:) ! Soil type
  REAL,    ALLOCATABLE :: stypfrct(:,:,:) ! Soil type
  INTEGER, ALLOCATABLE :: vegtyp(:,:)     ! Vegetation type
  REAL, ALLOCATABLE :: lai    (:,:)   ! Leaf Area Index
  REAL, ALLOCATABLE :: roufns (:,:)   ! Surface roughness
  REAL, ALLOCATABLE :: veg    (:,:)   ! Vegetation fraction

  REAL, ALLOCATABLE :: tsoil  (:,:,:,:) ! soil temperature (K)
  REAL, ALLOCATABLE :: qsoil  (:,:,:,:) ! soil moisture
  REAL, ALLOCATABLE :: wetcanp(:,:,:) ! Canopy water amount
  REAL, ALLOCATABLE :: snowdpth(:,:)  ! Snow depth (m)

  REAL, ALLOCATABLE :: rain (:,:)     ! Total rainfall
  REAL, ALLOCATABLE :: raing(:,:)     ! Grid supersaturation rain
  REAL, ALLOCATABLE :: rainc(:,:)     ! Cumulus convective rain
  REAL, ALLOCATABLE :: prcrate(:,:,:) ! precipitation rate (kg/(m**2*s))
                                      ! prcrate(1,1,1) = total precip. rate
                                      ! prcrate(1,1,2) = grid scale precip. rate
                                      ! prcrate(1,1,3) = cumulus precip. rate
                                      ! prcrate(1,1,4) = microphysics precip. rate

  REAL, ALLOCATABLE :: radfrc(:,:,:)  ! Radiation forcing (K/s)
  REAL, ALLOCATABLE :: radsw (:,:)    ! Solar radiation reaching the surface
  REAL, ALLOCATABLE :: rnflx (:,:)    ! Net radiation flux absorbed by surface
  REAL, ALLOCATABLE :: radswnet(:,:)  ! Net shortwave radiation
  REAL, ALLOCATABLE :: radlwin(:,:)   ! Incoming longwave radiation

  REAL, ALLOCATABLE :: usflx (:,:)    ! Surface flux of u-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: vsflx (:,:)    ! Surface flux of v-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: ptsflx(:,:)    ! Surface heat flux (K*kg/(m*s**2))
  REAL, ALLOCATABLE :: qvsflx(:,:)    ! Surface moisture flux (kg/(m**2*s))

!
! new stuff
!
  REAL, ALLOCATABLE :: e_mb   (:,:,:)    ! vapor pressure in mb
  REAL, ALLOCATABLE :: mix    (:,:,:)    ! total mixing ratio (kg/kg)
  REAL, ALLOCATABLE :: esat_mb(:,:,:)    ! saturation vapor pressure in mb
  REAL, ALLOCATABLE :: rh     (:,:,:)    ! Relative humidity in %
  REAL, ALLOCATABLE :: t_dew  (:,:,:)    ! dewpoint temp. in degrees K
!
!-----------------------------------------------------------------------
!
!  2-D stability index arrays
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: lcl(:,:)   ! lifting condensation level
  REAL, ALLOCATABLE :: lfc(:,:)   ! level of free convection
  REAL, ALLOCATABLE :: el(:,:)    ! equilibrium level
  REAL, ALLOCATABLE :: twdf(:,:)  ! max. wet bulb pot. temp. difference
  REAL, ALLOCATABLE :: li(:,:)    ! lifted index
  REAL, ALLOCATABLE :: pbe(:,:)   ! CAPE
  REAL, ALLOCATABLE :: mbe(:,:)   ! Moist CAPE
  REAL, ALLOCATABLE :: nbe(:,:)   ! CIN
  REAL, ALLOCATABLE :: tcap(:,:)  ! Cap Strength
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
  REAL, ALLOCATABLE :: tem2d(:,:)

  REAL, ALLOCATABLE :: wrk1(:),wrk2(:),wrk3(:),wrk4(:),wrk5(:),wrk6(:)
  REAL, ALLOCATABLE :: wrk7(:),wrk8(:),wrk9(:),wrk10(:),wrk11(:),wrk12(:)

  REAL, ALLOCATABLE :: xs(:)
  REAL, ALLOCATABLE :: ys(:)
  REAL, ALLOCATABLE :: zps(:,:,:)
!
!  height in km needed for arps_be
!
  REAL, ALLOCATABLE :: zps_km(:,:,:)
!
!-----------------------------------------------------------------------
!
!  Misc ARPS variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: nch
  REAL :: time
  REAL :: latnot(2)
!
!-----------------------------------------------------------------------
!
!  GEMPAK Output Levels...User Adjustable Parameters
!  to be put into the arps2gem.inc file and a new namelist
!
!-----------------------------------------------------------------------
!
!
!  update for GRIB packing of GEMPAK data
!
  integer nbits, ipktyp
  DATA nbits /16/
  DATA ipktyp /1/


!-----------------------------------------------------------------------
!
!  GEMPAK variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: gdfile
  CHARACTER (LEN=72)  :: gdarea
  REAL :: rnvblk (llnnav), anlblk (llnanl)
  INTEGER :: imxgrd
  PARAMETER (imxgrd=9999)
  INTEGER :: maxproj
  PARAMETER (maxproj=3)
  CHARACTER (LEN=3) :: cproj(0:maxproj)
  DATA cproj /'CED','STR','LCC','MER'/
  CHARACTER (LEN=72) :: chproj
  REAL :: zres
  PARAMETER(zres=10.)

  INTEGER, ALLOCATABLE :: zgem(:)

  REAL :: deltan
  REAL :: deltax,deltay
  REAL :: latll,lonll,latur,lonur
  REAL :: angle1,angle2,angle3
  LOGICAL :: angflg
  PARAMETER ( deltan= 1.,                                               &
              deltax= -9999.,                                           &
              deltay= -9999.,                                           &
              angflg=.true.)
  REAL :: gbnds(4)
  INTEGER :: ivsfc,ivprs,ivtheta,ivhgt
  PARAMETER (ivsfc=0,ivprs=1,ivtheta=2,ivhgt=3)
  INTEGER :: level(2)
  INTEGER :: ighdr(2)
  CHARACTER (LEN=20) :: gemtime(2)
  CHARACTER (LEN=12) :: parm
!
!-----------------------------------------------------------------------
!
!  Namelists
!
!-----------------------------------------------------------------------
!
  INTEGER :: igempr,igemz,kintvl
  NAMELIST /output/  gdfile, mstout,iceout,sfcout,igempr,               &
                     igemz,kintvl
!
!-----------------------------------------------------------------------
!
!  Physical constants
!  Normally these would be defined through include 'phycst.inc'
!  but there are some conflicts with the GEMPAK include files.
!
!-----------------------------------------------------------------------
!
  REAL :: rd        ! Gas constant for dry air  (m**2/(s**2*K))
  PARAMETER( rd     = 287.0 )
  REAL :: cp        ! Specific heat of dry air at constant pressure
                    ! (m**2/(s**2*K)).
  PARAMETER( cp     = 1004.0 )
  REAL :: rddcp
  PARAMETER( rddcp  = rd/cp )
  REAL :: p0        ! Surface reference pressure, is 100000 Pascal.
  PARAMETER( p0     = 1.0E5 )
!
!-----------------------------------------------------------------------
!
!  External functions
!
!-----------------------------------------------------------------------
!
  REAL :: wmr2td,oe,dpt
  EXTERNAL wmr2td,oe,dpt
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,nq,iret,igdfil,kk,klev,nzgem,ireturn, lens,navsz,ianlsz
  INTEGER :: is  ! loop index for soil type
  INTEGER :: iyr,ifhr,ifmin,ifile,ihd
  REAL :: rlnpgem,denom,prmb,tdew
  REAL :: ctrx,ctry,swx,swy,zsum,zmin,zmax,rzgem
  REAL :: gamma,ex2,rln700,p00

  REAL, ALLOCATABLE :: p_mb(:,:,:)

  LOGICAL :: wrtflg,rewrite,notopn
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Initialize GEMPAK sans TAE
!
!-----------------------------------------------------------------------
!
  CALL in_bdta ( iret )
!
  WRITE(6,'(6(/a))')                                                 &
      '###############################################################',&
      '#                                                             #',&
      '# Welcome to ARPS2GEM, a program that reads in history files  #',&
      '# generated by ARPS and produces a GEMPAK format file.        #',&
      '#                                                             #',&
      '###############################################################'

  101  CONTINUE
!
!-----------------------------------------------------------------------
!
!  Get the names of the input data files.
!
!-----------------------------------------------------------------------
!
  CALL get_input_file_names(5,hinfmt,grdbasfn,hisfile,nhisfile)

  lengbf = len_trim(grdbasfn)

  CALL get_dims_from_data(hinfmt,TRIM(hisfile(1)),                    &
       nx,ny,nz,nzsoil,nstyps, ireturn)

  IF (nstyps <= 0) nstyps = 1
  nstyp = nstyps ! Copy to global variable

  IF( ireturn /= 0 ) THEN
    PRINT*,'Problem occured when trying to get dimensions from data.'
    PRINT*,'Program stopped.'
    STOP
  END IF

  WRITE(6,'(4(a,i5))') 'nx =',nx,', ny=',ny,', nz=',nz,', nzsoil=',nzsoil

  READ(5,output,END=900)

  WRITE(6,'(a,i4)')'Option to store water data (1 or 0)   :',mstout
  WRITE(6,'(a,i4)')'Option to store ice data (1 or 0)     :',iceout
  WRITE(6,'(a,i4)')'Option to store sfc data (1 or 0)     :',sfcout
  WRITE(6,'(a,i4)')'Option to write pressure data (1 or 0):',igempr
  WRITE(6,'(a,i4)')'Option to write height data (1 or 0)  :',igemz
  WRITE(6,'(a,i4)')'Option to skip height levels (1 or 0) :',kintvl

  ALLOCATE(x      (nx))
  ALLOCATE(y      (ny))
  ALLOCATE(z      (nz))
  ALLOCATE(zp     (nx,ny,nz))
  ALLOCATE(zpsoil  (nx,ny,nzsoil))

  ALLOCATE(uprt   (nx,ny,nz))
  ALLOCATE(vprt   (nx,ny,nz))
  ALLOCATE(wprt   (nx,ny,nz))
  ALLOCATE(ptprt  (nx,ny,nz))
  ALLOCATE(pprt   (nx,ny,nz))
  ALLOCATE(qvprt  (nx,ny,nz))
  ALLOCATE(qc     (nx,ny,nz))
  ALLOCATE(qr     (nx,ny,nz))
  ALLOCATE(qi     (nx,ny,nz))
  ALLOCATE(qscalar(nx,ny,nz,nscalar))
  ALLOCATE(tke    (nx,ny,nz))
  ALLOCATE(kmh    (nx,ny,nz))
  ALLOCATE(kmv    (nx,ny,nz))
  ALLOCATE(ubar   (nx,ny,nz))
  ALLOCATE(vbar   (nx,ny,nz))
  ALLOCATE(wbar   (nx,ny,nz))
  ALLOCATE(ptbar  (nx,ny,nz))
  ALLOCATE(pbar   (nx,ny,nz))
  ALLOCATE(rhobar (nx,ny,nz))
  ALLOCATE(qvbar  (nx,ny,nz))
  ALLOCATE(u      (nx,ny,nz))
  ALLOCATE(v      (nx,ny,nz))
  ALLOCATE(w      (nx,ny,nz))
  ALLOCATE(qv     (nx,ny,nz))

  ALLOCATE(soiltyp (nx,ny,nstyps))
  ALLOCATE(stypfrct(nx,ny,nstyps))
  ALLOCATE(vegtyp (nx,ny))
  ALLOCATE(lai    (nx,ny))
  ALLOCATE(roufns (nx,ny))
  ALLOCATE(veg    (nx,ny))

  ALLOCATE(tsoil  (nx,ny,nzsoil,0:nstyps))
  ALLOCATE(qsoil  (nx,ny,nzsoil,0:nstyps))
  ALLOCATE(wetcanp(nx,ny,0:nstyps))
  ALLOCATE(snowdpth(nx,ny))

  ALLOCATE(rain (nx,ny))
  ALLOCATE(raing(nx,ny))
  ALLOCATE(rainc(nx,ny))
  ALLOCATE(prcrate(nx,ny,4))

  ALLOCATE(radfrc(nx,ny,nz))
  ALLOCATE(radsw (nx,ny))
  ALLOCATE(rnflx (nx,ny))
  ALLOCATE(radswnet(nx,ny))
  ALLOCATE(radlwin(nx,ny))

  ALLOCATE(usflx (nx,ny))
  ALLOCATE(vsflx (nx,ny))
  ALLOCATE(ptsflx(nx,ny))
  ALLOCATE(qvsflx(nx,ny))

  ALLOCATE(e_mb   (nx,ny,nz))
  ALLOCATE(mix    (nx,ny,nz))
  ALLOCATE(esat_mb(nx,ny,nz))
  ALLOCATE(rh     (nx,ny,nz))
  ALLOCATE(t_dew  (nx,ny,nz))

  ALLOCATE(lcl(nx,ny))
  ALLOCATE(lfc(nx,ny))
  ALLOCATE(el(nx,ny))
  ALLOCATE(twdf(nx,ny))
  ALLOCATE(li(nx,ny))
  ALLOCATE(pbe(nx,ny))
  ALLOCATE(mbe(nx,ny))
  ALLOCATE(nbe(nx,ny))
  ALLOCATE(tcap(nx,ny))

  ALLOCATE(tem1(nx,ny,nz))
  ALLOCATE(tem2(nx,ny,nz))
  ALLOCATE(tem3(nx,ny,nz))
  ALLOCATE(tem2d(nx,ny))

  ALLOCATE(wrk1(nz),wrk2(nz),wrk3(nz),wrk4(nz),wrk5(nz),wrk6(nz))
  ALLOCATE(wrk7(nz),wrk8(nz),wrk9(nz),wrk10(nz),wrk11(nz),wrk12(nz))

  ALLOCATE(xs(nx))
  ALLOCATE(ys(ny))
  ALLOCATE(zps(nx,ny,nz))
  ALLOCATE(zps_km(nx,ny,nz))

  ALLOCATE(zgem(nz))
  ALLOCATE(p_mb(nx,ny,nz))


  x      =0.0
  y      =0.0
  z      =0.0
  zp     =0.0

  uprt   =0.0
  vprt   =0.0
  wprt   =0.0
  ptprt  =0.0
  pprt   =0.0
  qvprt  =0.0
  qc     =0.0
  qr     =0.0
  qi     =0.0
  qscalar=0.0
  tke    =0.0
  kmh    =0.0
  kmv    =0.0
  ubar   =0.0
  vbar   =0.0
  wbar   =0.0
  ptbar  =0.0
  pbar   =0.0
  rhobar =0.0
  qvbar  =0.0
  u      =0.0
  v      =0.0
  w      =0.0
  qv     =0.0

  soiltyp =0.0
  stypfrct=0.0
  vegtyp =0.0
  lai    =0.0
  roufns =0.0
  veg    =0.0

  tsoil  =0.0
  qsoil  =0.0
  wetcanp=0.0
  snowdpth=0.0

  rain =0.0
  raing=0.0
  rainc=0.0
  prcrate=0.0

  radfrc=0.0
  radsw =0.0
  rnflx =0.0
  radswnet =0.0
  radlwin =0.0

  usflx =0.0
  vsflx =0.0
  ptsflx=0.0
  qvsflx=0.0

  e_mb    =0.0
  mix     =0.0
  esat_mb =0.0
  rh      =0.0
  t_dew   =0.0

  lcl =0.0
  lfc =0.0
  el  =0.0
  twdf=0.0
  li  =0.0
  pbe =0.0
  mbe =0.0
  nbe =0.0
  tcap=0.0

  tem1 =0.0
  tem2 =0.0
  tem3 =0.0
  tem2d=0.0

  wrk1 = 0.0
  wrk2 = 0.0
  wrk3 = 0.0
  wrk4 = 0.0
  wrk5 = 0.0
  wrk6 = 0.0
  wrk7 = 0.0
  wrk8 = 0.0
  wrk9 = 0.0
  wrk10= 0.0
  wrk11= 0.0
  wrk12= 0.0

  xs=0.0
  ys=0.0
  zps=0.0
  zps_km=0.0

  zgem = 0
  p_mb =0.0


  notopn=.true.
  DO ifile=1,nhisfile
    lenfil = LEN_trim(hisfile(ifile))
    WRITE(6,'(/a,/1x,a)')' The data set name is ',                      &
        hisfile(ifile)(1:lenfil)

!-----------------------------------------------------------------------
!
!  Read all input data arrays
!
!-----------------------------------------------------------------------
!
    CALL dtaread(nx,ny,nz,nzsoil,nstyps,hinfmt,nch,grdbasfn(1:lengbf),  &
               lengbf,                                                  &
               hisfile(ifile)(1:lenfil),lenfil,time,                    &
               x,y,z,zp,zpsoil,uprt,vprt,wprt,ptprt,pprt,               &
               qvprt, qscalar, tke, kmh, kmv,                           &
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
!
      curtim=time
!
!  MODIFICATION
!
      DO i=1,nx
        DO j=1,ny
          rain(i,j)=raing(i,j)+rainc(i,j)
        END DO
      END DO

!
!  MODIFICATION
!

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            pprt(i,j,k)=pprt(i,j,k)+pbar(i,j,k)
            ptprt(i,j,k)=ptprt(i,j,k)+ptbar(i,j,k)
            qvprt(i,j,k)=qvprt(i,j,k)+qvbar(i,j,k)
            tem1(i,j,k)=0.5*(uprt(i,j,k)+ubar(i,j,k)+                   &
                        uprt(i+1,j,k)+ubar(i+1,j,k))
            tem2(i,j,k)=0.5*(vprt(i,j,k)+vbar(i,j,k)+                   &
                        vprt(i,j+1,k)+vbar(i,j+1,k))
            tem3(i,j,k)=0.5*(wprt(i,j,k)+wbar(i,j,k)+                   &
                        wprt(i,j,k+1)+wbar(i,j,k+1))
!            qi(i,j,k)=qi(i,j,k)+qs(i,j,k)+qh(i,j,k)
            qi(i,j,k) = 0.0
            IF(P_QI > 0) THEN
              qi(i,j,k) = qi(i,j,k) + qscalar(i,j,k,P_QI)
            END IF
            IF(P_QS > 0) THEN
              qi(i,j,k) = qi(i,j,k) + qscalar(i,j,k,P_QS)
            END IF
            IF(P_QG > 0) THEN
              qi(i,j,k) = qi(i,j,k) + qscalar(i,j,k,P_QG)
            END IF
            IF(P_QH > 0) THEN
              qi(i,j,k) = qi(i,j,k) + qscalar(i,j,k,P_QH)
            END IF

            qc(i,j,k) = 0.0
            IF(P_QC > 0) THEN
              qc(i,j,k) = qscalar(i,j,k,P_QC)
            END IF

            qr(i,j,k) = 0.0
            IF(P_QR > 0) THEN
              qr(i,j,k) = qscalar(i,j,k,P_QR)
            END IF

          END DO
        END DO
      END DO

!
!  Swap wind data back into wind arrays
!  Change qv to mixing ratio
!
! MODIFICATION (LEAVE AS SPECIFIC HUMIDITY)

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            uprt(i,j,k)=tem1(i,j,k)
            vprt(i,j,k)=tem2(i,j,k)
            wprt(i,j,k)=tem3(i,j,k)
          END DO
        END DO
      END DO

      CALL edgfill(pprt, 1,nx,1,nx-1,1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL edgfill(ptprt,1,nx,1,nx-1,1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL edgfill(qvprt,1,nx,1,nx-1,1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL edgfill(qr,   1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)
      CALL edgfill(qi,   1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)
      CALL edgfill(uprt, 1,nx,1,nx-1,1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL edgfill(vprt, 1,nx,1,nx-1,1,ny,1,ny-1, 1,nz,1,nz-1)
      CALL edgfill(wprt, 1,nx,1,nx-1,1,ny,1,ny-1, 1,nz,1,nz-1)
      DO is=0,nstyps
        CALL edgfill(tsoil(1,1,1,is),1,nx,1,nx-1,1,ny,1,ny-1,          &
                     1,nzsoil,1,nzsoil)
        CALL edgfill(qsoil(1,1,1,is),1,nx,1,nx-1,1,ny,1,ny-1,          &
                     1,nzsoil,1,nzsoil)
      END DO
      CALL edgfill(wetcanp,1,nx,1,nx-1,1,ny,1,ny-1, 1,1,1,1)
!
!-----------------------------------------------------------------------
!
!  Build the  GEMPAK grid time string
!  It has format yymodd/hhmnFHHH
!  yy: year      mo: month  dd: GMT day
!  hh: GMT hour  mn: minute
!  F: seperation charcter
!  HHH: forecast hour (000 = analysis)
!  example  time(1)='950126/1200F000'
!
!-----------------------------------------------------------------------
!
      gemtime(1)='                    '
      gemtime(2)='                    '
      iyr=MOD(year,100)
      ifhr=INT(curtim/3600.)
      ifmin=nint((curtim-(ifhr*3600.))/60.)
      IF(curtim == 0) THEN
        WRITE(gemtime(1),                                               &
            '(i2.2,i2.2,i2.2,a1,i2.2,i2.2,a4)')                         &
            iyr,month,day,'/',hour,minute,'F000'
      ELSE
        WRITE(gemtime(1),                                               &
            '(i2.2,i2.2,i2.2,a1,i2.2,i2.2,a1,i3.3,i2.2)')               &
            iyr,month,day,'/',hour,minute,'F',ifhr,ifmin
      END IF
      WRITE(6,'(a,a)') '  GEMPAK time string ',gemtime(1)
!
!-----------------------------------------------------------------------
!
!  Initialize header, grid area coordinates and analysis block
!
!-----------------------------------------------------------------------
!
      IF(notopn) THEN
        ihd = 2
        ighdr(1)=0
        ighdr(2)=0
!
        DO i=1,llnanl
          anlblk(i) = 0.
        END DO
!
!-----------------------------------------------------------------------
!
!  Establish coordinate for scalar fields.
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
!  Assign constant height surfaces for GEMPAK output.
!  All output heights are rounded to nearest zres (normally 10 m).
!  The first level is the minimum height of scalar level 2.
!  The last level is the maximum height of scalar level nz-1.
!
!-----------------------------------------------------------------------
!
        denom=FLOAT((nx-1)*(ny-1))
        zmin=zps(1,1,2)
        zsum=0.
        DO j=1,ny-1
          DO i=1,nx-1
            zmin=AMIN1(zmin,zps(i,j,2))
            zsum=zsum+zps(i,j,2)
          END DO
        END DO
        zgem(1)=nint(nint(zmin/zres)*zres)
        zsum=zsum/denom
        zgem(2)=nint(nint(zsum/zres)*zres)
        IF(zgem(2) > zgem(1)) THEN
          kk=2
        ELSE
          kk=1
        END IF

        DO klev=3,nz-2
          kk=kk+1
          zsum=0.
          DO j=1,ny-1
            DO i=1,nx-1
              zsum=zsum+zps(i,j,klev)
            END DO
          END DO
          zsum=zsum/denom
          zgem(kk)=nint(nint(zsum/zres)*zres)
        END DO

        kk=kk+1
        zmax=zps(1,1,nz-1)
        zsum=0.
        DO j=1,ny-1
          DO i=1,nx-1
            zmax=AMAX1(zmax,zps(i,j,nz-1))
            zsum=zsum+zps(i,j,nz-1)
          END DO
        END DO

        zsum=zsum/denom
        zgem(kk)=nint(nint(zsum/zres)*zres)
        kk=kk+1
        zgem(kk)=nint(nint(zmax/zres)*zres)
        IF(zgem(kk) > zgem(kk-1)) THEN
          nzgem=kk
        ELSE
          nzgem=kk-1
        END IF
!
!-----------------------------------------------------------------------
!
!  Build navigation block
!
!-----------------------------------------------------------------------

        latnot(1)=trulat1
        latnot(2)=trulat2
        CALL setmapr( mapproj, 1.0 , latnot , trulon)

        CALL lltoxy( 1,1, ctrlat,ctrlon, ctrx, ctry )
        dx=x(3)-x(2)
        dy=y(3)-y(2)
        swx = ctrx - (FLOAT(nx-3)/2.) * dx
        swy = ctry - (FLOAT(ny-3)/2.) * dy
        CALL setorig( 1, swx, swy)
!
        CALL xytoll(1,1,xs(1),ys(1),latll,lonll)
        CALL xytoll(1,1,xs(nx),ys(ny),latur,lonur)
!
        chproj=cproj(mapproj)
        angle2=trulon
!
!  modify projection if polar stereographic
!
        IF (chproj == 'STR') THEN
          angle1=90.0
          angle3=0.0
          WRITE (*,*) '**********************************'
          WRITE (*,*) '**********************************'
          WRITE (*,*) ' '
          WRITE (*,*) 'SETTING ANGLE1=90. FOR GEMPAK FILE'
          WRITE (*,*) 'SETTING ANGLE3=0.0 FOR GEMPAK FILE'
          WRITE (*,*) '(GEMPAK origin is at North Pole)'
          WRITE (*,*) ' '
          WRITE (*,*) '**********************************'
          WRITE (*,*) '**********************************'
        ELSE
          angle1=trulat1
          angle3=trulat2
          WRITE (*,*) '**********************************'
          WRITE (*,*) '**********************************'
          WRITE (*,*) ' '
          WRITE (*,*) '      SETTING ANGLE1=TRULAT1      '
          WRITE (*,*) '      SETTING ANGLE3=TRULAT2      '
          WRITE (*,*) ' '
          WRITE (*,*) '**********************************'
          WRITE (*,*) '**********************************'
        END IF
!
        gbnds(1)=latll
        gbnds(2)=lonll
        gbnds(3)=latur
        gbnds(4)=lonur
!
        WRITE(gdarea,'(f8.4,a1,f9.4,a1,f8.4,a1,f9.4)')                  &
            latll,';',lonll,';',latur,';',lonur
!
        CALL gr_mnav  ( chproj, nx, ny, latll, lonll, latur, lonur,     &
                        angle1, angle2, angle3, angflg,                 &
                        rnvblk, iret )

        IF( iret /= 0 )  GO TO 950
!
!-----------------------------------------------------------------------
!
!  Build analysis block
!
!-----------------------------------------------------------------------
!
        CALL gr_mban  ( deltan, deltax, deltay,                         &
                          gbnds, gbnds, gbnds, anlblk, iret )

        IF( iret /= 0 )  GO TO 950
!
!-----------------------------------------------------------------------
!
!  Open/Create the grid file
!
!-----------------------------------------------------------------------
!
        CALL st_lstr  ( gdfile, lens, iret )
        WRITE (6, *) 'Opening GEMPAK file ',gdfile(1:lens)
        CALL gd_opnf  ( gdfile, wrtflg, igdfil, navsz, rnvblk,          &
                       ianlsz, anlblk, ihd, imxgrd, iret )
        IF ( iret /= 0 ) THEN
          WRITE (6, *) 'Error opening existing file',gdfile(1:lens)
          WRITE (6, *) 'Creating file ',gdfile(1:lens)
          CALL gd_cref  ( gdfile, llnnav, rnvblk, llnanl,               &
                         anlblk, ihd, imxgrd, igdfil, iret )

          IF( iret /= 0 )  GO TO 950
        END IF
!
!-----------------------------------------------------------------------
!
!  Set write flag to true
!
!-----------------------------------------------------------------------
!
        wrtflg=.true.
        CALL gd_swrt( igdfil, wrtflg, iret )
        IF( iret /= 0 )  GO TO 950
        notopn=.false.
!
!-----------------------------------------------------------------------
!
!  Output terrain data.
!
!-----------------------------------------------------------------------
!
        level(1)=0
        level(2)=-1
        PRINT *, ' Writing terrain data.'
        parm='ESFC'
        CALL gd_wpgd( igdfil, zp(1,1,2),                                &
                      nx, ny, ighdr,                                    &
                      gemtime, level, ivsfc, parm,                      &
                      rewrite, ipktyp, nbits, iret)
      END IF  ! notopn
!
!-----------------------------------------------------------------------
!
!  Put temperature into tem2 and -ln(p) into tem3
!
!-----------------------------------------------------------------------
!
      DO k=1,nz-1
        DO j=1,ny
          DO i=1,nx
            tem2(i,j,k)=ptprt(i,j,k)*                                   &
                          (pprt(i,j,k)/p0)**rddcp
            tem3(i,j,k)=-ALOG(pprt(i,j,k))
          END DO
        END DO
      END DO

      CALL edgfill(tem2,1,nx,1,nx-1,1,ny,1,ny-1,1,nz,1,nz-1)

!
! ------- modification ------------
!

      DO i=1,nx
        DO j=1,ny
          DO k=1,nz
            p_mb(i,j,k)=0.01*pprt(i,j,k)
            mix(i,j,k)=1000.0*(qvprt(i,j,k)/(1.-qvprt(i,j,k)))
            e_mb(i,j,k)=qvprt(i,j,k)*p_mb(i,j,k)/(qvprt(i,j,k)+         &
                        287.0/461.0)
            esat_mb(i,j,k)=6.1078*EXP((2.500780E6/461.0)*((1.0/273.15)- &
                          (1.0/tem2(i,j,k))))
            t_dew(i,j,k)=wmr2td(p_mb(i,j,k),mix(i,j,k))                 &
                       + 273.15
            rh(i,j,k)=100.0*(e_mb(i,j,k)/esat_mb(i,j,k))
            IF (rh(i,j,k) < 0) THEN
              rh(i,j,k)=0.0
            END IF
          END DO
        END DO
      END DO

!
!-----------------------------------------------------------------------
!
!  Calculate stability indices.
!  Use level k=2 as the "surface".
!
!-----------------------------------------------------------------------
!
!
!  convert height (zps) to km units, since that it what
!  it appears thermo3d.f uses in the arps_be routine.
!
      DO i=1,nx
        DO j=1,ny
          DO k=1,nz-1
            zps_km(i,j,k)=zps(i,j,k)/1000.0
          END DO
        END DO
      END DO

      WRITE(*,*) 'About to enter the dreaded arps_be'
      WRITE(*,*) 'subroutine.'
      WRITE(*,*) 'This will take a while.'
      WRITE(*,*) 'Please be patient'
      WRITE(*,*) ' '
      WRITE(*,*) 'Processing.............'

      CALL arps_be(nx,ny,nz,                                            &
           pprt,zps_km,tem2,qvprt,                                      &
           lcl,lfc,el,twdf,li,pbe,mbe,nbe,tcap,                         &
           wrk1,wrk2,wrk3,wrk4,wrk5,wrk6,                               &
           wrk7,wrk8,wrk9,wrk10,wrk11,wrk12,tem2d)


      WRITE(*,*) ' '
      WRITE(*,*) ' '
      WRITE(*,*) ' '
      WRITE(*,*) ' '
      WRITE(*,*) 'Now done with stability calculations.'
      WRITE(*,*) ' '

      CALL edgfill(lcl, 1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
      CALL edgfill(lfc, 1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
      CALL edgfill(el,  1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
      CALL edgfill(twdf,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
      CALL edgfill(li,  1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
      CALL edgfill(pbe, 1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
      CALL edgfill(mbe, 1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
      CALL edgfill(nbe, 1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)
      CALL edgfill(tcap,1,nx,1,nx-1,1,ny,1,ny-1,1,1,1,1)

      PRINT *, ' Sample stability values: '
      PRINT *, '   lcl,lfc : ',lcl(1,1),lfc(1,1)
      PRINT *, '   el, twdf: ',el(1,1), twdf(1,1)
      PRINT *, '   li, pbe : ',li(1,1), pbe(1,1)
      PRINT *, '   mbe, nbe: ',mbe(1,1),nbe(1,1)
      PRINT *, '   tcap    : ',tcap(1,1)
!
!  Store k=2 theta-e and dewpt in tem1,
!  level 1 and 2 respectively.
!
      DO j=1,ny
        DO i=1,nx
          prmb=0.01*pprt(i,j,2)
          tdew=wmr2td(prmb,(1000.*qvprt(i,j,2)))
          tem1(i,j,1)=oe((tem2(i,j,2)-273.15),tdew,prmb) + 273.15
          tem1(i,j,2)=tdew+273.15
          tem1(i,j,3)=prmb
          tem1(i,j,4)=raing(i,j)
          tem1(i,j,5)=rainc(i,j)
          tem1(i,j,6)=rain(i,j)
        END DO
      END DO
!
!-----------------------------------------------------------------------
!
!  Output near-sfc data.
!  Simularity theory or something could be applied here
!  to make these data be valid at sfc instrument height,
!  for now we output level 2.
!
!-----------------------------------------------------------------------
!
      level(1)=0
      level(2)=-1
      PRINT *, ' Writing near-sfc pressure'
      parm='PRES'
      CALL gd_wpgd( igdfil, tem1(1,1,3), nx, ny, ighdr,                 &
                    gemtime, level, ivsfc, parm,                        &
                    rewrite, ipktyp, nbits, iret)
!
! MODIFICATION
!
      PRINT *,' Writing grid-scale rainfall'
      parm='RAIN_G'
      CALL gd_wpgd( igdfil, tem1(1,1,4),                                &
                    nx, ny, ighdr,                                      &
                    gemtime, level, ivsfc, parm,                        &
                    rewrite, ipktyp, nbits, iret)

      PRINT *,' Writing convective rainfall'
      parm='RAIN_C'
      CALL gd_wpgd( igdfil, tem1(1,1,5),                                &
                    nx, ny, ighdr,                                      &
                    gemtime, level, ivsfc, parm,                        &
                    rewrite, ipktyp, nbits, iret)

      PRINT *,' Writing total accumulated rainfall'
      parm='RAIN'
      CALL gd_wpgd( igdfil, tem1(1,1,6),                                &
                    nx, ny, ighdr,                                      &
                    gemtime, level, ivsfc, parm,                        &
                    rewrite, ipktyp, nbits, iret)


      PRINT *, ' Writing near-sfc temperature'
      parm='TMPK'
      CALL gd_wpgd( igdfil, tem2(1,1,2),                                &
                    nx, ny, ighdr,                                      &
                    gemtime, level, ivsfc, parm,                        &
                    rewrite, ipktyp, nbits, iret)
!
      PRINT *, ' Writing near-sfc dew point temperature'
      parm='DWPK'
      CALL gd_wpgd( igdfil, tem1(1,1,2),                                &
                    nx, ny, ighdr,                                      &
                    gemtime, level, ivsfc, parm,                        &
                    rewrite, ipktyp, nbits, iret)
!
      PRINT *,'Writing near-sfc specific humidity'
      parm='QV'
      CALL gd_wpgd( igdfil, qvprt(1,1,2),                               &
                    nx, ny, ighdr,                                      &
                    gemtime, level, ivsfc, parm,                        &
                    rewrite, ipktyp, nbits, iret)
!
      PRINT *, ' Writing near-sfc u velocity'
      parm='UREL'
      CALL gd_wpgd( igdfil, uprt(1,1,2),                                &
                    nx, ny, ighdr,                                      &
                    gemtime, level, ivsfc, parm,                        &
                    rewrite, ipktyp, nbits, iret)
!
      PRINT *, ' Writing near-sfc v velocity'
      parm='VREL'
      CALL gd_wpgd( igdfil, vprt(1,1,2),                                &
                    nx, ny, ighdr,                                      &
                    gemtime, level, ivsfc, parm,                        &
                    rewrite, ipktyp, nbits, iret)
!
      PRINT *, ' Writing near-sfc theta-e'
      parm='THTE'
      CALL gd_wpgd( igdfil, tem1(1,1,1),                                &
                    nx, ny, ighdr,                                      &
                    gemtime, level, ivsfc, parm,                        &
                    rewrite, ipktyp, nbits, iret)

      PRINT *, ' Writing near-sfc LI'
      parm='LIFT'
      CALL gd_wpgd( igdfil, li,                                         &
                    nx, ny, ighdr,                                      &
                    gemtime, level, ivsfc, parm,                        &
                    rewrite, ipktyp, nbits, iret)

      PRINT *, ' Writing near-sfc CAPE'
      parm='CAPE'
      CALL gd_wpgd( igdfil, pbe,                                        &
                    nx, ny, ighdr,                                      &
                    gemtime, level, ivsfc, parm,                        &
                    rewrite, ipktyp, nbits, iret)

      PRINT *, ' Writing near-sfc moist CAPE'
      parm='MPBE'
      CALL gd_wpgd( igdfil, mbe,                                        &
                    nx, ny, ighdr,                                      &
                    gemtime, level, ivsfc, parm,                        &
                    rewrite, ipktyp, nbits, iret)

      PRINT *, ' Writing near-sfc CIN'
      parm='SCIN'
      CALL gd_wpgd( igdfil, nbe,                                        &
                    nx, ny, ighdr,                                      &
                    gemtime, level, ivsfc, parm,                        &
                    rewrite, ipktyp, nbits, iret)

      PRINT *, ' Writing near-sfc LCL'
      parm='HLCL'
      CALL gd_wpgd( igdfil, lcl,                                        &
                    nx, ny, ighdr,                                      &
                    gemtime, level, ivsfc, parm,                        &
                    rewrite, ipktyp, nbits, iret)
!
      PRINT *, ' Writing near-sfc LFC'
      parm='HLFC'
      CALL gd_wpgd( igdfil, lfc,                                        &
                    nx, ny, ighdr,                                      &
                    gemtime, level, ivsfc, parm,                        &
                    rewrite, ipktyp, nbits, iret)
!
      PRINT *, ' Writing near-sfc EL'
      parm='HSEL'
      CALL gd_wpgd( igdfil, lcl,                                        &
                    nx, ny, ighdr,                                      &
                    gemtime, level, ivsfc, parm,                        &
                    rewrite, ipktyp, nbits, iret)
!
      PRINT *, ' Writing wet bulb temp diff'
      parm='TWDF'
      CALL gd_wpgd( igdfil, twdf,                                       &
                    nx, ny, ighdr,                                      &
                    gemtime, level, ivsfc, parm,                        &
                    rewrite, ipktyp, nbits, iret)
!
      PRINT *, ' Writing Cap Strength'
      parm='TCAP'
      CALL gd_wpgd( igdfil, tcap,                                       &
                    nx, ny, ighdr,                                      &
                    gemtime, level, ivsfc, parm,                        &
                    rewrite, ipktyp, nbits, iret)
!
!-----------------------------------------------------------------------
!
!  Output single-level data.
!
!-----------------------------------------------------------------------
!
      IF (sfcout == 1) THEN
        PRINT *, ' Writing skin temperature data.'
        parm='SKTK'
        CALL gd_wpgd( igdfil,tsoil(1,1,1,0),                            &
                    nx, ny, ighdr,                                      &
                    gemtime, level, ivsfc, parm,                        &
                    rewrite, ipktyp, nbits, iret)
        PRINT *, ' Writing soil temp data.'
        parm='SLTK'
        CALL gd_wpgd( igdfil,tsoil(1,1,2,0),                            &
                      nx, ny, ighdr,                                    &
                      gemtime, level, ivsfc, parm,                      &
                      rewrite, ipktyp, nbits, iret)
        PRINT *, ' Writing soil wetness data.'
        parm='SLWT'
        CALL gd_wpgd( igdfil,qsoil(1,1,1,0),                            &
                      nx, ny, ighdr,                                    &
                      gemtime, level, ivsfc, parm,                      &
                      rewrite, ipktyp, nbits, iret)
        PRINT *, ' Writing deep soil wetness data.'
        parm='DPWT'
        CALL gd_wpgd( igdfil,qsoil(1,1,2,0),                            &
                      nx, ny, ighdr,                                    &
                      gemtime, level, ivsfc, parm,                      &
                      rewrite, ipktyp, nbits, iret)
        PRINT *, ' Writing canopy wetness data.'
        parm='CNWT'
        CALL gd_wpgd( igdfil, wetcanp,                                  &
                      nx, ny, ighdr,                                    &
                      gemtime, level, ivsfc, parm,                      &
                      rewrite, ipktyp, nbits, iret)
      END IF
!
!-----------------------------------------------------------------------
!
!  Output constant pressure level data
!
!-----------------------------------------------------------------------
!
      IF( igempr == 1 ) THEN
        rewrite=.false.
!
!-----------------------------------------------------------------------
!
!  Put temperature into tem2 and -ln(p) into tem3.
!
!-----------------------------------------------------------------------
!
        DO k=1,nz-1
          DO j=1,ny
            DO i=1,nx
              tem2(i,j,k)=ptprt(i,j,k)*(pprt(i,j,k)/p0)**rddcp
              tem3(i,j,k)=-ALOG(pprt(i,j,k))
            END DO
          END DO
        END DO
!
!
!-----------------------------------------------------------------------
!
!    Calculate temperature (K) at ARPS grid points
!   and 700mb pressure level
!
!-----------------------------------------------------------------------
!
        rln700=-ALOG(70000.0)
        CALL v2dinta(nx,ny,nz,1,nx,1,ny,1,nz-1,tem2,tem3,rln700,tem1)
!
!-----------------------------------------------------------------------
!
!    Calculate sea level pressure (mb)
!    Reduction method: Benjamin and Miller: 1990, MWR, vol.118, No.10,
!                   Page: 2100-2101
!
!-----------------------------------------------------------------------
!
        gamma=.0065      ! std lapse rate per meter
        ex2=5.2558774
        DO j=1,ny
          DO i=1,nx
            p00 = 0.01*(pprt(i,j,2))
            tem1(i,j,1)=p00*(1.0+gamma*zps(i,j,2)/tem1(i,j,1))**ex2
          END DO
        END DO
        level(1)=0
        level(2)=-1
        PRINT *, ' Writing MSL Pressure '
        parm='PMSL'
        CALL gd_wpgd( igdfil, tem1,                                     &
                      nx, ny, ighdr,                                    &
                      gemtime, level, ivsfc, parm,                      &
                      rewrite, ipktyp, nbits, iret)
!
!-----------------------------------------------------------------------
!
!    Calculate stability variables.
!    This is from the Oklahoma LAPS (O'LAPS) surface analysis
!    software.
!
!-----------------------------------------------------------------------
!
        DO klev=1,nprgem
          level(1)=iprgem(klev)
          level(2)=-1
          rlnpgem=-ALOG(100.*FLOAT(iprgem(klev)))
          PRINT *, ' Writing GEMPAK data at pr= ',iprgem(klev)
          parm='HGHT'
          CALL v2dinta(nx,ny,nz,1,nx,1,ny,1,nz-1,                       &
                       zps,tem3,rlnpgem,tem1)
          CALL extrph(nx,ny,nz,zps,tem2,pprt,                           &
                      iprgem(klev),tem1)
          CALL gd_wpgd( igdfil, tem1,                                   &
                        nx, ny, ighdr,                                  &
                        gemtime, level, ivprs, parm,                    &
                        rewrite, ipktyp, nbits, iret)
          parm='TMPK'
          CALL v2dinta(nx,ny,nz,1,nx,1,ny,1,nz-1,                       &
                        tem2,tem3,rlnpgem,tem1)
          CALL extrpt(nx,ny,nz,tem2,pprt,zps,                           &
                      iprgem(klev),tem1)
          CALL gd_wpgd( igdfil, tem1,                                   &
                        nx, ny, ighdr,                                  &
                        gemtime, level, ivprs, parm,                    &
                        rewrite, ipktyp, nbits, iret)
          parm='DWPK'
          CALL v2dinta(nx,ny,nz,1,nx,1,ny,1,nz-1,                       &
                        t_dew,tem3,rlnpgem,tem1)
          CALL extrpt(nx,ny,nz,t_dew,pprt,zps,                          &
                      iprgem(klev),tem1)
          CALL gd_wpgd( igdfil, tem1,                                   &
                        nx, ny, ighdr,                                  &
                        gemtime, level, ivprs, parm,                    &
                        rewrite, ipktyp, nbits, iret)
          parm='QV'
          CALL v2dinta(nx,ny,nz,1,nx,1,ny,1,nz-1,                       &
                       qvprt,tem3,rlnpgem,tem1)
          CALL extrpq(nx,ny,nz,qvprt,pprt,iprgem(klev),tem1)
          CALL gd_wpgd( igdfil, tem1,                                   &
                        nx, ny, ighdr,                                  &
                        gemtime, level, ivprs, parm,                    &
                        rewrite, ipktyp, nbits, iret)
          parm='RELH'
          CALL v2dinta(nx,ny,nz,1,nx,1,ny,1,nz-1,                       &
                       rh,tem3,rlnpgem,tem1)
          CALL extrpq(nx,ny,nz,rh,pprt,iprgem(klev),tem1)
          CALL gd_wpgd( igdfil, tem1,                                   &
                        nx, ny, ighdr,                                  &
                        gemtime, level, ivprs, parm,                    &
                        rewrite, ipktyp, nbits, iret)
          CALL v2dinta(nx,ny,nz,1,nx,1,ny,1,nz-1,                       &
                      uprt,tem3,rlnpgem,tem1(1,1,1))
          CALL v2dinta(nx,ny,nz,1,nx,1,ny,1,nz-1,                       &
                      vprt,tem3,rlnpgem,tem1(1,1,2))
          CALL extrpuv(nx,ny,nz,uprt,vprt,pprt,zps,                     &
              iprgem(klev),tem1(1,1,1),tem1(1,1,2))
          parm='UREL'
          CALL gd_wpgd( igdfil, tem1(1,1,1),                            &
                        nx, ny, ighdr,                                  &
                        gemtime, level, ivprs, parm,                    &
                        rewrite, ipktyp, nbits, iret)
          parm='VREL'
          CALL gd_wpgd( igdfil, tem1(1,1,2),                            &
                        nx, ny, ighdr,                                  &
                        gemtime, level, ivprs, parm,                    &
                        rewrite, ipktyp, nbits, iret)
          CALL v2dinta(nx,ny,nz,1,nx,1,ny,1,nz-1,                       &
                       wprt,tem3,rlnpgem,tem1)
          CALL extrpq(nx,ny,nz,wprt,pprt,iprgem(klev),tem1)
          parm='WWND'
          CALL gd_wpgd( igdfil, tem1,                                   &
                        nx, ny, ighdr,                                  &
                        gemtime, level, ivprs, parm,                    &
                        rewrite, ipktyp, nbits, iret)
          IF (mstout == 1) THEN
            parm='QCLD'
            CALL v2dinta(nx,ny,nz,1,nx,1,ny,1,nz-1,                     &
                                 qc,tem3,rlnpgem,tem1)
            CALL extrpq(nx,ny,nz,qc,pprt,iprgem(klev),tem1)
            CALL gd_wpgd( igdfil, tem1,                                 &
                          nx, ny, ighdr,                                &
                          gemtime, level, ivprs, parm,                  &
                          rewrite, ipktyp, nbits, iret)
            parm='QRAIN'
            CALL v2dinta(nx,ny,nz,1,nx,1,ny,1,nz-1,                     &
                                 qr,tem3,rlnpgem,tem1)
            CALL extrpq(nx,ny,nz,qr,pprt,iprgem(klev),tem1)
            CALL gd_wpgd( igdfil, tem1,                                 &
                          nx, ny, ighdr,                                &
                          gemtime, level, ivprs, parm,                  &
                          rewrite, ipktyp, nbits, iret)
            IF (iceout == 1 .and. P_QI > 0) THEN
              parm='QICE'
              CALL v2dinta(nx,ny,nz,1,nx,1,ny,1,nz-1,                   &
                                   qi,tem3,rlnpgem,tem1)
              CALL extrpq(nx,ny,nz,qi,pprt,iprgem(klev),tem1)
              CALL gd_wpgd( igdfil, tem1,                               &
                            nx, ny, ighdr,                              &
                            gemtime, level, ivprs, parm,                &
                            rewrite, ipktyp, nbits, iret)
            END IF ! iceout
          END IF   ! mstout
        END DO
      END IF ! igempr
!
!-----------------------------------------------------------------------
!
!  Output constant height level data
!
!-----------------------------------------------------------------------
!
      IF( igemz == 1 ) THEN
        rewrite=.false.
!
!  Put temperature into tem2
!
        DO k=1,nz-1
          DO j=1,ny
            DO i=1,nx
              tem2(i,j,k)=ptprt(i,j,k)*                                 &
                          (pprt(i,j,k)/p0)**rddcp
            END DO
          END DO
        END DO
!
        DO klev=2,nzgem,kintvl
          level(1)=zgem(klev)
          level(2)=-1
          rzgem=FLOAT(zgem(klev))
          PRINT *, ' Writing GEMPAK data at z= ',zgem(klev)
          parm='PRES'
          CALL v2dinta(nx,ny,nz,1,nx,1,ny,1,nz-1,                       &
                               pprt,zps,rzgem,tem1)
          CALL gd_wpgd( igdfil, tem1,                                   &
                        nx, ny, ighdr,                                  &
                        gemtime, level, ivhgt, parm,                    &
                        rewrite, ipktyp, nbits, iret)
          parm='TMPK'
          CALL v2dinta(nx,ny,nz,1,nx,1,ny,1,nz-1,                       &
                               tem2,zps,rzgem,tem1)
          CALL gd_wpgd( igdfil, tem1,                                   &
                        nx, ny, ighdr,                                  &
                        gemtime, level, ivhgt, parm,                    &
                        rewrite, ipktyp, nbits, iret)
          parm='DWPK'
          CALL v2dinta(nx,ny,nz,1,nx,1,ny,1,nz-1,                       &
                               t_dew,zps,rzgem,tem1)
          CALL gd_wpgd( igdfil, tem1,                                   &
                        nx, ny, ighdr,                                  &
                        gemtime, level, ivhgt, parm,                    &
                        rewrite, ipktyp, nbits, iret)
          parm='QV'
          CALL v2dinta(nx,ny,nz,1,nx,1,ny,1,nz-1,                       &
                               qvprt,zps,rzgem,tem1)
          CALL gd_wpgd( igdfil, tem1,                                   &
                        nx, ny, ighdr,                                  &
                        gemtime, level, ivhgt, parm,                    &
                        rewrite, ipktyp, nbits, iret)
          parm='RELH'
          CALL v2dinta(nx,ny,nz,1,nx,1,ny,1,nz-1,                       &
                       rh,tem3,rlnpgem,tem1)
          CALL gd_wpgd( igdfil, tem1,                                   &
                        nx, ny, ighdr,                                  &
                        gemtime, level, ivprs, parm,                    &
                        rewrite, ipktyp, nbits, iret)
          parm='UREL'
          CALL v2dinta(nx,ny,nz,1,nx,1,ny,1,nz-1,                       &
                               uprt,zps,rzgem,tem1)
          CALL gd_wpgd( igdfil, tem1,                                   &
                        nx, ny, ighdr,                                  &
                        gemtime, level, ivhgt, parm,                    &
                        rewrite, ipktyp, nbits, iret)
          parm='VREL'
          CALL v2dinta(nx,ny,nz,1,nx,1,ny,1,nz-1,                       &
                               vprt,zps,rzgem,tem1)
          CALL gd_wpgd( igdfil, tem1,                                   &
                        nx, ny, ighdr,                                  &
                        gemtime, level, ivhgt, parm,                    &
                        rewrite, ipktyp, nbits, iret)
          parm='WWND'
          CALL v2dinta(nx,ny,nz,1,nx,1,ny,1,nz-1,                       &
                               wprt,zps,rzgem,tem1)
          CALL gd_wpgd( igdfil, tem1,                                   &
                        nx, ny, ighdr,                                  &
                        gemtime, level, ivprs, parm,                    &
                        rewrite, ipktyp, nbits, iret)
          IF (mstout == 1) THEN
            parm='QCLD'
            CALL v2dinta(nx,ny,nz,1,nx,1,ny,1,nz-1,                     &
                                 qc,zps,rzgem,tem1)
            CALL gd_wpgd( igdfil, tem1,                                 &
                          nx, ny, ighdr,                                &
                          gemtime, level, ivhgt, parm,                  &
                          rewrite, ipktyp, nbits, iret)
            parm='QRAIN'
            CALL v2dinta(nx,ny,nz,1,nx,1,ny,1,nz-1,                     &
                                 qr,zps,rzgem,tem1)
            CALL gd_wpgd( igdfil, tem1,                                 &
                            nx, ny, ighdr,                              &
                          gemtime, level, ivhgt, parm,                  &
                          rewrite, ipktyp, nbits, iret)
            IF (iceout == 1 .and. P_QI > 0) THEN
              parm='QICE'
              CALL v2dinta(nx,ny,nz,1,nx,1,ny,1,nz-1,                   &
                                   qi,zps,rzgem,tem1)
              CALL gd_wpgd( igdfil, tem1,                               &
                            nx, ny, ighdr,                              &
                            gemtime, level, ivhgt, parm,                &
                            rewrite, ipktyp, nbits, iret)
            END IF ! iceout
          END IF ! mstout
        END DO
      END IF ! Pressure level output
    END IF   ! Good return from read
  END DO
  STOP
  900 CONTINUE
  WRITE(6,'(a)') ' Error reading input data'
  950 CONTINUE
  WRITE(6,'(a)') ' Error setting up GEMPAK file'
  STOP
END PROGRAM arps2gem

SUBROUTINE mvtmout(nx,ny,nxout,nyout,arrin,arrout)
  IMPLICIT NONE
  INTEGER :: nx,ny,nxout,nyout
  REAL :: arrin(nx,ny)
  REAL :: arrout(nxout,nyout)
  INTEGER :: i,j
!
  DO j=1,nyout
    DO i=1,nxout
      arrout(i,j)=arrin(i,j)
    END DO
  END DO

  RETURN
END SUBROUTINE mvtmout
!

SUBROUTINE extrph(nx,ny,nz,zps,t,pr,iprgem,hgtgem)
!
!  Extrapolate height by using a std atmos lapse rate
!  below the last physical level.  Above the domain,
!  assume a constant temperature above 300 mb, otherwise
!  use the std atmos lapse rate.
!

  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: zps(nx,ny,nz)
  REAL :: t(nx,ny,nz)
  REAL :: pr(nx,ny,nz)
  INTEGER :: iprgem
  REAL :: hgtgem(nx,ny)
!
  INCLUDE 'phycst.inc'
!
  REAL :: gamma,rddg,const
  PARAMETER ( gamma = 0.0065,    & ! degrees/m  lapse rate
          rddg  = (rd/g),                                               &
              const = (rd*gamma/g) )
!
  INTEGER :: i,j
  REAL :: prgem
!
  prgem=100.*FLOAT(iprgem)
  DO j=1,ny
    DO i=1,nx
      IF(prgem < pr(i,j,nz-1)) THEN
        IF(pr(i,j,nz-1) <= 30000.) THEN
          hgtgem(i,j)=zps(i,j,nz-1) +                                   &
              rddg*t(i,j,nz-1)*ALOG(pr(i,j,nz-1)/prgem)
        ELSE
          hgtgem(i,j)=zps(i,j,nz-1) + (t(i,j,nz-1)/gamma)*              &
                      (1.-(prgem/pr(i,j,nz-1))**const)
        END IF
      ELSE IF(prgem >= pr(i,j,2)) THEN
        hgtgem(i,j)=zps(i,j,2) + (t(i,j,2)/gamma)*                      &
               (1.-(prgem/pr(i,j,2))**const)

!       hgtgem(i,j)=0.0
      END IF
    END DO
  END DO
  RETURN
END SUBROUTINE extrph
!

SUBROUTINE extrpt(nx,ny,nz,t,pr,zps,iprgem,tgem)
!
!  Extrapolate temperature by using a std atmos lapse rate
!  below the last physical level.  Above the domain,
!  assume a constant temperature above 300 mb, otherwise
!  use the std atmos lapse rate.
!

  INTEGER :: nx,ny,nz
  REAL :: t(nx,ny,nz)
  REAL :: pr(nx,ny,nz)
  REAL :: zps(nx,ny,nz)
  INTEGER :: iprgem
  REAL :: tgem(nx,ny)
!
  INCLUDE 'phycst.inc'
!
  REAL :: gamma,const
  PARAMETER ( gamma = 0.0065,    & ! degrees/m  lapse rate
          const = (rd*gamma/g) )
!
  INTEGER :: i,j
  REAL :: prgem
!
  prgem=100.*FLOAT(iprgem)
  DO j=1,ny
    DO i=1,nx
      IF(prgem <= pr(i,j,nz-1)) THEN
        IF(pr(i,j,nz-1) <= 30000.) THEN
          tgem(i,j)=t(i,j,nz-1)
        ELSE
          tgem(i,j)=t(i,j,nz-1)*                                        &
                    ((prgem/pr(i,j,nz-1))**const)
        END IF
      ELSE IF(prgem >= pr(i,j,2)) THEN
        tgem(i,j)=t(i,j,2)*                                             &
                    ((prgem/pr(i,j,2))**const)
!
!  missing data flag
!       tgem(i,j)=-99.0
      END IF
    END DO
  END DO
  RETURN
END SUBROUTINE extrpt
!

SUBROUTINE extrpq(nx,ny,nz,q,pr,iprgem,qgem)
!
!  assign 0.0 missing value for below ground data
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: q(nx,ny,nz)
  REAL :: pr(nx,ny,nz)
  INTEGER :: iprgem
  REAL :: qgem(nx,ny)
!
  INTEGER :: i,j
  REAL :: prgem
!
  prgem=100.*FLOAT(iprgem)
  DO j=1,ny
    DO i=1,nx
      IF(prgem <= pr(i,j,nz-1)) THEN
        qgem(i,j)=q(i,j,nz-1)
      ELSE IF(prgem >= pr(i,j,2)) THEN
        qgem(i,j)=0.0
      END IF
    END DO
  END DO
  RETURN
END SUBROUTINE extrpq

SUBROUTINE extrpuv(nx,ny,nz,us,vs,pr,zps,                               &
           iprgem,ugem,vgem)
!
!  assign a "0" value for wind components for underground values.
!
  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: us(nx,ny,nz)
  REAL :: vs(nx,ny,nz)
  REAL :: pr(nx,ny,nz)
  REAL :: zps(nx,ny,nz)
  INTEGER :: iprgem
  REAL :: ugem(nx,ny)
  REAL :: vgem(nx,ny)
!
  INTEGER :: i,j
  REAL :: prgem
!
  prgem=100.*FLOAT(iprgem)
  DO j=1,ny
    DO i=1,nx
      IF(prgem <= pr(i,j,nz-1)) THEN
        ugem(i,j)=us(i,j,nz-1)
        vgem(i,j)=vs(i,j,nz-1)
      ELSE IF(prgem >= pr(i,j,2)) THEN
!       ugem(i,j)=us(i,j,2)
        ugem(i,j)=0.0
!       vgem(i,j)=vs(i,j,2)
        vgem(i,j)=0.0
      END IF
    END DO
  END DO
  RETURN
END SUBROUTINE extrpuv
