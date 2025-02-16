     PROGRAM rademul
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads an ARPS history file and creates data emulating that observed
!  by a radar.  Output is in polar (az-ran) coordinates.
!
!  The initial version uses the "nearest neighbor" in the horizontal and
!  linear vertical interpolation without Gaussian beam smoothing, a
!  feature to be added later.
!
!  Output are NetCDF files which are readable by the NSSL WDSS-II system.
!
!  Radar parameters (tilt angles, beamwidth, etc) are set in the input file.
!
!  AUTHOR:
!  Keith Brewster, Center for Analysis and Prediction of Storms
!  02/10/2004
!
!  MODIFICATIONS
!  10/30/2004 Version with beam-weight filtering.  KB, CAPS
!
!  11/15/2004 Added processing for effective horizontal beam width as a
!             function of rotation rate, PRT and number of samples.
!             Also added an option for X-band attenuation.
!             Keith Brewster, CAPS
!
!  12/05/2004 Added Gaussian random noise for velocity and reflectivity
!             Sometime this could be updated to parameterize the noise
!             Gaussian width as a function of SNR.
!             Keith Brewster, CAPS
!
!  02/06/2005 Added processing for reading individual variables from
!             ARPS hdf history files in order to save memory when reading
!             very large datasets.
!
!  03/18/2005 Added loop to process mulitple files (model output times)
!             in a single run.
!
!  05/10/2005  Keith Brewster, CAPS
!             Added vertical vorticity variable for
!             verification of vortex detection algorithms.
!
!  05/24/2005  Keith Brewster, CAPS
!             Added unattenuated variables and variable write flags.
!             Added smoothing option for vorticity variable.
!
!  06/14/2005 Keith Brewster, CAPS
!             Modified logic to compute weighted average of data re-
!             interpolated to points in a regular 7x7x7 az-ran-elv array.
!
!  12/14/2005 Keith Brewster, CAPS
!             Various updates to make computing on array more flexible, add
!             dual-pulse mode, add forced consistency between PRT, number
!             of samples, rotation rate and output radial spacing.
!
!  04/20/2006 Keith Brewster, CAPS
!             Added beam blockage.
!
!  07/21/2008 Keith Brewster, CAPS
!             Added T-matrix calculations.  Added RHI option.
!
!-----------------------------------------------------------------------
!
  USE DUALPARA
  IMPLICIT NONE

  INCLUDE 'grid.inc'
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'

  INTEGER :: nx, ny, nz
  INTEGER :: nzsoil
  INTEGER :: nstyps
!
!-----------------------------------------------------------------------
!
!  Arrays to be read in:
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: x(:)       ! The x-coord. of the physical and
                                  ! computational grid.
  REAL, ALLOCATABLE :: y(:)       ! The y-coord. of the physical and
                                  ! computational grid.
  REAL, ALLOCATABLE :: z(:)       ! The z-coord. of the computational
                                  ! grid. Defined at w-point.
  REAL, ALLOCATABLE :: zp(:,:,:)  ! The height of the terrain.
  REAL, ALLOCATABLE :: zpsoil(:,:,:)  ! The height for soil model

  REAL, ALLOCATABLE :: uprt   (:,:,:) ! Perturbation u-velocity (m/s)
  REAL, ALLOCATABLE :: vprt   (:,:,:) ! Perturbation v-velocity (m/s)
  REAL, ALLOCATABLE :: wprt   (:,:,:) ! Perturbation w-velocity (m/s)
  REAL, ALLOCATABLE :: ptprt  (:,:,:) ! Perturbation potential temperature (K)
  REAL, ALLOCATABLE :: pprt   (:,:,:) ! Perturbation pressure (Pascal)
  REAL, ALLOCATABLE :: qvprt  (:,:,:) ! Perturbation water vapor specific
                                      ! humidity (kg/kg)
  REAL, ALLOCATABLE, TARGET :: qscalar(:,:,:,:)

  REAL, POINTER :: qr(:,:,:)      ! for use within this program for HDF format  only
  REAL, POINTER :: qs(:,:,:)
  REAL, POINTER :: qh(:,:,:)

  REAL, ALLOCATABLE :: ubar   (:,:,:) ! Base state u-velocity (m/s)
  REAL, ALLOCATABLE :: vbar   (:,:,:) ! Base state v-velocity (m/s)
  REAL, ALLOCATABLE :: wbar   (:,:,:) ! Base state w-velocity (m/s)
  REAL, ALLOCATABLE :: ptbar  (:,:,:) ! Base state potential temperature (K)
  REAL, ALLOCATABLE :: pbar   (:,:,:) ! Base state pressure (Pascal)
  REAL, ALLOCATABLE :: rhobar (:,:,:) ! Base state air density (kg/m**3)
  REAL, ALLOCATABLE :: qvbar  (:,:,:) ! Base state water vapor specific
                                      ! humidity (kg/kg)

  INTEGER, ALLOCATABLE :: soiltyp (:,:,:) ! Soil type
  REAL, ALLOCATABLE ::  stypfrct(:,:,:)   ! Fraction of soil types
  INTEGER, ALLOCATABLE :: vegtyp(:,:)     ! Vegetation type

  REAL, ALLOCATABLE ::   tsoil  (:,:,:,:)   ! Deep soil temperature (K)
  REAL, ALLOCATABLE ::   qsoil  (:,:,:,:)   ! Deep soil temperature (K)
  REAL, ALLOCATABLE ::   wetcanp(:,:,:)   ! Canopy water amount

  REAL, ALLOCATABLE ::prcrate(:,:,:)      ! precipitation rate (kg/(m**2*s))
                                          ! prcrate(1,1,1) = total precip. rate
                                          ! prcrate(1,1,2) = grid scale precip. rate
                                          ! prcrate(1,1,3) = cumulus precip. rate
                                          ! prcrate(1,1,4) = microphysics precip. rate

  REAL, ALLOCATABLE ::tem2d(:,:)         ! 2D work array
  REAL, ALLOCATABLE ::t(:,:,:)           ! temperature (K)
  REAL, ALLOCATABLE ::tem1(:,:,:)        ! Work arrays
  REAL, ALLOCATABLE ::tem2(:,:,:)        ! Work arrays
  REAL, ALLOCATABLE ::tem3(:,:,:)        ! Work arrays
  REAL, ALLOCATABLE ::tem4(:,:,:)        ! Work arrays

  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE :: itmp(:,:,:) ! Temporary array
  REAL, ALLOCATABLE :: hmax(:)
  REAL, ALLOCATABLE :: hmin(:)

  REAL, ALLOCATABLE :: xsc(:)
  REAL, ALLOCATABLE :: ysc(:)
  REAL, ALLOCATABLE :: latsc(:,:)
  REAL, ALLOCATABLE :: lonsc(:,:)
  REAL, ALLOCATABLE :: azmsc(:,:)
  REAL, ALLOCATABLE :: sfcr(:,:)
  REAL, ALLOCATABLE :: cmpref(:,:)

  REAL, ALLOCATABLE :: usc(:,:,:)
  REAL, ALLOCATABLE :: vsc(:,:,:)
  REAL, ALLOCATABLE :: usm(:,:,:)
  REAL, ALLOCATABLE :: vsm(:,:,:)
  REAL, ALLOCATABLE :: wsc(:,:,:)
  REAL, ALLOCATABLE :: zps(:,:,:)
  REAL, ALLOCATABLE :: ref(:,:,:)
  REAL, ALLOCATABLE :: refz(:,:,:)
  REAL, ALLOCATABLE :: refh(:,:,:)
  REAL, ALLOCATABLE :: refv(:,:,:)
  REAL, ALLOCATABLE :: refhv(:,:,:)
  REAL, ALLOCATABLE :: zdr(:,:,:)
  REAL, ALLOCATABLE :: rhohv(:,:,:)
  REAL, ALLOCATABLE :: kdp(:,:,:)
  REAL, ALLOCATABLE :: elvsc(:,:,:)
  REAL, ALLOCATABLE :: rngsc(:,:,:)
  REAL, ALLOCATABLE :: vort(:,:,:)
  REAL, ALLOCATABLE :: vrsc(:,:,:)
!
!-----------------------------------------------------------------------
!
!  NEXRAD Radar variables
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: maxtilt=14
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
!  Grid adjustment variables -- adjustments to model grid specifications
!
!-----------------------------------------------------------------------
!
  INTEGER :: adjhgt,adjctr,adjmove
  REAL :: ctrlatem,ctrlonem
  REAL :: umovein,vmovein
  REAL :: hgtoffset
  NAMELIST /grid_adj/ adjhgt,adjctr,adjmove,ctrlatem,ctrlonem,          &
                      umovein,vmovein,hgtoffset
!
!-----------------------------------------------------------------------
!
!  I/O Options
!
!-----------------------------------------------------------------------
!
  INTEGER :: ifmt,creidx,ipktyp,nbits,nsmvort
  INTEGER :: wrtuaref,wrtuavel,wrtvort,wrtdualp,wrtqx,wrtuvwt
  CHARACTER (LEN=80) :: outdir
  NAMELIST /output_opts/ ifmt,creidx,ipktyp,nbits,outdir,              &
            nsmvort,wrtuaref,wrtuavel,wrtvort,wrtdualp,wrtqx,wrtuvwt
!-----------------------------------------------------------------------
!
!  Radar variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=4) :: radname
  CHARACTER (LEN=256) :: mobilefile
  INTEGER :: locopt,mobileopt
  INTEGER :: dualpol
  REAL :: xrad,yrad,radlat,radlon,radelv
  REAL :: wavelen,beamwid
  REAL :: pwrxmt,gaindb,lossdb,noisedbm
  NAMELIST /radar_specs/ radname,locopt,mobileopt,mobilefile,            &
           xrad,yrad,radlat,radlon,radelv,dualpol,wavelen,beamwid,       &
           pwrxmt,gaindb,lossdb,noisedbm
!
!-----------------------------------------------------------------------
!
!  Radar Operation Mode Variables
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: maxelv=200
  INTEGER :: ntilt,vcp,ngate,rhimode
  INTEGER :: unifpuls,dualprf,pulseopt,rotropt,rotdir
  REAL :: gatesp,rotrate
  REAL :: delazm,azim1,azim2
  REAL :: delelv,elev1,elev2
  REAL :: elev(maxelv)
  REAL :: pulselen1(maxelv),pulselen2(maxelv)
  REAL :: pulsewid1(maxelv),pulsewid2(maxelv)
  REAL :: prt1(maxelv),prt2(maxelv)
  INTEGER :: npulse1(maxelv),npulse2(maxelv)

  NAMELIST /radar_opmode/ ntilt,vcp,elev,ngate,gatesp,                 &
            rhimode,rotdir,rotropt,rotrate,                            &
            delazm,azim1,azim2,delelv,elev1,elev2,                     &
            unifpuls,dualprf,pulseopt,                                 &
            pulselen1,pulselen2,pulsewid1,pulsewid2,                   &
            prt1,prt2,npulse1,npulse2
!
!-----------------------------------------------------------------------
!
!  Radar emulator options
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: refchk = -99.
  REAL, PARAMETER :: refmax = 80.
  REAL, PARAMETER :: epsilon = 1.0E-10
  INTEGER :: nptsazm,nptselv,nptsrng
  REAL :: evalwid,refmin,rngmin
  REAL :: rmisval,rngfval
  REAL :: hblkmin,hblkmax,gndrefl,gndkdp
  INTEGER :: int_method,blockopt,kntrmin,kntvmin
  INTEGER :: dualpopt,attenopt,senstvopt
  INTEGER :: rferropt,vrerropt,nyqstopt
  INTEGER :: tmadvopt
  REAL :: sigmarf,sigmavr
  REAL :: samploi,stdvrmul
  REAL :: vnyquist(maxelv)
  REAL :: timeincr
  CHARACTER (LEN=256) :: tmatrix_dir
  NAMELIST /emul_specs/ int_method,nptsazm,nptselv,nptsrng,            &
           evalwid,refmin,rngmin,kntrmin,kntvmin,                      &
           blockopt,hblkmin,hblkmax,gndrefl,gndkdp,                    &
           tmadvopt,dualpopt,tmatrix_dir,attenopt,senstvopt,           &
           rferropt,sigmarf,vrerropt,sigmavr,samploi,stdvrmul,         &
           timeincr,rmisval,rngfval,nyqstopt,vnyquist

!-----------------------------------------------------------------------
!
!  Mobile radar location variables
!
!-----------------------------------------------------------------------
!
  CHARACTER(LEN=80) :: aline
  INTEGER, PARAMETER :: maxmob=200
  INTEGER, ALLOCATABLE :: mobtime(:)
  REAL, ALLOCATABLE :: mobilex(:)
  REAL, ALLOCATABLE :: mobiley(:)
  REAL, ALLOCATABLE :: moblat(:)
  REAL, ALLOCATABLE :: moblon(:)
  REAL, ALLOCATABLE :: mobhgt(:)
  REAL, ALLOCATABLE :: mobhead(:)
  REAL, ALLOCATABLE :: mobpitch(:)

  INTEGER :: imobf,jtime,ntime
  INTEGER :: myear,mmon,mday,mhour,mmin,msec
  REAL :: rpitch,rhead,tw0,tw1

!-----------------------------------------------------------------------
!
!  Radar tilt variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: maxazim
  INTEGER :: itilt,jazim,igate,iigate
  INTEGER :: iloc,jloc,iloc1,jloc1,ilc,jlc
  REAL :: deg2rad,rad2deg
  REAL :: irloc,jrloc
  REAL :: elvang(maxelv)
  REAL, ALLOCATABLE :: azim(:)
  REAL, ALLOCATABLE :: beamw(:)
  REAL, ALLOCATABLE :: gtspc(:)
  REAL, ALLOCATABLE :: vnyq(:)
  REAL, ALLOCATABLE :: radv(:,:)
  REAL, ALLOCATABLE :: attradv(:,:)
  REAL, ALLOCATABLE :: refl(:,:)
  REAL, ALLOCATABLE :: attrefl(:,:)
  REAL, ALLOCATABLE :: rzdr(:,:)
  REAL, ALLOCATABLE :: rrhohv(:,:)
  REAL, ALLOCATABLE :: rkdp(:,:)
  REAL, ALLOCATABLE :: stdvr(:,:)
  REAL, ALLOCATABLE :: vvort(:,:)
  REAL, ALLOCATABLE :: qrrad(:,:)
  REAL, ALLOCATABLE :: qsrad(:,:)
  REAL, ALLOCATABLE :: qhrad(:,:)
  REAL, ALLOCATABLE :: urad(:,:)
  REAL, ALLOCATABLE :: vrad(:,:)
  REAL, ALLOCATABLE :: wrad(:,:)
  REAL, ALLOCATABLE :: tkrad(:,:)
  REAL, ALLOCATABLE :: sens(:)
!
!-----------------------------------------------------------------------
!
!  Radar volume variables
!
!-----------------------------------------------------------------------
!
  INTEGER, ALLOCATABLE :: itimvol(:)
  REAL, ALLOCATABLE :: vnyqvol(:)
  REAL, ALLOCATABLE :: rngvol(:,:)
  REAL, ALLOCATABLE :: azmvol(:,:)
  REAL, ALLOCATABLE :: elvvol(:,:)
  REAL, ALLOCATABLE :: refvol(:,:,:)
  REAL, ALLOCATABLE :: uarefvol(:,:,:)
  REAL, ALLOCATABLE :: zdrvol(:,:,:)
  REAL, ALLOCATABLE :: kdpvol(:,:,:)
  REAL, ALLOCATABLE :: rhohvvol(:,:,:)
  REAL, ALLOCATABLE :: velvol(:,:,:)
  REAL, ALLOCATABLE :: uavelvol(:,:,:)
  REAL, ALLOCATABLE :: vorvol(:,:,:)
  REAL, ALLOCATABLE :: qrvol(:,:,:)
  REAL, ALLOCATABLE :: qsvol(:,:,:)
  REAL, ALLOCATABLE :: qhvol(:,:,:)
  REAL, ALLOCATABLE :: uvol(:,:,:)
  REAL, ALLOCATABLE :: vvol(:,:,:)
  REAL, ALLOCATABLE :: wvol(:,:,:)
  REAL, ALLOCATABLE :: tkvol(:,:,:)
!
!-----------------------------------------------------------------------
!
!  Misc. internal variables
!
!-----------------------------------------------------------------------
!
  REAL :: alatnot(2)
  REAL, PARAMETER :: c = 2.9979e08    ! speed of light m/s
  REAL, ALLOCATABLE :: delrng(:)
  REAL, ALLOCATABLE :: dazim(:)
  REAL, ALLOCATABLE :: delev(:)
  REAL, ALLOCATABLE :: wgtpt(:,:,:)
  REAL, ALLOCATABLE :: blockfct(:,:)
  INTEGER, PARAMETER :: nhisfile_max=1000
  INTEGER, PARAMETER :: dbz_out = 0
  CHARACTER (LEN=256) :: fname
  CHARACTER (LEN=256) :: fnamevol
  CHARACTER (LEN=256) :: idxfname
  CHARACTER (LEN=256) :: grdbasfn
  CHARACTER (LEN=256) :: hisfile(nhisfile_max)
  CHARACTER (LEN=80) :: vcptxt
  CHARACTER (LEN=40) :: varname
  CHARACTER (LEN=5) :: cmprext

  TYPE(T_obs_dual) :: obs_dual

  INTEGER, PARAMETER :: itim1970=315619200
  INTEGER, PARAMETER :: drysnowopt=4

  INTEGER :: hinfmt,nchin,nhisfile,idxunit,istatus
  INTEGER :: lengbf,lenfil,itime,itimcdf,npulsmn
  INTEGER :: itimtilt,iyear,imon,iday,ihour,imin,isec
  INTEGER :: ifile,initime,ifsecs,ireturn
  INTEGER :: i,j,k,ibgn,iend,jbgn,jend,ipt,jpt,kpt,iptmid
  INTEGER :: ipmid,jpmid,kpmid,nzm2
  INTEGER :: nelev,nazim,ielv,felv,ifold,irot,krot,ismth
  INTEGER :: iradius,jradius,kntr,kntv
  INTEGER :: kntgate,kntintr,kntintv,kntfewr,kntfewv
  INTEGER :: kntavgr,kntavgv
  INTEGER :: kntclr,knthob,kntvob,kntbob
  REAL :: twdx2inv,twdxinv,twdy2inv,twdyinv
  REAL :: dr,de,da,p0inv,azmsect,azmbgn,denom,ttime
  INTEGER :: itltbgn,itltend
  REAL :: time,time0,rfrgate,frtime,radlata,radlona,hlfgtsp
  REAL :: elv,azimuth,azrad,srange,hgt,hgtmsl,sfcrng
  REAL :: ctrx,ctry,xsw,ysw,gtlat,gtlon,gtx,gty,xpt,ypt
  REAL :: delx,dely,xrat,yrat,dhdr,dsdr,twovnyq,twovninv
  REAL :: ucmp,vcmp,upt,vpt,wpt,tkpt,vt,vr,vvmin,vvmax
  REAL :: xrada,yrada,ipos,jpos
  REAL :: w1i,w2i,w1j,w2j,w11,w12,w21,w22
  REAL :: r11h,r12h,r21h,r22h,rmax
  REAL :: r11l,r12l,r21l,r22l
  REAL :: refhpt,refvpt,refhvpt,kdppt
  REAL :: delv,dazm,drng,depth,elradius
  REAL :: whigh,wlow,a,b,fjp1,fj0,fjm1
  REAL :: flow,fhigh,reflz,refdbz,reflh,reflv,reflhv,reflzgnd
  REAL :: cinv,b6,hlftau,hlfctau,aconst,rngg
  REAL :: efbmwid,bmwrad,bmwm,sradius,vnyqstl
  REAL :: fourln4,sigr2,sr2inv,wesqinv,wasqinv
  REAL :: wgtg,wgtr,wgtv
  REAL :: sumwr,sumwv,sumv,sumr,sumrh,sumrv,sumrhv,sumkdp
  REAL :: sumv1,sumv2,flkntv,varvr
  REAL :: wgti,wgtj,wgtk
  REAL :: prtmn,pulsewmn,pulselmn
  REAL :: rngwid,srgpt,azmpt,azmrpt,elvpt,sfcrpt,sfcrngt,sfcrngb
  REAL :: elvtop,elvbot,trnpt,trngt,topgt,tophgt,topmsl,bothgt,botmsl
  REAL :: hgtpt,hgtagl,xmit,dhblkinv,blkcst,gndatten
  REAL :: qrmin,qsmin,qhmin,qrmax,qsmax,qhmax
  REAL :: reflmin,reflmax,tstmin,tstmax,err
  LOGICAL :: echo
!
!-----------------------------------------------------------------------
!
! External Functions
!
!-----------------------------------------------------------------------
!
  REAL :: randnor
  REAL :: effbmw
  REAL :: calculate_kdp
  REAL :: erf
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
!  Initializations of constants and namelist variables
!
!-----------------------------------------------------------------------
!
  itime=0
  frtime=0.
  time0=0.
  deg2rad=pi/180.
  rad2deg=1./deg2rad
  cinv=1./c
  p0inv=1./p0

  radname='NULL'
  mobileopt=0
  mobilefile='radlocs.txt'
  locopt=2
  xrad=0.
  yrad=0.
  radlat=0.
  radlon=0.
  radelv=0.
  rpitch=0.
  rhead=0.
  wavelen=10.5
  beamwid=0.89
  pwrxmt=475000.
  gaindb=44.5
  lossdb=5.
  noisedbm=-113.
  dualpol=0
  tmatrix_dir='./data/radaremul'

  int_method=2
  nptsazm=7
  nptselv=7
  nptsrng=7
  evalwid=2.0
  refmin=-20.
  rngmin=100.
  kntrmin=9
  kntvmin=9
  blockopt=2
  hblkmin=1.0
  hblkmax=10.0
  gndrefl=50.0
  gndkdp=2.0
  tmadvopt=1
  dualpopt=1
  attenopt=0
  senstvopt=0
  rferropt=1
  sigmarf=5.0
  vrerropt=1
  sigmavr=1.0
  samploi=0.5
  stdvrmul=1.5
  timeincr=1.
  rmisval=-99900.
  rngfval=-99901.
  nyqstopt=1
  vnyquist=100.

  nelev=1
  vcp=0
  elev=0.5
  ngate=100
  gatesp=1000.
  rhimode=0
  rotdir=1
  rotropt=1
  rotrate=18.0
  delelv=1.0
  elev1=-20.
  elev2=20.
  delazm=1.0
  azim1=0.
  azim2=360.
  unifpuls=1
  dualprf=0
  pulseopt=1
  pulselen1=1.57E-06
  pulselen2=1.57E-06
  pulsewid1=441.
  pulsewid2=441.
  prt1=1.06e-03
  prt2=1.06e-03
  npulse1=54
  npulse2=54

  ifmt=1
  creidx=0
  ipktyp=2
  nbits=16
  nsmvort=2
  wrtuaref=0
  wrtuavel=0
  wrtvort=0
  wrtdualp=0
  wrtqx=0
  wrtuvwt=0
  outdir='./'

  adjhgt=0
  adjctr=0
  adjmove=0
  ctrlatem=45.
  ctrlonem=-90.
  umovein=0.
  vmovein=0.
  hgtoffset=0.

  irloc=45.
  jrloc=45.
  cmprext='.gz'

!
!-----------------------------------------------------------------------
!
!  Get the names of the input data files.
!
!-----------------------------------------------------------------------
!
  CALL get_input_file_names(5,hinfmt,grdbasfn,hisfile,nhisfile)
  WRITE(6,'(a,i4,a/)') ' Processing ',nhisfile,' model files.'

  lengbf = len_trim(grdbasfn)
!
!-----------------------------------------------------------------------
!
!  Read other input namelists
!
!-----------------------------------------------------------------------
!
  WRITE(6,'(a)') ' Reading grid_adj namelist'
  READ(5,grid_adj)
  WRITE(6,'(a)') ' Reading output_opts namelist'
  READ(5,output_opts)
  WRITE(6,'(a)') ' Reading radar_specs namelist'
  READ(5,radar_specs)
  WRITE(6,'(a)') ' Reading radar_opmode namelist'
  READ(5,radar_opmode)
  WRITE(6,'(a)') ' Reading emul_specs namelist'
  READ(5,emul_specs)

  IF(rhimode < 1) THEN
    IF(vcp == 0) THEN
      nelev=ntilt
      write(vcptxt,'(a,i4,a,f4.1,a,f4.1,a)') &
        'Custom vcp',nelev,' tilts',elev(1),'-',elev(nelev),' deg'
      DO itilt=1,nelev
        elvang(itilt)=elev(itilt)
      END DO
    ELSE IF(vcp == 11) THEN
      nelev=14
      vcptxt='Storm mode  14 tilts 0.5-19.5 deg'
      DO itilt=1,nelev
        elvang(itilt)=ang11(itilt)
      END DO
    ELSE IF (vcp == 12) THEN
      nelev=9
      vcptxt='Storm mode   9 tilts 0.5-19.5 deg'
      DO itilt=1,nelev
        elvang(itilt)=ang11(itilt)
      END DO
    ELSE IF (vcp == 31) THEN
      nelev=5
      vcptxt='Clear-air    5 tilts 0.5- 4.5 deg'
      DO itilt=1,nelev
        elvang(itilt)=ang11(itilt)
      END DO
    ELSE IF (vcp == 32) THEN
      nelev=5
      vcptxt='Clear-air    5 tilts 0.5- 4.5 deg'
      DO itilt=1,nelev
        elvang(itilt)=ang11(itilt)
      END DO
    ELSE
      WRITE(6,*) vcp,' is not a recognized VCP number'
      STOP
    END IF
  ELSE  ! RHI mode
    vcp=-99
    nelev=INT(abs((elev2-elev1)/delelv))+1
    nelev=min(nelev,maxelv)
    WRITE(6,'(a,i4)') ' RHI mode, nelev=',nelev
    DO itilt=1,nelev
      elvang(itilt)=elev1+(itilt-1)*delelv
      WRITE(6,'(a,i4,a,f10.2)') ' elev(',itilt,') = ',elvang(itilt)
    END DO
  END IF

  IF(mobileopt > 0) THEN
    imobf=31
    OPEN(imobf,file=mobilefile,status='old',iostat=istatus)
    IF(istatus == 0) THEN
      ALLOCATE(mobtime(maxmob))
      ALLOCATE(mobhgt(maxmob))
      ALLOCATE(mobhead(maxmob))
      ALLOCATE(mobpitch(maxmob))
      ALLOCATE(mobilex(maxmob))
      ALLOCATE(mobiley(maxmob))
      ALLOCATE(moblat(maxmob))
      ALLOCATE(moblon(maxmob))
      IF(locopt == 1) THEN
        jtime=0
        DO
          READ(imobf,'(a80)',iostat=istatus) aline
          IF(istatus /= 0) EXIT
          IF(aline(1:1) == '#' .OR. aline(1:1) == '!') CYCLE
          jtime=jtime+1
          READ(aline,'(1x,i4,5(1x,i2),f11.4,f11.4,f9.1,f7.1,f6.1)', &
              iostat=istatus) &
              myear,mmon,mday,mhour,mmin,msec,moblat(jtime),moblon(jtime), &
              mobhgt(jtime),mobhead(jtime),mobpitch(jtime)
          print *, ' mobile lat, lon: ',moblat(jtime),moblon(jtime)
          CALL ctim2abss(myear,mmon,mday,mhour,mmin,msec,mobtime(jtime))
          IF(istatus /= 0) EXIT
        END DO
        ntime=jtime
        DO jtime=1,ntime-1
          IF(mobhead(jtime) < 0.) THEN
            CALL disthead(moblat(jtime),moblon(jtime), &
                          moblat(jtime+1),moblon(jtime+1),  &
                          mobhead(jtime),sfcrng)
            print *, ' heading(',jtime,') =',mobhead(jtime)
          END IF
        END DO
        IF(mobhead(ntime) < 0.) mobhead(ntime)=mobhead(ntime-1)
      ELSE
        jtime=0
        DO
          READ(imobf,'(a80)',iostat=istatus) aline
          IF(istatus /= 0) EXIT
          IF(aline(1:1) == '#' .OR. aline(1:1) == '!') CYCLE
          jtime=jtime+1
          READ(aline,'(1x,i4,5(1x,i2),f11.4,f11.4,f9.1,f7.1,f6.1)', &
              iostat=istatus) &
              myear,mmon,mday,mhour,mmin,msec,mobilex(jtime),mobiley(jtime), &
              mobhgt(jtime),mobhead(jtime),mobpitch(jtime)
          print *, ' mobile x, y: ',mobilex(jtime),mobiley(jtime)
          CALL ctim2abss(myear,mmon,mday,mhour,mmin,msec,mobtime(jtime))
          IF(istatus /= 0) EXIT
        END DO
        CALL xytoll(ntime,1,mobilex,mobiley,moblat,moblon)
        DO jtime=1,ntime-1
          IF(mobhead(jtime) < 0.) THEN
            CALL disthead(moblat(jtime),moblon(jtime), &
                          moblat(jtime+1),moblon(jtime+1),  &
                          mobhead(jtime),sfcrng)
            print *, ' heading(',jtime,') =',mobhead(jtime)
          END IF
        END DO
        IF(mobhead(ntime) < 0.) mobhead(ntime)=mobhead(ntime-1)
        ntime=jtime
      END IF
      WRITE(6,'(a,i8,a)') ' Successfully read ',jtime,' mobile locations'
      CLOSE(imobf)
    ELSE
      WRITE(6,'(a)') ' Error opening mobile locations file: ', &
        TRIM(mobilefile)
      STOP
    END IF
  END IF
!
!-----------------------------------------------------------------------
!
! Print a warning about blocking/clutter using integration method 1.
!
!-----------------------------------------------------------------------
!
  IF(int_method == 1 .AND. blockopt > 0) THEN
    WRITE(6,'(//a,/a,i4,/a,i4,/a//)') &
   ' ***********WARNING*********************', &
   ' Terrain blocking calculation not possible using integration method ', &
    int_method, ' Select integration method 2 with blockopt=',blockopt,&
    'Continuing using integration method 2...'
   int_method=2
  END IF
!
!-----------------------------------------------------------------------
!
! Calculate constants for antenna gain weighting.
!
!-----------------------------------------------------------------------
!
  kntvmin=max(kntvmin,2)
  fourln4=4.0*alog(4.0)
  wesqinv=fourln4/(beamwid*beamwid)
  reflzgnd=10.**(gndrefl*0.10)

  IF(azim2 > azim1) THEN
    azmsect=azim2-azim1
  ELSE
    azmsect=(360.+azim2)-azim1
  END IF

  IF(unifpuls == 1) THEN
    DO itilt=2,nelev
      prt1(itilt)=prt1(1)
      prt2(itilt)=prt2(1)
      npulse1(itilt)=npulse1(1)
      npulse2(itilt)=npulse2(1)
      pulselen1(itilt)=pulselen1(1)
      pulselen2(itilt)=pulselen2(1)
      pulsewid1(itilt)=pulsewid1(1)
      pulsewid2(itilt)=pulsewid2(1)
    END DO
  END IF

  IF(rhimode < 1) THEN
    IF(rotropt > 1) THEN
      nazim=NINT(azmsect/delazm)+1
      maxazim=max(nazim,1)
    ELSE
      maxazim=1
      DO itilt=1,nelev
        IF(dualprf > 0) THEN
          delazm=rotrate*((prt1(itilt)*npulse1(itilt))+                  &
                          (prt2(itilt)*npulse2(itilt)))
        ELSE
          delazm=rotrate*prt1(itilt)*npulse1(itilt)
        END IF
        nazim=NINT(azmsect/delazm)+1
        maxazim=max(nazim,maxazim)
      END DO
    END IF
  ELSE
    nazim=1
    maxazim=1
  END IF
  WRITE(6,'(a,i6)') ' Maximum azimuths: ',maxazim

  rfrgate=gatesp
  hlfgtsp=0.5*gatesp

  IF(blockopt == 0) THEN
    hblkmin=0.
    hblkmax=1.0
    dhblkinv=1.0
  ELSE
    hblkmin=max(hblkmin,0.0)
    hblkmax=max(hblkmax,(hblkmin+1.0))
    dhblkinv=1./(hblkmax-hblkmin)
  END IF

  nptsrng=max(nptsrng,5)
  ALLOCATE (delrng(nptsrng))
  dr=gatesp/float(nptsrng-1)
  WRITE(6,'(a,f10.2,a)') ' Range evaluation points (gatesp=',gatesp,'):'
  DO ipt=1,nptsrng
    delrng(ipt)=(-0.5*gatesp)+((ipt-1)*dr)
    WRITE(6,'(a,i4,a,f12.4)') ' delrng(',ipt,') = ',delrng(ipt)
  END DO
  iptmid=(nptsrng/2)+1
  blkcst=(0.001*gatesp)/float(nptsrng)

  nptsazm=max(nptsazm,5)
  ALLOCATE (dazim(nptsazm))

  nptselv=max(nptselv,5)
  ALLOCATE (delev(nptselv))
  elradius=evalwid*beamwid
  de=elradius/float(nptselv-1)
  WRITE(6,'(a,f10.2,a,f10.2)') &
      ' Beam width ',beamwid,' Eval region width:',(evalwid*beamwid)
  WRITE(6,'(a)') &
      ' Vertical cross-beam evaluation pts (elev deg above center):'
  DO kpt=1,nptselv
    delev(kpt)=(-0.5*evalwid*beamwid)+((kpt-1)*de)
    WRITE(6,'(a,i4,a,f12.4)') ' delev(',kpt,') = ',delev(kpt)
  END DO

  ipmid=(nptsrng/2)+1
  jpmid=(nptsazm/2)+1
  kpmid=(nptselv/2)+1

  ALLOCATE (wgtpt(nptsrng,nptsazm,nptselv))
  ALLOCATE (blockfct(nptsazm,nptselv))

!
!-----------------------------------------------------------------------
!
!  Obtain the grid dimensions from input data.
!
!-----------------------------------------------------------------------
!
  CALL get_dims_from_data(hinfmt,hisfile(1),                    &
       nx,ny,nz,nzsoil,nstyps,ireturn)

  nstyp = nstyps ! Copy to global variable

  IF( ireturn /= 0 ) THEN
    PRINT*,'Problem occured when trying to get dimensions from data.'
    PRINT*,'Program stopped.'
    STOP
  END IF

  WRITE(6,'(3(a,i5))') 'nx =',nx,', ny=',ny,', nz=',nz,', nzsoil',nzsoil
  nzm2=nz-2
!-----------------------------------------------------------------------
!
!  Allocate grid position arrays.
!
!-----------------------------------------------------------------------
!
  allocate(x(nx),stat=istatus)
  CALL check_alloc_status(istatus, "radaremul:x")
  x=0.
  allocate(y(ny),stat=istatus)
  CALL check_alloc_status(istatus, "radaremul:y")
  y=0.
  allocate(z(nz),stat=istatus)
  CALL check_alloc_status(istatus, "radaremul:z")
  z=0.
  allocate(zp(nx,ny,nz),stat=istatus)
  CALL check_alloc_status(istatus, "radaremul:zp")
  zp=0.

  IF(rhimode > 0) THEN
    krot=-1
  ELSE
    krot=1
  END IF

  DO ifile=1,nhisfile

    WRITE(6,'(a,a)') '  Processing file: ',TRIM(hisfile(ifile))

    IF( hinfmt .NE. 3) THEN

!-----------------------------------------------------------------------
!
!  Allocate arrays needed for dtaread
!  These are deallocated after being read and used so they need
!  to be reallocated with each new file.
!
!-----------------------------------------------------------------------
!
      allocate(zpsoil (nx,ny,nzsoil),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:zpsoil")
      zpsoil=0.

      allocate(uprt   (nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:uprt")
      uprt=0.
      allocate(vprt   (nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:vprt")
      vprt=0.
      allocate(wprt   (nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:wprt")
      wprt=0.
      allocate(ptprt  (nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:ptprt")
      ptprt=0.
      allocate(pprt   (nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:pprt")
      pprt=0.
      allocate(qvprt  (nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:qvprt")
      qvprt=0.

      allocate(ubar   (nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:ubar")
      ubar=0.
      allocate(vbar   (nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:vbar")
      vbar=0.
      allocate(wbar   (nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:wbar")
      wbar=0.
      allocate(ptbar  (nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:ptbar")
      ptbar=0.
      allocate(pbar   (nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:pbar")
      pbar=0.
      allocate(rhobar (nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:rhobar")
      rhobar=0.
      allocate(qvbar  (nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:qvbar")
      qvbar=0.

      ALLOCATE(qscalar(nx,ny,nz,nscalar),STAT=istatus)
      CALL check_alloc_status(istatus, "radaremul:qscalar")
      qscalar = 0.

      allocate(soiltyp (nx,ny,nstyps),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:soiltyp")
      soiltyp=0
      allocate(stypfrct(nx,ny,nstyps),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:stypfrct")
      stypfrct=0.
      allocate(vegtyp(nx,ny),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:vegtyp")
      vegtyp=0
      allocate(tsoil  (nx,ny,nzsoil,0:nstyps),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:tsoil")
      tsoil=0.
      allocate(qsoil  (nx,ny,nzsoil,0:nstyps),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:qsoil")
      qsoil=0.
      allocate(wetcanp(nx,ny,0:nstyps),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:wetcanp")
      wetcanp=0.
      allocate(prcrate(nx,ny,4),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:prcrate")
      prcrate=0.

      allocate(tem2d(nx,ny),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:tem2d")
      tem2d=0.
      allocate(tem1(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:tem1")
      tem1=0.
      allocate(tem2(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:tem2")
      tem2=0.
      allocate(tem3(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:tem3")
      tem3=0.
      allocate(tem4(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:tem4")
      tem4=0.
      allocate(t(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:t")
      t=0.

      lenfil = len_trim(hisfile(ifile))
      CALL dtaread(nx,ny,nz,nzsoil,nstyps,                                &
         hinfmt, nchin,grdbasfn(1:lengbf),lengbf,                         &
         hisfile(ifile),lenfil,time,                                      &
         x,y,z,zp,zpsoil,uprt,vprt,wprt,ptprt,pprt,                       &
         qvprt, qscalar, tem4,tem4,tem4,                                &
         ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,                    &
         soiltyp,stypfrct,vegtyp,tem2d,tem2d,tem2d,                       &
         tsoil,qsoil,wetcanp,tem2d,                                       &
         tem2d,tem2d,prcrate,                                             &
         tem4,tem2d,tem2d,tem2d,tem2d,                                    &
         tem2d,tem2d,tem2d,tem2d,                                         &
         ireturn, tem1,tem2,tem3)
!
!-----------------------------------------------------------------------
!
!  Deallocate the no-longer-needed dtaread arrays.
!
!-----------------------------------------------------------------------
!
      deallocate(qvprt)
      deallocate(qvbar)
      deallocate(wbar)
!      deallocate(qscalar)
      deallocate(zpsoil)
      deallocate(soiltyp)
      deallocate(stypfrct)
      deallocate(vegtyp)
      deallocate(tsoil)
      deallocate(qsoil)
      deallocate(wetcanp)
      deallocate(prcrate)
      deallocate(tem2d)
      deallocate(tem2)
      deallocate(tem3)
      deallocate(tem4)

      qr => qscalar(:,:,:,P_QR)
      qs => qscalar(:,:,:,P_QS)
      qh => qscalar(:,:,:,P_QH)


      allocate(usc(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:usc")
      usc=0.
      allocate(vsc(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:vsc")
      vsc=0.
      allocate(wsc(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:wsc")
      wsc=0.

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            usc(i,j,k)=0.5*((uprt(  i,j,k)+ubar(  i,j,k))+                 &
                            (uprt(i+1,j,k)+ubar(i+1,j,k)))+umove
            vsc(i,j,k)=0.5*((vprt(i,  j,k)+vbar(i,  j,k))+                 &
                            (vprt(i,j+1,k)+vbar(i,j+1,k)))+vmove
            wsc(i,j,k)=0.5*(wprt(i,j,k)+wprt(i,j,k+1))
          END DO
        END DO
      END DO

      deallocate(ubar)
      deallocate(uprt)
      deallocate(vbar)
      deallocate(vprt)
      deallocate(wprt)

!
!-----------------------------------------------------------------------
!
!   Calculate reflectivity at all scalar grid points.
!
!-----------------------------------------------------------------------
!
      t=0.
      tem2=0.
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            t(i,j,k)=(ptprt(i,j,k)+ptbar(i,j,k))*                       &
                    (((pprt(i,j,k) + pbar(i,j,k))*p0inv) ** rddcp)
          END DO
        END DO
      END DO

      deallocate(pbar)
      deallocate(pprt)
      deallocate(ptbar)
      deallocate(ptprt)

      allocate(ref(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:ref")
      ref=0.
      allocate(refz(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:refz")
      refz=0.

      allocate(refh(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:refh")
      refh=0.
      allocate(refv(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:refv")
      refv=0.
      allocate(refhv(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:refhv")
      refhv=0.
      allocate(kdp(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:kdp")
      kdp=0.

      IF(dualpol > 0 .AND. wrtdualp > 0) THEN
        allocate(zdr(nx,ny,nz),stat=istatus)
        CALL check_alloc_status(istatus, "radaremul:zdr")
        zdr=0.
        allocate(rhohv(nx,ny,nz),stat=istatus)
        CALL check_alloc_status(istatus, "radaremul:rhohv")
        rhohv=0.
      END IF

!     CALL reflec_ferrier(nx,ny,nz, tem2, qr, qs, qh, tem1, ref)
      CALL reflec(nx,ny,nz, rhobar, qscalar(:,:,:,P_QR), qscalar(:,:,:,P_QS), qscalar(:,:,:,P_QH), ref)

      IF(dualpol > 0 .AND. wrtdualp > 0) THEN
        IF(dualpopt == 1) THEN
          print *, ' Using S-band Rayleigh Approx for dual-pol variables'
          CALL init_dsd()
          CALL calcConstants()
          DO k=1,nz-1
            DO j=1,ny-1
              DO i=1,nx-1
                obs_dual=calculate_obs(rhobar(i,j,k),qr(i,j,k),          &
                           qs(i,j,k),qh(i,j,k),t(i,j,k),2,drysnowopt)
                refh(i,j,k) = obs_dual%T_sum_ref_h
                refv(i,j,k) = obs_dual%T_sum_ref_v
                ref(i,j,k) = refh(i,j,k)
                obs_dual=calculate_obs(rhobar(i,j,k),qr(i,j,k),          &
                           qs(i,j,k),qh(i,j,k),t(i,j,k),3,drysnowopt)
                refhv(i,j,k) = obs_dual%T_sum_ref_hv
                kdp(i,j,k) = calculate_kdp(rhobar(i,j,k),qr(i,j,k),      &
                         qs(i,j,k),qh(i,j,k),drysnowopt)
              END DO
            END DO
          END DO
        ELSE
          print *, ' Using T-matrix X-band for dual-pol variables 1'
          CALL init_dsd()
          CALL calcConstants()
          CALL tmdualref(nx,ny,nz,tmatrix_dir,dbz_out,rhobar,qr,qs,qh,   &
                       refh,refv,refhv,zdr,rhohv,kdp)
          print *, ' Finished Calculating dual pol variables 1'
        END IF
        deallocate(zdr)
        deallocate(rhohv)
      END IF

      deallocate(tem1)
      deallocate(rhobar)

      IF(dualpol > 0) THEN
        DO k=1,nz-1
          DO j=1,ny-1
            DO i=1,nx-1
              refh(i,j,k)=max(refh(i,j,k),0.)
              refz(i,j,k)=refh(i,j,k)
              IF(refz(i,j,k) > epsilon) THEN
                ref(i,j,k)=10.*LOG10(refz(i,j,k))
              ELSE
                ref(i,j,k)=0.
              END IF
              refv(i,j,k)=max(refv(i,j,k),0.)
              refhv(i,j,k)=max(refhv(i,j,k),0.)
              kdp(i,j,k)=min(max(kdp(i,j,k),-30.),30.)
            END DO
          END DO
        END DO
      ELSE
        DO k=1,nz-1
          DO j=1,ny-1
            DO i=1,nx-1
              IF(ref(i,j,k) > 0.) THEN
                refz(i,j,k)=10.0**(0.1*ref(i,j,k))
              ELSE
                refz(i,j,k)=0.0
              END IF
            END DO
          END DO
        END DO
      END IF

    ELSE  ! hdf format

      allocate(itmp(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:hdf:itmp")
      itmp=0
      allocate(hmax(nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:hdf:hmax")
      hmax=0.
      allocate(hmin(nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:hdf:hmin")
      hmin=0.

      IF(ifile == 1) THEN

        CALL get_gridxyzp_hdf(nx,ny,nz,grdbasfn(1:lengbf),             &
                          x,y,zp,                                      &
                          itmp,hmax,hmin,ireturn)

        IF( ireturn /= 0 ) THEN
          WRITE(6,'(a,a)') ' Problem getting x,y,zp from hdf file ',   &
                         grdbasfn(1:lengbf)
          WRITE(6,'(a)') ' Program stopped.'
          STOP
        END IF

        print *,' HDF info: x,y,z of input data read in.'
        print *,' x(1 )=',x(1)
        print *,' x(nx)=',x(nx)
        print *,' y(1 )=',y(1)
        print *,' y(ny)=',y(ny)
      END IF
!
!-----------------------------------------------------------------------
!
!   Read in the wind data and average to the scalar points
!   Here the tem arrays are used to store the total wind components
!   Note also: umove and vmove are returned via common
!
!-----------------------------------------------------------------------
!
      allocate(tem1(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:hdf:tem1")
      tem1=0.
      allocate(tem2(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:hdf:tem2")
      tem2=0.
      allocate(tem3(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:hdf:tem3")
      tem3=0.
      allocate(usc(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:hdf:usc")
      usc=0.
      allocate(vsc(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:hdf:vsc")
      vsc=0.
      allocate(wsc(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:hdf:wsc")
      wsc=0.
      allocate(rhobar(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:hdf:rhobar")
      rhobar=0.

      CALL hdfreaduvw(nx,ny,nz,TRIM(hisfile(ifile)),                   &
                      time,tem1,tem2,tem3,                             &
                      itmp,hmax,hmin,ireturn)
      IF(ireturn < 0) EXIT

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            usc(i,j,k)=(0.5*(tem1(i,j,k)+tem1(i+1,j,k)))+umove
            vsc(i,j,k)=(0.5*(tem2(i,j,k)+tem2(i+1,j,k)))+vmove
            wsc(i,j,k)=0.5*(tem3(i,j,k)+tem3(i,j,k+1))
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!
!   Get pressure and potential temperature to find temperature (K)
!   Here the tem1 = total pressure and tem2 is total potential temp
!   Use the rhobar array to store rho.
!
!-----------------------------------------------------------------------
!
      tem1=0.
      tem2=0.
      CALL hdfreadppt(nx,ny,nz,TRIM(hisfile(ifile)),                    &
                      time,tem1,tem2,                                   &
                      itmp,hmax,hmin,ireturn)
      IF(ireturn < 0) EXIT

      allocate(t(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:hdf:t")
      t=0.
      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            t(i,j,k)=(tem2(i,j,k))*((tem1(i,j,k)*p0inv) ** rddcp)
            rhobar(i,j,k)=tem1(i,j,k)/(rd*t(i,j,k))
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!
!   Get the hydrometeors qr,qs, and qh
!
!-----------------------------------------------------------------------
!
      deallocate(tem1)
      deallocate(tem2)
      deallocate(tem3)

      allocate(qr(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:hdf:qr")
      qr=0.
      allocate(qs(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:hdf:qs")
      qs=0.
      allocate(qh(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:hdf:qh")
      qh=0.
      CALL hdfreadqrsh(nx,ny,nz,TRIM(hisfile(ifile)),                   &
                       time,qr,qs,qh,                                   &
                       itmp,hmax,hmin,ireturn)
      IF(ireturn < 0) EXIT
      print *, ' Back from hdfreadqrsh, ireturn =',ireturn
!
!-----------------------------------------------------------------------
!
!   Free the hdf io arrays
!
!-----------------------------------------------------------------------
!
      deallocate(itmp)
      deallocate(hmax)
      deallocate(hmin)
!
!-----------------------------------------------------------------------
!
!   Compute the reflectivity
!
!-----------------------------------------------------------------------
!
      allocate(ref(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:hdf:ref")
      ref=0.
      allocate(refz(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:hdf:refz")
      refz=0.

      allocate(refh(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:refh")
      refh=0.
      allocate(refv(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:refv")
      refv=0.
      allocate(refhv(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:refhv")
      refhv=0.
      allocate(kdp(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:kdp")
      kdp=0.

      IF(dualpol > 0 .AND. wrtdualp > 0) THEN
        allocate(zdr(nx,ny,nz),stat=istatus)
        CALL check_alloc_status(istatus, "radaremul:zdr")
        zdr=0.
        allocate(rhohv(nx,ny,nz),stat=istatus)
        CALL check_alloc_status(istatus, "radaremul:rhohv")
        rhohv=0.
      END IF

      print *, ' Allocation of grid-based radar variables complete'

      print *, ' Calling reflect '
!     CALL reflec_ferrier(nx,ny,nz, rhobar, tem1, tem2, tem3, t, ref)
      CALL reflec(nx,ny,nz, rhobar, qr, qs, qh, ref)
      print *, ' Back from reflect'

      IF(dualpol > 0 .AND. wrtdualp > 0) THEN
        CALL init_dsd()
        CALL calcConstants()
        IF(dualpopt == 1) THEN
          DO k=1,nz-1
            DO j=1,ny-1
              DO i=1,nx-1
                obs_dual=calculate_obs(rhobar(i,j,k),qr(i,j,k),          &
                         qs(i,j,k),qh(i,j,k),t(i,j,k),2,drysnowopt)
                refh(i,j,k) = obs_dual%T_sum_ref_h
                refv(i,j,k) = obs_dual%T_sum_ref_v
                ref(i,j,k) = refh(i,j,k)
                obs_dual=calculate_obs(rhobar(i,j,k),qr(i,j,k),          &
                         qs(i,j,k),qh(i,j,k),t(i,j,k),3,drysnowopt)
                refhv(i,j,k) = obs_dual%T_sum_ref_hv
                kdp(i,j,k) = calculate_kdp(rhobar(i,j,k),qr(i,j,k),      &
                           qs(i,j,k),qh(i,j,k),drysnowopt)
              END DO
            END DO
          END DO
        ELSE IF(dualpol > 0 .AND. wrtdualp > 0 .AND. dualpopt == 2) THEN
          print *, ' Using T-matrix X-band for dual-pol variables 2'
          CALL tmdualref(nx,ny,nz,tmatrix_dir,dbz_out,rhobar,qr,qs,qh,   &
                       refh,refv,refhv,zdr,rhohv,kdp)
          print *, ' Finished Calculating dual pol variables 2'
        END IF

        deallocate(zdr)
        deallocate(rhohv)

      END IF
!
!  Assign refz
!
      IF(dualpol > 0 .AND. wrtdualp > 0) THEN
        DO k=1,nz-1
          DO j=1,ny-1
            DO i=1,nx-1
              refh(i,j,k)=max(refh(i,j,k),0.)
              refz(i,j,k)=refh(i,j,k)
              IF(refz(i,j,k) > epsilon) THEN
                ref(i,j,k)=10.*LOG10(refz(i,j,k))
              ELSE
                ref(i,j,k)=0.
              END IF
              refv(i,j,k)=max(refv(i,j,k),0.)
              refhv(i,j,k)=max(refhv(i,j,k),0.)
              kdp(i,j,k)=min(max(kdp(i,j,k),-30.),30.)
            END DO
          END DO
        END DO
      ELSE
        DO k=1,nz-1
          DO j=1,ny-1
            DO i=1,nx-1
              IF(ref(i,j,k) > 0.) THEN
                refz(i,j,k)=10.0**(0.1*ref(i,j,k))
              ELSE
                refz(i,j,k)=0.0
              END IF
            END DO
          END DO
        END DO
      END IF
!
!-----------------------------------------------------------------------
!
!   Free the temporary i/o arrays
!
!-----------------------------------------------------------------------
!
      deallocate(rhobar)

    END IF  ! hdf format
!
!-----------------------------------------------------------------------
!
!  Set up fake reflectivity for testing terrain blocking.
!
!-----------------------------------------------------------------------

!      CALL fake_reflec(nx,ny,nz,ref)
!      CALL fake_vel(nx,ny,nz,usc,vsc)
!
!-----------------------------------------------------------------------
!
!  Make a few calculations now that we have read the datafile.
!
!-----------------------------------------------------------------------
!
    dx=x(2)-x(1)
    dy=y(2)-y(1)
    dxinv=1./dx
    dyinv=1./dy
    twdx2inv=1./(2.*dx*dx)
    twdxinv=1./(2.*dx)
    twdy2inv=1./(2.*dy*dy)
    twdyinv=1./(2.*dy)

    curtim = time
    IF(ifile == 1 .AND. adjctr == 2) time0=curtim

    CALL ctim2abss(year,month,day,hour,minute,second,itime)
    ifsecs=NINT(curtim)
    itime=itime+ifsecs
    itimtilt=itime
    CALL gtlfnkey(runname,lfnkey)

    print *, ' runname:',runname(1:lfnkey)
    print *, ' yr,mo,day,hr,min,sec:',year,month,day
    print *, ' hr,min,sec:',hour,minute,second
    print *, ' ifsecs: ',ifsecs
    print *, ' curtim: ',curtim
    print *, ' mapproj,ctrlat,ctrlon: ',mapproj,ctrlat,ctrlon

    CALL abss2ctim(itime,year,month,day,hour,minute,second)
!
!-----------------------------------------------------------------------
!
!   Average spatial locations to scalar points
!   Need only do this once, for the first file.
!
!-----------------------------------------------------------------------
!
    IF(ifile == 1) THEN
      allocate(xsc(nx),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:xsc")
      xsc=0.
      allocate(ysc(ny),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:ysc")
      ysc=0.
      allocate(zps(nx,ny,nz),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:zps")
      zps=0.

      DO i=1,nx-1
        xsc(i)=0.5*(x(i)+x(i+1))
      END DO
      xsc(nx)=xsc(nx-1)+dx

      DO j=1,ny-1
        ysc(j)=0.5*(y(j)+y(j+1))
      END DO
      ysc(ny)=ysc(ny-1)+dy
      write(6,'(a,f10.2,a,f10.2,a)') ' Coordinate of x(2): ',(x(2)*0.001),    &
                                   ' km  y(2): ',(y(2)*0.001),' km'
      IF( adjhgt > 0) THEN
        write(6,'(a,f10.2,a)')                                             &
        ' Adding height offset to grid:',hgtoffset,' m.'
      ELSE
        hgtoffset=0.
      END IF

      DO k=1,nz-1
        DO j=1,ny-1
          DO i=1,nx-1
            zps(i,j,k)=(0.5*(zp(i,j,k)+zp(i,j,k+1)))+hgtoffset
          END DO
        END DO
      END DO
      DO j=1,ny-1
        DO i=1,nx-1
          zps(i,j,nz)=(2.0*zps(i,j,nz-1))-zps(i,j,nz-2)
        END DO
      END DO
    END IF ! first file
!
!-----------------------------------------------------------------------
!
! Allocate additional grid location and geometry arrays
! These speed calculations down the road.
!
!-----------------------------------------------------------------------
!
    allocate(latsc(nx,ny),stat=istatus)
    CALL check_alloc_status(istatus, "radaremul:latsc")
    latsc=0.
    allocate(lonsc(nx,ny),stat=istatus)
    CALL check_alloc_status(istatus, "radaremul:lonsc")
    lonsc=0.
    allocate(azmsc(nx,ny),stat=istatus)
    CALL check_alloc_status(istatus, "radaremul:azmsc")
    azmsc=0.
    allocate(sfcr(nx,ny),stat=istatus)
    CALL check_alloc_status(istatus, "radaremul:sfcr")
    sfcr=0.
    allocate(elvsc(nx,ny,nz),stat=istatus)
    CALL check_alloc_status(istatus, "radaremul:elvsc")
    elvsc=0.
    allocate(rngsc(nx,ny,nz),stat=istatus)
    CALL check_alloc_status(istatus, "radaremul:rngsc")
    rngsc=0.
    allocate(vrsc(nx,ny,nz),stat=istatus)
    CALL check_alloc_status(istatus, "radaremul:vrsc")
    vrsc=0.
    allocate(usm(nx,ny,nz),stat=istatus)
    CALL check_alloc_status(istatus, "radaremul:usm")
    usm=0.
    allocate(vsm(nx,ny,nz),stat=istatus)
    CALL check_alloc_status(istatus, "radaremul:vsm")
    vsm=0.
    allocate(vort(nx,ny,nz),stat=istatus)
    CALL check_alloc_status(istatus, "radaremul:vort")
    vort=0.
    allocate(cmpref(nx,ny),stat=istatus)
    CALL check_alloc_status(istatus, "radaremul:cmpref")
    cmpref=0.
!
!   For mobile radar interpolate the location from the
!   location arrays using the present time, itime
!
    IF(mobileopt > 0) THEN
      IF(locopt == 1) &
        CALL lltoxy(ntime,1,moblat,moblon,mobilex,mobiley)
      IF(ntime > 1) THEN
        DO jtime=1,(ntime-1)
          IF(mobtime(jtime) > itime) EXIT
        END DO
        tw1=float(itime-mobtime(jtime-1))/ &
            float(mobtime(jtime)-mobtime(jtime-1))
        print *, ' Data time: ',itime,'  jtime: ',jtime
        print *, ' mobile times:',mobtime(jtime-1),&
                                  mobtime(jtime)
        print *, ' Prelim time weight, tw1: ',tw1
        tw1=min(max(0.,tw1),1.0)
        tw0=1.0-tw1
        print *, ' Time weights, tw0: ',tw0,'  tw1: ',tw1
        xrad=tw0*mobilex(jtime-1)+tw1*mobilex(jtime)
        yrad=tw0*mobiley(jtime-1)+tw1*mobiley(jtime)
        radelv=tw0*mobhgt(jtime-1)+tw1*mobhgt(jtime)
        rhead=tw0*mobhead(jtime-1)+tw1*mobhead(jtime)
        rpitch=tw0*mobpitch(jtime-1)+tw1*mobpitch(jtime)
        IF(locopt == 1) CALL xytoll(1,1,xrad,yrad,radlat,radlon)
      ELSE
        IF(locopt == 1) THEN
          radlat=moblat(1)
          radlon=moblat(1)
        ELSE
          xrad=mobilex(1)
          yrad=mobiley(1)
        END IF
      END IF
    END IF
!
!-----------------------------------------------------------------------
!
! Define output data size and find location of radar in the grid.
! Account for umove and vmove by moving the radar relative to the grid.
!
! If no map projection has been defined, for the purposes of emulation,
! set the projection to Lambert conformal and center the grid on
! ctrlat, ctrlon.
!
!-----------------------------------------------------------------------
!
    IF ( adjmove > 0 ) THEN
      umove=umovein
      vmove=vmovein
    END IF
    IF ( ifile == 1 ) THEN
      IF ( mapproj == 0 ) THEN
        IF ( adjctr > 0 ) THEN
          ctrlat = ctrlatem
          ctrlon = ctrlonem
        END IF
        print *, ' ctrlat, ctrlon: ',ctrlat,ctrlon
        trulat1 = ctrlat
        trulat2 = ctrlat
        trulon  = ctrlon
        alatnot(1)=trulat1
        alatnot(2)=trulat2
        CALL setmapr(2,1.0,alatnot,trulon)
        CALL lltoxy(1,1,ctrlat,ctrlon,ctrx,ctry)
        xsw=ctrx-0.5*((nx-3)*dx)
        ysw=ctry-0.5*((ny-3)*dy)
        CALL setorig(1,xsw,ysw)
      ELSE IF ( adjctr > 0 ) THEN
        ctrlat = ctrlatem
        ctrlon = ctrlonem
        alatnot(1)=trulat1
        alatnot(2)=trulat2
        CALL setmapr(mapproj,1.0,alatnot,trulon)
        CALL lltoxy(1,1,ctrlat,ctrlon,ctrx,ctry)
        xsw=ctrx-0.5*((nx-3)*dx)
        ysw=ctry-0.5*((ny-3)*dy)
        CALL setorig(1,xsw,ysw)
      ELSE
        alatnot(1)=trulat1
        alatnot(2)=trulat2
        CALL setmapr(mapproj,1.0,alatnot,trulon)
        CALL lltoxy(1,1,ctrlat,ctrlon,ctrx,ctry)
        xsw=ctrx-0.5*((nx-3)*dx)
        ysw=ctry-0.5*((ny-3)*dy)
        CALL setorig(1,xsw,ysw)
      END IF
    END IF

    IF(locopt == 1) THEN
      print *, ' radlat,radlon=',radlat,radlon
      CALL lltoxy(1,1,radlat,radlon,xrad,yrad)
      xrada=xrad-umove*(curtim-time0)
      yrada=yrad-vmove*(curtim-time0)
      CALL xytoll(1,1,xrada,yrada,radlata,radlona)
      irloc=((xrada-xsc(1))*dxinv)+1.0
      jrloc=((yrada-ysc(1))*dyinv)+1.0
      WRITE(6,'(a,/a,f10.4,a,f10.4,/a,f10.4,a,f10.4,/a,f10.4,a,f10.4)') &
      '  Radar position (adjusted for umove, vmove)',                   &
      '  Radar latitude:',radlata,'  longitude: ',radlona,              &
      '  x-coord (km):',(0.001*xrada),'  y-coord (km):',(0.001*yrada),  &
      '  i-location: ',irloc,'  j-location:',jrloc
    ELSE
      xrada=xrad*1000.
      yrada=yrad*1000.
      CALL xytoll(1,1,xrada,yrada,radlat,radlon)
      radlata=radlat
      radlona=radlon
      irloc=((xrada-xsc(1))*dxinv)+1.0
      jrloc=((yrada-ysc(1))*dyinv)+1.0
      WRITE(6,'(a,f10.4,a,f10.4,/a,f10.4,a,f10.4,/a,f10.4,a,f10.4)')    &
        '  Radar latitude:',radlata,'  longitude: ',radlona,            &
        '  x-coord (km):',(0.001*xrada),'  y-coord (km):',(0.001*yrada),&
        '  i-location: ',irloc,'  j-location:',jrloc
    END IF
    iloc=INT(irloc)
    jloc=INT(jrloc)

    trngt=0.
    IF(iloc > 0 .AND. iloc < nx-2 .AND.                        &
       jloc > 0 .AND. jloc < ny-2 ) THEN
      w2i=irloc-iloc
      w1i=1.0-w2i
      w2j=jrloc-jloc
      w1j=1.0-w2j
      w11=w1i*w1j
      w21=w2i*w1j
      w12=w1i*w2j
      w22=w2i*w2j
      trngt=w11*zp(iloc,jloc,2) +          &
            w21*zp(iloc+1,jloc,2) +        &
            w12*zp(iloc,jloc+1,2) +        &
            w22*zp(iloc+1,jloc+1,2)
    END IF
    WRITE(6,'(a,f8.1,a,f8.1//)') &
        '  Radar Elevation: ',radelv,' Terrain Height= ',trngt
!
!-----------------------------------------------------------------------
!
! Calculate radar parameters at each scalar grid point.
!
!-----------------------------------------------------------------------
    print *, ' calculating radar parameters'
!
    IF( mapproj > 0 ) THEN
      CALL xytoll(nx,ny,xsc,ysc,latsc,lonsc)

      DO j=1,ny-1
        DO i=1,nx-1
          CALL disthead(radlata,radlona,latsc(i,j),lonsc(i,j),          &
                    azmsc(i,j),sfcr(i,j))
          IF(sfcr(i,j) > rngmin) THEN
            DO k=1,nz-1
              CALL beamelv((zps(i,j,k)-radelv),sfcr(i,j),               &
                            elvsc(i,j,k),rngsc(i,j,k))
            END DO
          ELSE
            rngsc(i,j,k)=sfcr(i,j)
            elvsc(i,j,k)=-99.
          END IF
        END DO
      END DO
    ELSE
      DO j=1,ny-1
        DO i=1,nx-1
          delx=xsc(i)-xrada
          dely=ysc(j)-yrada
          sfcr(i,j)=sqrt(delx*delx + dely*dely)
          IF(dely > 0.) THEN
            azmsc(i,j)=rad2deg*atan(delx/dely)
          ELSE IF (dely < 0.) THEN
            azmsc(i,j)=180.+rad2deg*atan(delx/dely)
          ELSE IF (delx > 0.) THEN
            azmsc(i,j)=90.
          ELSE
            azmsc(i,j)=270.
          END IF
          IF(azmsc(i,j) < 0.) azmsc(i,j)=azmsc(i,j)+360.
          IF(sfcr(i,j) > rngmin) THEN
            DO k=1,nz-1
              CALL beamelv((zps(i,j,k)-radelv),sfcr(i,j),              &
                            elvsc(i,j,k),rngsc(i,j,k))
            END DO
          ELSE
            rngsc(i,j,k)=sfcr(i,j)
            elvsc(i,j,k)=-99.
          END IF
        END DO
      END DO
    END IF

!
!-----------------------------------------------------------------------
!
! Calculate radial velocity, including terminal velocity, and reflectivity
! factor (in mm6/m3) at all scalar grid points.
!
! The termvel subroutine can be updated sometime to use actual freezing
! level and/or the model species variables to calculate vt.
!
!-----------------------------------------------------------------------
!
    print *, ' Calculating radial velocity'
    DO j=1,ny-1
      DO i=1,nx-1
        ucmp=sin(deg2rad*azmsc(i,j))
        vcmp=cos(deg2rad*azmsc(i,j))
        cmpref(i,j)=-99.
        DO k=2,nz-1
          vt=0.
          IF(ref(i,j,k) > epsilon) THEN
            CALL termvel(ref(i,j,k),zps(i,j,k),vt)
          END IF
          CALL dhdrange(elvsc(i,j,k),rngsc(i,j,k),dhdr)
          dsdr=sqrt(1.-dhdr*dhdr)
          vrsc(i,j,k)=dsdr*(ucmp*usc(i,j,k) + vcmp*vsc(i,j,k)) +        &
                      dhdr*(wsc(i,j,k)-vt)
          cmpref(i,j)=max(cmpref(i,j),ref(i,j,k))
        END DO
      END DO
    END DO

    usm=usc
    vsm=usc

    print *, ' Smoothing velocities, nsmvort=',nsmvort
    IF(nsmvort > 0) THEN
      DO ismth=1,nsmvort
        CALL smooth9p( usm, nx,ny, 1, nx-1, 1, ny-1, 0, vort )
        CALL smooth9p( vsm, nx,ny, 1, nx-1, 1, ny-1, 0, vort )
      END DO
      vort=0.
    END IF
!
!  Calculate vertical vorticity at scalar points
!
    WRITE(6,'(a,f10.8,a,f10.8)') ' Calculating vorticity: dxinv=',      &
         dxinv,' dyinv: ',dyinv
    DO k=1,nz
      DO j=2,ny-2
        DO i=2,nx-2
          vort(i,j,k)=0.5*((vsm(i+1,j,k)-vsm(i-1,j,k))*dxinv -          &
                           (usm(i,j+1,k)-usm(i,j-1,k))*dyinv)
        END DO
      END DO
    END DO
    DO k=1,nz
      DO j=2,ny-2
        vort(1,j,k)=vort(2,j,k)
        vort(nx-1,j,k)=vort(nx-2,j,k)
      END DO
    END DO
    DO k=1,nz
      DO i=1,nx-1
        vort(i,1,k)=vort(1,2,k)
        vort(i,ny-1,k)=vort(i,ny-2,k)
      END DO
    END DO
    vvmax=maxval(vort)
    vvmin=minval(vort)
    WRITE(6,'(a,g16.9,a,g16.9)') ' Grid min vert vort=',vvmin, &
      '  Max vert vort = ',vvmax
!
    deallocate(usm)
    deallocate(vsm)
!
!   Allocate output radar arrays
!
    allocate(azim(maxazim),stat=istatus)
    CALL check_alloc_status(istatus, "radaremul:azim")
    azim=0.
    allocate(beamw(maxazim),stat=istatus)
    CALL check_alloc_status(istatus, "radaremul:beamw")
    beamw=0.
    allocate(gtspc(maxazim),stat=istatus)
    CALL check_alloc_status(istatus, "radaremul:gtspc")
    gtspc=0.
    allocate(vnyq(maxazim),stat=istatus)
    CALL check_alloc_status(istatus, "radaremul:vnyq")
    vnyq=0.
    allocate(refl(ngate,maxazim),stat=istatus)
    CALL check_alloc_status(istatus, "radaremul:refl")
    refl=0.
    allocate(attrefl(ngate,maxazim),stat=istatus)
    CALL check_alloc_status(istatus, "radaremul:attrefl")
    attrefl=0.
    allocate(rzdr(ngate,maxazim),stat=istatus)
    CALL check_alloc_status(istatus, "radaremul:rzdr")
    rzdr=0.
    allocate(rkdp(ngate,maxazim),stat=istatus)
    CALL check_alloc_status(istatus, "radaremul:rkdp")
    rkdp=0.
    allocate(rrhohv(ngate,maxazim),stat=istatus)
    CALL check_alloc_status(istatus, "radaremul:rrhohv")
    rrhohv=0.
    allocate(radv(ngate,maxazim),stat=istatus)
    CALL check_alloc_status(istatus, "radaremul:radv")
    radv=0.
    allocate(attradv(ngate,maxazim),stat=istatus)
    CALL check_alloc_status(istatus, "radaremul:attradv")
    attradv=0.
    allocate(stdvr(ngate,maxazim),stat=istatus)
    CALL check_alloc_status(istatus, "radaremul:stdvr")
    stdvr=0.
    allocate(vvort(ngate,maxazim),stat=istatus)
    CALL check_alloc_status(istatus, "radaremul:vvort")
    vvort=0.
    IF(wrtqx > 0) THEN
      allocate(qrrad(ngate,maxazim),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:qrrad")
      qrrad=0.
      allocate(qsrad(ngate,maxazim),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:qsrad")
      qsrad=0.
      allocate(qhrad(ngate,maxazim),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:qhrad")
      qhrad=0.
    END IF
    IF(wrtuvwt > 0) THEN
      allocate(urad(ngate,maxazim),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:urad")
      urad=0.
      allocate(vrad(ngate,maxazim),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:vrad")
      vrad=0.
      allocate(wrad(ngate,maxazim),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:wrad")
      wrad=0.
      allocate(tkrad(ngate,maxazim),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:tkrad")
      tkrad=0.
    END IF
    allocate(sens(ngate),stat=istatus)
    CALL check_alloc_status(istatus, "radaremul:sens")
    sens=0.

    IF( ifmt == 1 ) THEN
      allocate(itimvol(nelev),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:itimvol")
      itimvol=0
      allocate(vnyqvol(nelev),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:vnyqvol")
      vnyqvol=0.
      allocate(rngvol(ngate,nelev),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:rngvol")
      rngvol=0.
      allocate(azmvol(maxazim,nelev),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:azmvol")
      azmvol=0.
      allocate(elvvol(maxazim,nelev),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:elvvol")
      elvvol=0.
      allocate(refvol(ngate,maxazim,nelev),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:refvol")
      refvol=0.
      allocate(uarefvol(ngate,maxazim,nelev),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:uarefvol")
      uarefvol=0.
      allocate(velvol(ngate,maxazim,nelev),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:velvol")
      velvol=0.
      allocate(uavelvol(ngate,maxazim,nelev),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:uavelvol")
      uavelvol=0.
      allocate(vorvol(ngate,maxazim,nelev),stat=istatus)
      CALL check_alloc_status(istatus, "radaremul:vorvol")
      vorvol=0.
      IF(dualpol > 0) THEN
        allocate(zdrvol(ngate,maxazim,nelev),stat=istatus)
        CALL check_alloc_status(istatus, "radaremul:zdrvol")
        zdrvol=0.
        allocate(kdpvol(ngate,maxazim,nelev),stat=istatus)
        CALL check_alloc_status(istatus, "radaremul:kdpvol")
        kdpvol=0.
        allocate(rhohvvol(ngate,maxazim,nelev),stat=istatus)
        CALL check_alloc_status(istatus, "radaremul:rhohvvol")
        rhohvvol=0.
      ELSE
        wrtdualp=0
        allocate(zdrvol(1,1,1),stat=istatus)
        CALL check_alloc_status(istatus, "radaremul:zdrvol")
        zdrvol=0.
        allocate(kdpvol(1,1,1),stat=istatus)
        CALL check_alloc_status(istatus, "radaremul:kdpvol")
        kdpvol=0.
        allocate(rhohvvol(1,1,1),stat=istatus)
        CALL check_alloc_status(istatus, "radaremul:rhohvvol")
        rhohvvol=0.
      END IF
      IF(wrtqx > 0) THEN
        allocate(qrvol(ngate,maxazim,nelev),stat=istatus)
        CALL check_alloc_status(istatus, "radaremul:qrvol")
        qrvol=0.
        allocate(qsvol(ngate,maxazim,nelev),stat=istatus)
        CALL check_alloc_status(istatus, "radaremul:qsvol")
        qsvol=0.
        allocate(qhvol(ngate,maxazim,nelev),stat=istatus)
        CALL check_alloc_status(istatus, "radaremul:qhvol")
        qhvol=0.
      ELSE
        allocate(qrvol(1,1,1),stat=istatus)
        CALL check_alloc_status(istatus, "radaremul:qrvol")
        qrvol=0.
        allocate(qsvol(1,1,1),stat=istatus)
        CALL check_alloc_status(istatus, "radaremul:qsvol")
        qsvol=0.
        allocate(qhvol(1,1,1),stat=istatus)
        CALL check_alloc_status(istatus, "radaremul:qhvol")
        qhvol=0.
      END IF
      IF(wrtuvwt > 0) THEN
        allocate(uvol(ngate,maxazim,nelev),stat=istatus)
        CALL check_alloc_status(istatus, "radaremul:uvol")
        uvol=0.
        allocate(vvol(ngate,maxazim,nelev),stat=istatus)
        CALL check_alloc_status(istatus, "radaremul:vvol")
        vvol=0.
        allocate(wvol(ngate,maxazim,nelev),stat=istatus)
        CALL check_alloc_status(istatus, "radaremul:wvol")
        wvol=0.
        allocate(tkvol(ngate,maxazim,nelev),stat=istatus)
        CALL check_alloc_status(istatus, "radaremul:tkvol")
        tkvol=0.
      ELSE
        allocate(uvol(1,1,1),stat=istatus)
        CALL check_alloc_status(istatus, "radaremul:uvol")
        uvol=0.
        allocate(vvol(1,1,1),stat=istatus)
        CALL check_alloc_status(istatus, "radaremul:vvol")
        vvol=0.
        allocate(wvol(1,1,1),stat=istatus)
        CALL check_alloc_status(istatus, "radaremul:wvol")
        wvol=0.
        allocate(tkvol(1,1,1),stat=istatus)
        CALL check_alloc_status(istatus, "radaremul:tkvol")
        tkvol=0.
      END IF
    END IF
!
!-----------------------------------------------------------------------
!
! Open log file for index records
!
!-----------------------------------------------------------------------
!
    IF( creidx > 0 .AND. ifmt == 2) THEN
      write(idxfname,'(a,a,a,a,a,a,i6.6,a)') &
        TRIM(outdir),'/',radname,'.',runname(1:lfnkey),'.',ifsecs,'.xml'
      CALL getunit(idxunit)
      write(6,'(a,a)') ' Writing index to:',TRIM(idxfname)
      OPEN(idxunit,file=idxfname,form='formatted',status='unknown')
      write(idxunit,'(a)')                                             &
                      '<?xml version="1.0" encoding="iso-8859-1" ?>'
      write(idxunit,'(a,a,a)') '<codeindex type="netcdf" dataset="',   &
         TRIM(outdir),'">'
    END IF
!
!-----------------------------------------------------------------------
!
! Create data from Cartesian grid.
!
!-----------------------------------------------------------------------
!
    irot=1
    IF(rotdir == 2) irot=-1
    IF(rhimode > 0) krot=-1*krot
    IF(krot > 0) THEN
      itltbgn=1
      itltend=nelev
    ELSE
      itltbgn=nelev
      itltend=1
    END IF

    WRITE(6,'(a,i6)') 'Entering tilt loop, nelev=',nelev
    DO itilt=itltbgn,itltend,krot
      kntgate=0
      kntintr=0
      kntintv=0
      kntfewr=0
      kntfewv=0
      kntclr=0
      knthob=0
      kntvob=0
      kntbob=0

      write(6,'(/a,i3,a,f6.2,a)')                                      &
        ' Processing tilt ',itilt,'  elev=',elvang(itilt),' deg'

      azim=rmisval
      beamw=rmisval
      gtspc=rmisval
      vnyq=rmisval
      refl=rmisval
      attrefl=rmisval
      radv=rmisval
      attradv=rmisval

      IF(wrtqx > 0) THEN
        qrrad=0.
        qsrad=0.
        qhrad=0.
      END IF

      IF(wrtuvwt > 0) THEN
        urad=rmisval
        vrad=rmisval
        wrad=rmisval
        tkrad=rmisval
      END IF

      vvort=0.
      vvmax=-999.
      vvmin=999.

      elv=elvang(itilt)+rpitch
      ielv=INT(elv)
      felv=NINT(100.*(elv-ielv))
      elvtop=elv+delev(nptselv)
      elvbot=elv+delev(1)
      write(6,'(a,f6.2,a,f6.2)')                                      &
        ' Elv ang beam bot=',elvbot,'  Elv ang beam top=',elvtop

      IF(tmadvopt == 1) THEN
        itimtilt=itime+NINT((itilt-1)*timeincr)
      END IF
      print *, ' itilt, itimtilt, tmadvopt, timeincr: ', &
                   itilt,itimtilt,tmadvopt,timeincr

      IF(rhimode < 1) THEN
        IF(rotropt > 1) THEN
          nazim=NINT(azmsect/delazm)+1
          IF(dualprf > 0) THEN
            ttime=nazim*((prt1(itilt)*npulse1(itilt))+                   &
                         (prt2(itilt)*npulse2(itilt)))
          ELSE
            ttime=nazim*prt1(itilt)*npulse1(itilt)
          END IF
          rotrate=azmsect/ttime
        ELSE
          IF(dualprf > 0) THEN
            delazm=rotrate*((prt1(itilt)*npulse1(itilt))+                 &
                           (prt2(itilt)*npulse2(itilt)))
          ELSE
            delazm=rotrate*prt1(itilt)*npulse1(itilt)
          END IF
          nazim=NINT(azmsect/delazm)+1
        END IF
      ELSE ! RHI
        nazim=1
      END IF
      WRITE(6,'(a,f10.2,a,/a,f10.2,a,/a,i6)') &
        ' Rotation rate:',rotrate,' deg/s', &
        ' Azimuth sampling increment:',delazm,' deg', &
        ' Number of azimuth samples:',nazim
!
!-----------------------------------------------------------------------
!
!  Calculate constants for range weighting.
!  Assumed here is a square wave transmitted with a matched filter
!  in the receiver resulting in a received wave envelope that is
!  approximately Gaussian.
!  See Eqs 11.118 and 5.76 in Doviak and Zrnic', 2nd Ed, 1993.
!
!  Note that pulselen has units of time, pulsewid units of meters
!
!-----------------------------------------------------------------------
!
      IF(dualprf == 0) THEN
        IF(pulseopt == 1) THEN
          pulselmn=pulselen1(itilt)
          pulsewmn=c*pulselmn
        ELSE
          pulsewmn=pulsewid1(itilt)
          pulselmn=cinv*pulsewmn
        END IF
        prtmn=prt1(itilt)
        npulsmn=npulse1(itilt)
        IF(nyqstopt == 1) THEN
          vnyqstl=vnyquist(itilt)
        ELSE
          vnyqstl=0.25*(0.01*wavelen)/prt1(itilt)
          vnyquist(itilt)=vnyqstl
        END IF
      ELSE
        IF(pulseopt == 1) THEN
          pulselmn=0.5*(pulselen1(itilt)*pulselen2(itilt))
          pulsewmn=c*pulselmn
        ELSE
          pulsewmn=0.5*(pulsewid1(itilt)*pulsewid2(itilt))
          pulselmn=cinv*pulsewmn
        END IF
        prtmn=0.5*(prt1(itilt)+prt2(itilt))
        npulsmn=NINT(0.5*(npulse1(itilt)+npulse2(itilt)))
        IF(nyqstopt == 1) THEN
          vnyqstl=vnyquist(itilt)
        ELSE
          vnyqstl=abs(0.25*(0.01*wavelen)/(prt1(itilt)-prt2(itilt)))
          vnyquist(itilt)=vnyqstl
        END IF
      END IF
      hlftau=0.5*pulselmn
      aconst=pi/(2.*sqrt(log(2.0)))
      b6=1.04/pulselmn
      hlfctau=0.5*pulsewmn
      sigr2=(0.35*hlfctau)*(0.35*hlfctau)
      sr2inv=2.0/(4.0*sigr2)
      rngwid=sqrt(sigr2)

      twovnyq=2.*vnyqstl
      twovninv=1./twovnyq

      WRITE(6,'(a,f10.2,a)') ' Mean PRT: ',(prtmn*1000.),' millisec'
      WRITE(6,'(a,f10.2,a)') ' Mean pulse length: ',(pulselmn*1.0E06), &
              ' microsec'
      WRITE(6,'(a,f10.2,a)') ' Mean pulse width: ',pulsewmn,' m'
      WRITE(6,'(a,f10.2,a)') &
        ' Gaussian filtered pulse half-width: ',rngwid,' m'
      WRITE(6,'(a,f10.2,a)') ' Unambiguous velocity: ',vnyqstl,' m/s'

      print *, 'rotrate,prtmn,npulsmn: ',rotrate,prtmn,npulsmn
      IF( abs(rotrate) > 0.) THEN
        efbmwid=EFFBMW(beamwid,rotrate,prtmn,npulsmn)
      ELSE
        efbmwid=beamwid
      END IF
      WRITE(6,'(a,f8.2,a)') ' Static beamwidth:',beamwid,' degrees.'
      WRITE(6,'(a,f8.2,a)') ' Effective beamwidth:',efbmwid,' degrees.'

      da=(evalwid*efbmwid)/float(nptsazm-1)
      WRITE(6,'(a)') &
       ' Horiz cross-beam evaluation pts (azimuth degrees from center):'
      DO jpt=1,nptsazm
        dazim(jpt)=(-0.5*evalwid*efbmwid)+((jpt-1)*da)
        WRITE(6,'(a,i4,a,f12.4)') ' dazim(',jpt,') = ',dazim(jpt)
      END DO
!
!-----------------------------------------------------------------------
!
! Calculate sensitivity as a function of range.
!
!-----------------------------------------------------------------------
!
      IF(senstvopt == 1) CALL xbsens(ngate,gatesp,pulselmn,beamwid,sens)

      wasqinv=fourln4/(efbmwid*efbmwid)
      bmwrad=deg2rad*efbmwid
      DO kpt=1,nptselv
      DO jpt=1,nptsazm
      DO ipt=1,nptsrng
        wgtpt(ipt,jpt,kpt)=exp(-((delrng(ipt)*delrng(ipt))*sr2inv  +    &
                                 (dazim(jpt)*dazim(jpt))*wasqinv +    &
                                 (delev(kpt)*delev(kpt))*wesqinv))
      END DO
      END DO
      END DO

      IF(rhimode < 1) THEN
        IF(rotdir == 1) THEN
          irot=1
        ELSE IF(rotdir == -1) THEN
          irot=-1
        ELSE IF (rotdir == 2 .OR. rotdir == -2) THEN
          irot=-1*irot
        END IF
        IF(irot > 0) THEN
          azmbgn=azim1+rhead
        ELSE
          azmbgn=azim2+rhead
        END IF
        IF(azmbgn > 360.) azmbgn=azmbgn-360.
      ELSE   ! RHI
        irot=1
        azmbgn=azim1
      END IF

      WRITE(6,'(a,f10.2,a,i4,a,f10.2)') &
        ' Elev:',elv,'  irot:',irot,'  azmbgn:',azmbgn

      IF(int_method == 1) THEN

      WRITE(6,'(a,i3)') 'Using integration method ',int_method

      kntavgr=0
      kntavgv=0
      DO jazim=1,nazim
        azimuth=azmbgn+irot*((jazim-1)*delazm)
        IF(azimuth > 360.) azimuth = azimuth-360.
        IF(azimuth < 0.) azimuth = azimuth+360.
!       IF(MOD((jazim-1),10) == 0) &
!         WRITE(6,'(a,f10.2)') ' Processing azimuth: ', azimuth
        azrad=deg2rad*azimuth
        azim(jazim)=azimuth
        beamw(jazim)=beamwid
        gtspc(jazim)=gatesp
        vnyq(jazim)=vnyquist(itilt)
        DO igate=1, ngate
          kntgate=kntgate+1
          srange=igate*gatesp
          CALL beamhgt(elv,srange,hgt,sfcrng)
          IF(mapproj > 0 ) THEN
            CALL gcircle(radlata,radlona,azimuth,sfcrng,gtlat,gtlon)
            CALL lltoxy(1,1,gtlat,gtlon,gtx,gty)
            delx=gtx-xrada
            dely=gty-yrada
            xrat=delx/srange
            yrat=dely/srange
            ipos=(gtx*dxinv)+1.5
            jpos=(gty*dyinv)+1.5
          ELSE
            delx=sin(azrad)*sfcrng
            dely=cos(azrad)*sfcrng
            xrat=delx/srange
            yrat=dely/srange
            ipos=irloc+(delx*dxinv)
            jpos=jrloc+(dely*dyinv)
          END IF
          iloc1=NINT(ipos)
          jloc1=NINT(jpos)
          iloc=INT(ipos)
          jloc=INT(jpos)

          IF(iloc > 0 .AND. iloc1 < nx-1 .AND.                         &
             jloc > 0 .AND. jloc1 < ny-1 ) THEN

            hgtmsl=hgt+radelv
            IF( hgtmsl > zps(iloc1,jloc1,2) .AND.                       &
                hgtmsl < zps(iloc1,jloc1,nz-1) ) THEN
              DO k=3,nz-2
                IF(zps(iloc1,jloc1,k) > hgtmsl) EXIT
              END DO

              bmwm=bmwrad*evalwid*srange
              sradius=max(bmwm,pulsewmn)
              iradius=max((1+INT(sradius/dx)),2)
              jradius=max((1+INT(sradius/dy)),2)
!             print *, 'bmwm,pulsewid = ',bmwm,pulsewmn
!             print *, 'iradius,jradius = ',iradius,jradius
!
!  First calculate provisional values using tri-linear interpolation
!  These will be used if there is insufficient density of data for a
!  weighted-average calculation
!
              rmax=0.
              r11l=max(refz(  iloc,  jloc,k-1),0.)
              rmax=max(rmax,r11l)
              r21l=max(refz(iloc+1,  jloc,k-1),0.)
              rmax=max(rmax,r21l)
              r12l=max(refz(  iloc,jloc+1,k-1),0.)
              rmax=max(rmax,r12l)
              r22l=max(refz(iloc+1,jloc+1,k-1),0.)
              rmax=max(rmax,r22l)
              r11h=max(refz(  iloc,  jloc,  k),0.)
              rmax=max(rmax,r11h)
              r21h=max(refz(iloc+1,  jloc,  k),0.)
              rmax=max(rmax,r21h)
              r12h=max(refz(  iloc,jloc+1,  k),0.)
              rmax=max(rmax,r12h)
              r22h=max(refz(iloc+1,jloc+1,  k),0.)
              rmax=max(rmax,r22h)
              whigh=(hgtmsl-zps(iloc1,jloc1,k-1))/                   &
                  (zps(iloc1,jloc1,k)-zps(iloc1,jloc1,k-1))
              wlow=1.-whigh
              w2i=ipos-iloc
              w1i=1.0-w2i
              w2j=jpos-jloc
              w1j=1.0-w2j
              w11=w1i*w1j
              w21=w2i*w1j
              w12=w1i*w2j
              w22=w2i*w2j
              vvort(igate,jazim)=whigh*(w11*vort(iloc,jloc,k) +        &
                                        w21*vort(iloc+1,jloc,k) +      &
                                        w12*vort(iloc,jloc+1,k) +      &
                                        w22*vort(iloc+1,jloc+1,k)) +   &
                                  wlow*(w11*vort(iloc,jloc,k-1) +      &
                                        w21*vort(iloc+1,jloc,k-1) +    &
                                        w12*vort(iloc,jloc+1,k-1) +    &
                                        w22*vort(iloc+1,jloc+1,k-1))
              vvmin=min(vvmin,vvort(igate,jazim))
              vvmax=max(vvmax,vvort(igate,jazim))
              IF(wrtqx > 0) THEN
                qrrad(igate,jazim)=whigh*(w11*qr(iloc,jloc,k) +        &
                                          w21*qr(iloc+1,jloc,k) +      &
                                          w12*qr(iloc,jloc+1,k) +      &
                                          w22*qr(iloc+1,jloc+1,k)) +   &
                                    wlow*(w11*qr(iloc,jloc,k-1) +      &
                                          w21*qr(iloc+1,jloc,k-1) +    &
                                          w12*qr(iloc,jloc+1,k-1) +    &
                                          w22*qr(iloc+1,jloc+1,k-1))
                qrrad(igate,jazim)=max(0.,qrrad(igate,jazim))
                qsrad(igate,jazim)=whigh*(w11*qs(iloc,jloc,k) +        &
                                          w21*qs(iloc+1,jloc,k) +      &
                                          w12*qs(iloc,jloc+1,k) +      &
                                          w22*qs(iloc+1,jloc+1,k)) +   &
                                    wlow*(w11*qs(iloc,jloc,k-1) +      &
                                          w21*qs(iloc+1,jloc,k-1) +    &
                                          w12*qs(iloc,jloc+1,k-1) +    &
                                          w22*qs(iloc+1,jloc+1,k-1))
                qsrad(igate,jazim)=max(0.,qsrad(igate,jazim))
                qhrad(igate,jazim)=whigh*(w11*qh(iloc,jloc,k) +        &
                                          w21*qh(iloc+1,jloc,k) +      &
                                          w12*qh(iloc,jloc+1,k) +      &
                                          w22*qh(iloc+1,jloc+1,k)) +   &
                                    wlow*(w11*qh(iloc,jloc,k-1) +      &
                                          w21*qh(iloc+1,jloc,k-1) +    &
                                          w12*qh(iloc,jloc+1,k-1) +    &
                                          w22*qh(iloc+1,jloc+1,k-1))
                qhrad(igate,jazim)=max(0.,qhrad(igate,jazim))
              END IF

              IF( rmax > 0. ) THEN
                CALL dhdrange(elv,srange,dhdr)
                dsdr=sqrt(1.-dhdr*dhdr)
                kntintv=kntintv+1
                kntintr=kntintr+1
                reflz=whigh*(w11*r11h+w21*r21h +                       &
                             w12*r12h+w22*r22h) +                      &
                       wlow*(w11*r11l+w21*r21l +                       &
                             w12*r12l+w22*r22l)
                refdbz=10.*alog10(reflz)
                refl(igate,jazim)=max(refdbz,refmin)
                r11l=max(refh(  iloc,  jloc,k-1),0.)
                r21l=max(refh(iloc+1,  jloc,k-1),0.)
                r12l=max(refh(  iloc,jloc+1,k-1),0.)
                r22l=max(refh(iloc+1,jloc+1,k-1),0.)
                r11h=max(refh(  iloc,  jloc,  k),0.)
                r21h=max(refh(iloc+1,  jloc,  k),0.)
                r12h=max(refh(  iloc,jloc+1,  k),0.)
                r22h=max(refh(iloc+1,jloc+1,  k),0.)
                reflh=whigh*(w11*r11h+w21*r21h +                       &
                             w12*r12h+w22*r22h) +                      &
                       wlow*(w11*r11l+w21*r21l +                       &
                             w12*r12l+w22*r22l)
                r11l=max(refv(  iloc,  jloc,k-1),0.)
                r21l=max(refv(iloc+1,  jloc,k-1),0.)
                r12l=max(refv(  iloc,jloc+1,k-1),0.)
                r22l=max(refv(iloc+1,jloc+1,k-1),0.)
                r11h=max(refv(  iloc,  jloc,  k),0.)
                r21h=max(refv(iloc+1,  jloc,  k),0.)
                r12h=max(refv(  iloc,jloc+1,  k),0.)
                r22h=max(refv(iloc+1,jloc+1,  k),0.)
                reflv=whigh*(w11*r11h+w21*r21h +                       &
                             w12*r12h+w22*r22h) +                      &
                       wlow*(w11*r11l+w21*r21l +                       &
                             w12*r12l+w22*r22l)
                r11l=max(refhv(  iloc,  jloc,k-1),0.)
                r21l=max(refhv(iloc+1,  jloc,k-1),0.)
                r12l=max(refhv(  iloc,jloc+1,k-1),0.)
                r22l=max(refhv(iloc+1,jloc+1,k-1),0.)
                r11h=max(refhv(  iloc,  jloc,  k),0.)
                r21h=max(refhv(iloc+1,  jloc,  k),0.)
                r12h=max(refhv(  iloc,jloc+1,  k),0.)
                r22h=max(refhv(iloc+1,jloc+1,  k),0.)
                reflhv=whigh*(w11*r11h+w21*r21h +                     &
                              w12*r12h+w22*r22h) +                    &
                        wlow*(w11*r11l+w21*r21l +                     &
                              w12*r12l+w22*r22l)
                IF(reflv > epsilon .AND. reflh > epsilon) THEN
                  rzdr(igate,jazim)=10.0*LOG10(reflh/reflv)
                  rzdr(igate,jazim)=min(max(rzdr(igate,jazim),-10.0),20.0)
                  rrhohv(igate,jazim)=reflhv/SQRT(reflh*reflv)
                ELSE
                  rzdr(igate,jazim)=0.
                  rrhohv(igate,jazim)=1.
                END IF
                r11l=kdp(  iloc,  jloc,k-1)
                r21l=kdp(iloc+1,  jloc,k-1)
                r12l=kdp(  iloc,jloc+1,k-1)
                r22l=kdp(iloc+1,jloc+1,k-1)
                r11h=kdp(  iloc,  jloc,  k)
                r21h=kdp(iloc+1,  jloc,  k)
                r12h=kdp(  iloc,jloc+1,  k)
                r22h=kdp(iloc+1,jloc+1,  k)
                rkdp(igate,jazim)=whigh*(w11*r11h+w21*r21h +           &
                                         w12*r12h+w22*r22h) +          &
                                   wlow*(w11*r11l+w21*r21l +           &
                                         w12*r12l+w22*r22l)
                upt=whigh*(w11*usc(  iloc,  jloc,k) +                  &
                              w21*usc(iloc+1,  jloc,k) +               &
                              w12*usc(  iloc,jloc+1,k) +               &
                              w22*usc(iloc+1,jloc+1,k))+               &
                        wlow*(w11*usc(  iloc,  jloc,k-1) +             &
                              w21*usc(iloc+1,  jloc,k-1) +             &
                              w12*usc(  iloc,jloc+1,k-1) +             &
                              w22*usc(iloc+1,jloc+1,k-1))
                vpt=whigh*(w11*vsc(  iloc,  jloc,k) +                  &
                              w21*vsc(iloc+1,  jloc,k) +               &
                              w12*vsc(  iloc,jloc+1,k) +               &
                              w22*vsc(iloc+1,jloc+1,k))+               &
                        wlow*(w11*vsc(  iloc,  jloc,k-1) +             &
                              w21*vsc(iloc+1,  jloc,k-1) +             &
                              w12*vsc(  iloc,jloc+1,k-1) +             &
                              w22*vsc(iloc+1,jloc+1,k-1))
                wpt=whigh*(w11*wsc(  iloc,  jloc,k) +                  &
                              w21*wsc(iloc+1,  jloc,k) +               &
                              w12*wsc(  iloc,jloc+1,k) +               &
                              w22*wsc(iloc+1,jloc+1,k))+               &
                        wlow*(w11*wsc(  iloc,  jloc,k-1) +             &
                              w21*wsc(iloc+1,  jloc,k-1) +             &
                              w12*wsc(  iloc,jloc+1,k-1) +             &
                              w22*wsc(iloc+1,jloc+1,k-1))
                IF(wrtuvwt > 0) THEN
                  urad(igate,jazim)=upt
                  vrad(igate,jazim)=vpt
                  wrad(igate,jazim)=wpt
                  tkrad(igate,jazim)=whigh*(w11*t(  iloc,  jloc,k) + &
                              w21*t(iloc+1,  jloc,k) +               &
                              w12*t(  iloc,jloc+1,k) +               &
                              w22*t(iloc+1,jloc+1,k))+               &
                        wlow*(w11*t(  iloc,  jloc,k-1) +             &
                              w21*t(iloc+1,  jloc,k-1) +             &
                              w12*t(  iloc,jloc+1,k-1) +             &
                              w22*t(iloc+1,jloc+1,k-1))
                END IF
!
                IF(refdbz > refmin ) THEN
!  Estimate terminal velocity (vt is positive down)
!
                  CALL termvel(refl(igate,jazim),hgtmsl,vt)
!
!  Calculate radial component
!
                  vr=dsdr*(xrat*upt+yrat*vpt)+dhdr*(wpt-vt)
!
!  Apply Doppler aliasing, if any
!
                  ifold=NINT(vr*twovninv)
                  radv(igate,jazim)=vr-ifold*twovnyq
                END IF
              END IF
!
!  Calculate weighted average values using antenna pattern gain.
!
              ibgn=max(1,(iloc1-iradius))
              iend=min((nx-1),(iloc1+iradius))
              jbgn=max(1,(jloc1-jradius))
              jend=min((ny-1),(jloc1+jradius))
!
!  First check to see of the entire observation volume area is clear.
!
              echo=.false.
              DO jpt=jbgn,jend
                DO ipt=ibgn,iend
                  IF(cmpref(ipt,jpt) > 0.) echo=.true.
                END DO
              END DO
!
!  If echo was found in 2d composite, proceed with processing 3d data.
!
              IF(echo) THEN
                kntr=0
                kntv=0
                sumwr=0.
                sumwv=0.
                sumv=0.
                sumr=0.
                sumrh=0.
                sumrv=0.
                sumrhv=0.
                sumkdp=0.
                DO jpt=jbgn,jend
                  DO ipt=ibgn,iend
                    DO kpt=2,nz-2
                      delv=elvsc(ipt,jpt,kpt)-elv
                      IF( abs(delv) < elradius .AND.                      &
                          refz(ipt,jpt,kpt) > 0. ) THEN
                        depth=0.5*(zps(ipt,jpt,kpt+1)-zps(ipt,jpt,kpt-1))
                        drng=rngsc(ipt,jpt,kpt)-srange
                        dazm=azmsc(ipt,jpt)-azimuth
                        IF(dazm > 180.) THEN
                          dazm=dazm-360.
                        ELSE IF( dazm < -180.) THEN
                          dazm=dazm+360.
                        END IF
                        wgtg=(drng*drng)*sr2inv  +                      &
                             (dazm*dazm)*wasqinv +                      &
                             (delv*delv)*wesqinv
                        wgtr=depth*exp(-wgtg)
                        kntr=kntr+1
                        sumwr=sumwr+wgtr
                        sumr  =  sumr+wgtr*refz(ipt,jpt,kpt)
                        sumrh = sumrh+wgtr*refh(ipt,jpt,kpt)
                        sumrv = sumrv+wgtr*refv(ipt,jpt,kpt)
                        sumrhv=sumrhv+wgtr*refhv(ipt,jpt,kpt)
                        sumkdp=sumkdp+wgtr*kdp(ipt,jpt,kpt)
                        IF( ref(ipt,jpt,kpt) > refmin ) THEN
                          wgtv=refz(ipt,jpt,kpt)*wgtr
                          kntv=kntv+1
                          sumwv=sumwv+wgtv
                          sumv=sumv+wgtv*vrsc(ipt,jpt,kpt)
                        END IF
                      END IF
                    END DO
                  END DO
                END DO
!
!  Solve for mean reflectivity
!
!               print *, ' igate,jazim,sumwr, kntr = ',               &
!                          igate,jazim,sumwr,kntr
                IF(sumwr > 0. .AND. kntr >= kntrmin ) THEN
                  kntavgr=kntavgr+1
                  reflz=sumr/sumwr
                  refdbz=10.*LOG10(reflz)
                  refl(igate,jazim)=max(refdbz,refmin)
                  reflh=sumrh/sumwr
                  reflv=sumrv/sumwr
                  IF(reflv > epsilon .AND. reflh > epsilon) THEN
                    rzdr(igate,jazim)=10.0*LOG10(reflh/reflv)
                    rzdr(igate,jazim)=min(max(rzdr(igate,jazim),-10.0),20.0)
                    rrhohv(igate,jazim)=reflhv/SQRT(reflh*reflv)
                  ELSE
                    rzdr(igate,jazim)=0.
                    rrhohv(igate,jazim)=1.
                  END IF
                  rkdp(igate,jazim)=sumkdp/sumwr
                END IF
!
!  Solve for mean radial velocity
!
!               print *, ' igate,jazim,sumwv, kntv = ',               &
!                          igate,jazim,sumwv,kntv
                IF(sumwv > 0. .AND. kntv >= kntvmin ) THEN
                  kntavgv=kntavgv+1
                  vr=sumv/sumwv
!
!  Apply Doppler aliasing, if any
!
                  ifold=NINT(vr*twovninv)
!                 radv(igate,jazim)=vr-ifold*twovnyq
                END IF
              ELSE ! no echo found in area
                kntclr=kntclr+1
                refl(igate,jazim)=0.
                rzdr(igate,jazim)=0.
                rrhohv(igate,jazim)=1.
                radv(igate,jazim)=rmisval
                IF(wrtuvwt > 0) THEN
                  upt=whigh*(w11*usc(  iloc,  jloc,k) +              &
                             w21*usc(iloc+1,  jloc,k) +              &
                             w12*usc(  iloc,jloc+1,k) +              &
                             w22*usc(iloc+1,jloc+1,k))+              &
                       wlow*(w11*usc(  iloc,  jloc,k-1) +            &
                             w21*usc(iloc+1,  jloc,k-1) +            &
                             w12*usc(  iloc,jloc+1,k-1) +            &
                             w22*usc(iloc+1,jloc+1,k-1))
                  vpt=whigh*(w11*vsc(  iloc,  jloc,k) +              &
                             w21*vsc(iloc+1,  jloc,k) +              &
                             w12*vsc(  iloc,jloc+1,k) +              &
                             w22*vsc(iloc+1,jloc+1,k))+              &
                       wlow*(w11*vsc(  iloc,  jloc,k-1) +            &
                             w21*vsc(iloc+1,  jloc,k-1) +            &
                             w12*vsc(  iloc,jloc+1,k-1) +            &
                             w22*vsc(iloc+1,jloc+1,k-1))
                  wpt=whigh*(w11*wsc(  iloc,  jloc,k) +              &
                             w21*wsc(iloc+1,  jloc,k) +              &
                             w12*wsc(  iloc,jloc+1,k) +              &
                             w22*wsc(iloc+1,jloc+1,k))+              &
                       wlow*(w11*wsc(  iloc,  jloc,k-1) +            &
                             w21*wsc(iloc+1,  jloc,k-1) +            &
                             w12*wsc(  iloc,jloc+1,k-1) +            &
                             w22*wsc(iloc+1,jloc+1,k-1))
                  urad(igate,jazim)=upt
                  vrad(igate,jazim)=vpt
                  wrad(igate,jazim)=wpt
                  tkrad(igate,jazim)=whigh*(w11*t(  iloc,  jloc,k) + &
                             w21*t(iloc+1,  jloc,k) +               &
                             w12*t(  iloc,jloc+1,k) +               &
                             w22*t(iloc+1,jloc+1,k))+               &
                       wlow*(w11*t(  iloc,  jloc,k-1) +             &
                             w21*t(iloc+1,  jloc,k-1) +             &
                             w12*t(  iloc,jloc+1,k-1) +             &
                             w22*t(iloc+1,jloc+1,k-1))
                END IF
              END IF
            ELSE ! outside vertical limits of grid
              kntvob=kntvob+1
            END IF
          ELSE ! outside grid domain
            knthob=knthob+1
          END IF
!
!  Impose physically realistic limits on dual-pol variables
!
          rzdr(igate,jazim)=min(max(rzdr(igate,jazim),-10.0),20.0)
          rkdp(igate,jazim)=min(max(rkdp(igate,jazim),-10.0),30.0)
          rrhohv(igate,jazim)=min(max(rrhohv(igate,jazim),0.0),1.0)

        END DO  ! igate
      END DO ! jazim
!
      ELSE  ! Integration method 2

      WRITE(6,'(a,i3)') 'Using integration method ',int_method

       DO jazim=1,nazim
        azimuth=azmbgn+irot*((jazim-1)*delazm)
        IF(azimuth > 360.) azimuth = azimuth-360.
        IF(azimuth < 0.) azimuth = azimuth+360.
!       IF(MOD((jazim-1),20) == 0) &
!         WRITE(6,'(a,f10.2)') ' Processing azimuth: ', azimuth
        azrad=deg2rad*azimuth
        azim(jazim)=azimuth
        beamw(jazim)=beamwid
        gtspc(jazim)=gatesp
        vnyq(jazim)=vnyqstl
        DO kpt=1,nptselv
        DO jpt=1,nptsazm
          blockfct(jpt,kpt)=1.0
        END DO
        END DO
        DO igate=1, ngate
          kntgate=kntgate+1
          srange=igate*gatesp
          CALL beamhgt(elv,srange,hgt,sfcrng)
          IF(mapproj > 0 ) THEN
            CALL gcircle(radlata,radlona,azimuth,sfcrng,gtlat,gtlon)
            CALL lltoxy(1,1,gtlat,gtlon,gtx,gty)
            delx=gtx-xrada
            dely=gty-yrada
            xrat=delx/srange
            yrat=dely/srange
            ipos=(gtx*dxinv)+1.5
            jpos=(gty*dyinv)+1.5
          ELSE
            delx=sin(azrad)*sfcrng
            dely=cos(azrad)*sfcrng
            xrat=delx/srange
            yrat=dely/srange
            ipos=irloc+(delx*dxinv)
            jpos=jrloc+(dely*dyinv)
          END IF
          iloc1=NINT(ipos)
          jloc1=NINT(jpos)
          iloc=INT(ipos)
          jloc=INT(jpos)

          IF(iloc > 0 .AND. iloc1 < nx-1 .AND.                        &
             jloc > 0 .AND. jloc1 < ny-1 ) THEN
            hgtmsl=hgt+radelv
            w2i=ipos-iloc
            w1i=1.0-w2i
            w2j=jpos-jloc
            w1j=1.0-w2j
            w11=w1i*w1j
            w21=w2i*w1j
            w12=w1i*w2j
            w22=w2i*w2j
            trngt=w11*zp(iloc,jloc,2) +          &
                  w21*zp(iloc+1,jloc,2) +        &
                  w12*zp(iloc,jloc+1,2) +        &
                  w22*zp(iloc+1,jloc+1,2)
            topgt=w11*zps(iloc,jloc,nzm2) +          &
                  w21*zps(iloc+1,jloc,nzm2) +        &
                  w12*zps(iloc,jloc+1,nzm2) +        &
                  w22*zps(iloc+1,jloc+1,nzm2)
            CALL beamhgt(elvtop,srange,tophgt,sfcrngt)
            CALL beamhgt(elvbot,srange,bothgt,sfcrngb)
            topmsl=tophgt+radelv
            botmsl=bothgt+radelv

            IF(botmsl < topgt .AND. topmsl > (trngt+hblkmin)) THEN

              DO k=3,nz-2
                IF(zps(iloc1,jloc1,k) > hgtmsl) EXIT
              END DO
!
!  Calculate vorticity by doing a linear interpolation directly to
!  the gate center location.
!
              whigh=(hgtmsl-zps(iloc1,jloc1,k-1))/                     &
                  (zps(iloc1,jloc1,k)-zps(iloc1,jloc1,k-1))
              wlow=1.-whigh
              vvort(igate,jazim)=whigh*(w11*vort(iloc,jloc,k) +      &
                                        w21*vort(iloc+1,jloc,k) +    &
                                        w12*vort(iloc,jloc+1,k) +    &
                                        w22*vort(iloc+1,jloc+1,k)) + &
                                  wlow*(w11*vort(iloc,jloc,k-1) +    &
                                        w21*vort(iloc+1,jloc,k-1) +  &
                                        w12*vort(iloc,jloc+1,k-1) +  &
                                        w22*vort(iloc+1,jloc+1,k-1))
              vvmin=min(vvmin,vvort(igate,jazim))
              vvmax=max(vvmax,vvort(igate,jazim))
              IF(wrtqx > 0) THEN
                qrrad(igate,jazim)=whigh*(w11*qr(iloc,jloc,k) +        &
                                          w21*qr(iloc+1,jloc,k) +      &
                                          w12*qr(iloc,jloc+1,k) +      &
                                          w22*qr(iloc+1,jloc+1,k)) +   &
                                    wlow*(w11*qr(iloc,jloc,k-1) +      &
                                          w21*qr(iloc+1,jloc,k-1) +    &
                                          w12*qr(iloc,jloc+1,k-1) +    &
                                          w22*qr(iloc+1,jloc+1,k-1))
                qrrad(igate,jazim)=max(0.,qrrad(igate,jazim))
                qsrad(igate,jazim)=whigh*(w11*qs(iloc,jloc,k) +        &
                                          w21*qs(iloc+1,jloc,k) +      &
                                          w12*qs(iloc,jloc+1,k) +      &
                                          w22*qs(iloc+1,jloc+1,k)) +   &
                                    wlow*(w11*qs(iloc,jloc,k-1) +      &
                                          w21*qs(iloc+1,jloc,k-1) +    &
                                          w12*qs(iloc,jloc+1,k-1) +    &
                                          w22*qs(iloc+1,jloc+1,k-1))
                qsrad(igate,jazim)=max(0.,qsrad(igate,jazim))
                qhrad(igate,jazim)=whigh*(w11*qh(iloc,jloc,k) +        &
                                          w21*qh(iloc+1,jloc,k) +      &
                                          w12*qh(iloc,jloc+1,k) +      &
                                          w22*qh(iloc+1,jloc+1,k)) +   &
                                    wlow*(w11*qh(iloc,jloc,k-1) +      &
                                          w21*qh(iloc+1,jloc,k-1) +    &
                                          w12*qh(iloc,jloc+1,k-1) +    &
                                          w22*qh(iloc+1,jloc+1,k-1))
                qhrad(igate,jazim)=max(0.,qhrad(igate,jazim))
              END IF
!
              IF(wrtuvwt > 0) THEN
                upt=whigh*(w11*usc(  iloc,  jloc,k) +              &
                           w21*usc(iloc+1,  jloc,k) +              &
                           w12*usc(  iloc,jloc+1,k) +              &
                           w22*usc(iloc+1,jloc+1,k))+              &
                     wlow*(w11*usc(  iloc,  jloc,k-1) +            &
                           w21*usc(iloc+1,  jloc,k-1) +            &
                           w12*usc(  iloc,jloc+1,k-1) +            &
                           w22*usc(iloc+1,jloc+1,k-1))
                urad(igate,jazim)=upt
                vpt=whigh*(w11*vsc(  iloc,  jloc,k) +              &
                           w21*vsc(iloc+1,  jloc,k) +              &
                           w12*vsc(  iloc,jloc+1,k) +              &
                           w22*vsc(iloc+1,jloc+1,k))+              &
                     wlow*(w11*vsc(  iloc,  jloc,k-1) +            &
                           w21*vsc(iloc+1,  jloc,k-1) +            &
                           w12*vsc(  iloc,jloc+1,k-1) +            &
                           w22*vsc(iloc+1,jloc+1,k-1))
                vrad(igate,jazim)=vpt
                wpt=whigh*(w11*wsc(  iloc,  jloc,k) +              &
                           w21*wsc(iloc+1,  jloc,k) +              &
                           w12*wsc(  iloc,jloc+1,k) +              &
                           w22*wsc(iloc+1,jloc+1,k))+              &
                     wlow*(w11*wsc(  iloc,  jloc,k-1) +            &
                           w21*wsc(iloc+1,  jloc,k-1) +            &
                           w12*wsc(  iloc,jloc+1,k-1) +            &
                           w22*wsc(iloc+1,jloc+1,k-1))
                wrad(igate,jazim)=wpt
                tkrad(igate,jazim)=whigh*(w11*t(  iloc,  jloc,k) + &
                           w21*t(iloc+1,  jloc,k) +               &
                           w12*t(  iloc,jloc+1,k) +               &
                           w22*t(iloc+1,jloc+1,k))+               &
                     wlow*(w11*t(  iloc,  jloc,k-1) +             &
                           w21*t(iloc+1,  jloc,k-1) +             &
                           w12*t(  iloc,jloc+1,k-1) +             &
                           w22*t(iloc+1,jloc+1,k-1))
              END IF
!  Check to see of the entire observation volume area is clear.
!
              bmwm=bmwrad*evalwid*srange
              sradius=max(bmwm,pulsewmn)
              iradius=1+INT(sradius/dx)
              jradius=1+INT(sradius/dy)
              ibgn=max(1,(iloc1-iradius))
              iend=min((nx-1),(iloc1+iradius))
              jbgn=max(1,(jloc1-jradius))
              jend=min((ny-1),(jloc1+jradius))
              echo=.false.
              DO jpt=jbgn,jend
                DO ipt=ibgn,iend
                  IF(cmpref(ipt,jpt) > 0.) THEN
                    echo=.true.
                    EXIT
                  END IF
                END DO
              END DO

              IF(echo .OR. (botmsl < trngt)) THEN
!
!  Calculate beam-weighted sums for an uniform array of points around this
!  range-gate using bi-quadratic horizontal and linear vertical interpolation
!
                kntr=0
                sumwr=0.
                sumr=0.
                sumrh=0.
                sumrv=0.
                sumrhv=0.
                sumkdp=0.
                kntv=0.
                sumv=0.
                sumv1=0.
                sumv2=0.
                sumwv=0.
                DO kpt=1,nptselv
                  elvpt=elv+delev(kpt)
!                  IF(igate == 100) print *, 'kpt, hgtmsl: ',kpt,hgtmsl
                DO jpt=1,nptsazm
                  azmpt=azimuth+dazim(jpt)
                  azmrpt=deg2rad*azmpt
                DO ipt=1,nptsrng
                  srgpt=(igate*gatesp)+delrng(ipt)
                  IF(srgpt <= 0.0) CYCLE
                  CALL beamhgt(elvpt,srgpt,hgtpt,sfcrpt)
                  hgtmsl=hgtpt+radelv
                  IF(mapproj > 0 ) THEN
                    CALL gcircle(radlata,radlona,azmpt,sfcrpt,gtlat,gtlon)
                    CALL lltoxy(1,1,gtlat,gtlon,gtx,gty)
                    delx=gtx-xrada
                    dely=gty-yrada
                    xrat=delx/sfcrpt
                    yrat=dely/sfcrpt
                    ipos=(gtx*dxinv)+1.5
                    jpos=(gty*dyinv)+1.5
                  ELSE
                    delx=sin(azmrpt)*sfcrpt
                    dely=cos(azmrpt)*sfcrpt
                    xrat=delx/sfcrpt
                    yrat=dely/sfcrpt
                    ipos=irloc+(delx*dxinv)
                    jpos=jrloc+(dely*dyinv)
                  END IF
                  iloc1=NINT(ipos)
                  jloc1=NINT(jpos)
                  iloc=INT(ipos)
                  jloc=INT(jpos)

                  IF(iloc > 0 .AND. iloc1 < nx-1 .AND.                &
                     jloc > 0 .AND. jloc1 < ny-1 ) THEN

                    ilc=max(min(iloc1,nx-2),2)
                    jlc=max(min(jloc1,ny-2),2)
                    xpt=gtx-xsc(ilc)
                    ypt=gty-ysc(jlc)

                    w2i=ipos-iloc
                    w1i=1.0-w2i
                    w2j=jpos-jloc
                    w1j=1.0-w2j
                    w11=w1i*w1j
                    w21=w2i*w1j
                    w12=w1i*w2j
                    w22=w2i*w2j
                    trnpt=w11*zp(iloc,jloc,2) +                        &
                            w21*zp(iloc+1,jloc,2) +                    &
                            w12*zp(iloc,jloc+1,2) +                    &
                            w22*zp(iloc+1,jloc+1,2)
                    hgtagl=hgtmsl-trnpt

                    IF(blockopt > 0 .AND. hgtagl < hblkmax) THEN
                      xmit=max(0.0,min(1.0,((hgtagl-hblkmin)*dhblkinv)))
                      gndatten=max(0.0,min(1.0,(blkcst*(1.-xmit))))
                      blockfct(jpt,kpt)=                               &
                        max(0.0,min(1.0,(blockfct(jpt,kpt)-gndatten)))
                    END IF

                    IF( hgtagl > hblkmin  .AND.                        &
                        hgtmsl < zps(iloc1,jloc1,nz-1) ) THEN
                      DO k=3,nz-2
                        IF(zps(iloc1,jloc1,k) > hgtmsl) EXIT
                      END DO
                      kntr=kntr+1
                      sumwr=sumwr+wgtpt(ipt,jpt,kpt)

                      whigh=(hgtmsl-zps(iloc1,jloc1,k-1))/             &
                          (zps(iloc1,jloc1,k)-zps(iloc1,jloc1,k-1))
                      wlow=1.-whigh

!                     CALL QUADINT(nx,ny,nz,refz,ilc,jlc,k,dx,dy,  &
!                                  xpt,ypt,whigh,wlow,reflz)
                      reflz=whigh*(w11*refz(  iloc,  jloc,  k)  +    &
                                   w21*refz(iloc+1,  jloc,  k)  +    &
                                   w12*refz(  iloc,jloc+1,  k)  +    &
                                   w22*refz(iloc+1,jloc+1,  k)) +    &
                             wlow*(w11*refz(  iloc,  jloc,k-1)  +    &
                                   w21*refz(iloc+1,  jloc,k-1)  +    &
                                   w12*refz(  iloc,jloc+1,k-1)  +    &
                                   w22*refz(iloc+1,jloc+1,k-1))
                      sumr=sumr+reflz*wgtpt(ipt,jpt,kpt)*blockfct(jpt,kpt)

!                     CALL QUADINT(nx,ny,nz,refh,ilc,jlc,k,dx,dy,  &
!                                  xpt,ypt,whigh,wlow,refhpt)
                      refhpt=whigh*(w11*refh(  iloc,  jloc,  k)  +    &
                                    w21*refh(iloc+1,  jloc,  k)  +    &
                                    w12*refh(  iloc,jloc+1,  k)  +    &
                                    w22*refh(iloc+1,jloc+1,  k)) +    &
                              wlow*(w11*refh(  iloc,  jloc,k-1)  +    &
                                    w21*refh(iloc+1,  jloc,k-1)  +    &
                                    w12*refh(  iloc,jloc+1,k-1)  +    &
                                    w22*refh(iloc+1,jloc+1,k-1))
                      sumrh=sumrh+refhpt*wgtpt(ipt,jpt,kpt)*blockfct(jpt,kpt)

!                     CALL QUADINT(nx,ny,nz,refv,ilc,jlc,k,dx,dy,  &
!                                  xpt,ypt,whigh,wlow,refvpt)
                      refvpt=whigh*(w11*refv(  iloc,  jloc,  k)  +    &
                                    w21*refv(iloc+1,  jloc,  k)  +    &
                                    w12*refv(  iloc,jloc+1,  k)  +    &
                                    w22*refv(iloc+1,jloc+1,  k)) +    &
                              wlow*(w11*refv(  iloc,  jloc,k-1)  +    &
                                    w21*refv(iloc+1,  jloc,k-1)  +    &
                                    w12*refv(  iloc,jloc+1,k-1)  +    &
                                    w22*refv(iloc+1,jloc+1,k-1))
                      sumrv=sumrv+refvpt*wgtpt(ipt,jpt,kpt)*blockfct(jpt,kpt)

!                     CALL QUADINT(nx,ny,nz,refhv,ilc,jlc,k,dx,dy,  &
!                                  xpt,ypt,whigh,wlow,refhvpt)
                      refhvpt=whigh*(w11*refhv(  iloc,  jloc,  k)  +    &
                                     w21*refhv(iloc+1,  jloc,  k)  +    &
                                     w12*refhv(  iloc,jloc+1,  k)  +    &
                                     w22*refhv(iloc+1,jloc+1,  k)) +    &
                               wlow*(w11*refhv(  iloc,  jloc,k-1)  +    &
                                     w21*refhv(iloc+1,  jloc,k-1)  +    &
                                     w12*refhv(  iloc,jloc+1,k-1)  +    &
                                     w22*refhv(iloc+1,jloc+1,k-1))
                      sumrhv=sumrhv+refhvpt*wgtpt(ipt,jpt,kpt)*blockfct(jpt,kpt)

!                     CALL QUADINT(nx,ny,nz,kdp,ilc,jlc,k,dx,dy,  &
!                                  xpt,ypt,whigh,wlow,kdppt)
                      kdppt=whigh*(w11*kdp(  iloc,  jloc,  k)  +    &
                                   w21*kdp(iloc+1,  jloc,  k)  +    &
                                   w12*kdp(  iloc,jloc+1,  k)  +    &
                                   w22*kdp(iloc+1,jloc+1,  k)) +    &
                             wlow*(w11*kdp(  iloc,  jloc,k-1)  +    &
                                   w21*kdp(iloc+1,  jloc,k-1)  +    &
                                   w12*kdp(  iloc,jloc+1,k-1)  +    &
                                   w22*kdp(iloc+1,jloc+1,k-1))
                      sumkdp=sumkdp+kdppt*wgtpt(ipt,jpt,kpt)*blockfct(jpt,kpt)

                      IF( reflz > 0. ) THEN
                        CALL dhdrange(elvpt,srgpt,dhdr)
                        dsdr=sqrt(1.-dhdr*dhdr)
                        refdbz=10.*alog10(reflz)
                        IF(refdbz > refmin ) THEN
!                         CALL QUADINT(nx,ny,nz,usc,ilc,jlc,k,dx,dy,  &
!                                    xpt,ypt,whigh,wlow,urad)
                          upt =whigh*(w11*usc(  iloc,  jloc,  k)  +    &
                                      w21*usc(iloc+1,  jloc,  k)  +    &
                                      w12*usc(  iloc,jloc+1,  k)  +    &
                                      w22*usc(iloc+1,jloc+1,  k)) +    &
                                wlow*(w11*usc(  iloc,  jloc,k-1)  +    &
                                      w21*usc(iloc+1,  jloc,k-1)  +    &
                                      w12*usc(  iloc,jloc+1,k-1)  +    &
                                      w22*usc(iloc+1,jloc+1,k-1))
!                         CALL QUADINT(nx,ny,nz,usc,ilc,jlc,k,dx,dy,  &
!                                    xpt,ypt,whigh,wlow,urad)
                          vpt =whigh*(w11*vsc(  iloc,  jloc,  k)  +    &
                                      w21*vsc(iloc+1,  jloc,  k)  +    &
                                      w12*vsc(  iloc,jloc+1,  k)  +    &
                                      w22*vsc(iloc+1,jloc+1,  k)) +    &
                                wlow*(w11*vsc(  iloc,  jloc,k-1)  +    &
                                      w21*vsc(iloc+1,  jloc,k-1)  +    &
                                      w12*vsc(  iloc,jloc+1,k-1)  +    &
                                      w22*vsc(iloc+1,jloc+1,k-1))
!                         CALL QUADINT(nx,ny,nz,usc,ilc,jlc,k,dx,dy,  &
!                                    xpt,ypt,whigh,wlow,urad)
                          wpt =whigh*(w11*wsc(  iloc,  jloc,  k)  +    &
                                      w21*wsc(iloc+1,  jloc,  k)  +    &
                                      w12*wsc(  iloc,jloc+1,  k)  +    &
                                      w22*wsc(iloc+1,jloc+1,  k)) +    &
                                wlow*(w11*wsc(  iloc,  jloc,k-1)  +    &
                                      w21*wsc(iloc+1,  jloc,k-1)  +    &
                                      w12*wsc(  iloc,jloc+1,k-1)  +    &
                                      w22*wsc(iloc+1,jloc+1,k-1))
!
!  Estimate terminal velocity (vt is positive down)
!
                          CALL termvel(refdbz,hgtmsl,vt)
!
!  Calculate radial component
!
                          vr=dsdr*(xrat*upt+yrat*vpt)+dhdr*(wpt-vt)
                          kntv=kntv+1
                          sumv=sumv+vr*wgtpt(ipt,jpt,kpt)* &
                                       blockfct(jpt,kpt)*reflz
                          sumwv=sumwv+wgtpt(ipt,jpt,kpt)*reflz
                          sumv1=sumv1+vr
                          sumv2=sumv2+vr*vr
                        END IF ! refdbz > refmin
                      END IF ! reflz > 0
                    ELSE IF(hgtagl <= hblkmin .AND. blockopt > 1) THEN
                      vr=0.
                      reflz=reflzgnd
                      kntr=kntr+1
                      sumr=sumr+reflz*wgtpt(ipt,jpt,kpt)*blockfct(jpt,kpt)
                      sumrh=sumrh+reflz*wgtpt(ipt,jpt,kpt)*blockfct(jpt,kpt)
                      sumrv=sumrv+reflz*wgtpt(ipt,jpt,kpt)*blockfct(jpt,kpt)
                      sumrhv=sumrhv+reflz*wgtpt(ipt,jpt,kpt)*blockfct(jpt,kpt)
                      sumkdp=sumkdp+gndkdp*wgtpt(ipt,jpt,kpt)*blockfct(jpt,kpt)
                      sumwr=sumwr+wgtpt(ipt,jpt,kpt)
                      kntv=kntv+1
                      sumv=sumv+vr*wgtpt(ipt,jpt,kpt)* &
                                   blockfct(jpt,kpt)*reflz
                      sumwv=sumwv+wgtpt(ipt,jpt,kpt)*reflz
                      sumv1=sumv1+vr
                      sumv2=sumv2+vr*vr
                    END IF ! hgt in grid
                  END IF ! iloc in grid
                END DO ! ipt
                END DO ! jpt
                END DO ! kpt
!
!  Calculate reflectivity and velocity from the uniform
!  array center on this range gate.
!
                IF(sumwr > epsilon .AND. kntr > kntrmin) THEN
                  kntintr=kntintr+1
                  reflz=sumr/sumwr
                  IF(reflz > epsilon) THEN
                    refdbz=10.*LOG10(reflz)
                    refl(igate,jazim)=min(max(refdbz,refmin),refmax)
                  ELSE
                    refdbz=refmin
                  END IF

                  reflh=sumrh/sumwr
                  reflv=sumrv/sumwr
                  IF(reflh > epsilon .AND. reflv > epsilon) THEN
                    rzdr(igate,jazim)=10.0*LOG10(reflh/reflv)
                    rrhohv(igate,jazim)=(sumrhv/sumwr)/SQRT(reflh*reflv)
                  ELSE
                    rzdr(igate,jazim)=0.
                    rrhohv(igate,jazim)=1.
                  END IF
                  rkdp(igate,jazim)=sumkdp/sumwr
                ELSE
                  kntfewr=kntfewr+1
                  refl(igate,jazim)=0.
                  rzdr(igate,jazim)=0.
                  rrhohv(igate,jazim)=1.0
                  rkdp(igate,jazim)=0.
                END IF
                IF(sumwv > epsilon .AND. kntv > kntvmin) THEN
                  kntintv=kntintv+1
                  vr=sumv/sumwv
                  ifold=NINT(vr*twovninv)
                  radv(igate,jazim)=vr-ifold*twovnyq
                  flkntv=float(kntv)
                  varvr=max(0.01,                                      &
                        ((sumv2-(sumv1*sumv1/flkntv))/(flkntv-1.)))
                  stdvr(igate,jazim)=stdvrmul*SQRT(varvr)
                ELSE
                  kntfewv=kntfewv+1
                  radv(igate,jazim)=rmisval
                END IF
              ELSE ! no echo
                kntclr=kntclr+1
                refl(igate,jazim)=0.
                rzdr(igate,jazim)=0.
                rrhohv(igate,jazim)=1.0
                rkdp(igate,jazim)=0.
                radv(igate,jazim)=rmisval
              END IF
            ELSE IF ( botmsl > topgt) THEN ! beam completely above grid
              kntvob=kntvob+(ngate-igate)+1
              kntgate=kntgate+(ngate-igate)
              EXIT
            ELSE IF ( topmsl < (trnpt+hblkmin) ) THEN ! beam completely below ground
              print *, ' Entire beam below ground: ',azimuth,igate
              IF(blockopt > 1) THEN
                sumr=0.
                sumkdp=0.
                sumwr=0.
                DO kpt=1,nptselv
                DO jpt=1,nptsazm
                  reflz=reflzgnd
                  sumr=sumr+reflz*wgtpt(ipt,jpt,kpt)*blockfct(jpt,kpt)
                  sumkdp=sumkdp+gndkdp*wgtpt(ipt,jpt,kpt)*blockfct(jpt,kpt)
                  sumwr=sumwr+wgtpt(ipt,jpt,kpt)
                END DO
                END DO
                IF(sumwr > 0.) THEN
                  reflz=sumr/sumwr
                  refdbz=10.*alog10(reflz)
                  refl(igate,jazim)=min(max(refdbz,refmin),refmax)
                  rzdr(igate,jazim)=0.
                  rrhohv(igate,jazim)=1.0
                  rkdp(igate,jazim)=sumkdp/sumwr
                ELSE
                  refl(igate,jazim)=0.
                  rzdr(igate,jazim)=0.
                  rrhohv(igate,jazim)=1.0
                  rkdp(igate,jazim)=0.
                END IF
              ELSE
                refl(igate,jazim)=0.
                rzdr(igate,jazim)=0.
                rrhohv(igate,jazim)=1.0
                rkdp(igate,jazim)=0.
              END IF
              radv(igate,jazim)=0.
              kntbob=kntbob+(ngate-igate)+1
              kntgate=kntgate+(ngate-igate)
              DO iigate=igate+1,ngate
                radv(iigate,jazim)=0.
                refl(iigate,jazim)=0.
                rzdr(igate,jazim)=0.
                rrhohv(igate,jazim)=1.0
                rkdp(igate,jazim)=0.
              END DO
              EXIT
            END IF
          ELSE ! outside grid domain
            knthob=knthob+1
          END IF
          rzdr(igate,jazim)=min(max(rzdr(igate,jazim),-10.0),20.0)
          rkdp(igate,jazim)=min(max(rkdp(igate,jazim),-30.0),30.0)
          rrhohv(igate,jazim)=min(max(rrhohv(igate,jazim),0.0),1.0)
        END DO  ! igate
        IF(mod(jazim,10) == 0) THEN
          reflmin=999.
          reflmax=-999.
          DO igate=1,ngate
            IF(refl(igate,jazim) > refchk) THEN
              reflmin=min(reflmin,refl(igate,jazim))
              reflmax=max(reflmax,refl(igate,jazim))
            END IF
          END DO
          print *, ' Azimuth,reflmin,reflmax=',azimuth,reflmin,reflmax
        END IF
      END DO ! jazim

      print *, ' Done with jazim loop'

      END IF ! integration method
!
!-----------------------------------------------------------------------
!
!   Calculate attenuation for this tilt
!
!-----------------------------------------------------------------------
!
      WRITE(6,'(a,i6)') ' Attenopt block, attenopt= ',attenopt
      IF (attenopt > 0) THEN
        CALL xbatten(ngate,maxazim,nazim,                              &
                     refl,gatesp,rmisval,refmin,attrefl)
      ELSE
        DO jazim=1,nazim
          DO igate=1,ngate
            attrefl(igate,jazim)=refl(igate,jazim)
          END DO
        END DO
      END IF
!
!-----------------------------------------------------------------------
!
!   Apply noise to the reflectivity and velocity according the
!   indicated option.
!
!-----------------------------------------------------------------------
!
      WRITE(6,'(a,i4)') ' Radial vel error, vrerropt= ',vrerropt
      IF(vrerropt > 0 )                                                &
        CALL radverr(ngate,maxazim,nazim,radv,attrefl,stdvr,gatesp,    &
                   vrerropt,sigmavr,samploi,vnyqstl,                   &
                   wavelen,pwrxmt,gaindb,lossdb,noisedbm,beamwid,      &
                   pulselmn,npulsmn,rmisval)
      WRITE(6,'(a,i4)') ' Reflectivity error, rferropt= ',rferropt
      IF(rferropt > 0 ) THEN
        WRITE(6,'(a,f10.2)') ' Applying sigmarf: ',sigmarf
        CALL reflerr(ngate,maxazim,nazim,refl,attrefl,gatesp,          &
                     rferropt,sigmarf,refmin,samploi,                  &
                     wavelen,pwrxmt,gaindb,lossdb,noisedbm,beamwid,    &
                     pulselmn,npulsmn,rmisval)
        DO jazim=1,nazim
          DO igate=1,ngate
            err=sigmarf*randnor()
            rzdr(i,j)=min(max((rzdr(i,j)+err),-10.),20.)
          END DO
        END DO
      END IF
!
!-----------------------------------------------------------------------
!
!   Apply X-band radar sensitivity thresholding to the
!   attenuated velocity using the attenuated reflectivity.
!
!-----------------------------------------------------------------------
!
      IF(senstvopt > 0) THEN
        DO jazim=1,nazim
          DO igate=1,ngate
            IF(attrefl(igate,jazim) > sens(igate)) THEN
              attradv(igate,jazim)=radv(igate,jazim)
            ELSE
              attrefl(igate,jazim)=rmisval
              attradv(igate,jazim)=rmisval
            END IF
          END DO
        END DO
      ELSE
        DO jazim=1,nazim
          DO igate=1,ngate
            attradv(igate,jazim)=radv(igate,jazim)
          END DO
        END DO
      END IF
!
!-----------------------------------------------------------------------
!
!   Report statistics
!
!-----------------------------------------------------------------------
!
      denom=1./float(kntgate)
      write(6,'(a,i9)')         ' Total gates:           ',kntgate
      write(6,'(a,i9,f7.1,a)') ' Out-of-bounds Horiz:    ',knthob,     &
                                (100.*float(knthob)*denom),' percent'
      write(6,'(a,i9,f7.1,a)') ' Above grid points:      ',kntvob,     &
                                (100.*float(kntvob)*denom),' percent'
      write(6,'(a,i9,f7.1,a)') ' Ground-blocked points:  ',kntbob,     &
                                (100.*float(kntbob)*denom),' percent'
      write(6,'(a,i9,f7.1,a)') ' Gates no echo:          ',kntclr,     &
                                (100.*float(kntclr)*denom),' percent'
      write(6,'(a,i9,f7.1,a)') ' Gates interp refl:      ',kntintr,    &
                                (100.*float(kntintr)*denom),' percent'
      write(6,'(a,i9,f7.1,a)') ' Gates interp rvel:      ',kntintv,    &
                                (100.*float(kntintv)*denom),' percent'
      write(6,'(a,i9,f7.1,a)') ' Gates too few pts refl: ',kntfewr,    &
                                (100.*float(kntfewr)*denom),' percent'
      write(6,'(a,i9,f7.1,a)') ' Gates too few pts rvel: ',kntfewv,    &
                                (100.*float(kntfewv)*denom),' percent'
      write(6,'(a,f16.8,a,f16.8)')                                     &
      ' Tilt min vert vort (*10**3)=',vvmin,'  Max vert vort = ',vvmax

      print *, ' Stats by azimuth, refl:'
      DO jazim=1,nazim,40
        reflmin=999.
        reflmax=-999.
        DO igate=1,ngate
          IF(refl(igate,jazim) > refchk) THEN
            reflmin=min(reflmin,refl(igate,jazim))
            reflmax=max(reflmax,refl(igate,jazim))
          END IF
        END DO
        print *, ' Azimuth,reflmin,reflmax=',azim(jazim),reflmin,reflmax
      END DO

      print *, ' Stats by azimuth, Attenuted refl:'
      DO jazim=1,nazim,40
        reflmin=999.
        reflmax=-999.
        DO igate=1,ngate
          reflmin=min(reflmin,attrefl(igate,jazim))
          reflmax=max(reflmax,attrefl(igate,jazim))
        END DO
        print *, ' Azimuth,reflmin,reflmax=',azim(jazim),reflmin,reflmax
      END DO

      IF(dualpol > 0 .AND. wrtdualp > 0) THEN

        print *, ' Stats by azimuth, zdr:'
        DO jazim=1,nazim,40
          tstmin=9999.
          tstmax=-9999.
          DO igate=1,ngate
            IF(rzdr(igate,jazim) > -90.) THEN
              tstmin=min(tstmin,rzdr(igate,jazim))
              tstmax=max(tstmax,rzdr(igate,jazim))
            END IF
          END DO
          print *, ' Azimuth,zdrmin,zdrmax=',azim(jazim),tstmin,tstmax
        END DO

        print *, ' Stats by azimuth, kdp:'
        DO jazim=1,nazim,40
          tstmin=9999.
          tstmax=-9999.
          DO igate=1,ngate
            IF(rkdp(igate,jazim) > -90.) THEN
              tstmin=min(tstmin,rkdp(igate,jazim))
              tstmax=max(tstmax,rkdp(igate,jazim))
            END IF
          END DO
          print *, ' Azimuth,kdpmin,kdpmax=',azim(jazim),tstmin,tstmax
        END DO

        print *, ' Stats by azimuth, rhohv:'
        DO jazim=1,nazim,40
          tstmin=9999.
          tstmax=-9999.
          DO igate=1,ngate
            IF(rrhohv(igate,jazim) > -1.) THEN
              tstmin=min(tstmin,rrhohv(igate,jazim))
              tstmax=max(tstmax,rrhohv(igate,jazim))
            END IF
          END DO
          print *, ' Azimuth,rhohvmin,rhohvmax=',azim(jazim),tstmin,tstmax
        END DO

      END IF


      CALL abss2ctim(itimtilt,iyear,imon,iday,ihour,imin,isec)
      itimcdf=itimtilt-itim1970
!
!-----------------------------------------------------------------------
!
!   Write tilt netCDF file - Reflectivity
!
!-----------------------------------------------------------------------
!
      IF( ifmt == 1) THEN
        itimvol(itilt)=itimtilt
        vnyqvol(itilt)=vnyqstl
        DO igate=1,ngate
          rngvol(igate,itilt)=igate*gatesp
        END DO
        DO jazim=1,nazim
          elvvol(jazim,itilt)=elvang(itilt)
          azmvol(jazim,itilt)=azim(jazim)
        END DO
        DO jazim=1,nazim
          DO igate=1,ngate
            refvol(igate,jazim,itilt)=attrefl(igate,jazim)
            velvol(igate,jazim,itilt)=attradv(igate,jazim)
            uarefvol(igate,jazim,itilt)=refl(igate,jazim)
            uavelvol(igate,jazim,itilt)=radv(igate,jazim)
            vorvol(igate,jazim,itilt)=vvort(igate,jazim)
          END DO
        END DO
        IF(wrtqx > 0) THEN
          qrmin=999.
          qrmax=-999.
          qsmin=999.
          qsmax=-999.
          qhmin=999.
          qhmax=-999.
          DO jazim=1,nazim
            DO igate=1,ngate
              qrvol(igate,jazim,itilt)=qrrad(igate,jazim)
              qsvol(igate,jazim,itilt)=qsrad(igate,jazim)
              qhvol(igate,jazim,itilt)=qhrad(igate,jazim)
              qrmin=min(qrmin,qrrad(igate,jazim))
              qsmin=min(qsmin,qsrad(igate,jazim))
              qhmin=min(qhmin,qhrad(igate,jazim))
              qrmax=max(qrmax,qrrad(igate,jazim))
              qsmax=max(qsmax,qsrad(igate,jazim))
              qhmax=max(qhmax,qhrad(igate,jazim))
            END DO
          END DO
          print *, ' qr range (g/kg): ',(1000.0*qrmin),' to ',(1000.0*qrmax)
          print *, ' qs range (g/kg): ',(1000.0*qsmin),' to ',(1000.0*qsmax)
          print *, ' qh range (g/kg): ',(1000.0*qhmin),' to ',(1000.0*qhmax)
        END IF
        IF(dualpol > 0) THEN
          DO jazim=1,nazim
            DO igate=1,ngate
              zdrvol(igate,jazim,itilt)=rzdr(igate,jazim)
              kdpvol(igate,jazim,itilt)=rkdp(igate,jazim)
              rhohvvol(igate,jazim,itilt)=rrhohv(igate,jazim)
            END DO
          END DO
        END IF
        IF(wrtuvwt > 0) THEN
          DO jazim=1,nazim
            DO igate=1,ngate
              uvol(igate,jazim,itilt)=urad(igate,jazim)
              vvol(igate,jazim,itilt)=vrad(igate,jazim)
              wvol(igate,jazim,itilt)=wrad(igate,jazim)
              tkvol(igate,jazim,itilt)=tkrad(igate,jazim)
            END DO
          END DO
        END IF
      ELSE IF( ifmt == 2 ) THEN
!
!-----------------------------------------------------------------------
!
!  Write Reflectivity
!
!-----------------------------------------------------------------------
!
        initime=itimcdf-ifsecs
        write(fname,'(a,i2.2,a,i2.2,a,i4.4,2(i2.2),a,3(i2.2),a)')      &
         'Reflectivity_',ielv,'.',felv,'_',iyear,imon,iday,'-',        &
                     ihour,imin,isec,'.netcdf'
        varname='Reflectivity'
        CALL wtrftiltcdf(ngate,maxazim,nazim,fname,outdir,varname,     &
                       radname,radlata,radlona,radelv,vcp,elv,         &
                       rmisval,rngfval,itimtilt,frtime,initime,        &
                       vnyqstl,rfrgate,                                &
                       azim,beamw,gtspc,attrefl)
!
!-----------------------------------------------------------------------
!
!  Write index record
!
!-----------------------------------------------------------------------
!
        IF( creidx > 0 ) THEN
         write(idxunit,'(a)') '<item>'
         write(idxunit,'(a,f8.6,a,i10,a)')                             &
            '<time fractional="',frtime,'">',itimcdf,'</time>'
         write(idxunit,'(a,a,a,a)') '<params>netcdf {indexlocation} ', &
             TRIM(fname),TRIM(cmprext),'</params>'
         write(idxunit,'(a,i4.4,2i2.2,a,3i2.2,1x,a,1x,i2.2,a1,i2.2,a)')&
            '<selections>',iyear,imon,iday,'-',                        &
            ihour,imin,isec,'Reflectivity',ielv,'.',felv,              &
            '</selections>'
         write(idxunit,'(a)') '</item>'
        END IF ! creidx
!
!-----------------------------------------------------------------------
!
!  Write tilt netCDF file - Velocity
!
!-----------------------------------------------------------------------
!
        write(fname,'(a,i2.2,a,i2.2,a,i4.4,2(i2.2),a,3(i2.2),a)')      &
         'Velocity_',ielv,'.',felv,'_',iyear,imon,iday,'-',            &
                     ihour,imin,isec,'.netcdf'
        varname='RadialVelocity'
        CALL wtvrtiltcdf(ngate,maxazim,nazim,fname,outdir,varname,     &
                         radname,radlat,radlon,radelv,vcp,elv,         &
                         rmisval,rngfval,itimcdf,frtime,initime,       &
                         vnyqstl,rfrgate,                              &
                         azim,beamw,gtspc,vnyq,attradv)
!
!-----------------------------------------------------------------------
!
!  Write index record
!
!-----------------------------------------------------------------------
!
        IF( creidx > 0 ) THEN
         write(idxunit,'(a)') '<item>'
         write(idxunit,'(a,f8.6,a,i10,a)')                             &
           '<time fractional="',frtime,'">',itimcdf,'</time>'
         write(idxunit,'(a,i4.4,2i2.2,a,3i2.2,1x,a,1x,i2.2,a1,i2.2,a)')&
            '<selections>',iyear,imon,iday,'-',                        &
            ihour,imin,isec,'Velocity',ielv,'.',felv,                  &
            '</selections>'
         write(idxunit,'(a)') '</item>'
        END IF ! creidx
!
        write(fname,'(a,i2.2,a,i2.2,a,i4.4,2(i2.2),a,3(i2.2),a)')      &
         'ReflectUnatten_',ielv,'.',felv,'_',iyear,imon,iday,'-',      &
                     ihour,imin,isec,'.netcdf'
        varname='UnattenuatedReflectivity'
        CALL wtrftiltcdf(ngate,maxazim,nazim,fname,outdir,varname,     &
                         radname,radlata,radlona,radelv,vcp,elv,       &
                         rmisval,rngfval,itimcdf,frtime,initime,       &
                         vnyqstl,rfrgate,                              &
                         azim,beamw,gtspc,refl)
!
!-----------------------------------------------------------------------
!
!  Write index record
!
!-----------------------------------------------------------------------
!
        IF( creidx > 0 ) THEN
         write(idxunit,'(a)') '<item>'
         write(idxunit,'(a,f8.6,a,i10,a)')                             &
            '<time fractional="',frtime,'">',itimcdf,'</time>'
         write(idxunit,'(a,a,a,a)') '<params>netcdf {indexlocation} ', &
             TRIM(fname),TRIM(cmprext),'</params>'
         write(idxunit,'(a,i4.4,2i2.2,a,3i2.2,1x,a,1x,i2.2,a1,i2.2,a)')&
            '<selections>',iyear,imon,iday,'-',                        &
            ihour,imin,isec,'Reflectivity',ielv,'.',felv,              &
            '</selections>'
         write(idxunit,'(a)') '</item>'
        END IF ! creidx
!
!-----------------------------------------------------------------------
!
!  Write tilt netCDF file - Velocity
!
!-----------------------------------------------------------------------
!
        write(fname,'(a,i2.2,a,i2.2,a,i4.4,2(i2.2),a,3(i2.2),a)')      &
         'VelocityUnatten_',ielv,'.',felv,'_',iyear,imon,iday,'-',     &
                     ihour,imin,isec,'.netcdf'
        varname='UnattenuatedVelocity'
        CALL wtvrtiltcdf(ngate,maxazim,nazim,fname,outdir,varname,     &
                         radname,radlat,radlon,radelv,vcp,elv,         &
                         rmisval,rngfval,itimcdf,frtime,initime,       &
                         vnyqstl,rfrgate,                              &
                         azim,beamw,gtspc,vnyq,radv)
!
!-----------------------------------------------------------------------
!
!  Write index record
!
!-----------------------------------------------------------------------
!
        IF( creidx > 0 ) THEN
         write(idxunit,'(a)') '<item>'
         write(idxunit,'(a,f8.6,a,i10,a)') &
           '<time fractional="',frtime,'">',itimcdf,'</time>'
         write(idxunit,'(a,i4.4,2i2.2,a,3i2.2,1x,a,1x,i2.2,a1,i2.2,a)')&
            '<selections>',iyear,imon,iday,'-',                        &
            ihour,imin,isec,'Velocity',ielv,'.',felv,                  &
            '</selections>'
         write(idxunit,'(a)') '</item>'
        END IF ! creidx
!-----------------------------------------------------------------------
!
!  Write tilt netCDF file - Vorticity
!
!-----------------------------------------------------------------------
!
        write(fname,'(a,i2.2,a,i2.2,a,i4.4,2(i2.2),a,3(i2.2),a)')      &
         'Vorticity_',ielv,'.',felv,'_',iyear,imon,iday,'-',           &
                     ihour,imin,isec,'.netcdf'
        varname='VerticalVorticity'
        CALL wtvvtiltcdf(ngate,maxazim,nazim,fname,outdir,varname,     &
                         radname,radlat,radlon,radelv,vcp,elv,         &
                         rmisval,rngfval,itimcdf,frtime,initime,       &
                         vnyqstl,rfrgate,                              &
                         azim,beamw,gtspc,vvort)
!
!-----------------------------------------------------------------------
!
!  Write index record
!
!-----------------------------------------------------------------------
!
        IF( creidx > 0 ) THEN
         write(idxunit,'(a)') '<item>'
         write(idxunit,'(a,f8.6,a,i10,a)') &
           '<time fractional="',frtime,'">',itimcdf,'</time>'
         write(idxunit,'(a,i4.4,2i2.2,a,3i2.2,1x,a,1x,i2.2,a1,i2.2,a)')&
            '<selections>',iyear,imon,iday,'-',                        &
            ihour,imin,isec,'Vorticity',ielv,'.',felv,                 &
            '</selections>'
         write(idxunit,'(a)') '</item>'
        END IF ! creidx

        IF(dualpol > 0) THEN
!
!-----------------------------------------------------------------------
!
!  Write Reflectivity
!
!-----------------------------------------------------------------------
!
        write(fname,'(a,i2.2,a,i2.2,a,i4.4,2(i2.2),a,3(i2.2),a)')      &
         'Reflectivity_H_',ielv,'.',felv,'_',iyear,imon,iday,'-',        &
                     ihour,imin,isec,'.netcdf'
        varname='Reflectivity_h'
        CALL wtrftiltcdf(ngate,maxazim,nazim,fname,outdir,varname,     &
                       radname,radlata,radlona,radelv,vcp,elv,         &
                       rmisval,rngfval,itimtilt,frtime,initime,        &
                       vnyqstl,rfrgate,                                &
                       azim,beamw,gtspc,reflh)
!
!-----------------------------------------------------------------------
!
!  Write index record
!
!-----------------------------------------------------------------------
!
        IF( creidx > 0 ) THEN
         write(idxunit,'(a)') '<item>'
         write(idxunit,'(a,f8.6,a,i10,a)')                             &
            '<time fractional="',frtime,'">',itimcdf,'</time>'
         write(idxunit,'(a,a,a,a)') '<params>netcdf {indexlocation} ', &
             TRIM(fname),TRIM(cmprext),'</params>'
         write(idxunit,'(a,i4.4,2i2.2,a,3i2.2,1x,a,1x,i2.2,a1,i2.2,a)')&
            '<selections>',iyear,imon,iday,'-',                        &
            ihour,imin,isec,'Reflectivity_h',ielv,'.',felv,              &
            '</selections>'
         write(idxunit,'(a)') '</item>'
        END IF ! creidx
        END IF
      END IF ! ifmt=2

      IF(tmadvopt == 2) itimtilt=itimtilt+NINT(azmsect/rotrate)
!
    END DO  ! itilt

    print *, ' Out of tilt loop'
    print *, ' ifmt =',ifmt

    IF( ifmt == 1 ) THEN
      write(fnamevol,'(4a,i4.4,i2.2,i2.2,a,i2.2,i2.2,i2.2,a)')         &
         TRIM(outdir),'/',TRIM(radname),'_',year,month,day,'_',        &
                     hour,minute,second,'.vol'
      print *, ' fnamevol =',TRIM(fnamevol)
!
      print *, ' calling wtradvol'
      print *, ' ngate,maxazim,nelev =',ngate,maxazim,nelev
      CALL wtradvol(fnamevol,ngate,ngate,maxazim,nelev,                &
                rngvol,azmvol,elvvol,velvol,uavelvol,                  &
                vorvol,qrvol,qsvol,qhvol,uvol,vvol,wvol,tkvol,         &
                vnyqvol,itimvol,rngvol,azmvol,elvvol,refvol,uarefvol,  &
                zdrvol,kdpvol,rhohvvol,                                &
                wrtuaref,wrtuavel,wrtvort,wrtdualp,wrtqx,wrtuvwt,      &
                runname,radname,vcp,curtim,beamwid,                    &
                rmisval,rngfval,radelv,radlat,radlon)
    END IF
!
!-----------------------------------------------------------------------
!
!   Release memory for gridded data
!
!-----------------------------------------------------------------------
    WRITE(6,'(a)') ' Deallocating grid memory'
!
    deallocate(usc)
    deallocate(vsc)
    deallocate(wsc)

    deallocate(qr)
    deallocate(qs)
    deallocate(qh)
    deallocate(t)

    deallocate(latsc)
    deallocate(lonsc)
    deallocate(azmsc)
    deallocate(sfcr)
    deallocate(elvsc)
    deallocate(rngsc)
    deallocate(vort)
    deallocate(vrsc)

    deallocate(ref)
    deallocate(refz)
    deallocate(cmpref)

    deallocate(refh)
    deallocate(refv)
    deallocate(refhv)
    deallocate(kdp)
!
!   Release memory for polar coordinate radar data
!
    WRITE(6,'(a)') ' Deallocating radar tilt memory'

    deallocate(azim)
    deallocate(beamw)
    deallocate(gtspc)
    deallocate(vnyq)
    deallocate(radv)
    deallocate(attradv)
    deallocate(refl)
    deallocate(attrefl)
    deallocate(stdvr)
    deallocate(vvort)
    deallocate(sens)
    deallocate(rzdr)
    deallocate(rrhohv)
    deallocate(rkdp)
    IF(wrtqx > 0) THEN
      deallocate(qrrad)
      deallocate(qsrad)
      deallocate(qhrad)
    END IF
    IF(wrtuvwt > 0) THEN
      deallocate(urad)
      deallocate(vrad)
      deallocate(wrad)
      deallocate(tkrad)
    END IF

    WRITE(6,'(a)') ' Deallocating radar volume memory'

    IF( ifmt == 1 ) THEN
      deallocate(itimvol)
      deallocate(vnyqvol)
      deallocate(rngvol)
      deallocate(azmvol)
      deallocate(elvvol)
      deallocate(refvol)
      deallocate(uarefvol)
      deallocate(velvol)
      deallocate(uavelvol)
      deallocate(vorvol)
      deallocate(zdrvol)
      deallocate(kdpvol)
      deallocate(rhohvvol)
      deallocate(qrvol)
      deallocate(qsvol)
      deallocate(qhvol)
      deallocate(uvol)
      deallocate(vvol)
      deallocate(wvol)
      deallocate(tkvol)
    END IF

    WRITE(6,'(a)') ' End of file loop'
!
  END DO ! ifile loop
!
! Finish index file and close it out
!
  IF( creidx > 0  .AND. ifmt == 2 ) THEN
    write(idxunit,'(a)') '</codeindex>'
    CLOSE(idxunit)
    CALL retunit(idxunit)
  END IF

  STOP

END PROGRAM rademul

SUBROUTINE termvel(refl,hgtmsl,vt)
!
!-----------------------------------------------------------------------
!
!  Calculates the precipitation terminal velocity from
!  the observed radial velocity.  From Zeigler, 1978.
!
!  An empirical formula based on the reflectivity is
!  used to estimate the terminal velocity.
!
!  For simplicity, atmospheric density is approxmated
!  by a standard atmosphere, rho(z)=rhonot exp(-z/h0)
!
!-----------------------------------------------------------------------
!
!
  IMPLICIT NONE
  REAL :: refl
  REAL :: hgtmsl
  REAL :: vt
!
  REAL, PARAMETER :: zfrez=3000.
  REAL, PARAMETER :: zice=8000.
  REAL, PARAMETER :: rho0=1.2250
  REAL, PARAMETER :: h0=7000.
  REAL, PARAMETER :: denom=(1./(zice-zfrez))
  REAL, PARAMETER :: dbzmin=10.
  REAL, PARAMETER :: dbzmax=100.
!
  REAL :: refz,rhofact,s1,s2
!
  refz=10.**(0.1*refl)
  rhofact=EXP(0.4*hgtmsl/h0)
  IF(hgtmsl < zfrez) THEN
    vt=2.6*(refz**0.107)*rhofact
  ELSE IF(hgtmsl < zice) THEN
    s1=(zice-hgtmsl)*denom
    s2=2.*(hgtmsl-zfrez)*denom
    vt=s1*2.6*(refz**0.107)*rhofact + s2
  ELSE
    vt=2.0
  END IF
  RETURN
END SUBROUTINE termvel
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE WTRADVOL                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
  SUBROUTINE wtradvol(fnamevol,maxgatevel,maxgateref,maxazim,maxelev,  &
               rngvvol,azmvvol,elvvvol,velvol,uavelvol,                &
               vorvol,qrvol,qsvol,qhvol,uvol,vvol,wvol,tkvol,          &
               vnyqvol,itimvol,rngrvol,azmrvol,elvrvol,refvol,uarefvol,&
               zdrvol,kdpvol,rhohvvol,                                 &
               wrtuaref,wrtuavel,wrtvort,wrtdualp,wrtqx,wrtuvwt,       &
               runname,radname,vcp,curtim,beamwid,                     &
               rmisval,rngfval,radar_alt,radar_lat,radar_lon)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write out radar observations in radar coordinate
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Hu, CAPS
!
!  MODIFICATION HISTORY:
!    Converted to subroutine WTRADVOL, Keith Brewster
!
!    11/13/2004  Keith Brewster, CAPS
!    Updated to write more variables for use by radbin2cdf
!
!    05/24/2005  Keith Brewster, CAPS
!    Added unattenuated variables and variable write switches
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    maxgate   Maximum gates in a radial
!    maxazim   Maximum radials per tilt
!    maxelev   Maximum number of tilts
!
  IMPLICIT NONE
!
  CHARACTER (LEN=256), INTENT(IN) :: fnamevol

  INTEGER, INTENT(IN) :: maxgatevel
  INTEGER, INTENT(IN) :: maxgateref
  INTEGER, INTENT(IN) :: maxazim
  INTEGER, INTENT(IN) :: maxelev

  INTEGER, INTENT(IN) :: itimvol(maxelev)
  REAL, INTENT(IN)    :: rngvvol(maxgatevel,maxelev)
  REAL, INTENT(IN)    :: azmvvol(maxazim,maxelev)
  REAL, INTENT(IN)    :: elvvvol(maxazim,maxelev)
  REAL, INTENT(IN)    :: velvol(maxgatevel,maxazim,maxelev)
  REAL, INTENT(IN)    :: uavelvol(maxgatevel,maxazim,maxelev)
  REAL, INTENT(IN)    :: vorvol(maxgatevel,maxazim,maxelev)
  REAL, INTENT(IN)    :: qrvol(maxgateref,maxazim,maxelev)
  REAL, INTENT(IN)    :: qsvol(maxgateref,maxazim,maxelev)
  REAL, INTENT(IN)    :: qhvol(maxgateref,maxazim,maxelev)
  REAL, INTENT(IN)    :: uvol(maxgateref,maxazim,maxelev)
  REAL, INTENT(IN)    :: vvol(maxgateref,maxazim,maxelev)
  REAL, INTENT(IN)    :: wvol(maxgateref,maxazim,maxelev)
  REAL, INTENT(IN)    :: tkvol(maxgateref,maxazim,maxelev)
  REAL, INTENT(IN)    :: vnyqvol(maxelev)

  REAL, INTENT(IN)    :: rngrvol(maxgateref,maxelev)
  REAL, INTENT(IN)    :: azmrvol(maxazim,maxelev)
  REAL, INTENT(IN)    :: elvrvol(maxazim,maxelev)
  REAL, INTENT(IN)    :: refvol(maxgateref,maxazim,maxelev)
  REAL, INTENT(IN)    :: uarefvol(maxgateref,maxazim,maxelev)
  REAL, INTENT(IN)    :: zdrvol(maxgateref,maxazim,maxelev)
  REAL, INTENT(IN)    :: kdpvol(maxgateref,maxazim,maxelev)
  REAL, INTENT(IN)    :: rhohvvol(maxgateref,maxazim,maxelev)
  INTEGER, INTENT(IN) :: wrtuaref,wrtuavel
  INTEGER, INTENT(IN) :: wrtvort,wrtdualp,wrtqx,wrtuvwt
  CHARACTER (LEN=80 ), INTENT(IN) :: runname
  CHARACTER (LEN=4  ), INTENT(IN) :: radname

  REAL, INTENT(IN)    :: beamwid,rmisval,rngfval

  INTEGER, INTENT(IN) :: vcp
  REAL, INTENT(IN)    :: curtim
  REAL, INTENT(IN)    :: radar_alt
  REAL, INTENT(IN)    :: radar_lat
  REAL, INTENT(IN)    :: radar_lon
!
!-----------------------------------------------------------------------
!
! Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iyear,imon,iday,ihour,imin,isec,idummy
  INTEGER :: igate,jazim
  REAL :: qrmin,qrmax,qsmin,qsmax,qhmin,qhmax
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  idummy=0
  CALL abss2ctim(itimvol(1),iyear,imon,iday,ihour,imin,isec)

  qrmin=999.
  qrmax=-999.
  qsmin=999.
  qsmax=-999.
  qhmin=999.
  qhmax=-999.
  DO jazim=1,maxazim
    DO igate=1,maxgateref
      qrmin=min(qrmin,qrvol(igate,jazim,1))
      qsmin=min(qsmin,qsvol(igate,jazim,1))
      qhmin=min(qhmin,qhvol(igate,jazim,1))
      qrmax=max(qrmax,qrvol(igate,jazim,1))
      qsmax=max(qsmax,qsvol(igate,jazim,1))
      qhmax=max(qhmax,qhvol(igate,jazim,1))
    END DO
  END DO
  print *, ' qr range (g/kg): ',(1000.0*qrmin),' to ',(1000.0*qrmax)
  print *, ' qs range (g/kg): ',(1000.0*qsmin),' to ',(1000.0*qsmax)
  print *, ' qh range (g/kg): ',(1000.0*qhmin),' to ',(1000.0*qhmax)
!
!-----------------------------------------------------------------------
!
!  Open file for output
!
!-----------------------------------------------------------------------
!
  WRITE(6,'(a,a)') 'Opening: ',TRIM(fnamevol)
  OPEN(13,FILE=trim(fnamevol),STATUS='unknown',                        &
       FORM='unformatted')

  print *, ' runname=',runname
  WRITE(13) runname
  WRITE(13) radname
  WRITE(13) maxgatevel,maxgateref,maxazim,maxelev
  WRITE(13) iyear,imon,iday,ihour,imin,isec
  WRITE(13) radar_alt,radar_lat,radar_lon
  WRITE(13) vcp,beamwid,rmisval,rngfval,curtim
  WRITE(13) wrtuaref,wrtuavel,wrtvort,wrtdualp,                        &
            wrtqx,wrtuvwt,idummy,idummy

  WRITE(13) itimvol
  WRITE(13) rngrvol
  WRITE(13) azmrvol
  WRITE(13) elvrvol
  WRITE(13) refvol
  IF(wrtuaref /= 0) WRITE(13) uarefvol
  IF(wrtdualp /= 0) THEN
    WRITE(13) zdrvol
    WRITE(13) kdpvol
    WRITE(13) rhohvvol
  END IF

  WRITE(13) vnyqvol
  WRITE(13) rngvvol
  WRITE(13) azmvvol
  WRITE(13) elvvvol
  WRITE(13) velvol
  IF(wrtuavel /= 0) WRITE(13) uavelvol
  IF(wrtvort  /= 0) WRITE(13) vorvol
  IF(wrtqx  /= 0) THEN
    WRITE(13) qrvol
    WRITE(13) qsvol
    WRITE(13) qhvol
  END IF
  IF(wrtuvwt  /= 0) THEN
    WRITE(13) uvol
    WRITE(13) vvol
    WRITE(13) wvol
    WRITE(13) tkvol
  END IF

  close(13)

  WRITE(6,'(a,a)') ' Writing complete for file ',TRIM(fnamevol)

  RETURN
END SUBROUTINE wtradvol
!
SUBROUTINE xbatten(ngate,maxazim,nazim,                                &
                   refl,gatesp,rmisval,refmin,attrefl)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate X-band (3.2 cm) attenuated reflectivity given unattenuated
!  reflectivty (from S-band data or from simulation model hydrometeors)
!  for a single tilt of data.
!
!  Uses attenuation coefficients from Burrows and Atwood (1949) data as
!  reported in Doviak and Zrnic', Doppler Radar and Weather Observations,
!  2nd Ed., Eq. 3.16c.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, CAPS  and Erin Fay, OU/SoM
!          November, 2004
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  ngate:     Number of gates in each radial
!  nazim:     Number of radials of data
!  refl:      Unattenuated reflectivity factor dBZ
!  gatesp:    Gate spacing (m)
!  rmisval:   Missing data value, assumed below threshold
!  refmin:    Threshold minimum reflectivity (dBZ)
!
!
!  OUTPUT:
!  attrefl:    Attenuated reflectivity.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER, INTENT (IN)  :: ngate      ! Number of gates
  INTEGER, INTENT (IN)  :: maxazim     ! Maximum number of radials
  INTEGER, INTENT (IN)  :: nazim       ! Number of radials
  REAL, INTENT (IN):: refl(ngate,maxazim)     ! Reflectivity (dBZ)
  REAL, INTENT (IN)  :: gatesp         ! Gate spacing (m)
  REAL, INTENT (IN)  :: rmisval
  REAL, INTENT (IN)  :: refmin
  REAL, INTENT (OUT) :: attrefl(ngate,maxazim) ! Attenuated reflectivities

  INTEGER :: igate,jazim ! Index.
  REAL :: att_tally      ! Running tally of attenuation values.
  REAL :: att_const1     ! Attenuation constant 1.
  REAL :: att_const2     ! Attenuation constant 2.
  REAL :: att_exp        ! Attenuation exponent
  REAL :: atten,atten_last,hlfgatsp

  REAL, PARAMETER :: zrconst=400.
  REAL, PARAMETER :: zrexp=1.4
  REAL, PARAMETER :: krconst=0.01
  REAL, PARAMETER :: krexp=1.21
  REAL, PARAMETER :: mperkm=0.001
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  att_const1 = 2.0*krconst*mperkm ! leading constant (two-way)
  att_const2 = 1.0/zrconst
  att_exp = krexp/zrexp
  hlfgatsp = 0.5*gatesp

  DO jazim=1, nazim
    IF(refl(1,jazim) .gt. rmisval) THEN
      atten = &
        att_const1*((att_const2*(10.0**(0.1*refl(1,jazim))))**att_exp)
      att_tally = atten*hlfgatsp
      attrefl(1,jazim) = refl(1,jazim)-att_tally
      attrefl(1,jazim) = max(attrefl(1,jazim),refmin)
    ELSE
      atten=0.
      att_tally=0.
      attrefl(1,jazim)=refl(1,jazim)
    END IF
    atten_last=atten
    DO igate = 2, ngate
      IF(refl(igate,jazim) .gt. rmisval) THEN
        atten = &
        att_const1*((att_const2*(10.0**(0.1*refl(igate,jazim))))**att_exp)
        att_tally = att_tally + hlfgatsp*atten_last + hlfgatsp*atten
        attrefl(igate,jazim) = refl(igate,jazim)-att_tally
        attrefl(igate,jazim) = max(attrefl(igate,jazim),refmin)
        atten_last=atten
      ELSE
        attrefl(igate,jazim) = refl(igate,jazim)
        attrefl(igate,jazim) = max(attrefl(igate,jazim),refmin)
        atten_last=0.
      END IF
    END DO
  END DO
  RETURN
END SUBROUTINE xbatten

SUBROUTINE xbsens(ngate,gatesp,pulselen,beamwid,sens)
!
!-----------------------------------------------------------------------
! Subroutine xbsens
! Calculate X-band radar sensitivity (dBZ) as a function of range.
!
! Based on information provided by Francesc Junyent
!
! Keith Brewster and Erin Fay, CAPS
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ngate
  REAL, INTENT(IN) :: gatesp
  REAL, INTENT(IN) :: pulselen
  REAL, INTENT(IN) :: beamwid
  REAL, INTENT (OUT):: sens(ngate)     ! Sensitivity values for all dist from radar.
!
! Parameters for CASA radar
!
  REAL, PARAMETER :: lmbda = 3.0E-02 ! Wavelength (m)
  REAL, PARAMETER :: pt = 12500      ! Peak power (Watts)
  REAL, PARAMETER :: G  = 38.0       ! Antenna gain (dB)
  REAL, PARAMETER :: F  = 5.5        ! Noise figure (dB)
  REAL, PARAMETER :: B = 2.5         ! Bandwidth (MHz)
  REAL, PARAMETER :: lm = 1.0        ! Receiver mis-match loss(dB)
  REAL, PARAMETER :: Kw2 = 0.91      ! Refractive Index of water
  REAL, PARAMETER :: T0 = 300.       ! Rx temperature (K)
  REAL, PARAMETER :: Ta = 200.       ! Antenna temperature (K)
!
! Physical constants
!
  REAL, PARAMETER :: K = 1.38E-23    ! Boltsmann's Constant (J/K)
  REAL, PARAMETER :: c = 2.99792E08  ! Speed of light (m/s)
  REAL, PARAMETER :: rconst =1.0E18  ! m6/m3 to mm6/m3
!
! Misc internal variables
!
  INTEGER :: igate
  REAL :: pulsel,pulsew
  REAL :: ln2,pi,pi3,four3,bwrad,rnoise,BHz,Ni,sconst,rlm,rG
  REAL :: range
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ln2=log(2.0)
  pi=acos(-1.)
  pi3=pi*pi*pi
  four3=4.*4.*4.
  bwrad=beamwid*pi/180.
  rnoise=10.**(0.1*F)
  rlm=10.**(0.1*lm)
  rG=10.**(0.1*G)
  BHz=1.0E06*b
  Ni=K*(Ta+(rnoise-1.0)*T0)*BHz
  sconst=rconst *(Ni/(pi3*Kw2)) * ((8.*ln2)/(bwrad*bwrad)) *         &
         (2./(c*pulselen)) * ((lmbda*lmbda*four3*rlm)/(Pt*rG*rG))
  print *, ' pulselen=',pulselen
  print *, ' beamwid=',beamwid
  print *, ' Sensitivity................'
  DO igate = 1, ngate
     range=igate*gatesp
     sens(igate) = 10.*LOG10(range*range*sconst)
     IF(mod(igate,10) == 0) print *, (0.001*range),sens(igate)
  END DO
END SUBROUTINE xbsens

FUNCTION EFFBMW(beamwid,rotrate,prt,npulse)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Estimate effective half-power beamwidth given the antenna half-power
!  beamwidth, the rotation rate, pulse repetition time and the number of
!  pulses used in estimating the moment data.
!
!  Based on solutions to Eq. 7.34 in Doviak and Zrnic',
!  Doppler Radar and Weather Observations, 2nd Ed. (1993)
!
!  AUTHOR: Keith Brewster, CAPS
!          November, 2004
!
!  INPUT:
!    beamwid    Half-power antenna beamwidth (degrees)
!    rotrate    Rotation rate (degrees per second)
!    prt        Pulse repetition time (seconds)
!    npulse     Number of pulses averaged to get moment data
!
!  OUTPUT:
!    effbmw     Effective half-power beamwidth (degrees)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  REAL :: effbmw
!
  REAL :: beamwid
  REAL :: rotrate
  REAL :: prt
  INTEGER :: npulse
!
!-----------------------------------------------------------------------
!
! Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: ntry=30
  REAL,PARAMETER :: eps = 1.0E-06
!
  REAL :: sqrtln4,twortln4,bmwinv,amts,erfcst
  REAL :: gtheta,theta,test,test0,test2,delth
  REAL :: theta1,theta2
  INTEGER :: itry
!
!-----------------------------------------------------------------------
!
! External functions
!
!-----------------------------------------------------------------------
!
  REAL :: erf
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  sqrtln4=sqrt(alog(4.))
  twortln4=2.0*sqrtln4
  bmwinv=1.0/beamwid
  amts=abs(rotrate)*npulse*prt
  erfcst=0.5*erf(sqrtln4*amts*bmwinv)
!
!-----------------------------------------------------------------------
!
! Iteratively solve for zeros.  First for positive solution angle.
!
!-----------------------------------------------------------------------
!
  gtheta=beamwid+amts
  theta=0.
  test0=erf(twortln4*theta*bmwinv)-erf(twortln4*(theta-amts)*bmwinv)- &
       erfcst
  delth=gtheta
  theta=gtheta
  test2=erf(twortln4*theta*bmwinv)-erf(twortln4*(theta-amts)*bmwinv)- &
       erfcst
  DO itry=1,ntry
    delth=0.5*sign(delth,(test0*test2))
    theta=theta+delth
    test2=erf(twortln4*theta*bmwinv)-erf(twortln4*(theta-amts)*bmwinv)- &
       erfcst
    IF(abs(test2) < eps) EXIT
  END DO
  theta2=theta
  WRITE(6,'(a,f10.4,a,g12.4,a,i3)')                                   &
    ' EFFBMW Positive sol, theta=',theta2,'  diffr=',test2,           &
    ' iter=',itry
!
!-----------------------------------------------------------------------
!
! Solve for negative solution angle.
!
!-----------------------------------------------------------------------
!
  theta=-gtheta
  test0=erf(twortln4*theta*bmwinv)-erf(twortln4*(theta-amts)*bmwinv)- &
       erfcst
  delth=gtheta
  theta=0.
  test2=erf(twortln4*theta*bmwinv)-erf(twortln4*(theta-amts)*bmwinv)- &
       erfcst
  DO itry=1,ntry
    delth=0.5*sign(delth,(test0*test2))
    theta=theta+delth
    test2=erf(twortln4*theta*bmwinv)-erf(twortln4*(theta-amts)*bmwinv)- &
       erfcst
    IF(abs(test2) < eps) EXIT
  END DO
  theta1=theta
  WRITE(6,'(a,f10.4,a,g12.4,a,i3)')                                   &
    ' EFFBMW Negative sol, theta=',theta1,'  diffr=',test2,           &
    ' iter=',itry
!
!-----------------------------------------------------------------------
!
! Effective beamwidth is difference between positive and negative
! solution angles.
!
!-----------------------------------------------------------------------
!
  effbmw=(theta2-theta1)
  WRITE(6,'(a,f12.4,a)')                                              &
    ' EFFBMW Effective half-power beamwidth =',effbmw,' degrees'

  RETURN
END FUNCTION EFFBMW
!
FUNCTION RANDNOR()
!
! Returns a random number with a normal (Gaussian) distribution
! with a mean of zero and standard deviation of one.
! A distribution with a mean of m and standard deviation of s would
! use the result of this function in the following way:
!    x = m + s*randnor()
!
! The algorithm used here to convert a uniform random distribution (0,1)
! to a Gaussian is the polar form of the Box-Muller transformation.  This
! code based on sample C code and a description of the technique given
! by Dr. Everett Carter of Taygeta Scientific on the following website:
! http://www.taygeta.com/random/gaussian.html
!
! I've added a test to avoid dividing by zero, rsq > epsilon, which doesn't
! seem to affect the statistics of the distribution.
!
! random_number is a fortran-90 intrinsic function.  Implementation on
! different machines and compilers may vary.
!
! Keith Brewster, CAPS
! August, 2004
!
! MODIFICATIONS
!
! 02/28/2005 Keith Brewster
! Replaced intrinsic function rand with random_number, that seems to be
! available on more machines than "rand".
!
!
  IMPLICIT NONE
  REAL, PARAMETER :: epsilon = 1.0E-12
  REAL :: randnor
  REAL :: xr,yr,rsq,r
  REAL :: myranf

  DO
    xr = 2.*myranf() - 1.
    yr = 2.*myranf() - 1.
    rsq = xr*xr + yr*yr
    IF(rsq < 1.0 .AND. rsq > epsilon) EXIT
  END DO
  r = sqrt( (-2.0 * log( rsq ) ) / rsq )
  randnor = xr * r
  RETURN
  END FUNCTION RANDNOR

  REAL FUNCTION myranf ()
!
!  A simple random number generator provided by Vijay Lakamraju UMASS, ECE
!  as ranf is not implemented on all compilers.
!
!  Generates a random number between 0 and 1.
!

  INTEGER     L, C, M
  PARAMETER ( L = 1029, C = 221591, M = 1048576 )

  INTEGER     SEED
  SAVE        SEED
  DATA        SEED / 0 /

  SEED = MOD ( SEED * L + C, M )
  myranf = REAL ( SEED ) / M

  RETURN
  END FUNCTION myranf
!
!##################################################################
!##################################################################
!######                                                      ######
!######                    FUNCTION ERF                      ######
!######                                                      ######
!##################################################################
!##################################################################
!
  FUNCTION erf(x)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  ERF computes the error function.
!
!
!  Definition:
!
!    ERF(X) = ( 2 / SQRT ( PI ) ) * Integral ( 0 <= T <= X ) EXP ( -T**2 ) dT
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  INPUT:
!
!    real x, the argument of the error function.
!
!  OUTPUT:
!    real erf, the value of the error function at X.
!
!
!  AUTHOR: John Burkardt, Pittsburgh Supercomputing Center
!          NMS Fortran Software Package
!
!  MODIFICATION HISTORY:
!  23-Oct-2004  Modified for ARPS Standard Fortran-90
!               Keith Brewster, CAPS Univ of Oklahoma
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  REAL :: x
  REAL :: erf
!
  REAL :: csevl
  REAL :: erfc
!
  REAL :: erfcs(13)
  INTEGER inits
  REAL :: r1mach
  REAL :: sqeps
  REAL :: sqrtpi
  REAL :: y
!
  REAL, SAVE :: xbig = 0.
  INTEGER, SAVE :: nterf = 0
!
  data erfcs( 1) /   -0.049046121234691808 /
  data erfcs( 2) /   -0.14226120510371364 /
  data erfcs( 3) /    0.010035582187599796 /
  data erfcs( 4) /   -0.000576876469976748 /
  data erfcs( 5) /    0.000027419931252196 /
  data erfcs( 6) /   -0.000001104317550734 /
  data erfcs( 7) /    0.000000038488755420 /
  data erfcs( 8) /   -0.000000001180858253 /
  data erfcs( 9) /    0.000000000032334215 /
  data erfcs(10) /   -0.000000000000799101 /
  data erfcs(11) /    0.000000000000017990 /
  data erfcs(12) /   -0.000000000000000371 /
  data erfcs(13) /    0.000000000000000007 /
  data sqeps / 0. /
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  sqrtpi=sqrt(acos(-1.))
!
  IF ( nterf == 0 ) THEN
    nterf = inits ( erfcs, 13, 0.1*r1mach(3) )
    xbig = sqrt ( - log ( sqrtpi * r1mach(3) ) )
    sqeps = sqrt ( 2. * epsilon ( sqeps ) )
  END IF

  y = abs ( x )

  IF ( y <= sqeps ) THEN
    erf = 2.0 * x / sqrtpi
  ELSE IF ( y <= 1.0 ) THEN
    erf = x * ( 1. + csevl ( (2.*x*x - 1.), erfcs, nterf ) )
  ELSE IF ( y <= xbig ) THEN
    erf = sign ( 1. - erfc ( y ), x )
  ELSE
    erf = sign ( 1., x )
  END IF

  RETURN
  END FUNCTION erf
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                    FUNCTION ERFC                     ######
!######                                                      ######
!##################################################################
!##################################################################
!
  FUNCTION erfc ( x )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!    erfc(x) calculates the single precision complementary error
!    function for single precision argument x.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!
! series for erf        on the interval  0.          to  1.00000d+00
!                                        with weighted error   7.10e-18
!                                         log weighted error  17.15
!                               significant figures required  16.31
!                                    decimal places required  17.71
!
! series for erfc       on the interval  0.          to  2.50000d-01
!                                        with weighted error   4.81e-17
!                                         log weighted error  16.32
!                        approx significant figures required  15.0E+00
!
!
! series for erc2       on the interval  2.50000d-01 to  1.00000d+00
!                                        with weighted error   5.22e-17
!                                         log weighted error  16.28
!                        approx significant figures required  15.0E+00
!                                    decimal places required  16.96
!
!
!  AUTHOR: John Burkardt, Pittsburgh Supercomputing Center
!          NMS Fortran Software Package
!
!  MODIFICATION HISTORY:
!  23-Oct-2004  Modified for ARPS Standard Fortran-90
!               Keith Brewster, CAPS Univ of Oklahoma
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  REAL :: erfc
  REAL :: x
!
  REAL :: csevl
  REAL :: r1mach
!
  REAL :: erfcs(13)
  REAL :: erfccs(24)
  REAL :: erc2cs(23)
  REAL :: eta
  INTEGER :: inits
  INTEGER :: nterc2
  INTEGER :: nterf
  INTEGER :: nterfc
  REAL :: sqeps
  REAL :: sqrtpi
  REAL :: xmax
  REAL :: xsml
  REAL :: y
!
  data erfcs( 1) /   -.049046121234691808e0 /
  data erfcs( 2) /   -.14226120510371364e0 /
  data erfcs( 3) /    .010035582187599796e0 /
  data erfcs( 4) /   -.000576876469976748e0 /
  data erfcs( 5) /    .000027419931252196e0 /
  data erfcs( 6) /   -.000001104317550734e0 /
  data erfcs( 7) /    .000000038488755420e0 /
  data erfcs( 8) /   -.000000001180858253e0 /
  data erfcs( 9) /    .000000000032334215e0 /
  data erfcs(10) /   -.000000000000799101e0 /
  data erfcs(11) /    .000000000000017990e0 /
  data erfcs(12) /   -.000000000000000371e0 /
  data erfcs(13) /    .000000000000000007e0 /
  data erc2cs( 1) /   -.069601346602309501e0 /
  data erc2cs( 2) /   -.041101339362620893e0 /
  data erc2cs( 3) /    .003914495866689626e0 /
  data erc2cs( 4) /   -.000490639565054897e0 /
  data erc2cs( 5) /    .000071574790013770e0 /
  data erc2cs( 6) /   -.000011530716341312e0 /
  data erc2cs( 7) /    .000001994670590201e0 /
  data erc2cs( 8) /   -.000000364266647159e0 /
  data erc2cs( 9) /    .000000069443726100e0 /
  data erc2cs(10) /   -.000000013712209021e0 /
  data erc2cs(11) /    .000000002788389661e0 /
  data erc2cs(12) /   -.000000000581416472e0 /
  data erc2cs(13) /    .000000000123892049e0 /
  data erc2cs(14) /   -.000000000026906391e0 /
  data erc2cs(15) /    .000000000005942614e0 /
  data erc2cs(16) /   -.000000000001332386e0 /
  data erc2cs(17) /    .000000000000302804e0 /
  data erc2cs(18) /   -.000000000000069666e0 /
  data erc2cs(19) /    .000000000000016208e0 /
  data erc2cs(20) /   -.000000000000003809e0 /
  data erc2cs(21) /    .000000000000000904e0 /
  data erc2cs(22) /   -.000000000000000216e0 /
  data erc2cs(23) /    .000000000000000052e0 /
  data erfccs( 1) /   0.0715179310202925e0 /
  data erfccs( 2) /   -.026532434337606719e0 /
  data erfccs( 3) /    .001711153977920853e0 /
  data erfccs( 4) /   -.000163751663458512e0 /
  data erfccs( 5) /    .000019871293500549e0 /
  data erfccs( 6) /   -.000002843712412769e0 /
  data erfccs( 7) /    .000000460616130901e0 /
  data erfccs( 8) /   -.000000082277530261e0 /
  data erfccs( 9) /    .000000015921418724e0 /
  data erfccs(10) /   -.000000003295071356e0 /
  data erfccs(11) /    .000000000722343973e0 /
  data erfccs(12) /   -.000000000166485584e0 /
  data erfccs(13) /    .000000000040103931e0 /
  data erfccs(14) /   -.000000000010048164e0 /
  data erfccs(15) /    .000000000002608272e0 /
  data erfccs(16) /   -.000000000000699105e0 /
  data erfccs(17) /    .000000000000192946e0 /
  data erfccs(18) /   -.000000000000054704e0 /
  data erfccs(19) /    .000000000000015901e0 /
  data erfccs(20) /   -.000000000000004729e0 /
  data erfccs(21) /    .000000000000001432e0 /
  data erfccs(22) /   -.000000000000000439e0 /
  data erfccs(23) /    .000000000000000138e0 /
  data erfccs(24) /   -.000000000000000048e0 /
  data nterf, nterfc, nterc2, xsml, xmax, sqeps /3*0, 3*0./
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  sqrtpi=sqrt(acos(-1.))
  IF ( nterf == 0 ) THEN

    eta = 0.1 * r1mach(3)
    nterf = inits ( erfcs, 13, eta )
    nterfc = inits ( erfccs, 24, eta )
    nterc2 = inits ( erc2cs, 23, eta )

    xsml = -sqrt ( - log ( sqrtpi * r1mach(3) ) )
    xmax = sqrt ( - log ( sqrtpi * r1mach(1) ) )
    xmax = xmax - 0.5 * log ( xmax ) / xmax - 0.01
    sqeps = sqrt ( 2. * r1mach(3) )

  END IF

  IF ( x <= xsml ) THEN
    erfc = 2.
    RETURN
  END IF

  IF ( x > xmax ) THEN
    WRITE(6,'(a)') 'erfc    x so big erfc underflows'
    erfc = 0.
    RETURN
  END IF

  y = abs(x)
  IF (y>1.0) THEN
!
! erfc(x) = 1.0E+00 - erf(x) for 1. < abs(x) <= xmax
!
    y = y*y

    IF (y<=4.) THEN
      erfc = exp(-y)/abs(x) * (0.5 + csevl ((8./y-5.)/3.,erc2cs, nterc2) )
    ELSE
      erfc = exp(-y)/abs(x) * (0.5 + csevl (8./y-1.,erfccs, nterfc) )
    END IF

    IF ( x < 0.0 ) erfc = 2.0 - erfc

    RETURN
  ELSE
!
!  erfc(x) = 1. - erf(x) for -1. <= x <= 1.
!
    IF ( y < sqeps ) THEN
      erfc = 1.0 - 2.0*x/sqrtpi
    ELSE IF (y>=sqeps) THEN
      erfc = 1.0 - x*(1.0 + csevl((2.0*x*x-1.0), erfcs, nterf) )
    END IF

    RETURN
  END IF
  END FUNCTION erfc
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  FUNCTION CSEVL                      ######
!######                                                      ######
!##################################################################
!##################################################################
!
  FUNCTION csevl ( x, cs, n )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  CSEVL evaluates an N term Chebyshev series.
!
!
!  Reference:
!
!    R Broucke,
!    algorithm 446, c.a.c.m.,
!    volume 16, page 254, 1973.
!
!    Fox and Parker,
!    chebyshev polynomials in numerical analysis,
!    oxford press, page 56.
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!  INPUT:
!
!    x    value at which the series is to be evaluated.
!    cs   array of n terms of a chebyshev series.  in eval-
!         uating cs, only half the first coefficient is summed.
!    n    number of terms in array cs.
!
!
!  AUTHOR: John Burkardt, Pittsburgh Supercomputing Center
!          NMS Fortran Software Package
!
!  MODIFICATION HISTORY:
!  23-Oct-2004  Modified for ARPS Standard Fortran-90
!               Keith Brewster, CAPS Univ of Oklahoma
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  REAL :: csevl
  INTEGER :: n
  REAL :: x
  REAL :: cs(n)
!
  REAL :: b0
  REAL :: b1
  REAL :: b2
  INTEGER :: i
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF ( n < 1 ) THEN
    WRITE ( *, * ) ' '
    WRITE ( *, * ) 'CSEVL - Fatal error!'
    WRITE ( *, * ) '  Number of terms N is less than 1.'
    STOP
  END IF

  IF ( n > 1000 ) THEN
    WRITE ( *, * ) ' '
    WRITE ( *, * ) 'CSEVL - Fatal error!'
    WRITE ( *, * ) '  The number of terms is more than 1000.'
    STOP
  END IF

  IF ( x < -1.0 .or. x > 1.0 ) THEN
    WRITE ( *, * ) ' '
    WRITE ( *, * ) 'CSEVL - Fatal error!'
    WRITE ( *, * ) '  The input argument X is outside the interval [-1,1].'
    STOP
  END IF

  b1 = 0.
  b0 = 0.

  DO i = n, 1, -1
    b2 = b1
    b1 = b0
    b0 = 2.*x*b1 - b2 + cs(i)
  END DO

  csevl = 0.5*(b0-b2)

  RETURN
  END FUNCTION csevl
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  FUNCTION INITS                      ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION inits ( os, nos, eta )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!
!  INITS estimates the order of an orthogonal series guaranteeing a given accuracy.
!
!
!  Reference:
!
!    David Kahaner, Clever Moler, Steven Nash,
!    Numerical Methods and Software,
!    Prentice Hall, 1988.
!
!
!  INPUT:
!
!    os     array of nos coefficients in an orthogonal series.
!    nos    number of coefficients in os.
!    eta    requested accuracy of series.
!           Ordinarily, eta will be chosen to be one-tenth machine precision.
!
!
!  AUTHOR: John Burkardt, Pittsburgh Supercomputing Center
!          NMS Fortran Software Package
!
!  MODIFICATION HISTORY:
!  23-Oct-2004  Modified for ARPS Standard Fortran-90
!               Keith Brewster, CAPS Univ of Oklahoma
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: inits
  INTEGER :: nos
  REAL :: eta
!
  INTEGER :: ii,i
  REAL :: err
  REAL :: os(nos)
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF ( nos < 1) THEN
    WRITE(6,'(a)') 'inits   number of coefficients lt 1'
    STOP
  END IF

  err = 0.
  DO ii=1,nos
    i = nos + 1 - ii
    err = err + abs(os(i))
    if (err>eta) EXIT
  END DO

  if (i==nos) WRITE(6,'(a)') 'inits   eta may be too small'
  inits = i

  RETURN
  END FUNCTION inits
!
!##################################################################
!##################################################################
!######                                                      ######
!######                    FUNCTION ERF                      ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
  FUNCTION r1mach ( i )
!
!*******************************************************************************
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  R1MACH returns single precision machine constants.
!
!
!  Assume that single precision numbers are stored with a mantissa of T digits
!  in base B, with an exponent whose value must lie between EMIN and EMAX.  Then
!  for values of I between 1 and 5, R1MACH will return the following values:
!
!    R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
!    R1MACH(2) = B**EMAX*(1-B**(-T)), the largest magnitude.
!    R1MACH(3) = B**(-T), the smallest relative spacing.
!    R1MACH(4) = B**(1-T), the largest relative spacing.
!    R1MACH(5) = log10(B)
!
!  To alter this function for a particular environment, the desired set of data
!  statements should be activated by removing the C from column 1.
!
!  On rare machines a STATIC statement may need to be added.  But probably more
!  systems prohibit it that require it.
!
!  For IEEE-arithmetic machines (binary standard), the first set of constants
!  below should be appropriate.
!
!  Where possible, octal or hexadecimal constants have been used to specify the
!  constants exactly which has in some cases required the use of EQUIVALENCED
!  integer arrays.  If your compiler uses half-word integers by default
!  (sometimes called INTEGER*2), you may need to change INTEGER to INTEGER*4 or
!  otherwise instruct your compiler to use full-word integers in the next 5
!  declarations.
!
!
!  AUTHOR: John Burkardt, Pittsburgh Supercomputing Center
!          NMS Fortran Software Package
!
!  MODIFICATION HISTORY:
!  23-Oct-2004  Modified for ARPS Standard Fortran-90
!               Keith Brewster, CAPS Univ of Oklahoma
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  REAL :: r1mach
  INTEGER :: i
!
  INTEGER :: diver(2)
  INTEGER :: large(2)
  INTEGER :: ilog10(2)
  INTEGER :: right(2)
  REAL :: rmach(5)
  INTEGER ::small(2)
!
  equivalence (rmach(1),small(1))
  equivalence (rmach(2),large(1))
  equivalence (rmach(3),right(1))
  equivalence (rmach(4),diver(1))
  equivalence (rmach(5),ilog10(1))
!
!  IEEE arithmetic machines, such as the ATT 3B series, Motorola 68000 based
!  machines such as the SUN 3 and ATT PC 7300, and 8087 based micros such as
!  the IBM PC and ATT 6300.
!
!  data small(1) /     8388608 /
!  data large(1) /  2139095039 /
!  data right(1) /   864026624 /
!  data diver(1) /   872415232 /
!  data ilog10(1) /  1050288283 /
!
!  ALLIANT FX/8 UNIX Fortran compiler with the -r8 command line option.  This
!  option causes all variables declared with 'REAL' to be of type 'REAL*8' or
!  DOUBLE PRECISION.  This option does not override the 'REAL*4' declarations.
!  These R1MACH numbers below and the coresponding I1MACH are simply the DOUBLE
!  PRECISION or 'REAL*8' numbers.  If you use the -r8 your whole code (and the
!  user libraries you link with, the system libraries are taken care of
!  automagicly) must be compiled with this option.
!
!      data rmach(1) / 2.22507385850721D-308 /
!      data rmach(2) / 1.79769313486231D+308 /
!      data rmach(3) / 1.1101827117665D-16 /
!      data rmach(4) / 2.2203654423533D-16 /
!      data rmach(5) / 3.01029995663981E-1 /
!
!  AMDAHL machines.
!
!      data small(1) /    1048576 /
!      data large(1) / 2147483647 /
!      data right(1) /  990904320 /
!      data diver(1) / 1007681536 /
!      data ilog10(1) / 1091781651 /
!
!  BURROUGHS 1700 system.
!
!      data rmach(1) / Z400800000 /
!      data rmach(2) / Z5FFFFFFFF /
!      data rmach(3) / Z4E9800000 /
!      data rmach(4) / Z4EA800000 /
!      data rmach(5) / Z500E730E8 /
!
!  BURROUGHS 5700/6700/7700 systems.
!
!      data rmach(1) / O1771000000000000 /
!      data rmach(2) / O0777777777777777 /
!      data rmach(3) / O1311000000000000 /
!      data rmach(4) / O1301000000000000 /
!      data rmach(5) / O1157163034761675 /
!
!  CDC CYBER 170/180 series using NOS
!
!      data rmach(1) / O"00014000000000000000" /
!      data rmach(2) / O"37767777777777777777" /
!      data rmach(3) / O"16404000000000000000" /
!      data rmach(4) / O"16414000000000000000" /
!      data rmach(5) / O"17164642023241175720" /
!
!  CDC CYBER 170/180 series using NOS/VE
!
!      data rmach(1) / Z"3001800000000000" /
!      data rmach(2) / Z"4FFEFFFFFFFFFFFE" /
!      data rmach(3) / Z"3FD2800000000000" /
!      data rmach(4) / Z"3FD3800000000000" /
!      data rmach(5) / Z"3FFF9A209A84FBCF" /
!
!  CDC CYBER 200 series
!
!      data rmach(1) / X'9000400000000000' /
!      data rmach(2) / X'6FFF7FFFFFFFFFFF' /
!      data rmach(3) / X'FFA3400000000000' /
!      data rmach(4) / X'FFA4400000000000' /
!      data rmach(5) / X'FFD04D104D427DE8' /
!
!  CDC 6000/7000 series using FTN4.
!
!      data rmach(1) / 00564000000000000000B /
!      data rmach(2) / 37767777777777777776B /
!      data rmach(3) / 16414000000000000000B /
!      data rmach(4) / 16424000000000000000B /
!      data rmach(5) / 17164642023241175720B /
!
!  CDC 6000/7000 series using FTN5.
!
!      data rmach(1) / O"00564000000000000000" /
!      data rmach(2) / O"37767777777777777776" /
!      data rmach(3) / O"16414000000000000000" /
!      data rmach(4) / O"16424000000000000000" /
!      data rmach(5) / O"17164642023241175720" /
!
!  CONVEX C-1.
!
!      data rmach(1) / '00800000'X /
!      data rmach(2) / '7FFFFFFF'X /
!      data rmach(3) / '34800000'X /
!      data rmach(4) / '35000000'X /
!      data rmach(5) / '3F9A209B'X /
!
!  CONVEX C-120 (native mode) without -R8 option
!
!      data rmach(1) / 2.9387360E-39 /
!      data rmach(2) / 1.7014117E+38 /
!      data rmach(3) / 5.9604645E-08 /
!      data rmach(4) / 1.1920929E-07 /
!      data rmach(5) / 3.0102999E-01 /
!
!  CONVEX C-120 (native mode) with -R8 option
!
!      data rmach(1) / 5.562684646268007D-309 /
!      data rmach(2) / 8.988465674311577D+307 /
!      data rmach(3) / 1.110223024625157D-016 /
!      data rmach(4) / 2.220446049250313D-016 /
!      data rmach(5) / 3.010299956639812D-001 /
!
!  CONVEX C-120 (IEEE mode) without -R8 option
!
!      data rmach(1) / 1.1754945E-38 /
!      data rmach(2) / 3.4028234E+38 /
!      data rmach(3) / 5.9604645E-08 /
!      data rmach(4) / 1.1920929E-07 /
!      data rmach(5) / 3.0102999E-01 /
!
!  CONVEX C-120 (IEEE mode) with -R8 option
!
!      data rmach(1) / 2.225073858507202D-308 /
!      data rmach(2) / 1.797693134862315D+308 /
!      data rmach(3) / 1.110223024625157D-016 /
!      data rmach(4) / 2.220446049250313D-016 /
!      data rmach(5) / 3.010299956639812D-001 /
!
!  CRAY 1, 2, XMP and YMP.
!
!      data rmach(1) / 200034000000000000000B /
!      data rmach(2) / 577767777777777777776B /
!      data rmach(3) / 377224000000000000000B /
!      data rmach(4) / 377234000000000000000B /
!      data rmach(5) / 377774642023241175720B /
!
!  DATA GENERAL ECLIPSE S/200.
!  Note - It may be appropriate to include the line: STATIC RMACH(5)
!
!      data small /20K,0/
!      data large /77777K,177777K/
!      data right /35420K,0/
!      data diver /36020K,0/
!      data ilog10 /40423K,42023K/
!
!  ELXSI 6400, assuming REAL*4 is the default real type.
!
!      data small(1) / '00800000'X /
!      data large(1) / '7F7FFFFF'X /
!      data right(1) / '33800000'X /
!      data diver(1) / '34000000'X /
!      data ilog10(1) / '3E9A209B'X /
!
!  HARRIS 220
!
!      data small(1),small(2) / '20000000, '00000201 /
!      data large(1),large(2) / '37777777, '00000177 /
!      data right(1),right(2) / '20000000, '00000352 /
!      data diver(1),diver(2) / '20000000, '00000353 /
!      data ilog10(1),ilog10(2) / '23210115, '00000377 /
!
!  HARRIS SLASH 6 and SLASH 7.
!
!      data small(1),small(2) / '20000000, '00000201 /
!      data large(1),large(2) / '37777777, '00000177 /
!      data right(1),right(2) / '20000000, '00000352 /
!      data diver(1),diver(2) / '20000000, '00000353 /
!      data ilog10(1),ilog10(2) / '23210115, '00000377 /
!
!  HONEYWELL DPS 8/70 and 600/6000 series.
!
!      data rmach(1) / O402400000000 /
!      data rmach(2) / O376777777777 /
!      data rmach(3) / O714400000000 /
!      data rmach(4) / O716400000000 /
!      data rmach(5) / O776464202324 /
!
!  HP 2100, 3 word double precision with FTN4
!
!      data small(1), small(2) / 40000B,       1 /
!      data large(1), large(2) / 77777B, 177776B /
!      data right(1), right(2) / 40000B,    325B /
!      data diver(1), diver(2) / 40000B,    327B /
!      data ilog10(1), ilog10(2) / 46420B,  46777B /
!
!  HP 2100, 4 word double precision with FTN4
!
!      data small(1), small(2) / 40000B,       1 /
!      data large91), large(2) / 77777B, 177776B /
!      data right(1), right(2) / 40000B,    325B /
!      data diver(1), diver(2) / 40000B,    327B /
!      data ilog10(1), ilog10(2) / 46420B,  46777B /
!
!  HP 9000
!
!      r1mach(1) = 1.17549435E-38
!      r1mach(2) = 1.70141163E+38
!      r1mach(3) = 5.960464478E-8
!      r1mach(4) = 1.119209290E-7
!      r1mach(5) = 3.01030010E-1
!
!      data small(1) / 00040000000B /
!      data large(1) / 17677777777B /
!      data right(1) / 06340000000B /
!      data diver(1) / 06400000000B /
!      data ilog10(1) / 07646420233B /
!
!  IBM 360/370 series, XEROX SIGMA 5/7/9, SEL systems 85/86, PERKIN ELMER 3230,
!  and PERKIN ELMER (INTERDATA) 3230.
!
!      data rmach(1) / Z00100000 /
!      data rmach(2) / Z7FFFFFFF /
!      data rmach(3) / Z3B100000 /
!      data rmach(4) / Z3C100000 /
!      data rmach(5) / Z41134413 /
!
!  IBM PC - Microsoft FORTRAN
!
!      data small(1) / #00800000 /
!      data large(1) / #7F7FFFFF /
!      data right(1) / #33800000 /
!      data diver(1) / #34000000 /
!      data ilog10(1) / #3E9A209A /
!
!  IBM PC - Professional FORTRAN and Lahey FORTRAN
!
!      data small(1)/ Z'00800000' /
!      data large(1)/ Z'7F7FFFFF' /
!      data right(1)/ Z'33800000' /
!      data diver(1)/ Z'34000000' /
!      data ilog10(1)/ Z'3E9A209A' /
!
!  INTERDATA 8/32 with the UNIX system FORTRAN 77 compiler.
!  For the INTERDATA FORTRAN VII compiler replace the Z'S specifying HEX
!  constants with Y'S.
!
!      data rmach(1) / Z'00100000' /
!      data rmach(2) / Z'7EFFFFFF' /
!      data rmach(3) / Z'3B100000' /
!      data rmach(4) / Z'3C100000' /
!      data rmach(5) / Z'41134413' /
!
!  PDP-10 (KA or KI processor).
!
!      data rmach(1) / "000400000000 /
!      data rmach(2) / "377777777777 /
!      data rmach(3) / "146400000000 /
!      data rmach(4) / "147400000000 /
!      data rmach(5) / "177464202324 /
!
!  PDP-11 FORTRANS supporting 32-bit integers (integer version).
!
!      data small(1) /    8388608 /
!      data large(1) / 2147483647 /
!      data right(1) /  880803840 /
!      data diver(1) /  889192448 /
!      data ilog10(1) / 1067065499 /
!
!  PDP-11 FORTRANS supporting 32-bit integers (octal version).
!
!      data rmach(1) / O00040000000 /
!      data rmach(2) / O17777777777 /
!      data rmach(3) / O06440000000 /
!      data rmach(4) / O06500000000 /
!      data rmach(5) / O07746420233 /
!
!  PDP-11 FORTRANS supporting 16-bit integers (integer version).
!
!      data small(1),small(2) /   128,     0 /
!      data large(1),large(2) / 32767,    -1 /
!      data right(1),right(2) / 13440,     0 /
!      data diver(1),diver(2) / 13568,     0 /
!      data ilog10(1),ilog10(2) / 16282,  8347 /
!
!  PDP-11 FORTRANS supporting 16-bit integers (octal version).
!
!      data small(1),small(2) / O000200, O000000 /
!      data large(1),large(2) / O077777, O177777 /
!      data right(1),right(2) / O032200, O000000 /
!      data diver(1),diver(2) / O032400, O000000 /
!      data ilog10(1),ilog10(2) / O037632, O020233 /
!
!  SEQUENT BALANCE 8000.
!
!      data small(1) / $00800000 /
!      data large(1) / $7F7FFFFF /
!      data right(1) / $33800000 /
!      data diver(1) / $34000000 /
!      data ilog10(1) / $3E9A209B /
!
!  SUN Microsystems UNIX F77 compiler.
!
       data rmach(1) / 1.17549435E-38 /
       data rmach(2) / 3.40282347E+38 /
       data rmach(3) / 5.96016605E-08 /
       data rmach(4) / 1.19203321E-07 /
       data rmach(5) / 3.01030010E-01 /
!
!  SUN 3 (68881 or FPA)
!
!      data small(1) / X'00800000' /
!      data large(1) / X'7F7FFFFF' /
!      data right(1) / X'33800000' /
!      data diver(1) / X'34000000' /
!      data ilog10(1) / X'3E9A209B' /
!
!  UNIVAC 1100 series.
!
!      data rmach(1) / O000400000000 /
!      data rmach(2) / O377777777777 /
!      data rmach(3) / O146400000000 /
!      data rmach(4) / O147400000000 /
!      data rmach(5) / O177464202324 /
!
!  VAX/ULTRIX F77 compiler.
!
!      data small(1) /       128 /
!      data large(1) /    -32769 /
!      data right(1) /     13440 /
!      data diver(1) /     13568 /
!      data ilog10(1) / 547045274 /
!
!  VAX-11 with FORTRAN IV-PLUS compiler.
!
!      data rmach(1) / Z00000080 /
!      data rmach(2) / ZFFFF7FFF /
!      data rmach(3) / Z00003480 /
!      data rmach(4) / Z00003500 /
!      data rmach(5) / Z209B3F9A /
!
!  VAX/VMS version 2.2.
!
!      data rmach(1) /       '80'X /
!      data rmach(2) / 'FFFF7FFF'X /
!      data rmach(3) /     '3480'X /
!      data rmach(4) /     '3500'X /
!      data rmach(5) / '209B3F9A'X /
!
!  VAX/VMS 11/780
!
!      data small(1) / Z00000080 /
!      data large(1) / ZFFFF7FFF /
!      data right(1) / Z00003480 /
!      data diver(1) / Z00003500 /
!      data ilog10(1) / Z209B3F9A /
!
!  Z80 microprocessor.
!
!      data small(1), small(2) /     0,    256 /
!      data large(1), large(2) /    -1,   -129 /
!      data right(1), right(2) /     0,  26880 /
!      data diver(1), diver(2) /     0,  27136 /
!      data ilog10(1), ilog10(2) /  8347,  32538 /
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF ( i > 0 .AND. i < 6 ) THEN
    r1mach = rmach(i)
  ELSE
    WRITE(6,'(/a)') 'R1MACH - Fatal error!'
    WRITE(6,'(a,i9)')'I is out of bounds=',i
    r1mach=0.
    STOP
  END IF
  RETURN
END FUNCTION r1mach
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE GET_GRIDXYZP_HDF              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE get_gridxyzp_hdf(nx,ny,nz,filename,                          &
                               x,y,zp,                                  &
                               itmp,hmax,hmin,ireturn)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Read in grid variables from base state/grid history data.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, after a similar program developed by Ming Xue
!  02/06/2005.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    filename HDF file name of grid/base file.
!
!  OUTPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    x        x-coordinate data
!    y        y-coordinate data
!    zp       zp-coordinate data
!
!  SCRATCH:
!    itmp     Temporary array for hdf compression
!    hmin     hmin temporary array
!    hmax     hmax temporary array
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: stat, sd_id
  CHARACTER (LEN=*) :: filename

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  REAL :: x(nx)
  REAL :: y(ny)
  REAL :: zp(nx,ny,nz)

  INTEGER (KIND=selected_int_kind(4)) :: itmp(nx,ny,nz) ! Temporary array
  REAL :: hmin(nz)
  REAL :: hmax(nz)

  INTEGER :: ireturn           ! Return status indicator

!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: istat
  REAL :: alatpro(2)
  REAL :: sclf,dxscl,dyscl,ctrx,ctry,swx,swy
  LOGICAL :: success
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'grid.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  CALL hdfopen(filename,1,sd_id)

  IF (sd_id < 0) THEN
    WRITE (6,*) "get_gridxyzprho_hdf: ERROR opening ",                 &
                 trim(filename)," for reading."
    ireturn=-1
    RETURN
  ELSE
    WRITE(6,*) 'File ',filename,' opened.'
  END IF

  success=.true.

  CALL hdfrdi(sd_id,"mapproj",mapproj,istat)
  CALL hdfrdr(sd_id,"trulat1",trulat1,istat)
  CALL hdfrdr(sd_id,"trulat2",trulat2,istat)
  CALL hdfrdr(sd_id,"trulon",trulon,istat)
  CALL hdfrdr(sd_id,"sclfct",sclfct,istat)
  CALL hdfrdr(sd_id,"ctrlat",ctrlat,istat)
  CALL hdfrdr(sd_id,"ctrlon",ctrlon,istat)

  CALL hdfrd1d(sd_id,"x",nx,x,istat)
  IF (istat /= 0) success=.false.

  CALL hdfrd1d(sd_id,"y",ny,y,istat)
  IF (istat /= 0) success=.false.

  CALL hdfrd3d(sd_id,"zp",nx,ny,nz,zp,istat,itmp,hmax,hmin)
  IF (istat /= 0) success=.false.

  IF (success) THEN
    ireturn = 0
  ELSE

!-----------------------------------------------------------------------
!
!  Error during read
!
!-----------------------------------------------------------------------

    WRITE(6,'(/a/)') ' Error reading data in GET_GRIDXYZPRHO_HDF.'
    ireturn=-2
  END IF

  CALL hdfclose(sd_id,stat)

  RETURN
END SUBROUTINE get_gridxyzp_hdf

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE HDFREADUVW                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfreaduvw(nx,ny,nz,filename,                               &
                      time,u,v,w,                                      &
                      itmp,hmax,hmin,ireturn)

!-----------------------------------------------------------------------
!  PURPOSE:
!  Read in ARPS wind data in the NCSA HDF4 format.
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  2000/04/15
!
!  MODIFICATION HISTORY:
!  Keith Brewster, CAPS
!  2005/02/06 Created hdfreaduvw, streamlined from hdfreadwind
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    filename  Character variable naming the input HDF file

!-----------------------------------------------------------------------
!  Variable Declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: nx,ny,nz

  CHARACTER (LEN=*) :: filename

  REAL :: time
  REAL :: u(nx,ny,nz)
  REAL :: v(nx,ny,nz)
  REAL :: w(nx,ny,nz)

  INTEGER (KIND=selected_int_kind(4)) :: itmp(nx,ny,nz) ! Temporary array
  REAL :: hmin(nz) ! Temporary array
  REAL :: hmax(nz) ! Temporary array

  INTEGER :: ireturn

!-----------------------------------------------------------------------
!  Parameters describing routine that wrote the gridded data
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
!-----------------------------------------------------------------------

  CHARACTER (LEN=40) :: fmtver410,fmtver500
  INTEGER  :: intver,intver410,intver500

  PARAMETER (fmtver410='004.10 HDF4 Coded Data',intver410=410)
  PARAMETER (fmtver500='005.00 HDF4 Coded Data',intver500=500)

  CHARACTER (LEN=40) :: fmtverin
  CHARACTER (LEN=10) :: tmunit

!-----------------------------------------------------------------------
!  Misc. local variables
!-----------------------------------------------------------------------

  INTEGER :: lchanl
  PARAMETER (lchanl=6)      ! Channel number for formatted printing.

  INTEGER :: i,j,k,is
  INTEGER :: nxin,nyin,nzin

  INTEGER :: istat, sd_id
  INTEGER :: varflg, istatus

  LOGICAL :: success

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'        ! Grid parameters
  INCLUDE 'indtflg.inc'


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  Beginning of executable code...
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  WRITE(*,*) 'HDFREADUVW: Reading HDF file: ', trim(filename)

!-----------------------------------------------------------------------
! Open file for reading
!-----------------------------------------------------------------------

  CALL hdfopen(filename,1,sd_id)
  IF (sd_id < 0) THEN
    WRITE (6,*) "HDFREADUVW: ERROR opening ", trim(filename)," for reading."
    ireturn=-1
    RETURN
  END IF

  fmtverin = fmtver500

  WRITE(6,'(/1x,a,a/)') 'Incoming data format, fmtverin=',fmtverin

  CALL hdfrdc(sd_id,80,"runname",runname,istat)
  CALL hdfrdi(sd_id,"nocmnt",nocmnt,istat)
  IF( nocmnt > 0 ) THEN
    CALL hdfrdc(sd_id,80*nocmnt,"cmnt",cmnt,istat)
  END IF

  WRITE(6,'(//''  THE NAME OF THIS RUN IS:  '',A//)') trim(runname)

  WRITE (6,*) "Comments:"
  IF( nocmnt > 0 ) THEN
    DO i=1,nocmnt
      WRITE(6,'(1x,a)') cmnt(i)
    END DO
  END IF

  WRITE (6,*) " "

  CALL hdfrdc(sd_id,10,"tmunit",tmunit,istat)
  CALL hdfrdr(sd_id,"time",time,istat)

!-----------------------------------------------------------------------
!  Get dimensions of data in binary file and check against
!  the dimensions passed to HDFREAD
!-----------------------------------------------------------------------

  CALL hdfrdi(sd_id,"nx",nxin,istat)
  CALL hdfrdi(sd_id,"ny",nyin,istat)
  CALL hdfrdi(sd_id,"nz",nzin,istat)

  IF ( nxin /= nx .OR. nyin /= ny .OR. nzin /= nz ) THEN
    WRITE(6,'(1x,a)') ' Dimensions in HDFREADUVW inconsistent with data.'
    WRITE(6,'(1x,a,3I15)') ' Read were: ', nxin, nyin, nzin
    WRITE(6,'(1x,a,3I15)') ' Expected:  ', nx, ny, nz
    WRITE(6,'(1x,a)') ' Program aborted in HDFREAD.'
    CALL arpsstop('arpsstop called from HDFREADUVW due to nxin...',1)
  END IF

  WRITE(lchanl,'(1x,a,f8.1,a,f8.3,a/)')'To read data for time:',      &
         time,' secs = ',(time/60.),' mins.'

  CALL hdfrdi(sd_id,"grdflg",grdin,istat)
  CALL hdfrdi(sd_id,"basflg",basin,istat)
  CALL hdfrdi(sd_id,"varflg",varin,istat)

  WRITE(6,'(a)') ' Done reading parameters'

  CALL hdfrdi(sd_id,"month",month,istat)
  CALL hdfrdi(sd_id,"day",day,istat)
  CALL hdfrdi(sd_id,"year",year,istat)
  CALL hdfrdi(sd_id,"hour",hour,istat)
  CALL hdfrdi(sd_id,"minute",minute,istat)
  CALL hdfrdi(sd_id,"second",second,istat)

  CALL hdfrdr(sd_id,"umove",umove,istat)
  CALL hdfrdr(sd_id,"vmove",vmove,istat)
  CALL hdfrdr(sd_id,"xgrdorg",xgrdorg,istat)
  CALL hdfrdr(sd_id,"ygrdorg",ygrdorg,istat)

  success=.true.
  IF( varin == 1 ) then

!-----------------------------------------------------------------------
!  Read in total values of variables from history dump
!-----------------------------------------------------------------------

    CALL hdfrd3d(sd_id,"u",nx,ny,nz,u,istat,itmp,hmax,hmin)
    IF (istat /= 0) success=.false.

    CALL hdfrd3d(sd_id,"v",nx,ny,nz,v,istat,itmp,hmax,hmin)
    IF (istat /= 0) success=.false.

    CALL hdfrd3d(sd_id,"w",nx,ny,nz,w,istat,itmp,hmax,hmin)
    IF (istat /= 0) success=.false.

!-----------------------------------------------------------------------
!
!  Friendly exit message
!
!-----------------------------------------------------------------------

    IF(success) THEN
      WRITE(6,'(/a,F8.1,a/)')                                           &
        ' Data at time=', time/60,' (min) were successfully read.'

      ireturn = 0
    ELSE
      WRITE(6,'(/a,F8.1,a/)')                                           &
        ' Error reading u,v,w data at time=', time/60,' (min)'
      ireturn = -3
    END IF
    CALL hdfclose(sd_id,istat)
  ELSE
    WRITE(6,'(/a/a,a/)') ' Error reading data in HDFREADUVW:',          &
      ' No time-dependent data in file',trim(filename)
    ireturn = -2
    CALL hdfclose(sd_id,istat)
  END IF

  RETURN
END SUBROUTINE hdfreaduvw
!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE HDFREADPPTRHO                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfreadppt(nx,ny,nz,filename,                               &
                      time,p,pt,                                       &
                      itmp,hmax,hmin,ireturn)

!-----------------------------------------------------------------------
!  PURPOSE:
!  Read in ARPS pressure and temperature data in the NCSA HDF4 format.
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster after hdfreaduvw
!  02/06/2005
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
!
!    filename  Character variable naming the input HDF file

!-----------------------------------------------------------------------
!  Variable Declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: nx,ny,nz

  CHARACTER (LEN=*) :: filename

  REAL :: time
  REAL :: p(nx,ny,nz)
  REAL :: pt(nx,ny,nz)

  INTEGER (KIND=selected_int_kind(4)) :: itmp(nx,ny,nz) ! Temporary array
  REAL :: hmin(nz) ! Temporary array
  REAL :: hmax(nz) ! Temporary array

  INTEGER :: ireturn

!-----------------------------------------------------------------------
!  Parameters describing routine that wrote the gridded data
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
!-----------------------------------------------------------------------

  CHARACTER (LEN=40) :: fmtver410,fmtver500
  INTEGER  :: intver,intver410,intver500

  PARAMETER (fmtver410='004.10 HDF4 Coded Data',intver410=410)
  PARAMETER (fmtver500='005.00 HDF4 Coded Data',intver500=500)

  CHARACTER (LEN=40) :: fmtverin
  CHARACTER (LEN=10) :: tmunit

!-----------------------------------------------------------------------
!  Misc. local variables
!-----------------------------------------------------------------------

  INTEGER :: lchanl
  PARAMETER (lchanl=6)      ! Channel number for formatted printing.

  INTEGER :: i,j,k,is
  INTEGER :: nxin,nyin,nzin

  INTEGER :: istat, sd_id
  INTEGER :: varflg, istatus

  LOGICAL :: success

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'        ! Grid parameters
  INCLUDE 'indtflg.inc'


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  Beginning of executable code...
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  WRITE(*,*) 'HDFREADPPTRHO: Reading HDF file: ', trim(filename)

!-----------------------------------------------------------------------
! Open file for reading
!-----------------------------------------------------------------------

  CALL hdfopen(filename,1,sd_id)
  IF (sd_id < 0) THEN
    WRITE (6,'(3a)') 'HDFREADPPTRHO: ERROR opening ', trim(filename),   &
                     ' for reading.'
    ireturn=-1
    RETURN
  END IF

  fmtverin = fmtver500

  WRITE(6,'(/1x,a,a/)') 'Incoming data format, fmtverin=',fmtverin

  CALL hdfrdc(sd_id,80,"runname",runname,istat)
  CALL hdfrdi(sd_id,"nocmnt",nocmnt,istat)
  IF( nocmnt > 0 ) THEN
    CALL hdfrdc(sd_id,80*nocmnt,"cmnt",cmnt,istat)
  END IF

  WRITE(6,'(//''  THE NAME OF THIS RUN IS:  '',A//)') trim(runname)

  WRITE (6,*) "Comments:"
  IF( nocmnt > 0 ) THEN
    DO i=1,nocmnt
      WRITE(6,'(1x,a)') cmnt(i)
    END DO
  END IF

  WRITE (6,*) " "

  CALL hdfrdc(sd_id,10,"tmunit",tmunit,istat)
  CALL hdfrdr(sd_id,"time",time,istat)

!-----------------------------------------------------------------------
!  Get dimensions of data in binary file and check against
!  the dimensions passed to HDFREAD
!-----------------------------------------------------------------------

  CALL hdfrdi(sd_id,"nx",nxin,istat)
  CALL hdfrdi(sd_id,"ny",nyin,istat)
  CALL hdfrdi(sd_id,"nz",nzin,istat)

  IF ( nxin /= nx .OR. nyin /= ny .OR. nzin /= nz ) THEN
    WRITE(6,'(1x,a)') ' Dimensions in HDFREADUVW inconsistent with data.'
    WRITE(6,'(1x,a,3I15)') ' Read were: ', nxin, nyin, nzin
    WRITE(6,'(1x,a,3I15)') ' Expected:  ', nx, ny, nz
    WRITE(6,'(1x,a)') ' Program aborted in HDFREAD.'
    CALL arpsstop('arpsstop called from HDFREADUVW due to nxin...',1)
  END IF

  WRITE(lchanl,'(1x,a,f8.1,a,f8.3,a/)')'To read data for time:',      &
         time,' secs = ',(time/60.),' mins.'

  CALL hdfrdi(sd_id,"grdflg",grdin,istat)
  CALL hdfrdi(sd_id,"basflg",basin,istat)
  CALL hdfrdi(sd_id,"varflg",varin,istat)

  WRITE(6,'(a)') ' Done reading parameters'

  CALL hdfrdi(sd_id,"month",month,istat)
  CALL hdfrdi(sd_id,"day",day,istat)
  CALL hdfrdi(sd_id,"year",year,istat)
  CALL hdfrdi(sd_id,"hour",hour,istat)
  CALL hdfrdi(sd_id,"minute",minute,istat)
  CALL hdfrdi(sd_id,"second",second,istat)

  success=.true.
  IF( varin == 1 ) THEN

!-----------------------------------------------------------------------
!  Read in total values of variables from history dump
!-----------------------------------------------------------------------

    CALL hdfrd3d(sd_id,"p",nx,ny,nz,p,istat,itmp,hmax,hmin)
    IF (istat /= 0) success=.false.

    CALL hdfrd3d(sd_id,"pt",nx,ny,nz,pt,istat,itmp,hmax,hmin)
    IF (istat /= 0) success=.false.

!-----------------------------------------------------------------------
!
!  Friendly exit message
!
!-----------------------------------------------------------------------

    IF(success) THEN
      WRITE(6,'(/a,F8.1,a/)')                                           &
        ' Data at time=', time/60,' (min) were successfully read.'

      ireturn = 0
    ELSE
      WRITE(6,'(/a,F8.1,a/)')                                           &
        ' Error reading u,v,w data at time=', time/60,' (min)'
      ireturn = -3
    END IF
    CALL hdfclose(sd_id,istat)
  ELSE
    WRITE(6,'(/a/a,a/)') ' Error reading data in HDFREADPPT:',          &
      ' No time-dependent data in file',trim(filename)
    ireturn = -2
    CALL hdfclose(sd_id,istat)
  END IF

  RETURN
END SUBROUTINE hdfreadppt
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE HDFREADQRSH                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE hdfreadqrsh(nx,ny,nz,filename,time,qr,qs,qh,                &
                       itmp,hmax,hmin,ireturn)

!-----------------------------------------------------------------------
!  PURPOSE:
!  Read in ARPS pressure and temperature data in the NCSA HDF4 format.
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster after hdfreaduvw
!  02/06/2005
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
!
!    filename  Character variable naming the input HDF file

!-----------------------------------------------------------------------
!  Variable Declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER :: nx,ny,nz

  CHARACTER (LEN=*) :: filename

  REAL :: time
  REAL :: qr(nx,ny,nz)
  REAL :: qs(nx,ny,nz)
  REAL :: qh(nx,ny,nz)

  INTEGER (KIND=selected_int_kind(4)) :: itmp(nx,ny,nz) ! Temporary array
  REAL :: hmin(nz) ! Temporary array
  REAL :: hmax(nz) ! Temporary array

  INTEGER :: ireturn

!-----------------------------------------------------------------------
!  Parameters describing routine that wrote the gridded data
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
!-----------------------------------------------------------------------

  CHARACTER (LEN=40) :: fmtver410,fmtver500
  INTEGER  :: intver,intver410,intver500

  PARAMETER (fmtver410='004.10 HDF4 Coded Data',intver410=410)
  PARAMETER (fmtver500='005.00 HDF4 Coded Data',intver500=500)

  CHARACTER (LEN=40) :: fmtverin
  CHARACTER (LEN=10) :: tmunit

!-----------------------------------------------------------------------
!  Misc. local variables
!-----------------------------------------------------------------------

  INTEGER :: lchanl
  PARAMETER (lchanl=6)      ! Channel number for formatted printing.

  INTEGER :: i,j,k,is
  INTEGER :: nxin,nyin,nzin

  INTEGER :: istat, sd_id
  INTEGER :: varflg, istatus

  LOGICAL :: success

!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'        ! Grid parameters
  INCLUDE 'indtflg.inc'


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  Beginning of executable code...
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  WRITE(*,*) 'HDFREADQRSH: Reading HDF file: ', trim(filename)

!-----------------------------------------------------------------------
! Open file for reading
!-----------------------------------------------------------------------

  CALL hdfopen(filename,1,sd_id)
  IF (sd_id < 0) THEN
    WRITE (6,*) "HDFREADQRSH: ERROR opening ", trim(filename)," for reading."
    ireturn=-1
    RETURN
  END IF

  fmtverin = fmtver500

  WRITE(6,'(/1x,a,a/)') 'Incoming data format, fmtverin=',fmtverin

  CALL hdfrdc(sd_id,80,"runname",runname,istat)
  CALL hdfrdi(sd_id,"nocmnt",nocmnt,istat)
  IF( nocmnt > 0 ) THEN
    CALL hdfrdc(sd_id,80*nocmnt,"cmnt",cmnt,istat)
  END IF

  WRITE(6,'(//''  THE NAME OF THIS RUN IS:  '',A//)') trim(runname)

  WRITE (6,*) "Comments:"
  IF( nocmnt > 0 ) THEN
    DO i=1,nocmnt
      WRITE(6,'(1x,a)') cmnt(i)
    END DO
  END IF

  WRITE (6,*) " "

  CALL hdfrdc(sd_id,10,"tmunit",tmunit,istat)
  CALL hdfrdr(sd_id,"time",time,istat)

!-----------------------------------------------------------------------
!  Get dimensions of data in binary file and check against
!  the dimensions passed to HDFREAD
!-----------------------------------------------------------------------

  CALL hdfrdi(sd_id,"nx",nxin,istat)
  CALL hdfrdi(sd_id,"ny",nyin,istat)
  CALL hdfrdi(sd_id,"nz",nzin,istat)

  IF ( nxin /= nx .OR. nyin /= ny .OR. nzin /= nz ) THEN
    WRITE(6,'(1x,a)') ' Dimensions in HDFREADUVW inconsistent with data.'
    WRITE(6,'(1x,a,3I15)') ' Read were: ', nxin, nyin, nzin
    WRITE(6,'(1x,a,3I15)') ' Expected:  ', nx, ny, nz
    WRITE(6,'(1x,a)') ' Program aborted in HDFREAD.'
    CALL arpsstop('arpsstop called from HDFREADUVW due to nxin...',1)
  END IF

  WRITE(lchanl,'(1x,a,f8.1,a,f8.3,a/)')'To read data for time:',      &
         time,' secs = ',(time/60.),' mins.'

  CALL hdfrdi(sd_id,"grdflg",grdin,istat)
  CALL hdfrdi(sd_id,"basflg",basin,istat)
  CALL hdfrdi(sd_id,"varflg",varin,istat)
  CALL hdfrdi(sd_id,"iceflg",icein,istat)

  WRITE(6,'(a)') ' Done reading parameters'

  CALL hdfrdi(sd_id,"month",month,istat)
  CALL hdfrdi(sd_id,"day",day,istat)
  CALL hdfrdi(sd_id,"year",year,istat)
  CALL hdfrdi(sd_id,"hour",hour,istat)
  CALL hdfrdi(sd_id,"minute",minute,istat)
  CALL hdfrdi(sd_id,"second",second,istat)

  success=.true.
  IF( varin == 1 ) THEN

!-----------------------------------------------------------------------
!  Read in total values of variables from history dump
!-----------------------------------------------------------------------

    CALL hdfrd3d(sd_id,"qr",nx,ny,nz,qr,istat,itmp,hmax,hmin)
    IF (istat /= 0) THEN
      success=.false.
      qr = 0.
    END IF

!    IF(icein == 1) THEN
      CALL hdfrd3d(sd_id,"qs",nx,ny,nz,qs,istat,itmp,hmax,hmin)
    IF (istat /= 0) THEN
      success=.false.
      qs = 0.
    END IF

      CALL hdfrd3d(sd_id,"qh",nx,ny,nz,qh,istat,itmp,hmax,hmin)
    IF (istat /= 0) THEN
      success=.false.
      qh=0.
    END IF
!    ELSE
!      qs=0.
!      qh=0.
!    END IF

!-----------------------------------------------------------------------
!
!  Friendly exit message
!
!-----------------------------------------------------------------------

    IF(success) THEN
      WRITE(6,'(/a,F8.1,a/)')                                           &
        ' Data at time=', time/60,' (min) were successfully read.'

      ireturn = 0
    ELSE
      WRITE(6,'(/a,F8.1,a/)')                                           &
        ' Error reading qr,qs,qh data at time=', time/60,' (min)'
      ireturn = -3
    END IF
    CALL hdfclose(sd_id,istat)
  ELSE
    WRITE(6,'(/a/a,a/)') ' Error reading data in HDFREADQRSH:',          &
      ' No time-dependent data in file',trim(filename)
    ireturn = -2
    CALL hdfclose(sd_id,istat)
  END IF

  print *, 'Exiting hdfreadqrsh'
  RETURN
END SUBROUTINE hdfreadqrsh

SUBROUTINE radverr(ngate,maxazim,nazim,radv,attrefl,stdvr,gatesp, &
                   vrerropt,sigmavr,samploi,vnyqstl,              &
                   wavelen,pwrxmt,gaindb,lossdb,noisedbm,beamwid, &
                   pulselen,npulse,rmisval)
!
!  Calculate an expected error for reflectivity based in the SNR
!  of the attenuated signal and apply to the relectivity arrays.
!
!  Calculation of expected radial velocity error based on SNR,
!  estimated spectrum width and number of samples is derived
!  from a fit to the modelled errors presented in Fig 6.5 of
!  Doviak and Zrnic', 1993.
!
!  Keith Brewster, CAPS
!  29 Sept 2005
!
  IMPLICIT NONE
  INTEGER :: ngate
  INTEGER :: maxazim
  INTEGER :: nazim
  REAL :: radv(ngate,maxazim)
  REAL :: attrefl(ngate,maxazim)
  REAL :: stdvr(ngate,maxazim)
  REAL :: gatesp
  INTEGER :: vrerropt
  REAL :: sigmavr
  REAL :: samploi
  REAL :: vnyqstl
  REAL :: wavelen
  REAL :: pwrxmt
  REAL :: gaindb
  REAL :: lossdb
  REAL :: noisedbm
  REAL :: beamwid
  REAL :: pulselen
  INTEGER :: npulse
  REAL :: rmisval
!
! Misc. local variables
!
  REAL, PARAMETER :: kw2=0.91
  REAL :: taums,gatspkm,nplsind,nplsinv,twovnyq,twovninv
  REAL :: pi,rkm,dBZ,ze,pwrmw,noisemw,snr,snrinv,pwrcst
  REAL :: gain,gain2,loss,wavln2,beamwid2,pi5,ln2,two14,loss2,err
  REAL :: stdvrnor,sdnor,sigmav,vr,sigmvcst,sigmvlim
  INTEGER :: i,j,mid,ifold
  REAL :: randnor

!
! Calculate constants, including power constant needed to
! compute SRM from reflectivity.   Change units of radar parameters
! to those used in Eq. 4.35 of Doviak and Zrnic', 1993.
!
  mid=(nazim/2)+1
  pi=acos(-1.)
  twovnyq=2.*vnyqstl
  twovninv=1./twovnyq
  sigmvcst=twovnyq/sqrt(npulse*samploi)
  sigmvlim=0.8*vnyqstl
  taums=1.0E06*pulselen
  gatspkm=0.001*gatesp
  nplsind=samploi*npulse
  nplsinv=1./nplsind
  noisemw=10.**(0.1*noisedbm)
  gain=10.**(0.1*gaindb)
  loss=10.**(0.1*lossdb)
  gain2=gain*gain
  wavln2=wavelen*wavelen
  beamwid2=beamwid*beamwid
  pi5=pi**5
  two14=2.0**14
  ln2=log(2.)
  pwrcst=(pi5*1.0E-17*pwrxmt*gain2*taums*beamwid2*kw2)/                &
         (6.75*two14*ln2*wavln2*loss)
  IF(vrerropt == 1) THEN
    DO j=1,nazim
      DO i=1,ngate
        IF(radv(i,j) > rmisval) THEN
          vr=radv(i,j)+sigmavr*randnor()
          ifold=NINT(vr*twovninv)
          radv(i,j)=vr-ifold*twovnyq
        END IF
      END DO
    END DO
  ELSE IF(vrerropt > 1) THEN
   WRITE(6,'(a,f12.2)') ' vnyqstl=',vnyqstl
   WRITE(6,'(a,f12.0,a,f12.3)') ' gatesp=',gatesp,'gatspkm',gatspkm
!  OPEN(51,file='snrp20.txt',status='unknown',form='formatted')
!  WRITE(51,'(6x,a)') &
!   'rkm   attrefl       snr        stdvr      sigmav         vr'
!  OPEN(52,file='snrp15.txt',status='unknown',form='formatted')
!  WRITE(52,'(6x,a)') &
!   'rkm   attrefl       snr        stdvr      sigmav         vr'
!  OPEN(53,file='snrp08.txt',status='unknown',form='formatted')
!  WRITE(53,'(6x,a)') &
!   'rkm   attrefl       snr        stdvr      sigmav         vr'
!  OPEN(54,file='snrp00.txt',status='unknown',form='formatted')
!  WRITE(54,'(6x,a)') &
!   'rkm   attrefl       snr        stdvr      sigmav         vr'
!  OPEN(55,file='snrm05.txt',status='unknown',form='formatted')
!  WRITE(55,'(6x,a)') &
!   'rkm   attrefl       snr        stdvr      sigmav         vr'
!  OPEN(56,file='snrm06.txt',status='unknown',form='formatted')
!  WRITE(56,'(6x,a)') &
!   'rkm   attrefl       snr        stdvr      sigmav         vr'
    DO j=1,nazim
      DO i=1,ngate
        IF(radv(i,j) > rmisval) THEN
          IF(attrefl(i,j) > rmisval ) THEN
            rkm=i*gatspkm
            ze=10.**(attrefl(i,j)*0.1)
            pwrmw=pwrcst*ze/(rkm*rkm)
            snr=pwrmw/noisemw
          ELSE
            rkm=i*gatspkm
            snr=0.01
          END IF
          snr=max(snr,0.01)
          snrinv=1.0/snr
          stdvrnor=stdvr(i,j)*twovninv
          stdvrnor=min(max(stdvrnor,0.01),0.4)
          sdnor=(0.0177*exp(12.0*stdvrnor)) +                             &
                ((0.159-0.755*stdvrnor+6.5*stdvrnor*stdvrnor)*snrinv)
          sigmav=sigmvcst*sdnor
          sigmav=min(max(sigmav,0.1),sigmvlim)
!         IF(snr > 100.) THEN   ! 20 dB
!           WRITE(51,'(f10.1,2x,g10.5,2x,g10.5,2x,g10.5,2x,g10.4,f10.2)') &
!             rkm,attrefl(i,j),snr,stdvrnor,sdnor,sigmav
!         ELSE IF(snr > 15.85) THEN   ! 12 dB
!           WRITE(52,'(f10.1,2x,g10.5,2x,g10.5,2x,g10.5,2x,g10.4,f10.2)') &
!             rkm,attrefl(i,j),snr,stdvrnor,sdnor,sigmav
!         ELSE IF(snr > 3.16) THEN  ! 5 dB
!           WRITE(53,'(f10.1,2x,g10.5,2x,g10.5,2x,g10.5,2x,g10.4,f10.2)') &
!             rkm,attrefl(i,j),snr,stdvrnor,sdnor,sigmav
!         ELSE IF(snr > 0.631) THEN  ! -2 dB
!           WRITE(54,'(f10.1,2x,g10.5,2x,g10.5,2x,g10.5,2x,g10.4,f10.2)') &
!             rkm,attrefl(i,j),snr,stdvrnor,sdnor,sigmav
!         ELSE IF(snr > 0.20) THEN  ! -7 dB
!           WRITE(55,'(f10.1,2x,g10.5,2x,g10.5,2x,g10.5,2x,g10.4,f10.2)') &
!             rkm,attrefl(i,j),snr,stdvrnor,sdnor,sigmav
!         ELSE
!           WRITE(56,'(f10.1,2x,g10.5,2x,g10.5,2x,g10.5,2x,g10.4,f10.2)') &
!             rkm,attrefl(i,j),snr,stdvr(i,j),sigmav,radv(i,j)
!         END IF
          vr=radv(i,j)+sigmav*randnor()
          ifold=NINT(vr*twovninv)
          radv(i,j)=vr-ifold*twovnyq
        END IF
      END DO
    END DO
  END IF
  RETURN
END SUBROUTINE radverr
SUBROUTINE reflerr(ngate,maxazim,nazim,refl,attrefl,gatesp,          &
                     rferropt,sigmarf,refmin,samploi,                  &
                     wavelen,pwrxmt,gaindb,lossdb,noisedbm,beamwid,    &
                     pulselen,npulse,rmisval)
!
!  Calculate an expected error for reflectivity based in the SNR
!  of the attenuated signal and apply to the relectivity arrays.
!
!  Keith Brewster, CAPS
!  29 Sept 2005
!
  IMPLICIT NONE
  INTEGER :: ngate
  INTEGER :: maxazim
  INTEGER :: nazim
  REAL :: refl(ngate,maxazim)
  REAL :: attrefl(ngate,maxazim)
  REAL :: gatesp
  INTEGER :: rferropt
  REAL :: sigmarf
  REAL :: refmin
  REAL :: samploi
  REAL :: wavelen
  REAL :: pwrxmt
  REAL :: gaindb
  REAL :: lossdb
  REAL :: noisedbm
  REAL :: beamwid
  REAL :: pulselen
  INTEGER :: npulse
  REAL :: rmisval
!
! Misc. local variables
!
  REAL, PARAMETER :: kw2=0.91
  REAL :: taums,gatspkm,nplsind,nplsinv,rminze
  REAL :: pi,rkm,rkm2,dBZ,ze,pwrmw,noisemw,randerr
  REAL :: snr,sigmaz,snrinv,pwrcst,pcstinv
  REAL :: gain,gain2,loss,wavln2,beamwid2,pi5,ln2,two14,loss2,err
  INTEGER :: i,j,mid
  REAL :: randnor

  WRITE(6,'(a,i4)') ' Calculating Reflectivity Error, Opt=',rferropt

  mid=(nazim/2)+1
  pi=acos(-1.)
  taums=1.0E06*pulselen
  gatspkm=0.001*gatesp
  nplsind=samploi*npulse
  nplsinv=1./nplsind
  noisemw=10.**(0.1*noisedbm)
  gain=10.**(0.1*gaindb)
  loss=10.**(0.1*lossdb)
  gain2=gain*gain
  wavln2=wavelen*wavelen
  beamwid2=beamwid*beamwid
  pi5=pi**5
  two14=2.0**14
  ln2=log(2.)
  pwrcst=(pi5*1.0E-17*pwrxmt*gain2*taums*beamwid2*kw2)/  &
         (6.75*two14*ln2*wavln2*loss)
  pcstinv=1.0/pwrcst

  IF(rferropt == 1) THEN
    DO j=1,nazim
      DO i=1,ngate
        IF(attrefl(i,j) > rmisval) THEN
          err=sigmarf*randnor()
          refl(i,j)=max(refmin,(refl(i,j)+err))
          attrefl(i,j)=max(refmin,(attrefl(i,j)+err))
        END IF
      END DO
    END DO
  ELSE IF(rferropt > 1) THEN
!    WRITE(6,'(6x,a)') 'rkm      attrefl       snr         sigmaz'
    DO j=1,nazim
      DO i=1,ngate
        IF(attrefl(i,j) > rmisval) THEN
          rkm=i*gatspkm
          rkm2=rkm*rkm
          ze=10.**(attrefl(i,j)*0.1)
          pwrmw=pwrcst*ze/rkm2
          snr=pwrmw/noisemw
          snrinv=1.0/snr
          sigmaz= &
            1.0+sqrt(nplsinv*((1.+snrinv)*(1.+snrinv)+(snrinv*snrinv)))
!          IF(mod(i,100) == 0) &
!            WRITE(6,'(/f10.1,2x,g10.4,2x,g10.4,2x,g10.4)')   &
!              rkm,attrefl(i,j),snr,sigmaz
          randerr=randnor()
          err=sigmaz*noisemw*randerr
          pwrmw=pwrmw+err
          IF (pwrmw > 0. ) THEN
            ze=rkm2*pcstinv*pwrmw
            attrefl(i,j)=10.*log10(ze)
            attrefl(i,j)=max(refmin,attrefl(i,j))
          ELSE
            attrefl(i,j)=refmin
          END IF
!          IF(mod(i,100) == 0) THEN
!            WRITE(6,'(12x,g10.4,2x,g10.4)')   &
!             attrefl(i,j),pwrmw
!            print *, 'randerr, err=',randerr,err
!          END IF
          ze=10.**(refl(i,j)*0.1)
          pwrmw=pwrcst*ze/(rkm*rkm)
          pwrmw=pwrmw+err
          IF (pwrmw > 0.) THEN
            ze=rkm2*pcstinv*pwrmw
            refl(i,j)=10.*log10(ze)
            refl(i,j)=max(refmin,refl(i,j))
          ELSE
            refl(i,j)=refmin
          END IF
        END IF
      END DO
    END DO
  END IF
  RETURN
END SUBROUTINE reflerr

SUBROUTINE fake_reflec(nx,ny,nz,ref)
  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: ref(nx,ny,nz)

  INTEGER :: i,j,k
  REAL, PARAMETER :: refcst=30.0

  WRITE(6,'(//a,/,/a,f10.2//)')  '*******WARNING***************', &
                  '  Test code: Setting all reflectivity in domain to:',refcst

  DO k=1,nz
  DO j=1,ny
  DO i=1,nx
    ref(i,j,k)=refcst
  END DO
  END DO
  END DO
  RETURN
END SUBROUTINE fake_reflec
SUBROUTINE fake_vel(nx,ny,nz,usc,vsc)
  IMPLICIT NONE
  INTEGER :: nx,ny,nz
  REAL :: usc(nx,ny,nz)
  REAL :: vsc(nx,ny,nz)

  INTEGER :: i,j,k
  REAL, PARAMETER :: ucst=50.0
  REAL, PARAMETER :: vcst=0.0

  WRITE(6,'(//a,/,/a,f10.2,/a,f10.2//)')  &
           '*******WARNING***************', &
           '  Test code: Setting all u velocity in domain to:',ucst, &
           '  Test code: Setting all v velocity in domain to:',vcst

  DO k=1,nz
  DO j=1,ny
  DO i=1,nx
    usc(i,j,k)=ucst
    vsc(i,j,k)=vcst
  END DO
  END DO
  END DO
  RETURN
END SUBROUTINE fake_vel
!
SUBROUTINE QUADINT(nx,ny,nz,var,ilc,jlc,k,dx,dy,xpt,ypt,whigh,wlow,varint)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nx,ny,nz
  REAL, INTENT(IN) :: var(nx,ny,nz)
  INTEGER, INTENT(IN) :: ilc,jlc,k
  REAL, INTENT(IN) :: dx,dy
  REAL, INTENT(IN) :: xpt,ypt
  REAL, INTENT(IN) :: whigh,wlow
  REAL, INTENT(OUT) :: varint

  REAL :: twdxinv,twdx2inv
  REAL :: twdyinv,twdy2inv

  REAL :: a,b,fjm1,fj0,fjp1,flow,fhigh

  twdx2inv=1./(2.*dx*dx)
  twdxinv=1./(2.*dx)
  twdy2inv=1./(2.*dy*dy)
  twdyinv=1./(2.*dy)

  a=(var(ilc-1,jlc-1,k-1)-2.*var(ilc,jlc-1,k-1)+ &
     var(ilc+1,jlc-1,k-1))*twdx2inv
  b=(var(ilc+1,jlc-1,k-1)-var(ilc-1,jlc-1,k-1))*twdxinv
  fjm1=a*xpt*xpt+b*xpt+var(ilc,jlc-1,k-1)
  a=(var(ilc-1,jlc,k-1)-2.*var(ilc,jlc,k-1)+     &
     var(ilc+1,jlc,k-1))*twdx2inv
  b=(var(ilc+1,jlc,k-1)-var(ilc-1,jlc,k-1))*twdxinv
  fj0=a*xpt*xpt+b*xpt+var(ilc,jlc,k-1)
  a=(var(ilc-1,jlc+1,k-1)-2.*var(ilc,jlc+1,k-1)+ &
     var(ilc+1,jlc+1,k-1))*twdx2inv
  b=(var(ilc+1,jlc+1,k-1)-var(ilc-1,jlc+1,k-1))*twdxinv
  fjp1=a*xpt*xpt+b*xpt+var(ilc,jlc+1,k-1)
  a=(fjm1-2.*fj0+fjp1)*twdy2inv
  b=(fjp1-fjm1)*twdyinv
  flow=a*ypt*ypt + b*ypt + fj0

  a=(var(ilc-1,jlc-1,k)-2.*var(ilc,jlc-1,k)+     &
     var(ilc+1,jlc-1,k))*twdx2inv
  b=(var(ilc+1,jlc-1,k)-var(ilc-1,jlc-1,k))*twdxinv
  fjm1=a*xpt*xpt+b*xpt+var(ilc,jlc-1,k)
  a=(var(ilc-1,jlc,k)-2.*var(ilc,jlc,k)+         &
     var(ilc+1,jlc,k))*twdx2inv
  b=(var(ilc+1,jlc,k)-var(ilc-1,jlc,k))*twdxinv
  fj0=a*xpt*xpt+b*xpt+var(ilc,jlc,k)
  a=(var(ilc-1,jlc+1,k)-2.*var(ilc,jlc+1,k)+     &
     var(ilc+1,jlc+1,k))*twdx2inv
  b=(var(ilc+1,jlc+1,k)-var(ilc-1,jlc+1,k))*twdxinv
  fjp1=a*xpt*xpt+b*xpt+var(ilc,jlc+1,k)
  a=(fjm1-2.*fj0+fjp1)*twdy2inv
  b=(fjp1-fjm1)*twdyinv
  fhigh=a*ypt*ypt + b*ypt + fj0

  varint=whigh*fhigh+wlow*flow
 RETURN
 END SUBROUTINE QUADINT
