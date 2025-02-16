PROGRAM attenuation
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Sample program to read history data files produced by ARPS.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!    12/12/2003: Rewrote based on arpscvt. Input file names will be
!    specified in arpsread.input.
!
!  MODIFICATION HISTORY: edited arpsread for attenuation program
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz,nzsoil   ! Grid dimensions.
  INTEGER :: nstyps            ! Maximum number of soil types.

  INTEGER :: hinfmt,nhisfile_max,nhisfile,lengbf,nf,lenfil
  PARAMETER (nhisfile_max=200)
  CHARACTER (LEN=256) :: grdbasfn,hisfile(nhisfile_max)
!
!-----------------------------------------------------------------------
!  Include files:
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
!
!-----------------------------------------------------------------------
!  Arrays to be read in:
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: x(:)   ! The x-coord. of the physical and
                              ! computational grid. Defined at u-point.
  REAL, ALLOCATABLE :: y(:)   ! The y-coord. of the physical and
                              ! computational grid. Defined at v-point.
  REAL, ALLOCATABLE :: z(:)   ! The z-coord. of the computational
                              ! grid. Defined at w-point.

  REAL, ALLOCATABLE :: zp(:,:,:)  ! The height of the terrain.
  REAL, ALLOCATABLE :: zpsoil(:,:,:)  ! The height of the terrain.

  REAL, ALLOCATABLE :: uprt   (:,:,:) ! Perturbation u-velocity (m/s)
  REAL, ALLOCATABLE :: vprt   (:,:,:) ! Perturbation v-velocity (m/s)
  REAL, ALLOCATABLE :: wprt   (:,:,:) ! Perturbation w-velocity (m/s)
  REAL, ALLOCATABLE :: ptprt  (:,:,:) ! Perturbation potential temperature (K)
  REAL, ALLOCATABLE :: pprt   (:,:,:) ! Perturbation pressure (Pascal)
  REAL, ALLOCATABLE :: qvprt  (:,:,:) ! Perturbation water vapor specific
                                      ! humidity (kg/kg)
  REAL, ALLOCATABLE :: qscalar(:,:,:,:)

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
  REAL, ALLOCATABLE :: stypfrct(:,:,:)    ! Soil type
  INTEGER, ALLOCATABLE :: vegtyp(:,:)     ! Vegetation type
  REAL, ALLOCATABLE :: lai    (:,:)   ! Leaf Area Index
  REAL, ALLOCATABLE :: roufns (:,:)   ! Surface roughness
  REAL, ALLOCATABLE :: veg    (:,:)   ! Vegetation fraction

  REAL, ALLOCATABLE :: tsoil (:,:,:,:) ! Soil temperature (K)
  REAL, ALLOCATABLE :: qsoil (:,:,:,:) ! Soil moisture (m**3/m**3)
  REAL, ALLOCATABLE :: wetcanp(:,:,:) ! Canopy water amount
  REAL, ALLOCATABLE :: snowdpth(:,:)  ! Snow depth (m)

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
  REAL, ALLOCATABLE :: radswnet (:,:) ! Net solar radiation, SWin - SWout
  REAL, ALLOCATABLE :: radlwin  (:,:) ! Incoming longwave radiation

  REAL, ALLOCATABLE :: usflx (:,:)    ! Surface flux of u-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: vsflx (:,:)    ! Surface flux of v-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: ptsflx(:,:)    ! Surface heat flux (K*kg/(m*s**2))
  REAL, ALLOCATABLE :: qvsflx(:,:)    ! Surface moisture flux (kg/(m**2*s))

  REAL, ALLOCATABLE :: xsc(:)
  REAL, ALLOCATABLE :: ysc(:)
  REAL, ALLOCATABLE :: zps(:,:,:)
  REAL, ALLOCATABLE :: latsc(:,:)
  REAL, ALLOCATABLE :: lonsc(:,:)
  REAL, ALLOCATABLE :: azmsc(:,:)
  REAL, ALLOCATABLE :: sfcr(:,:)
  REAL, ALLOCATABLE :: cmpref(:,:)
  REAL, ALLOCATABLE :: elvsc(:,:,:)
  REAL, ALLOCATABLE :: rngsc(:,:,:)

  REAL, ALLOCATABLE :: reflect(:,:,:)  ! reflectivity (dBZ)
  REAL, ALLOCATABLE :: refz(:,:,:)  ! reflectivity (mm6/m3)
!
!-----------------------------------------------------------------------
! Temporary working arrays:
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: tem1(:,:,:)
  REAL, ALLOCATABLE :: tem2(:,:,:)
  REAL, ALLOCATABLE :: tem3(:,:,:)
  REAL, ALLOCATABLE :: tem4(:,:,:)
!
!-----------------------------------------------------------------------
! Az-ran weighting arrays:
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: delrng(:)
  REAL, ALLOCATABLE :: delazm(:)
  REAL, ALLOCATABLE :: delelv(:)
  REAL, ALLOCATABLE :: wgtpt(:,:,:)
!
!-----------------------------------------------------------------------
! Radar radial arrays
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: refl(:)
  REAL, ALLOCATABLE :: refatt(:)
  REAL, ALLOCATABLE :: sens(:)
!
!-----------------------------------------------------------------------
! Hit-Miss bookkeeping arrays
!-----------------------------------------------------------------------
!
  INTEGER, ALLOCATABLE :: score(:,:,:)
  INTEGER, ALLOCATABLE :: netscore(:,:)
  REAL, ALLOCATABLE :: tornlat(:,:)
  REAL, ALLOCATABLE :: tornlon(:,:)
!
!-----------------------------------------------------------------------
! Original Observing Radar data variables
!-----------------------------------------------------------------------
!
  CHARACTER (len=4) :: radname
  REAL :: radlat,radlon,radelv,lowest_elv,bmwidth
  INTEGER :: fillopt

  NAMELIST /radar_src/ radname,radlat,radlon,radelv, &
                       lowest_elv,bmwidth,fillopt
!
!-----------------------------------------------------------------------
! Tornado location information
!-----------------------------------------------------------------------
!
  REAL :: torn_azim,torn_rngkm,torn_refl
  NAMELIST /tornado_loc/ torn_azim,torn_rngkm,torn_refl
!
!-----------------------------------------------------------------------
! Spatial offset variables
!-----------------------------------------------------------------------
!
  INTEGER :: adjctr,adjhgt
  REAL :: ctrlatnew,ctrlonnew,hgtoffset
  REAL :: deltax,deltay
  INTEGER :: ioffbgn,ioffend,joffbgn,joffend
  NAMELIST /offsetxy/ adjctr,ctrlatnew,ctrlonnew,                       &
                      adjhgt,hgtoffset,                                 &
                      deltax,deltay,ioffbgn,ioffend,joffbgn,joffend
!
!-----------------------------------------------------------------------
! Radar Network variables
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: mxradnet = 40
  INTEGER :: nradnet,ngate,pulseopt,dualprf
  INTEGER :: npulse1,npulse2
  REAL :: elvmin,gatesp
  REAL :: pulselen1,pulselen2,pulsewid1,pulsewid2,prt1,prt2
  REAL :: rotrate,beamwid
  REAL :: radnet_lat(mxradnet),radnet_lon(mxradnet)
  REAL :: radnet_elev(mxradnet)
  CHARACTER (len=4) :: radnet(mxradnet)

  NAMELIST /radar_net/ nradnet,ngate,pulseopt,dualprf,                 &
                       elvmin,gatesp,pulselen1,pulselen2,              &
                       pulsewid1,pulsewid2,prt1,prt2,                  &
                       npulse1,npulse2,beamwid,rotrate,                &
                       radnet_lat,radnet_lon,radnet_elev,radnet

  INTEGER :: nptsrng,nptsazm,nptselv
  REAL :: evalwid
  NAMELIST /emul_specs/ nptsrng,nptsazm,nptselv,evalwid
!
!-----------------------------------------------------------------------
! Misc. internal variables
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: kntrmin=5
  REAL, PARAMETER :: rngmin=1000.
  REAL, PARAMETER :: rmisval=-99.
  REAL, PARAMETER :: c = 2.9979e08    ! speed of light m/s
  REAL, PARAMETER :: sumweps=1.0E-04

  INTEGER :: nchin,idxdir
  INTEGER :: grdbas
  INTEGER :: i,j,k,ireturn
  INTEGER :: idx,igate,iradius,jradius,irad
  INTEGER :: ibgn,iend,jbgn,jend,ipt,jpt,ioff,joff
  INTEGER :: iloc1,jloc1,iloc,jloc,ilc,jlc,klc,itgate
  INTEGER :: kpt,npulsmn
  INTEGER :: kntclr,knthob,kntvob,kntfewr,kntr,kntpt,kntintr
  INTEGER :: length, ierr, istatus
  INTEGER :: nt,nhits,nattn,nrng,kntrad

  REAL :: alatnot(2)
  REAL :: rad2deg,deg2rad
  REAL :: twdx2inv,twdxinv,twdy2inv,twdyinv
  REAL :: prtmn,pulselmn,pulsewmn
  REAL :: r11l,r21l,r12l,r22l,r11h,r21h,r12h,r22h
  REAL :: time,hgtmsl,bmwrad,bmwm,elv,azimuth,azrad
  REAL :: depth,dsdr,srange,hgt,sradius,elradius
  REAL :: hlfgtsp,efbmwid,fourln4,wasqinv,wesqinv,da,de,dr
  REAL :: xsw,ysw,ctrx,ctry,xrad,yrad,xrada,yrada,radlata,radlona
  REAL :: srgpt,azmpt,elvpt,hgtpt,sfcrpt
  REAL :: rmax,dhdr,gatspkm,delx,dely
  REAL :: gtlat,gtlon,gtx,gty,xrat,yrat,xpt,ypt
  REAL :: ipos,jpos,irloc,jrloc
  REAL :: w1i,w2i,w1j,w2j,w11,w12,w21,w22
  REAL :: wgt0,wgt1,whigh,wlow,a,b,fjp1,fj0,fjm1
  REAL :: flow,fhigh,reflz,refdbz
  REAL :: reflzkm1,reflzk,reflz0,reflz1,reflz2
  REAL :: pi,cinv,hlftau,aconst,b6,hlfctau,sigr2,sr2inv
  REAL :: wgtg,wgtr,sumwr,sumr
  REAL :: torn_rng,torn0x,torn0y,tornx,torny
  REAL :: tornsfcr,tornhgt,tornelv,tornrngi,torngate
  REAL :: tornref,tornaref,tornsen
  REAL :: hpct,apct,rpct,hsum,asum

  CHARACTER(len=24) :: varunits
  CHARACTER(len=6) :: varid
  CHARACTER(len=24) :: varname
  CHARACTER(len=180), PARAMETER :: tornfile='tornscore.txt'
  CHARACTER(len=180), PARAMETER :: tornstat='tornstat.txt'

  LOGICAL :: echo
  LOGICAL :: hit,rng
  LOGICAL :: printit
!
!-----------------------------------------------------------------------
! Functions
!-----------------------------------------------------------------------
!
  REAL :: effbmw
!
!-----------------------------------------------------------------------
! Variables used in finding height of lowest elv + half beamwidth
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: refmiss = -90.

  INTEGER :: klow

  REAL :: radarx
  REAL :: radary
  REAL :: bmtop
  REAL :: sfcrng
  REAL :: bmhgt,bmhgtmsl
  REAL :: reflow,hgtlow

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  pi=acos(-1.)
  deg2rad=pi/180.
  rad2deg=1./deg2rad
!
!-----------------------------------------------------------------------
!  Get the names of the input data files.
!-----------------------------------------------------------------------
!
  CALL get_input_file_names(5,hinfmt,grdbasfn,hisfile,nhisfile)

  !lengbf = len_trim(grdbasfn)

  CALL get_dims_from_data(hinfmt,hisfile(1),                            &
                          nx,ny,nz,nzsoil,nstyps, ireturn)

  IF (nstyps <= 0) nstyps = 1
  nstyp = nstyps

  IF( ireturn /= 0 ) THEN
    PRINT*,'Problem occured when trying to get dimensions from data.'
    PRINT*,'Program stopped.'
    STOP
  END IF

  WRITE(6,'(4(a,i5))') 'nx =',nx,', ny=',ny,', nz=',nz ,'nzsoil=',nzsoil

  ALLOCATE(x      (nx))
  ALLOCATE(y      (ny))
  ALLOCATE(z      (nz))
  ALLOCATE(zp     (nx,ny,nz))
  ALLOCATE(zpsoil (nx,ny,nzsoil))

  ALLOCATE(uprt (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "uprt")

  ALLOCATE(vprt (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "vprt")

  ALLOCATE(wprt (nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "wprt")

  ALLOCATE(ptprt  (nx,ny,nz))
  ALLOCATE(pprt   (nx,ny,nz))
  ALLOCATE(qvprt  (nx,ny,nz))
  ALLOCATE(qscalar(nx,ny,nz,nscalar), STAT = istatus)

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

  ALLOCATE(tsoil (nx,ny,nzsoil,0:nstyps))
  ALLOCATE(qsoil (nx,ny,nzsoil,0:nstyps))
  ALLOCATE(wetcanp(nx,ny,0:nstyps))
  ALLOCATE(snowdpth(nx,ny))

  ALLOCATE(raing(nx,ny))
  ALLOCATE(rainc(nx,ny))
  ALLOCATE(prcrate(nx,ny,4))

  ALLOCATE(radsw (nx,ny))
  ALLOCATE(rnflx (nx,ny))
  ALLOCATE(radswnet (nx,ny))
  ALLOCATE(radlwin (nx,ny))

  ALLOCATE(radfrc(nx,ny,nz))

  ALLOCATE(usflx (nx,ny))
  ALLOCATE(vsflx (nx,ny))
  ALLOCATE(ptsflx(nx,ny))
  ALLOCATE(qvsflx(nx,ny))
  ALLOCATE(tem1(nx,ny,nz))
  ALLOCATE(tem2(nx,ny,nz))
  ALLOCATE(tem3(nx,ny,nz))
  ALLOCATE(tem4(nx,ny,nz))

  ALLOCATE(xsc   (nx))
  ALLOCATE(ysc   (ny))
  ALLOCATE(latsc(nx,ny))
  ALLOCATE(lonsc(nx,ny))
  ALLOCATE(azmsc(nx,ny))
  ALLOCATE(sfcr(nx,ny))
  ALLOCATE(cmpref(nx,ny))

  ALLOCATE(zps    (nx,ny,nz))
  ALLOCATE(reflect(nx,ny,nz))
  ALLOCATE(refz(nx,ny,nz))
  ALLOCATE(elvsc(nx,ny,nz))
  ALLOCATE(rngsc(nx,ny,nz))

  x      =0.0
  y      =0.0
  z      =0.0
  zp     =0.0
  zpsoil =0.0

  uprt   =0.0
  vprt   =0.0
  wprt   =0.0
  ptprt  =0.0
  pprt   =0.0
  qvprt  =0.0
  qscalar(:,:,:,:) = 0.0
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

  raing=0.0
  rainc=0.0
  prcrate=0.0

  radfrc=0.0
  radsw =0.0
  rnflx =0.0
  radswnet = 0.0
  radlwin = 0.0

  usflx =0.0
  vsflx =0.0
  ptsflx=0.0
  qvsflx=0.0
  tem1=0.0
  tem2=0.0
  tem4=0.0

  xsc=0.
  ysc=0.
  latsc=0.
  lonsc=0.
  azmsc=0.
  sfcr=0.
  cmpref=0.
  zps=0.
  reflect=0.
  refz=0.
  elvsc=0.
  rngsc=0.

  lengbf=len_trim(grdbasfn)
  WRITE(6,'(/a,a)')' The grid/base name is ', grdbasfn(1:lengbf)

  lenfil=len_trim(hisfile(1))
  WRITE(6,'(/a,a,a)')                                                   &
        ' Data set ', trim(hisfile(1)),' to be converted.'
!
!-----------------------------------------------------------------------
!  Read all input data arrays
!-----------------------------------------------------------------------

  CALL dtaread(nx,ny,nz,nzsoil,nstyps,                                  &
                 hinfmt, nchin,grdbasfn(1:lengbf),lengbf,               &
                 hisfile(1)(1:lenfil),lenfil,time,                      &
                 x,y,z,zp,zpsoil, uprt ,vprt ,wprt ,ptprt, pprt ,       &
                 qvprt, qscalar, tke,kmh,kmv,                           &
                 ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,          &
                 soiltyp,stypfrct,vegtyp,lai,roufns,veg,                &
                 tsoil, qsoil, wetcanp,snowdpth,                        &
                 raing,rainc,prcrate,                                   &
                 radfrc,radsw,rnflx,radswnet,radlwin,                   &
                 usflx,vsflx,ptsflx,qvsflx,                             &
                 ireturn, tem1, tem2, tem3)

  CALL gtlfnkey( runname, lfnkey )
  idxdir=index(hisfile(1)(1:lenfil),runname(1:lfnkey),.true.)
  IF(idxdir<=1) THEN
    dirname='.'
  ELSE
    IF(hisfile(1)(idxdir-1:idxdir-1) == '/') THEN
      dirname=hisfile(1)(1:idxdir-2)
    ELSE
      dirname=hisfile(1)(1:idxdir-1)
    END IF
  END IF
  print *, ' ATTEN dirname= ',dirname
  ldirnam=LEN(dirname)
  CALL strlnth( dirname , ldirnam)

  curtim = time
  dx=(x(2)-x(1))
  dy=(y(2)-y(1))
  dxinv=1./dx
  dyinv=1./dy

  twdx2inv=1./(2.*dx*dx)
  twdxinv=1./(2.*dx)
  twdy2inv=1./(2.*dy*dy)
  twdyinv=1./(2.*dy)

  alatnot(1)=trulat1
  alatnot(2)=trulat2
  CALL setmapr(mapproj,1.0,alatnot,trulon)
  CALL lltoxy(1,1,ctrlat,ctrlon,ctrx,ctry)
  xsw=ctrx-0.5*((nx-3)*dx)
  ysw=ctry-0.5*((ny-3)*dy)
  CALL setorig(1,xsw,ysw)

  if( ireturn /= 0 ) then
     WRITE(6,'(1x,a,i2,/1x,a)')                                       &
    'Data read was unsuccessful. ireturn =', ireturn,' Job stopped.'
     STOP
  endif
!
  radname='KTLX'
  radlat=35.33
  radlon=-97.28
  radelv=350.0
  lowest_elv=0.5

  WRITE(6,'(a)') ' Reading radar_src namelist...'
  READ(5,radar_src)
  WRITE(6,'(a/)') ' Namelist radar_src successfully read.'
!
!-----------------------------------------------------------------------
! Get tornado location relative to original observing radar
!-----------------------------------------------------------------------
!
  torn_azim=300.
  torn_rngkm=65.
  torn_refl=40.

  WRITE(6,'(a)') ' Reading tornado_loc namelist...'
  READ(5,tornado_loc)
  WRITE(6,'(a/)') ' Namelist tornado_loc successfully read.'
  torn_rng=1000.0*torn_rngkm
!
!-----------------------------------------------------------------------
! Get offset control variables
!-----------------------------------------------------------------------
!
  adjctr=0
  adjhgt=0
  ctrlatnew=35.
  ctrlonnew=-98.
  hgtoffset=0.
  deltax=3000.
  deltay=3000.
  ioffbgn=-10
  ioffend=10
  joffbgn=-10
  joffend=10

  WRITE(6,'(a)') ' Reading offsetxy namelist...'
  READ(5,offsetxy)
  WRITE(6,'(a/)') ' Namelist offsetxy successfully read.'
!
!-----------------------------------------------------------------------
! Get radar network observing locations
! and operating specs for network radar
!-----------------------------------------------------------------------
!
  nradnet=1
  ngate=150
  elvmin=1.0
  gatesp=200.
  dualprf=1
  pulseopt=1
  pulselen1=1.57E-06
  pulselen2=1.57E-06
  pulsewid1=441.
  pulsewid2=441.
  prt1=1.06e-03
  prt2=1.06e-03
  npulse1=54
  npulse2=54
  beamwid = 2.0
  rotrate=18.
  evalwid=2.0

  WRITE(6,'(a)') ' Reading radar_net namelist...'
  READ(5,radar_net)
  WRITE(6,'(a/)') ' Namelist radar_net successfully read.'
!
  gatspkm=0.001*gatesp
  hlfgtsp=0.5*gatesp
  ALLOCATE(refl(ngate))
  refl=0.
  ALLOCATE(refatt(ngate))
  refatt=0.
  ALLOCATE(sens(ngate))
  sens=0.
  ALLOCATE(score(ioffbgn:ioffend,joffbgn:joffend,nradnet))
  score=-99
  ALLOCATE(netscore(ioffbgn:ioffend,joffbgn:joffend))
  score=-99
  ALLOCATE(tornlat(ioffbgn:ioffend,joffbgn:joffend))
  tornlat=-999.
  ALLOCATE(tornlon(ioffbgn:ioffend,joffbgn:joffend))
  tornlon=-999.
!
  CALL radcoord(nx,ny,nz,                                              &
                     x,y,z,zp,xsc,ysc,zps,                             &
                     radlat,radlon,radarx,radary)
  varid='reflct'
  CALL readvar1(nx,ny,nz, reflect, varid, varname,                     &
                varunits, time, runname, dirname, istatus)

  elv=elvmin

  IF(fillopt > 0 ) THEN
!
! Fill in radar reflectivity below the lowest beam height using
! zero gradient assumption.
!
    bmtop=lowest_elv+0.5*bmwidth
    DO j=1,ny
      DO i=1,nx
!
!   find height of lowest elev + half beamwidth
!
        delx=xsc(i)-radarx
        dely=ysc(j)-radary
        sfcrng=sqrt(delx*delx+dely*dely)
        CALL bmhgtsfr(bmtop,sfcrng,bmhgt)
        bmhgtmsl=bmhgt+radelv
        klow=nz
        reflow=0.
        hgtlow=1.0E06
        DO k=nz-1,2,-1
          IF(reflect(i,j,k) > refmiss) THEN
            klow=k
            reflow=reflect(i,j,k)
            hgtlow=zps(i,j,k)
          END IF
        END DO
        IF(hgtlow < bmhgtmsl) THEN
          DO k=2,klow-1
            reflect(i,j,k)=reflow
          END DO
        END IF
      END DO
    END DO
  END IF

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        IF(reflect(i,j,k) < -30.) reflect(i,j,k)=-30.
      END DO
    END DO
  END DO

  IF(fillopt > 0 ) THEN
    varid='rflctx'
    CALL wrtvar1(nx,ny,nz,reflect, varid,varname,varunits,              &
           time,runname,dirname, istatus)
  END IF

  DO j=1,ny-1
    DO i=1,nx-1
      DO k=2,nz-1
        IF(reflect(i,j,k) > -30.) THEN
          refz(i,j,k)=10.**(0.1*reflect(i,j,k))
        ELSE
          refz(i,j,k)=0.
        END IF
        cmpref(i,j)=max(cmpref(i,j),reflect(i,j,k))
      END DO
    END DO
  END DO
!
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
      pulselmn=pulselen1
      pulsewmn=c*pulselmn
    ELSE
      pulsewmn=pulsewid1
      pulselmn=cinv*pulsewmn
    END IF
    prtmn=prt1
    npulsmn=npulse1
  ELSE
    IF(pulseopt == 1) THEN
      pulselmn=0.5*(pulselen1*pulselen2)
      pulsewmn=c*pulselmn
    ELSE
      pulsewmn=0.5*(pulsewid1*pulsewid2)
      pulselmn=cinv*pulsewmn
    END IF
    prtmn=0.5*(prt1+prt2)
    npulsmn=NINT(0.5*(npulse1+npulse2))
  END IF

  hlftau=0.5*pulselmn
  aconst=pi/(2.*sqrt(log(2.0)))
  b6=1.04/pulselmn
  hlfctau=0.5*pulsewmn
  sigr2=(0.35*hlfctau)*(0.35*hlfctau)
  sr2inv=2.0/(4.0*sigr2)

  WRITE(6,'(a,f10.2,a)') ' Mean pulse width: ',pulsewmn,' m'
  WRITE(6,'(a,f10.2,a)') &
    ' Mean pulse length: ',(pulselmn*1.0E06),' microsec'
  WRITE(6,'(a,f10.2,a)') &
    ' Gaussian filtered pulse half-width: ',sqrt(sigr2),' m'

!
!-----------------------------------------------------------------------
! Calculate constants for antenna gain weighting.
!-----------------------------------------------------------------------
!
  IF( abs(rotrate) > 0.) THEN
    efbmwid=EFFBMW(beamwid,rotrate,prtmn,npulsmn)
  ELSE
    efbmwid=beamwid
  END IF

  fourln4=4.0*alog(4.0)
  wasqinv=fourln4/(efbmwid*efbmwid)
  wesqinv=fourln4/(beamwid*beamwid)
  bmwrad=deg2rad*efbmwid
  elradius=evalwid*beamwid
  hlfgtsp=0.5*gatesp

  nptsrng=max(nptsrng,3)
  ALLOCATE (delrng(nptsrng))
  dr=gatesp/float(nptsrng-1)
  WRITE(6,'(a,f10.2,a)') ' Range evaluation points (gatesp=',gatesp,'):'
  DO ipt=1,nptsrng
    delrng(ipt)=(-0.5*gatesp)+((ipt-1)*dr)
    WRITE(6,'(a,i4,a,f12.4)') ' delrng(',ipt,') = ',delrng(ipt)
  END DO

  nptsazm=max(nptsazm,3)
  ALLOCATE (delazm(nptsazm))

  nptselv=max(nptselv,3)
  ALLOCATE (delelv(nptselv))
  de=(evalwid*beamwid)/float(nptselv-1)
  WRITE(6,'(a,f10.2,a,f10.2)') &
   ' Beam width ',beamwid,' Eval region width:',(evalwid*beamwid)
  WRITE(6,'(a)') &
   ' Cross-beam evaluation pts (elev deg from center):'
  DO kpt=1,nptselv
    delelv(kpt)=(-0.5*evalwid*beamwid)+((kpt-1)*de)
    WRITE(6,'(a,i4,a,f12.4)') ' delelv(',kpt,') = ',delelv(kpt)
  END DO

  WRITE(6,'(a,f8.2,a)') ' Static beamwidth:',beamwid,' degrees.'
  WRITE(6,'(a,f8.2,a)') ' Effective beamwidth:',efbmwid,' degrees.'

  da=(evalwid*efbmwid)/float(nptsazm-1)
  WRITE(6,'(a)') &
       ' Cross-beam evaluation pts (azimuth degrees from center):'
  DO jpt=1,nptsazm
    delazm(jpt)=(-0.5*evalwid*efbmwid)+((jpt-1)*da)
    WRITE(6,'(a,i4,a,f12.4)') ' delazm(',jpt,') = ',delazm(jpt)
  END DO

  wasqinv=fourln4/(efbmwid*efbmwid)
  bmwrad=deg2rad*efbmwid
  ALLOCATE (wgtpt(nptsrng,nptsazm,nptselv))

  DO kpt=1,nptselv
    DO jpt=1,nptsazm
      DO ipt=1,nptsrng
        wgtpt(ipt,jpt,kpt)=exp(-((delrng(ipt)*delrng(ipt))*sr2inv  +    &
                                 (delazm(jpt)*delazm(jpt))*wasqinv +    &
                                 (delelv(kpt)*delelv(kpt))*wesqinv))
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
! If no map projection has been defined, for the purposes of emulation,
! set the projection to Lambert conformal and center the grid on
! ctrlat, ctrlon.
!
!-----------------------------------------------------------------------
!
  IF ( mapproj == 0 ) THEN
    IF ( adjctr > 0 ) THEN
      ctrlat = ctrlatnew
      ctrlon = ctrlonnew
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
    ctrlat = ctrlatnew
    ctrlon = ctrlonnew
    alatnot(1)=trulat1
    alatnot(2)=trulat2
    CALL setmapr(mapproj,1.0,alatnot,trulon)
    CALL lltoxy(1,1,ctrlat,ctrlon,ctrx,ctry)
    xsw=ctrx-0.5*((nx-3)*dx)
    ysw=ctry-0.5*((ny-3)*dy)
    CALL setorig(1,xsw,ysw)
  END IF

  torn0x=radarx+torn_rng*sin(deg2rad*torn_azim)
  torn0y=radary+torn_rng*cos(deg2rad*torn_azim)
!
!-----------------------------------------------------------------------
! Adjust grid heights to be compatible with the netrad domain
!-----------------------------------------------------------------------
!
  IF(adjhgt > 0) THEN
    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          zps(i,j,k)=zps(i,j,k)+hgtoffset
        END DO
      END DO
    END DO
  END IF
!
!-----------------------------------------------------------------------
! Calculate sensitivity values for all distances from radar.
!-----------------------------------------------------------------------
!
  CALL xbsens(ngate,gatesp,sens)
!
!-----------------------------------------------------------------------
! Loop through all possible offsets
!-----------------------------------------------------------------------
!
  DO joff=joffbgn,joffend
  DO ioff=ioffbgn,ioffend
!
!-----------------------------------------------------------------------
!
! Calculate new tornado location for this storm offset.
!
!-----------------------------------------------------------------------
!
    printit=((mod(ioff,10) == 0) .AND. (mod(joff,10) == 0))
    tornx=torn0x+(ioff*deltax)
    torny=torn0y+(joff*deltay)
    CALL xytoll(1,1,tornx,torny,tornlat(ioff,joff),tornlon(ioff,joff))
    IF (printit) &
      WRITE(6,'(a,f12.4,a,f12.4)') ' Tornado lat: ',tornlat(ioff,joff),   &
                                   ' Tornado lon: ',tornlon(ioff,joff)
!
    DO irad=1,nradnet

      IF (printit) &
      WRITE(6,'(a,a)') ' Processing network radar: ',radnet(irad)
!
! Set up radar location info.
! For now, consider just one radar.
! ioff, joff refer to an offset of the grid and hence storm location
! This is accounted for by moving the radar locations by an equal
! distance in the opposite direction.
!
      CALL lltoxy(1,1,radnet_lat(irad),radnet_lon(irad),xrad,yrad)
      xrada=xrad-umove*curtim
      yrada=yrad-vmove*curtim
      xrada=xrada-(ioff*deltax)
      yrada=yrada-(joff*deltay)
      CALL xytoll(1,1,xrada,yrada,radlata,radlona)
      irloc=((xrada-xsc(1))*dxinv)+1.0
      jrloc=((yrada-ysc(1))*dyinv)+1.0
      IF (printit) &
        WRITE(6,'(a,i3,a,i3,/,a,f10.3,a,f10.3,/a,f10.3,a,f10.3)') &
          ' Ioffset ',ioff,'  Joffset ',joff, &
          ' x-coord (km) ',(0.001*xrada),'  y-coord (km) ',(0.001*yrada),&
          ' i-location: ',irloc,'  j-location:',jrloc
!
!-----------------------------------------------------------------------
!
! Calculate radar parameters at each scalar grid point.
!
!-----------------------------------------------------------------------
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
                CALL beamelv((zps(i,j,k)-radelv),sfcr(i,j),elvsc(i,j,k),rngsc(i,j,k))
              END DO
            ELSE
              rngsc(i,j,k)=sfcr(i,j)
              elvsc(i,j,k)=-99.
            END IF
          END DO
        END DO
      END IF
!
! Calculate beam heading to the position of the tornado
!
      delx=torn0x-xrada
      dely=torn0y-yrada
      IF(abs(dely) .gt. 1.0E-06) THEN
        IF(dely .gt. 0.) THEN
          azrad=atan(delx/dely)
        ELSE
          azrad=pi+atan(delx/dely)
        END IF
      ELSE IF(delx .gt. 0.) THEN
        azrad=0.5*pi
      ELSE
        azrad=1.5*pi
      END IF

      azimuth=rad2deg*azrad
      IF(azimuth < 0.) azimuth=azimuth+360.
      IF(azimuth > 360.) azimuth=azimuth-360.
      tornsfcr=sqrt((delx*delx)+(dely*dely))
      CALL bmhgtsfr(elv,tornsfcr,tornhgt)
      CALL beamelv(tornhgt,tornsfcr,tornelv,tornrngi)
      IF (printit) THEN
        WRITE(6,'(a,3(f10.2))') ' torn delx,dely,azimuth: ', &
              (0.001*delx),(0.001*dely),azimuth
        WRITE(6,'(a,3(f10.2))') ' torn elv in, elv out: ', &
              elv,tornelv
        WRITE(6,'(a,3(f10.2))') ' torn sfcrange, range: ',  &
              (0.001*tornsfcr),(0.001*tornrngi)
      END IF
      torngate=(tornrngi/gatesp)
      IF(torngate < float(ngate)) THEN
!
! For this radial, fill a 1d-array with the reflectivity data
!
        kntpt=0
        kntclr=0
        knthob=0
        kntvob=0
        kntintr=0
        kntfewr=0
        DO igate=1, ngate
          kntpt=kntpt+1
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

          IF(iloc1 > 0 .AND. iloc1 < nx-1 .AND.                         &
             jloc1 > 0 .AND. jloc1 < ny-1 ) THEN
            hgtmsl=hgt+radelv
            IF( hgtmsl > zps(iloc1,jloc1,2) .AND.                       &
                hgtmsl < zps(iloc1,jloc1,nz-1) ) THEN
              DO k=3,nz-2
                IF(zps(iloc1,jloc1,k) > hgtmsl) EXIT
              END DO
!
!  Check to see of the entire observation volume area is clear.
!
              bmwm=bmwrad*evalwid*srange
              sradius=max(bmwm,pulsewmn)
              iradius=1+INT(sradius/dx)
              jradius=1+INT(sradius/dy)
!             print *, 'bmwm,pulsewid = ',bmwm,pulsewmn
!             print *, 'iradius,jradius = ',iradius,jradius
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

              IF(echo) THEN
!
!  Calculate beam-weighted sums for an uniform array of points around this
!  range-gate using tri-linear interpolation
!
                kntr=0
                sumr=0.
                sumwr=0.
                DO kpt=1,nptselv
                DO jpt=1,nptsazm
                DO ipt=1,nptsrng
                  srgpt=(igate*gatesp)+delrng(ipt)
                  azmpt=azimuth+delazm(jpt)
                  elvpt=elv+delelv(kpt)
                  CALL beamhgt(elvpt,srgpt,hgtpt,sfcrpt)
                  IF(mapproj > 0 ) THEN
                    CALL gcircle(radlata,radlona,azmpt,sfcrpt,gtlat,gtlon)
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

                  IF(iloc1 > 0 .AND. iloc1 < nx-1 .AND.                &
                     jloc1 > 0 .AND. jloc1 < ny-1 ) THEN

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

                    hgtmsl=hgtpt+radelv
                    IF( hgtmsl > zps(iloc1,jloc1,2) .AND.              &
                        hgtmsl < zps(iloc1,jloc1,nz-1) ) THEN
                      DO k=3,nz-2
                        IF(zps(iloc1,jloc1,k) > hgtmsl) EXIT
                      END DO
                      whigh=(hgtmsl-zps(iloc1,jloc1,k-1))/             &
                          (zps(iloc1,jloc1,k)-zps(iloc1,jloc1,k-1))
                      wlow=1.-whigh

                      a=(refz(ilc-1,jlc-1,k-1)-2.*refz(ilc,jlc-1,k-1)+ &
                         refz(ilc+1,jlc-1,k-1))*twdx2inv
                      b=(refz(ilc+1,jlc-1,k-1)-refz(ilc-1,jlc-1,k-1))  &
                         *twdxinv
                      fjm1=a*xpt*xpt+b*xpt+refz(ilc,jlc-1,k-1)
                      a=(refz(ilc-1,jlc,k-1)-2.*refz(ilc,jlc,k-1)+     &
                         refz(ilc+1,jlc,k-1))*twdx2inv
                      b=(refz(ilc+1,jlc,k-1)-refz(ilc-1,jlc,k-1))      &
                         *twdxinv
                      fj0=a*xpt*xpt+b*xpt+refz(ilc,jlc,k-1)
                      a=(refz(ilc-1,jlc+1,k-1)-2.*refz(ilc,jlc+1,k-1)+ &
                         refz(ilc+1,jlc+1,k-1))*twdx2inv
                      b=(refz(ilc+1,jlc+1,k-1)-refz(ilc-1,jlc+1,k-1))  &
                         *twdxinv
                      fjp1=a*xpt*xpt+b*xpt+refz(ilc,jlc+1,k-1)

                      a=(fjm1-2.*fj0+fjp1)*twdy2inv
                      b=(fjp1-fjm1)*twdyinv
                      flow=a*ypt*ypt + b*ypt + fj0

                      a=(refz(ilc-1,jlc-1,k)-2.*refz(ilc,jlc-1,k)+     &
                         refz(ilc+1,jlc-1,k))*twdx2inv
                      b=(refz(ilc+1,jlc-1,k)-refz(ilc-1,jlc-1,k))      &
                         *twdxinv
                      fjm1=a*xpt*xpt+b*xpt+refz(ilc,jlc-1,k)
                      a=(refz(ilc-1,jlc,k)-2.*refz(ilc,jlc,k)+         &
                         refz(ilc+1,jlc,k))*twdx2inv
                      b=(refz(ilc+1,jlc,k)-refz(ilc-1,jlc,k))*twdxinv
                      fj0=a*xpt*xpt+b*xpt+refz(ilc,jlc,k)
                      a=(refz(ilc-1,jlc+1,k)-2.*refz(ilc,jlc+1,k)+     &
                         refz(ilc+1,jlc+1,k))*twdx2inv
                      b=(refz(ilc+1,jlc+1,k)-refz(ilc-1,jlc+1,k))      &
                         *twdxinv
                      fjp1=a*xpt*xpt+b*xpt+refz(ilc,jlc+1,k)

                      a=(fjm1-2.*fj0+fjp1)*twdy2inv
                      b=(fjp1-fjm1)*twdyinv
                      fhigh=a*ypt*ypt + b*ypt + fj0

                      reflz2=whigh*fhigh+wlow*flow
!
!   TEST CODE
!
                      reflzk=refz(iloc1,jloc1,k)
                      reflzkm1=refz(iloc1,jloc1,k-1)
                      reflz0=whigh*reflzk+wlow*reflzkm1
!
!   TEST CODE 2
!
                      reflz1=whigh*(w11*refz(iloc,jloc,k) +         &
                                   w21*refz(iloc+1,jloc,k) +        &
                                   w12*refz(iloc,jloc+1,k) +        &
                                   w22*refz(iloc+1,jloc+1,k)) +     &
                             wlow*(w11*refz(iloc,jloc,k-1) +        &
                                   w21*refz(iloc+1,jloc,k-1) +      &
                                   w12*refz(iloc,jloc+1,k-1) +      &
                                   w22*refz(iloc+1,jloc+1,k-1))

            WRITE(36,'(5f12.1)') reflzkm1,reflzk,reflz0,reflz1,reflz2
                      kntr=kntr+1
                      sumr=sumr+reflz0*wgtpt(ipt,jpt,kpt)
                      sumwr=sumwr+wgtpt(ipt,jpt,kpt)

                    END IF ! hgt in grid
                  END IF ! iloc in grid
                END DO
                END DO
                END DO
!
!  Calculate reflectivity and velocity from the uniform
!  array center on this range gate.
!
                IF(sumwr > 0. .AND. kntr > kntrmin) THEN
                  kntintr=kntintr+1
                  reflz=sumr/sumwr
                  refdbz=10.*alog10(reflz)
                  refl(igate)=max(refdbz,0.)
                ELSE
                  kntfewr=kntfewr+1
                  refl(igate)=0.
                END IF
              ELSE ! no echo found in area
                kntclr=kntclr+1
                refl(igate)=0.
              END IF
            ELSE ! outside vertical limits of grid
              kntvob=kntvob+1
            END IF
          ELSE ! outside grid domain
            knthob=knthob+1
          END IF
        END DO  ! igate

!-----------------------------------------------------------------------
! Calculate attenuation values for all distances from radar.
!-----------------------------------------------------------------------
!
        CALL xbatten(ngate,1,refl,gatesp,rmisval,refatt)

!-----------------------------------------------------------------------
! Output R(km), K'(dB), K(dB), Z(dB), Z'(dB), and Ze_min.
!-----------------------------------------------------------------------

        IF(ioff == 0 .AND. joff == 0) THEN
        WRITE(6,'(5(1x,a10))') &
                'R(km)','Z(dB)','Za(dB)','Za-Sens','Sens'
        WRITE(6,'(a)')                                         &
        '-----------------------------------------------------------------'

        DO idx = 5, ngate, 5
          WRITE(6,'(5(1x,f10.3))') &
            (idx*gatspkm),refl(idx),refatt(idx), &
            (refatt(idx)-sens(idx)),sens(idx)
        END DO
        END IF

!-----------------------------------------------------------------------
! Decide if the tornado is detectable from this radar.
!-----------------------------------------------------------------------
        itgate=INT(torngate)
        wgt1=torngate-float(itgate)
        wgt0=1.-wgt1
        tornref=wgt0*refl(itgate)+wgt1*refl(itgate+1)
        tornaref=wgt0*refatt(itgate)+wgt1*refatt(itgate+1)
        tornsen=wgt0*sens(itgate)+wgt1*sens(itgate+1)
        IF(printit) &
          WRITE(6,'(a,f10.2,/a,f10.2,a,f10.2,a,f10.2)') &
            ' At tornado location:  Refl=',tornref, &
            ' Att Ref=',tornaref,' Sens=',tornsen,' Att Ref-Sens=',&
            (tornaref-tornsen)
        IF(tornaref > tornsen) THEN
          score(ioff,joff,irad)=1      ! hit, detectable
        ELSE
          score(ioff,joff,irad)=0      ! attenuated, undetectable
        END IF

      ELSE
        score(ioff,joff,irad)=-1     ! out of range.undetectable
      END IF

    END DO  ! i-radar

  END DO  ! i-offset
  END DO  ! j-offset


!
! Compute summary scores of hits-n-misses
! First summarize by individual radar
!
  print *, ' writing tornstat file'
  open(42,file=tornstat,status='unknown',form='formatted')
  WRITE(42,'(a)') TRIM(runname)
  WRITE(42,'(i4.4,a,i2.2,a,i2.2,1x,i2.2,a,i2.2,a,i2.2,a)') &
    year,'-',month,'-',day,hour,':',minute,':',second,' UTC'
  WRITE(42,'(a,f8.1,a,f8.1,a,f8.1)')  &
    ' Tornado Azim:',torn_azim,' Range:',torn_rngkm,' Refl:',torn_refl
  WRITE(42,'(/a)') '  Individual Radar Statistics'
  WRITE(42,'(a)') '  Radar    Trials     Hits      Miss   Out-of-Range'
  kntrad=0
  hsum=0.
  asum=0.
  DO irad=1,nradnet
    nt=0
    nhits=0
    nattn=0
    nrng=0
    DO ioff=ioffbgn,ioffend
      DO joff=joffbgn,joffend
        IF(score(ioff,joff,irad) > -90) THEN
          nt=nt+1
          IF(score(ioff,joff,irad) == 1) THEN
            nhits=nhits+1
          ELSE IF (score(ioff,joff,irad) == 0) THEN
            nattn=nattn+1
          ELSE IF (score(ioff,joff,irad) == -1) THEN
            nrng=nrng+1
          END IF
        END IF
      END DO
    END DO

    IF((nt-nrng) > 0 ) THEN
      hpct=100.*float(nhits)/float(nt-nrng)
      apct=100.*float(nattn)/float(nt-nrng)
      hsum=hsum+hpct
      asum=asum+apct
      kntrad=kntrad+1
    ELSE
      hpct=0.
      apct=0.
    END IF
    rpct=100.*float(nrng)/float(nt)
    WRITE(42,'(a4,i6,3f12.1)') &
      radnet(irad),nt,hpct,apct,rpct
    IF(kntrad > 0) THEN
      hsum=hsum/float(kntrad)
      asum=hsum/float(kntrad)
    END IF
    WRITE(42,'(a4,i6,2f12.1)') &
      'AVG:',kntrad,hsum,asum
  END DO
!
! Compute hits and misses for the network as a whole
!
  WRITE(6,'(a)') ' Computing hits and misses'
  nt=0
  nhits=0
  nattn=0
  nrng=0
  DO joff=joffbgn,joffend
    DO ioff=ioffbgn,ioffend
      nt=nt+1
      hit=.false.
      rng=.false.
      DO irad=1,nradnet
        IF(score(ioff,joff,irad) == 1) hit=.true.
        IF(score(ioff,joff,irad) /= -1) rng=.true.
      END DO
      IF(hit) THEN
        nhits=nhits+1
        netscore(ioff,joff)=1      ! hit, detectable
      ELSE IF (.not.rng) THEN
        nrng=nrng+1
        netscore(ioff,joff)=-1      ! out-of-range
      ELSE
        nattn=nattn+1
        netscore(ioff,joff)=0      ! network miss
      END IF
    END DO
  END DO

  IF((nt-nrng) > 0) THEN
    hpct=100.*float(nhits)/float(nt-nrng)
    apct=100.*float(nattn)/float(nt-nrng)
  ELSE
    hpct=0.
    apct=0.
  END IF
  rpct=100.*float(nrng)/float(nt)

  print *, ' writing tornfile file'
  open(41,file=tornfile,status='unknown',form='formatted')
  WRITE(41,'(a)') TRIM(runname)
  WRITE(41,'(i4.4,a,i2.2,a,i2.2,1x,i2.2,a,i2.2,a,i2.2,a)') &
    year,'-',month,'-',day,hour,':',minute,':',second,' UTC'
  WRITE(41,'(a,f8.1,a,f8.1,a,f8.1)')  &
    ' Tornado Azim:',torn_azim,' Range:',torn_rngkm,' Refl:',torn_refl
  DO joff=joffbgn,joffend
    DO ioff=ioffbgn,ioffend
      WRITE(41,'(2i6,2f12.4,13i5)') &
      ioff,joff,tornlat(ioff,joff),tornlon(ioff,joff), &
      netscore(ioff,joff), &
      (score(ioff,joff,irad),irad=1,nradnet)
    END DO
  END DO
  CLOSE(41)

  WRITE(42,'(//a)') 'Network Statistics'
  WRITE(42,'(a)') 'Trials    Hits       Miss      Out-of-Range'
  WRITE(42,'(i6,3f12.1)') nt,hpct,apct,rpct
  CLOSE(42)

  CLOSE(9)

  STOP
END PROGRAM attenuation
!-----------------------------------------------------------------------
! Subroutine atten
!-----------------------------------------------------------------------
!
SUBROUTINE atten(all_rngs,all_refs,seg,att_values,att)
  IMPLICIT NONE

  INTEGER :: idx         ! Index.
  INTEGER :: seg         ! Number of distance segments.
  REAL :: att_tally      ! Running tally of attenuation values.
  REAL :: att_const1     ! Attenuation constant 1.
  REAL :: att_const2     ! Attenuation constant 2.
  REAL :: att_const3     ! Attenuation constant 3.

  REAL, DIMENSION(seg),INTENT (IN):: all_rngs      ! Range values for all dist from radar.
  REAL, DIMENSION(seg),INTENT (IN):: all_refs      ! Reflectivity values for all dist from radar.
  REAL, DIMENSION(seg),INTENT (OUT) :: att_values  ! Attenuation values for all dist from radar.
  REAL :: att(seg)                                 ! Attenuation value for a specific distance.
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  att_const1 = 2.0*0.01
  att_const2 = 1.0/400.0
  att_const3 = 1.21/1.4

  att_tally = 0

  DO idx = 1, seg
    att(idx) = &
      att_const1*((att_const2*(10.0**(0.1*all_refs(idx))))**att_const3)
    IF(idx == 1) THEN
      att_tally = att_tally + att(idx)
      att_values(idx) = att_tally
    ELSE IF(all_refs(idx) == all_refs(idx-1)) THEN
      att_tally = att_tally + att(idx)
      att_values(idx) = att_tally
    ELSE
      att_tally = att_tally + .75*att(idx-1) + &
        .25*att_const1*((att_const2*(10.0**(0.1*all_refs(idx))))**att_const3)
      att_values(idx) = att_tally
    END IF
  END DO
END SUBROUTINE atten

!-----------------------------------------------------------------------
! Subroutine sens
!-----------------------------------------------------------------------
!
!
SUBROUTINE xbsens(ngate,gatesp,sens)
!
!-----------------------------------------------------------------------
! Subroutine xbsens
! Calculate X-band radar sensitivity (dBZ) as a function of range.
!
! Based on information provided by Francesc Joyvount
!
! Keith Brewster and Erin Fay, CAPS
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ngate
  REAL, INTENT(IN) :: gatesp
  REAL, INTENT (OUT):: sens(ngate)     ! Sensitivity values for all dist from radar.
!
! Parameters for CASA radar
!
  REAL, PARAMETER :: lmbda = 3.0E-02 ! Wavelength (m)
  REAL, PARAMETER :: pt = 12500      ! Peak power (Watts)
  REAL, PARAMETER :: G  = 38.0       ! Antenna gain (dB)
  REAL, PARAMETER :: F  = 5.5        ! Noise figure (dB)
  REAL, PARAMETER :: tau = 0.67E-06  ! Radar pulse length (s)
  REAL, PARAMETER :: theta = 2.0     ! Antenna half-power beamwidth (deg)
  REAL, PARAMETER :: B = 2.5         ! Bandwidth (MHz)
  REAL, PARAMETER :: lm = 1.0        ! Receiver mis-match loss(dB)
  REAL, PARAMETER :: Kw2 = 0.91      ! Refractive Index of water squarred
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
  bwrad=theta*pi/180.
  rnoise=10.**(0.1*F)
  rlm=10.**(0.1*lm)
  rG=10.**(0.1*G)
  BHz=1.0E06*b
  Ni=K*(Ta+(rnoise-1.0)*T0)*BHz
  sconst=rconst *(Ni/(pi3*Kw2)) * ((8.*ln2)/(bwrad*bwrad)) *         &
         (2./(c*tau)) * ((lmbda*lmbda*four3*rlm)/(Pt*rG*rG))
  DO igate = 1, ngate
     range=igate*gatesp
     sens(igate) = 10.*LOG10(range*range*sconst)
  END DO
END SUBROUTINE xbsens
SUBROUTINE xbatten(ngate,nazim,refl,gatesp,rmisval,refatt)
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
!
!
!  OUTPUT:
!  refatt:    Attenuated reflectivity.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER, INTENT (IN)  :: ngate      ! Number of gates
  INTEGER, INTENT (IN)  :: nazim       ! Number of radials
  REAL, INTENT (IN):: refl(ngate,nazim)     ! Reflectivity (dBZ)
  REAL, INTENT (IN)  :: gatesp         ! Gate spacing (m)
  REAL, INTENT (IN)  :: rmisval
  REAL, INTENT (OUT) :: refatt(ngate,nazim) ! Attenuated reflectivities

  INTEGER :: igate,jazim ! Index.
  REAL :: att_tally      ! Running tally of attenuation values.
  REAL :: att_const      ! Attenuation constant
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
!
! Use Z-R relationship to develop attenuation constant for rain.
!   refl(dBZ)=10 log10(zrconst*R**zrexp)
!   So R=(10**(refl*.1))**(1/zrexp)
!   Kr(dbZ/km)=krconst*(R**krexp)
!   So
!   Kr=(krconst/zrconst)*(10**(refl*0.1)**(krexp/zrexp))
!
!   Calculating two-way attenuation requires this figure be doubled.
!

  att_exp=krexp/zrexp
  att_const=2.*krconst/(zrconst**att_exp)
  hlfgatsp=0.5*gatesp*mperkm

  DO jazim=1, nazim
    IF(refl(1,jazim) .gt. rmisval) THEN
      atten = att_const*((10.0**(0.1*refl(1,jazim)))**att_exp)
      att_tally = atten*hlfgatsp
      refatt(1,jazim) = refl(1,jazim)-att_tally
    ELSE
      atten=0.
      att_tally=0.
      refatt(1,jazim)=refl(1,jazim)
    END IF
    atten_last=atten
    DO igate = 2, ngate
      IF(refl(igate,jazim) .gt. rmisval) THEN
        atten = att_const*((10.0**(0.1*refl(igate,jazim)))**att_exp)
        att_tally = att_tally + hlfgatsp*atten_last + hlfgatsp*atten
        refatt(igate,jazim) = refl(igate,jazim)-att_tally
        atten_last=atten
      ELSE
        refatt(igate,jazim) = refl(igate,jazim)
        atten_last=0.
      END IF
    END DO
  END DO
  RETURN
END SUBROUTINE xbatten
!
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
