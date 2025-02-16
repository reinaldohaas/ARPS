PROGRAM arpsraindiff
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Sample program to read history data file produced by ARPS 4.0
!
!  Link arpsread using the following command on a IBM RISC/6000
!  system, assuming HDF, NetCDF and Savi3D libraries are not
!  available:
!
!  f90 arpsread.f read3d.o nohdfio3d.o nonetio3d.o nosviio3d.o \
!      gradsio3d.o pakio3d.o outlib3d.o ibmlib3d.o
!
!-----------------------------------------------------------------------
!
! Modified 1 June 2002 by Eric Kemp
! Soil variable updates.  Also changed name of program.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'globcst.inc'

!-----------------------------------------------------------------------
!
! Dimension declaration
!
!-----------------------------------------------------------------------

  INTEGER :: nx, ny, nz, nzsoil

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

  REAL, ALLOCATABLE :: zp     (:,:,:) ! The height of the terrain.
  REAL, ALLOCATABLE :: zpsoil (:,:,:) ! The depth of the soil.

  REAL, ALLOCATABLE :: hterain(:,:)   ! Terrain height.

  REAL, ALLOCATABLE :: j1     (:,:,:) ! Coordinate transformation
                                      ! Jacobian -d(zp)/d(x)
  REAL, ALLOCATABLE :: j2     (:,:,:) ! Coordinate transformation
                                      ! Jacobian -d(zp)/d(y)
  REAL, ALLOCATABLE :: j3     (:,:,:) ! Coordinate transformation
                                      ! Jacobian  d(zp)/d(z)
  REAL, ALLOCATABLE :: j3soil (:,:,:) ! Coordinate transformation
                                      ! Jacobian  d(zpsoil)/d(z)

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

  INTEGER :: nstyps
!  PARAMETER ( nstyps = 4 )

  INTEGER, ALLOCATABLE :: soiltyp (:,:,:)    ! Soil type
  REAL, ALLOCATABLE :: stypfrct(:,:,:)       ! Fraction of soil types
  INTEGER, ALLOCATABLE :: vegtyp(:,:)        ! Vegetation type
  REAL, ALLOCATABLE :: lai    (:,:)          ! Leaf Area Index
  REAL, ALLOCATABLE :: roufns (:,:)          ! Surface roughness
  REAL, ALLOCATABLE :: veg    (:,:)          ! Vegetation fraction

  REAL, ALLOCATABLE :: tsoil (:,:,:,:)  ! Soil temperature (K)
  REAL, ALLOCATABLE :: qsoil (:,:,:,:)  ! Soil moisture (m**3/m**3)
  REAL, ALLOCATABLE :: wetcanp(:,:,:)   ! Canopy water amount
  REAL, ALLOCATABLE :: snowdpth(:,:)    ! Snow depth (m)

  REAL, ALLOCATABLE :: raing(:,:)       ! Grid supersaturation rain
  REAL, ALLOCATABLE :: rainc(:,:)       ! Cumulus convective rain
  REAL, ALLOCATABLE :: prcrate(:,:,:)   ! precipitation rate (kg/(m**2*s))
                                        ! prcrate(1,1,1) = total precip. rate
                                        ! prcrate(1,1,2) = grid scale precip. rate
                                        ! prcrate(1,1,3) = cumulus precip. rate
                                        ! prcrate(1,1,4) = microphysics precip. rate

  REAL, ALLOCATABLE :: radfrc(:,:,:)    ! Radiation forcing (K/s)
  REAL, ALLOCATABLE :: radsw (:,:)      ! Solar radiation reaching the surface
  REAL, ALLOCATABLE :: rnflx (:,:)      ! Net radiation flux absorbed by surface
  REAL, ALLOCATABLE :: radswnet (:,:)   ! Net shortwave radiation
  REAL, ALLOCATABLE :: radlwin(:,:)     ! Incoming longwave radiation

  REAL, ALLOCATABLE :: usflx (:,:)      ! Surface flux of u-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: vsflx (:,:)      ! Surface flux of v-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: ptsflx(:,:)      ! Surface heat flux (K*kg/(m*s**2))
  REAL, ALLOCATABLE :: qvsflx(:,:)      ! Surface moisture flux (kg/(m**2*s))

  REAL, ALLOCATABLE :: tem1(:,:,:)      ! Work arrays
  REAL, ALLOCATABLE :: tem2(:,:,:)      ! Work arrays
  REAL, ALLOCATABLE :: tem3(:,:,:)      ! Work arrays

  REAL, ALLOCATABLE :: raing1(:,:)      ! Grid supersaturation rain
  REAL, ALLOCATABLE :: rainc1(:,:)      ! Cumulus convective rain
!
!-----------------------------------------------------------------------
!
!  Misc. internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: hinfmt, nchin
  INTEGER :: lengbf, lenfile1,lenfile2, ireturn
  CHARACTER (LEN=256) :: file1
  CHARACTER (LEN=256) :: grdbasfn
  CHARACTER (LEN=256) :: file2
  CHARACTER (LEN=256) :: basdmpfn
  INTEGER :: lbasdmpf, i,j, istatus

  REAL :: time
!-----------------------------------------------------------------------
!
! NAMELIST declaration
!
!-----------------------------------------------------------------------

!  NAMELIST /grid_dims/ nx, ny, nz, nzsoil
  NAMELIST /input_files/ hinfmt, grdbasfn, file1, file2
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!   READ(5,grid_dims,END=100)
!   WRITE(6,'(a)')'Namelist block grid_dims sucessfully read in.'
!
!-----------------------------------------------------------------------
!
!  Get the name of the input data set.
!
!-----------------------------------------------------------------------
!
  READ(5,input_files,END=100)
  WRITE(6,'(a)')'Namelist block input_files sucessfully read in.'

  lengbf = 256
  CALL strlnth( grdbasfn, lengbf)

  lenfile1 = 256
  CALL strlnth( file1, lenfile1 )

  lenfile2 = 256
  CALL strlnth( file2, lenfile2 )

  CALL get_dims_from_data(hinfmt,file1(1:lenfile1),                     &
                          nx,ny,nz,nzsoil,nstyps, ireturn)

  IF (nstyps <= 0) nstyps = 1
  nstyp = nstyps ! Copy to global variable


!-----------------------------------------------------------------------
!
!  Allocate the variables and initialize the them to zero
!
!-----------------------------------------------------------------------

  ALLOCATE(x(nx), STAT=istatus)
  x = 0
  ALLOCATE(y(ny), STAT=istatus)
  y = 0
  ALLOCATE(z(nz), STAT=istatus)
  z = 0
  ALLOCATE(zp(nx,ny,nz), STAT=istatus)
  zp = 0
  ALLOCATE(zpsoil(nx,ny,nzsoil), STAT=istatus)
  zpsoil = 0
  ALLOCATE(hterain(nx,ny), STAT=istatus)
  hterain = 0
  ALLOCATE(j1(nx,ny,nz), STAT=istatus)
  j1 = 0
  ALLOCATE(j2(nx,ny,nz), STAT=istatus)
  j2 = 0
  ALLOCATE(j3(nx,ny,nz), STAT=istatus)
  j3 = 0
  ALLOCATE(j3soil(nx,ny,nzsoil), STAT=istatus)
  j3soil = 0
  ALLOCATE(uprt(nx,ny,nz), STAT=istatus)
  uprt = 0
  ALLOCATE(vprt(nx,ny,nz), STAT=istatus)
  vprt = 0
  ALLOCATE(wprt(nx,ny,nz), STAT=istatus)
  wprt = 0
  ALLOCATE(ptprt(nx,ny,nz), STAT=istatus)
  ptprt = 0
  ALLOCATE(pprt(nx,ny,nz), STAT=istatus)
  pprt = 0
  ALLOCATE(qvprt(nx,ny,nz), STAT=istatus)
  qvprt = 0
  ALLOCATE(qscalar(nx,ny,nz,nscalar), STAT=istatus)
  qscalar = 0.0
  ALLOCATE(tke(nx,ny,nz), STAT=istatus)
  tke = 0
  ALLOCATE(kmh(nx,ny,nz), STAT=istatus)
  kmh = 0
  ALLOCATE(kmv(nx,ny,nz), STAT=istatus)
  kmv = 0
  ALLOCATE(ubar(nx,ny,nz), STAT=istatus)
  ubar = 0
  ALLOCATE(vbar(nx,ny,nz), STAT=istatus)
  vbar = 0
  ALLOCATE(wbar(nx,ny,nz), STAT=istatus)
  wbar = 0
  ALLOCATE(ptbar(nx,ny,nz), STAT=istatus)
  ptbar = 0
  ALLOCATE(pbar(nx,ny,nz), STAT=istatus)
  pbar = 0
  ALLOCATE(rhobar(nx,ny,nz), STAT=istatus)
  rhobar = 0
  ALLOCATE(qvbar(nx,ny,nz), STAT=istatus)
  qvbar = 0
  ALLOCATE(u(nx,ny,nz), STAT=istatus)
  u = 0
  ALLOCATE(v(nx,ny,nz), STAT=istatus)
  v = 0
  ALLOCATE(w(nx,ny,nz), STAT=istatus)
  w = 0
  ALLOCATE(qv(nx,ny,nz), STAT=istatus)
  qv = 0
  ALLOCATE(soiltyp(nx,ny,nstyps), STAT=istatus)
  soiltyp = 0
  ALLOCATE(stypfrct(nx,ny,nstyps), STAT=istatus)
  stypfrct = 0
  ALLOCATE(vegtyp(nx,ny), STAT=istatus)
  vegtyp = 0
  ALLOCATE(lai(nx,ny), STAT=istatus)
  lai = 0
  ALLOCATE(roufns(nx,ny), STAT=istatus)
  roufns = 0
  ALLOCATE(veg(nx,ny), STAT=istatus)
  veg = 0
  ALLOCATE(tsoil(nx,ny,nzsoil,0:nstyps), STAT=istatus)
  tsoil = 0
  ALLOCATE(qsoil(nx,ny,nzsoil,0:nstyps), STAT=istatus)
  qsoil = 0
  ALLOCATE(wetcanp(nx,ny,0:nstyps), STAT=istatus)
  wetcanp = 0
  ALLOCATE(snowdpth(nx,ny), STAT=istatus)
  snowdpth = 0
  ALLOCATE(raing(nx,ny), STAT=istatus)
  raing = 0
  ALLOCATE(rainc(nx,ny), STAT=istatus)
  rainc = 0
  ALLOCATE(prcrate(nx,ny,4), STAT=istatus)
  prcrate = 0
  ALLOCATE(radfrc(nx,ny,nz), STAT=istatus)
  radfrc = 0
  ALLOCATE(radsw(nx,ny), STAT=istatus)
  radsw = 0
  ALLOCATE(rnflx(nx,ny), STAT=istatus)
  rnflx = 0
  ALLOCATE(radswnet(nx,ny), STAT=istatus)
  radswnet = 0
  ALLOCATE(radlwin(nx,ny), STAT=istatus)
  radlwin = 0
  ALLOCATE(usflx(nx,ny), STAT=istatus)
  usflx = 0
  ALLOCATE(vsflx(nx,ny), STAT=istatus)
  vsflx = 0
  ALLOCATE(ptsflx(nx,ny), STAT=istatus)
  ptsflx = 0
  ALLOCATE(qvsflx(nx,ny), STAT=istatus)
  qvsflx = 0
  ALLOCATE(tem1(nx,ny,nz), STAT=istatus)
  tem1 = 0
  ALLOCATE(tem2(nx,ny,nz), STAT=istatus)
  tem2 = 0
  ALLOCATE(tem3(nx,ny,nz), STAT=istatus)
  tem3 = 0
  ALLOCATE(raing1(nx,ny), STAT=istatus)
  raing1 = 0
  ALLOCATE(rainc1(nx,ny), STAT=istatus)
  rainc1 = 0

!-----------------------------------------------------------------------
!
! Read in files.
!
!-----------------------------------------------------------------------

  nchin = 9

  CALL dtaread(nx,ny,nz,nzsoil,nstyps,                                  &
       hinfmt, nchin,grdbasfn(1:lengbf),lengbf,                         &
       file1(1:lenfile1),lenfile1,time,                                 &
       x,y,z,zp,zpsoil, uprt ,vprt ,wprt ,ptprt, pprt ,                 &
       qvprt, qscalar, tke,kmh,kmv,                                     &
       ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,                    &
       soiltyp,stypfrct,vegtyp,lai,roufns,veg,                          &
       tsoil, qsoil, wetcanp,snowdpth,                                  &
       raing,rainc,prcrate,                                             &
       radfrc,radsw,rnflx,radswnet,radlwin,                             &
       usflx,vsflx,ptsflx,qvsflx,                                       &
       ireturn, tem1,tem2,tem3)

  CALL dtaread(nx,ny,nz,nzsoil,nstyps,                                  &
       hinfmt, nchin,grdbasfn(1:lengbf),lengbf,                         &
       file2(1:lenfile2),lenfile2,time,                                 &
       x,y,z,zp,zpsoil,uprt,vprt ,wprt ,ptprt, pprt,                    &
       qvprt, qscalar, tke,kmh,kmv,                                     &
       ubar, vbar, wbar, ptbar, pbar, rhobar, qvbar,                    &
       soiltyp,stypfrct,vegtyp,lai,roufns,veg,                          &
       tsoil, qsoil, wetcanp,snowdpth,                                  &
       raing1,rainc1,prcrate,                                           &
       radfrc,radsw,rnflx,radswnet,radlwin,                             &
       usflx,vsflx,ptsflx,qvsflx,                                       &
       ireturn, tem1,tem2,tem3)

  curtim = time


  DO i=1,nx
    DO j=1,ny
      rainc1(i,j) = rainc1(i,j) - rainc(i,j)
      raing1(i,j) = raing1(i,j) - raing(i,j)
    END DO
  END DO

  dirname = './'

  CALL wrtvar(nx,ny,1, rainc1,'rainc',time,runname,dirname)
  CALL wrtvar(nx,ny,1, raing1,'raing',time,runname,dirname)

  DO i=1,nx
    DO j=1,ny
      raing1(i,j) = raing1(i,j) + rainc1(i,j)
    END DO
  END DO

  CALL wrtvar(nx,ny,1, raing1,'raint',time,runname,dirname)

  GOTO 101

  100 CONTINUE

  WRITE(6,'(a)')'Namelist block READ in error. Then program will terminated.'

  101 CONTINUE

  STOP

END PROGRAM arpsraindiff
