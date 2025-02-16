PROGRAM testret
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   PROGRAM TESTRET                    ######
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
!  Test writing of retrieval columns.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  DATA ARRAYS READ IN:
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (km)
!    zp       z coordinate of grid points in computational space (m)
!    zpsoil   z coordinate of soil model 
!
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
!    pt       Potential temperature (K)
!    qv       Water vapor mixing ratio (kg/kg)
!
!
!  WORK ARRAYS:
!
!    tem1     Temporary work array.
!    tem2     Temporary work array.
!    tem3     Temporary work array.
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
  INCLUDE 'indtflg.inc'
  INCLUDE 'grid.inc'

  INTEGER :: nx, ny, nz            ! Dimensions declaration
  INTEGER :: nzsoil 
!
!-----------------------------------------------------------------------
!
!  Arrays to be read in:
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: x     (:)   ! The x-coord. of the physical and
                                   ! computational grid. Defined at u-point.
  REAL, ALLOCATABLE :: y     (:)   ! The y-coord. of the physical and
                                   ! computational grid. Defined at v-point.
  REAL, ALLOCATABLE :: z     (:)   ! The z-coord. of the computational grid.
                                   ! Defined at w-point on the staggered grid.
  REAL, ALLOCATABLE :: zp  (:,:,:) ! The physical height coordinate defined at
                                   ! w-point of the staggered grid.
  REAL, ALLOCATABLE :: zpsoil(:,:,:) ! The physical height coordinate of 
                                   !  the soil model 
  REAL, ALLOCATABLE :: j1  (:,:,:) ! Coordinate transformation Jacobian defined
                                   ! as - d( zp )/d( x )
  REAL, ALLOCATABLE :: j2  (:,:,:) ! Coordinate transformation Jacobian defined
                                   ! as - d( zp )/d( y )
  REAL, ALLOCATABLE :: j3  (:,:,:) ! Coordinate transformation Jacobian defined
                                   ! as d( zp )/d( z )

  REAL, ALLOCATABLE :: u(:,:,:)    ! Total u-velocity (m/s)
  REAL, ALLOCATABLE :: v(:,:,:)    ! Total v-velocity (m/s)
  REAL, ALLOCATABLE :: w(:,:,:)    ! Total w-velocity (m/s)
  REAL, ALLOCATABLE :: pt(:,:,:)   ! Total potential temperature (K)

  REAL, ALLOCATABLE :: qv    (:,:,:)  ! Water vapor specific humidity (kg/kg)
  REAL, ALLOCATABLE :: qc    (:,:,:)  ! Cloud water mixing ratio (kg/kg)
  REAL, ALLOCATABLE :: qr    (:,:,:)  ! Rain water mixing ratio (kg/kg)
  REAL, ALLOCATABLE :: qi    (:,:,:)  ! Cloud ice mixing ratio (kg/kg)
  REAL, ALLOCATABLE :: qs    (:,:,:)  ! Snow mixing ratio (kg/kg)
  REAL, ALLOCATABLE :: qh    (:,:,:)   ! Hail mixing ratio (kg/kg)
  REAL, ALLOCATABLE :: tke   (:,:,:)  ! Turbulent Kinetic Energy ((m/s)**2)
  REAL, ALLOCATABLE :: kmh   (:,:,:)  ! Horizontal turb. mixing coef. for
                                      ! momentum. ( m**2/s )
  REAL, ALLOCATABLE :: kmv   (:,:,:)  ! Vertical turb. mixing coef. for
                                      ! momentum. ( m**2/s )
  REAL, ALLOCATABLE :: ubar  (:,:,:)  ! Base state u-velocity (m/s)
  REAL, ALLOCATABLE :: vbar  (:,:,:)  ! Base state v-velocity (m/s)
  REAL, ALLOCATABLE :: wbar  (:,:,:)  ! Base state w-velocity (m/s)
  REAL, ALLOCATABLE :: ptbar (:,:,:)  ! Base state potential temperature (K)
  REAL, ALLOCATABLE :: rhobar(:,:,:)  ! Base state air density (kg/m**3)
  REAL, ALLOCATABLE :: pbar  (:,:,:)  ! Base state pressure (Pascal)
  REAL, ALLOCATABLE :: qvbar (:,:,:)  ! Base state water vapor specific humidity
                                      ! (kg/kg)
  INTEGER :: nstyps
!  PARAMETER (nstyps=4)

  INTEGER, ALLOCATABLE:: soiltyp (:,:,:)  ! Soil type
  REAL, ALLOCATABLE :: stypfrct(:,:,:)    ! Fraction of soil types
  INTEGER, ALLOCATABLE :: vegtyp  (:,:)   ! Vegetation type
  REAL, ALLOCATABLE :: lai     (:,:)      ! Leaf Area Index
  REAL, ALLOCATABLE :: roufns  (:,:)      ! Surface roughness
  REAL, ALLOCATABLE :: veg     (:,:)      ! Vegetation fraction

  REAL, ALLOCATABLE :: tsoil  (:,:,:,:)     ! Soil Temperature (K)
  REAL, ALLOCATABLE :: qsoil  (:,:,:,:)     ! soil moisture
  REAL, ALLOCATABLE :: wetcanp(:,:,:)     ! Canopy water amount
  REAL, ALLOCATABLE :: snowdpth(:,:)      ! Snow depth (m)

  REAL, ALLOCATABLE :: raing  (:,:)       ! Grid supersaturation rain
  REAL, ALLOCATABLE :: rainc  (:,:)       ! Cumulus convective rain
  REAL, ALLOCATABLE :: prcrate(:,:,:)     ! precipitation rate (kg/(m**2*s))
                                          ! prcrate(1,1,1) = total precip. rate
                                          ! prcrate(1,1,2) = grid scale precip. rate
                                          ! prcrate(1,1,3) = cumulus precip. rate
                                          ! prcrate(1,1,4) = microphysics precip. rate

  REAL, ALLOCATABLE :: radfrc(:,:,:)      ! Radiation forcing (K/s)
  REAL , ALLOCATABLE:: radsw (:,:)        ! Solar radiation reaching the surface
  REAL, ALLOCATABLE :: rnflx (:,:)        ! Net radiation flux absorbed by surface
  REAL, ALLOCATABLE :: radswnet (:,:) ! Net solar radiation, SWin - SWout
  REAL, ALLOCATABLE :: radlwin  (:,:) ! Incoming longwave radiation

  REAL, ALLOCATABLE :: usflx (:,:)        ! Surface flux of u-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: vsflx (:,:)        ! Surface flux of v-momentum (kg/(m*s**2))
  REAL, ALLOCATABLE :: ptsflx(:,:)        ! Surface heat flux (K*kg/(m*s**2))
  REAL, ALLOCATABLE :: qvsflx(:,:)        ! Surface moisture flux (kg/(m**2*s))

  REAL, ALLOCATABLE :: uprt  (:,:,:)      ! Perturbation u-velocity (m/s)
  REAL, ALLOCATABLE :: vprt  (:,:,:)      ! Perturbation v-velocity (m/s)
  REAL, ALLOCATABLE :: wprt  (:,:,:)      ! Perturbation w-velocity (m/s)
  REAL, ALLOCATABLE :: ptprt (:,:,:)      ! Perturbation potential temperature (K)
  REAL, ALLOCATABLE :: pprt  (:,:,:)      ! Perturbation pressure (Pascal)
  REAL, ALLOCATABLE :: rhoprt(:,:,:)      ! Perturbation air density (kg/m**3)
  REAL, ALLOCATABLE :: qvprt (:,:,:)      ! Perturbation water vapor specific
                                          ! humidity (kg/kg)
!
!-----------------------------------------------------------------------
!
!  Other data variables
!
!-----------------------------------------------------------------------
!
  REAL :: time
!
!-----------------------------------------------------------------------
!
!  Temporary working arrays:
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: tem1(:,:,:)
  REAL, ALLOCATABLE :: tem2(:,:,:)
  REAL, ALLOCATABLE :: tem3(:,:,:)

  REAL, ALLOCATABLE :: tem1d1(:),tem1d2(:),tem1d3(:),tem1d4(:),    &
       tem1d5(:), tem1d6(:),tem1d7(:),tem1d8(:),tem1d9(:)
!
!-----------------------------------------------------------------------
!
!  "Fake" radar id stuff
!
!-----------------------------------------------------------------------
!
  INTEGER :: iretfmt
  CHARACTER (LEN=256) :: retfname
  PARAMETER(iretfmt=1)
  CHARACTER (LEN=4) :: radid
  REAL :: latrad,lonrad,elvrad
  PARAMETER(radid='KTLX',                                               &
            latrad=35.3331,                                             &
            lonrad=-97.2778,                                            &
            elvrad=389.4)
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: xsc(:)
  REAL, ALLOCATABLE :: ysc(:)
  REAL, ALLOCATABLE :: zpsc(:,:,:)
  REAL :: latnot(2)
!
  INTEGER :: i,j,k,ireturn
  INTEGER :: istride,jstride,kstride
  PARAMETER (istride=2,                                                 &
             jstride=2,                                                 &
             kstride=1)
  CHARACTER (LEN=256) :: filename
  CHARACTER (LEN=256) :: grdbasfn
  INTEGER :: ngchan,nchanl,lenfil,lengbf
  INTEGER :: nchin,iyr
  INTEGER :: istatus
!----------------------------------------------------------------------
!
!  NAMELIST declaration
!
!----------------------------------------------------------------------
   
!   NAMELIST /grid_dims/ nx, ny, nz

   NAMELIST /input_fn/ hdmpfmt, grdbasfn, filename 
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
     '###                  Welcome to ARPSPRT                     ###', &
     '###      This program reads in the history dump data        ###', &
     '###      sets generated by ARPS, and prints the array       ###', &
     '###      contents as 2-D arrays tables at user specified    ###', &
     '###      slices.                                            ###', &
     '###                                                         ###', &
     '###############################################################', &
     '###############################################################'

!  READ(5,grid_dims,END=100)
!  WRITE(6,'(/a)') 'Namelist block grid_dims successfully read.'

!-----------------------------------------------------------------------
!
!  Get the name of the input data set.
!
!-----------------------------------------------------------------------
!
  READ(5,input_fn,END=100)
  WRITE(6,'(/a,a)')' Namelist block input_fn successfully read in.'

  lengbf=LEN_trim(grdbasfn)
  lenfil=LEN_trim(filename)
  WRITE(6,'(/a,a)')' The data set name is ', filename(1:lenfil)

  CALL get_dims_from_data(hdmpfmt,grdbasfn(1:lenfil),nx,ny,nz,nzsoil,   &
                          nstyps,ireturn)
     
  nstyp = nstyps ! Copy to global variable

!-----------------------------------------------------------------------
!
!  Allocate all arrays and fill them with zero, which is the default 
!  value of the data arrays.
!
!-----------------------------------------------------------------------
!

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
  ALLOCATE(xsc(nx),STAT=istatus)
  xsc=0
  ALLOCATE(ysc(ny),STAT=istatus)
  ysc=0
  ALLOCATE(zpsc(nx,ny,nz),STAT=istatus)
  zpsc=0
  ALLOCATE(j1(nx,ny,nz),STAT=istatus)
  j1=0
  ALLOCATE(j2(nx,ny,nz),STAT=istatus)
  j2=0
  ALLOCATE(j3(nx,ny,nz),STAT=istatus)
  j3=0
  ALLOCATE(u(nx,ny,nz),STAT=istatus)
  u=0
  ALLOCATE(v(nx,ny,nz),STAT=istatus)
  v=0
  ALLOCATE(w(nx,ny,nz),STAT=istatus)
  w=0
  ALLOCATE(pt(nx,ny,nz),STAT=istatus)
  pt=0
  ALLOCATE(qv(nx,ny,nz),STAT=istatus)
  qv=0
  ALLOCATE(qc(nx,ny,nz),STAT=istatus)
  qc=0
  ALLOCATE(qr(nx,ny,nz),STAT=istatus)
  qr=0
  ALLOCATE(qi(nx,ny,nz),STAT=istatus)
  qi=0
  ALLOCATE(qs(nx,ny,nz),STAT=istatus)
  qs=0
  ALLOCATE(qh(nx,ny,nz),STAT=istatus)
  qh=0
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
  ALLOCATE(rhobar(nx,ny,nz),STAT=istatus)
  rhobar=0
  ALLOCATE(pbar(nx,ny,nz),STAT=istatus)
  pbar=0
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
  ALLOCATE(rhoprt(nx,ny,nz),STAT=istatus)
  rhoprt=0
  ALLOCATE(qvprt(nx,ny,nz),STAT=istatus)
  qvprt=0  

  ALLOCATE(tem1(nx,ny,nz),STAT=istatus)
  tem1=0  
  ALLOCATE(tem2(nx,ny,nz),STAT=istatus)
  tem2=0  
  ALLOCATE(tem3(nx,ny,nz),STAT=istatus)
  tem3=0  

  ALLOCATE(tem1d1(nz),STAT=istatus)
  tem1d1=0  
  ALLOCATE(tem1d2(nz),STAT=istatus)
  tem1d2=0  
  ALLOCATE(tem1d3(nz),STAT=istatus)
  tem1d3=0  
  ALLOCATE(tem1d4(nz),STAT=istatus)
  tem1d4=0  
  ALLOCATE(tem1d5(nz),STAT=istatus)
  tem1d5=0  
  ALLOCATE(tem1d6(nz),STAT=istatus)
  tem1d6=0  
  ALLOCATE(tem1d7(nz),STAT=istatus)
  tem1d7=0  
  ALLOCATE(tem1d8(nz),STAT=istatus)
  tem1d8=0  
  ALLOCATE(tem1d9(nz),STAT=istatus)
  tem1d9=0  
  
!
!-----------------------------------------------------------------------
!
!  Read all input data arrays
!
!-----------------------------------------------------------------------
!
  CALL dtaread(nx,ny,nz,nzsoil,nstyps,                                  &
               hdmpfmt,nchin,grdbasfn(1:lengbf),lengbf,                 &
               filename(1:lenfil),lenfil,time,                          &
               x,y,z,zp,zpsoil,uprt,vprt,wprt,ptprt,pprt,               &
               qvprt, qc, qr, qi, qs, qh, tke,kmh,kmv,                  &
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
  IF( ireturn == 0 ) THEN    ! successful read

    DO i=1,nx-1
      xsc(i)=0.5*(x(i)+x(i+1))
    END DO
    xsc(nx)=2.*xsc(nx-1)-xsc(nx-2)
    DO j=1,ny-1
      ysc(j)=0.5*(y(j)+y(j+1))
    END DO
    ysc(ny)=2.*ysc(ny-1)-ysc(ny-2)

    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          zpsc(i,j,k)=0.5*(zp(i,j,k)+zp(i,j,k+1))
          tem1(i,j,k)=0.5*(uprt(  i,j,k)+ubar(  i,j,k)                  &
                          +uprt(i+1,j,k)+ubar(i+1,j,k))
          tem2(i,j,k)=0.5*(vprt(i,  j,k)+vbar(i,  j,k)                  &
                          +vprt(1,j+1,k)+vbar(i,j+1,k))
          IF(MOD(i,istride) == 0 .AND. MOD(j,jstride) == 0 .AND.        &
                MOD(k,kstride) == 0 ) THEN
            tem3(i,j,k)=1.0
          ELSE
            tem3(i,j,k)=-999.
          END IF
        END DO
      END DO
    END DO

    latnot(1) = trulat1
    latnot(2) = trulat2
    CALL setmapr(mapproj,1.0,latnot,trulon)

    iyr=MOD(year,100)
    WRITE(retfname,'(a,a,i2.2,i2.2,i2.2,a,i2.2,i2.2)')                  &
        radid,'.',iyr,month,day,'.',hour,minute
    PRINT *, ' Writing data into ',retfname

    CALL wtretcol(nx,ny,nz,                                             &
                  2,nx-2,2,ny-2,2,nz-2,                                 &
                  iyr,month,day,hour,minute,second,                     &
                  iretfmt,retfname,radid,latrad,lonrad,elvrad,          &
                  xsc,ysc,zpsc,                                         &
                  tem1,tem2,ptprt,pprt,qvprt,qr,                        &
                  ptbar,pbar,qvbar,tem3,                                &
                  tem1d1,tem1d2,tem1d3,tem1d4,                          &
                  tem1d5,tem1d6,tem1d7,tem1d8,tem1d9)

  END IF                                       ! successful read

  GOTO 101
  
  100 CONTINUE
  
  WRITE(6,'(/a,a)') 'Error reading NAMELIST file. The program will abort.'
  
  101 CONTINUE

  STOP
END PROGRAM testret
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WTRETCOL                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE wtretcol(nx,ny,nz,                                           &
           ibeg,iend,jbeg,jend,kbeg,kend,                               &
           iyr,imon,iday,ihr,imin,isec,                                 &
           iretfmt,retfname,radid,latrad,lonrad,elvrad,                 &
           xsc,ysc,zpsc,                                                &
           us,vs,ptprt,pprt,qvprt,qr,                                   &
           ptbar,pbar,qvbar,retrflg,                                    &
           outk,outhgt,outu,outv,                                       &
           outpr,outpt,outqv,outqr,outret)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Test writing of retrieval columns.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  Writes gridded radar data to a file
!
!  INPUT
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (km)
!    zp       z coordinate of grid points in computational space (m)
!
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
!    qr       Rainwater mixing ratio (kg/kg)
!
!-----------------------------------------------------------------------
!
!
  IMPLICIT NONE
!
  INTEGER :: nx,ny,nz
  INTEGER :: ibeg,iend,jbeg,jend,kbeg,kend
!
  REAL :: xsc(nx)
  REAL :: ysc(ny)
  REAL :: zpsc(nx,ny,nz)
  REAL :: us(nx,ny,nz)         ! total u velocity component at scalar points
  REAL :: vs(nx,ny,nz)         ! total v velocity component at scalar points
  REAL :: ptprt (nx,ny,nz)     ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)     ! Perturbation pressure (Pascal)
  REAL :: qvprt (nx,ny,nz)     ! Perturbation water vapor specific
                               ! humidity (kg/kg)
  REAL :: qr    (nx,ny,nz)     ! Rain water mixing ratio (kg/kg)
  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal)
  REAL :: qvbar (nx,ny,nz)     ! Base state water vapor specific humidity
                               ! (kg/kg)
  REAL :: retrflg(nx,ny,nz)
!
  INTEGER :: iretfmt
  CHARACTER (LEN=*) :: retfname
  CHARACTER (LEN=4) :: radid
  REAL :: latrad
  REAL :: lonrad
  REAL :: elvrad
  INTEGER :: iyr,imon,iday,ihr,imin,isec
!
!-----------------------------------------------------------------------
!
!  Retrieval output variables
!
!-----------------------------------------------------------------------
!
  REAL :: outk(nz)
  REAL :: outhgt(nz)
  REAL :: outu(nz)
  REAL :: outv(nz)
  REAL :: outpr(nz)
  REAL :: outpt(nz)
  REAL :: outqv(nz)
  REAL :: outqr(nz)
  REAL :: outret(nz)
!
!-----------------------------------------------------------------------
!
!  Retrieved data threshold
!
!-----------------------------------------------------------------------
!
  REAL :: retrthr
  PARAMETER(retrthr=0.)
!
!-----------------------------------------------------------------------
!
!  Include file
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iunit,myr,itime
  INTEGER :: i,j,k,klev,kk,kntcol
  INTEGER :: nradvr,iradvr,ireftim,idummy
  REAL :: gridlat,gridlon,elev,rdummy
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  myr=1900+iyr
  IF(myr < 1960) myr=myr+100
  CALL ctim2abss(myr,imon,iday,ihr,imin,isec,itime)
!
  CALL getunit(iunit)
!
!  Open file for output
!
  OPEN(iunit,FILE=TRIM(retfname),STATUS='unknown',                      &
       FORM='unformatted')
!
!  Write retrieval description variables
!
  WRITE(iunit) radid
  WRITE(iunit) ireftim,itime,idummy,idummy,idummy,                      &
               idummy,idummy,idummy,idummy,idummy
!
!  Write grid description variables
!  This should provide enough info to verify that the
!  proper grid has been chosen.  To recreate the grid,
!  icluding elevation information,
!  the reading program should get a grid-base-file
!  named runname.grdbasfil
!
  idummy=0
  rdummy=0.
  WRITE(iunit) runname
  WRITE(iunit) hdmpfmt,strhopt,mapproj,idummy,idummy,                   &
               idummy,idummy,idummy,idummy,idummy
  WRITE(iunit) dx,dy,dz,dzmin,ctrlat,                                   &
               ctrlon,trulat1,trulat2,trulon,sclfct,                    &
               latrad,lonrad,elvrad,rdummy,rdummy
  WRITE(iunit) nradvr,iradvr
!
!  For each horizontal grid point form a column of remapped
!  data containing the non-missing grid points
!
  kntcol=0
  DO j=jbeg,jend
    DO i=ibeg,iend
      klev=0
      DO k=kbeg,kend
        IF(retrflg(i,j,k) > retrthr) THEN
          klev=klev+1
          outk(klev)=FLOAT(k)
          outhgt(klev)=zpsc(i,j,k)
          outu(klev)=us(i,j,k)
          outv(klev)=vs(i,j,k)
          outpr(klev)=pprt(i,j,k)+pbar(i,j,k)
          outpt(klev)=ptprt(i,j,k)+ptbar(i,j,k)
          outqv(klev)=qvprt(i,j,k)+qvbar(i,j,k)
          outqr(klev)=qr(i,j,k)
          outret(klev)=retrflg(i,j,k)
        END IF
      END DO
!
!  If there are data in this column, write them to the file.
!
      IF(klev > 0) THEN
        kntcol=kntcol+1
        CALL xytoll(1,1,xsc(i),ysc(j),gridlat,gridlon)
        elev=0.5*(zpsc(i,j,1)+zpsc(i,j,2))
        WRITE(iunit) i,j,xsc(i),ysc(j),                                 &
                     gridlat,gridlon,elev,klev
        WRITE(iunit) (outk(kk),kk=1,klev)
        WRITE(iunit) (outhgt(kk),kk=1,klev)
        WRITE(iunit) (outu(kk),kk=1,klev)
        WRITE(iunit) (outv(kk),kk=1,klev)
        WRITE(iunit) (outpr(kk),kk=1,klev)
        WRITE(iunit) (outpt(kk),kk=1,klev)
        WRITE(iunit) (outqv(kk),kk=1,klev)
        WRITE(iunit) (outqr(kk),kk=1,klev)
        WRITE(iunit) (outret(kk),kk=1,klev)
      END IF
    END DO
  END DO
!
  CLOSE(iunit)
  CALL retunit(iunit)
!
!  Report on what data were written
!
  WRITE(6,'(//a,i2.2,i2.2,i2.2,a1,i2.2,a1,i2.2)')                       &
                    ' Output statistics for time ',                     &
                      iyr,imon,iday,' ',ihr,':',imin
  WRITE(6,'(a,i6,a,/a,i6,a//)')                                         &
           ' There were ',kntcol,' columns written ',                   &
           ' of a total ',(nx*ny),' possible.'
!
  RETURN
END SUBROUTINE wtretcol
