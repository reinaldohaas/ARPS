!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETARPS                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getarps(nx_ext,ny_ext,nz_ext,nzsoil_ext,                     &
           dir_extd,extdname,extdopt,extdfmt,                           &
           extdinit,extdfcst,julfname,nstyps,                           &
           iproj_ext,scale_ext,                                         &
           trlon_ext,latnot_ext,x0_ext,y0_ext,                          &
           lat_ext,lon_ext,latu_ext,lonu_ext,latv_ext,lonv_ext,         &
           p_ext,hgt_ext,zp_ext,zpsoil_ext,t_ext,qv_ext,                &
           u_ext,vatu_ext,v_ext,uatv_ext,w_ext,                         &
           qscalar_ext,                                                 &
           soiltyp_ext,stypfrct_ext,vegtyp_ext,                         &
           lai_ext,roufns_ext,veg_ext,                                  &
           tsoil_ext,qsoil_ext,wetcanp_ext,                             &
           snowdpth_ext,ubar_ext,vbar_ext,wbar_ext,                     &
           ptbar_ext,pbar_ext,qvbar_ext,                                &
           tem1_ext,tem2_ext,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  ARPS version.
!
!  Reads an ARPS file for processing by ext2arps, a program
!  that converts external files to ARPS variables and format.
!  This version is useful when you want to use an ARPS file
!  with a different orientation or terrain than your final
!  ARPS product so arpsr2h does not work.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  November, 1995
!
!  MODIFICATION HISTORY:
!    01/16/1996 (Yuhe Liu)
!    Added three arguments to specify the "external" ARPS data file
!    names and format.
!
!  2000/08/14 (Gene Bassett)
!  Added multiple soil types, sfcdata variables and grid staggering
!  for use with arps history format of external model data.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    dir_extd      Directory name for external file
!    extdname      Prefix string of external file name
!    extdopt       Option of external data sources
!    extdfmt       Option of external data format
!
!    extdinit      Initialized time in mm-dd-yyyy:hh:mm:ss format
!    extdfcst      Forecast hour in HHH:MM:SS format
!    julfname      File name in yyjjjhhmm format
!
!  OUTPUT:
!
!    iproj_ext     Map projection number of external data
!    scale_ext     Scale factor of external data
!    trlon_ext     True longitude of external data (degrees E)
!    latnot_ext(2) True latitude(s) of external data (degrees N)
!    x0_ext        x coordinate of origin of external data
!    y0_ext        y coordinate of origin of external data
!    lat_ext       latitude of external data points (degrees N)
!    lon_ext       longitude of external data points (degrees E)
!    latu_ext      latitude of external u points (degrees N)
!    lonu_ext      longitude of external u points (degrees E)
!    latv_ext      latitude of external v points (degrees N)
!    lonv_ext      longitude of external v points (degrees E)
!    p_ext         pressure (Pascal)
!    hgt_ext       height (m)
!    zp_ext        height (m) (on arps grid)
!    zpsoil_ext    height (m) (soil model on arps grid)
!    t_ext         temperature (K)
!    qv_ext        specific humidity (kg/kg)
!    u_ext         u wind component (m/s) (staggered grid, true n-s component)
!    vatu_ext      v wind component at u location (m/s)
!    v_ext         v wind component (m/s) (staggered grid, true e-w component)
!    uatv_ext      u wind component at v location (m/s)
!    qc_ext        Cloud water mixing ratio (kg/kg)
!    qr_ext        Rain water mixing ratio (kg/kg)
!    qi_ext        Ice mixing ratio (kg/kg)
!    qs_ext        Snow mixing ratio (kg/kg)
!    qh_ext        Hail mixing ratio (kg/kg)
!
!    soiltyp_ext   Soil type
!    stypfrct_ext  Soil type fraction
!    vegtyp_ext    Vegetation type
!    lai_ext       Leaf Area Index
!    roufns_ext    Surface roughness
!    veg_ext       Vegetation fraction
!
!    tsoil_ext     Soil temperature (K)
!    qsoil_ext     Soil moisture
!    wetcanp_ext   Water content on canopy
!
!    snowdpth_ext  Snow depth (m)
!    ubar_ext      Base state u-velocity (m/s)
!    vbar_ext      Base state v-velocity (m/s)
!    wbar_ext      Base state w-velocity (m/s)
!
!    ptbar_ext     Base state potential temperature (K)
!    pbar_ext      Base state pressure (Pascal)
!    qvbar_ext     Base state water vapor specific humidity (kg/kg)
!
!    istatus       status indicator
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!

  CHARACTER (LEN=*) :: dir_extd
  CHARACTER (LEN=*) :: extdname

  INTEGER :: extdopt
  INTEGER :: extdfmt

  INTEGER :: nstyps

  CHARACTER (LEN=19) :: extdinit
  CHARACTER (LEN=9) :: extdfcst
  CHARACTER (LEN=9) :: julfname
!
!-----------------------------------------------------------------------
!
!  External grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj_ext
  REAL :: scale_ext,trlon_ext
  REAL :: latnot_ext(2)
  REAL :: x0_ext,y0_ext
!
!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx_ext,ny_ext,nz_ext,nzsoil_ext

  REAL :: lat_ext(nx_ext,ny_ext)
  REAL :: lon_ext(nx_ext,ny_ext)
  REAL :: latu_ext(nx_ext,ny_ext)
  REAL :: lonu_ext(nx_ext,ny_ext)
  REAL :: latv_ext(nx_ext,ny_ext)
  REAL :: lonv_ext(nx_ext,ny_ext)
  REAL :: p_ext(nx_ext,ny_ext,nz_ext)     ! Pressure (Pascals)
  REAL :: hgt_ext(nx_ext,ny_ext,nz_ext)   ! Height (m)
  REAL :: zp_ext(nx_ext,ny_ext,nz_ext)    ! Height (m) (on arps grid)
  REAL :: zpsoil_ext(nx_ext,ny_ext,nzsoil_ext) ! Height (m) (soil model)
  REAL :: t_ext(nx_ext,ny_ext,nz_ext)     ! Temperature (K)
  REAL :: qv_ext(nx_ext,ny_ext,nz_ext)    ! Specific humidity (kg/kg)
  REAL :: u_ext(nx_ext,ny_ext,nz_ext)     ! Eastward wind component (m/s)
  REAL :: v_ext(nx_ext,ny_ext,nz_ext)     ! Northward wind component (m/s)
  REAL :: uatv_ext(nx_ext,ny_ext,nz_ext)
  REAL :: vatu_ext(nx_ext,ny_ext,nz_ext)
  REAL :: w_ext(nx_ext,ny_ext,nz_ext)     ! Vertical velocity component (m/s)
  REAL :: qscalar_ext(nx_ext,ny_ext,nz_ext,nscalar)

  INTEGER soiltyp_ext (nx_ext,ny_ext,nstyps)    ! Soil type
  REAL    stypfrct_ext(nx_ext,ny_ext,nstyps)    ! Soil type
  INTEGER vegtyp_ext  (nx_ext,ny_ext)           ! Vegetation type
  REAL    lai_ext     (nx_ext,ny_ext)           ! Leaf Area Index
  REAL    roufns_ext  (nx_ext,ny_ext)           ! Surface roughness
  REAL    veg_ext     (nx_ext,ny_ext)           ! Vegetation fraction

  REAL :: tsoil_ext(nx_ext,ny_ext,nzsoil_ext,0:nstyps) ! Soil temperature(K)
  REAL :: qsoil_ext(nx_ext,ny_ext,nzsoil_ext,0:nstyps) ! Soil moisture (m3/m3)
  REAL :: wetcanp_ext(nx_ext,ny_ext,0:nstyps)  ! Canopy water amount
  REAL :: snowdpth_ext(nx_ext,ny_ext)          ! Snow depth (m)

  REAL :: ubar_ext(nx_ext,ny_ext,nz_ext)   ! Base state x velocity
                                           ! component (m/s)
  REAL :: vbar_ext(nx_ext,ny_ext,nz_ext)   ! Base state y velocity
                                           ! component (m/s)
  REAL :: wbar_ext(nx_ext,ny_ext,nz_ext)   ! Base state z velocity
                                           ! component (m/s)
  REAL :: ptbar_ext(nx_ext,ny_ext,nz_ext)  ! Base state potential
                                           ! temperature (K)
  REAL :: pbar_ext(nx_ext,ny_ext,nz_ext)   ! Base state pressure (Pascal)
  REAL :: qvbar_ext(nx_ext,ny_ext,nz_ext)  ! Base state water vapor
                                           ! mixing ratio (kg/kg)
  REAL :: tem1_ext(nx_ext,ny_ext,nz_ext)   ! Temporary work array
  REAL :: tem2_ext(nx_ext,ny_ext,nz_ext)   ! Temporary work array
!
!-----------------------------------------------------------------------
!
!  Other  external variable arrays
!
!-----------------------------------------------------------------------
!
  REAL :: x_ext(nx_ext)
  REAL :: y_ext(ny_ext)
  REAL :: z_ext(nz_ext)
  REAL :: xsc(nx_ext)
  REAL :: ysc(ny_ext)
!
!  REAL :: uprt_ext(nx_ext,ny_ext,nz_ext)   ! x velocity component (m/s)
!  REAL :: vprt_ext(nx_ext,ny_ext,nz_ext)   ! y velocity component (m/s)

!  real raing_ext  (nx_ext,ny_ext)               ! Grid supersaturation rain
!  real rainc_ext  (nx_ext,ny_ext)               ! Cumulus convective rain
!  real prcrate_ext(nx_ext,ny_ext,4)
                                ! precipitation rate (kg/(m**2*s))
! prcrate(1,1,1) = total precip. rate
! prcrate(1,1,2) = grid scale precip. rate
! prcrate(1,1,3) = cumulus precip. rate
! prcrate(1,1,4) = microphysics precip. rate

!  real radfrc_ext(nx,ny,nz) ! Radiation forcing (K/s)
!  real radsw_ext (nx,ny)    ! Solar radiation reaching the surface
!  real rnflx_ext (nx,ny)    ! Net radiation flux absorbed by surface

!  real usflx_ext (nx,ny)    ! Surface flux of u-momentum (kg/(m*s**2))
!  real vsflx_ext (nx,ny)    ! Surface flux of v-momentum (kg/(m*s**2))
!  real ptsflx_ext(nx,ny)    ! Surface heat flux (K*kg/(m*s**2))
!  real qvsflx_ext(nx,ny)    ! Surface moisture flux (kg/(m**2*s))

  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=3 ) :: fmtn
  CHARACTER (LEN=80) :: timsnd
  INTEGER :: lenrun, tmstrln
  INTEGER :: i,j,k,ldir,ireturn
  INTEGER :: ihr,imin,isec
  REAL :: govrd
  REAL :: xctr,yctr,dx_ext,dy_ext,tvbot,tvtop,tvbar

  CHARACTER (LEN=80) :: runname_org
  CHARACTER (LEN=80) :: cmnt_org(50)
  INTEGER :: nocmnt_org,mapproj_org
  INTEGER :: month_org,day_org,year_org
  INTEGER :: hour_org,minute_org,second_org

  REAL :: latnot(2)
  REAL :: umove_org,vmove_org,xgrdorg_org,ygrdorg_org
  REAL :: trulat1_org,trulat2_org,trulon_org,sclfct_org
  REAL :: latitud_org,ctrlat_org,ctrlon_org
  REAL :: dx_org,dy_org
  REAL :: lat_org,lon_org
  REAL :: dz_org,dzmin_org,zrefsfc_org,dlayer1_org,dlayer2_org
  REAL :: zflat_org,strhtune_org
  INTEGER :: strhopt_org

  CHARACTER (LEN=256) :: grdbasfn_ext
  CHARACTER (LEN=256) :: datafn_ext
  INTEGER :: nchanl_ext,lengbf_ext
  INTEGER :: lendtf_ext

  INTEGER :: mp_opt_save

  REAL :: time_ext

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
!  Build file names
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0 ) THEN

  IF ( extdfcst == '         ') extdfcst='000:00:00'

  lenrun=LEN(dir_extd)
  ldir=lenrun
  CALL strlnth( dir_extd, ldir )

  IF ( ldir == 0 .OR. dir_extd(1:ldir) == ' ' ) THEN
    dir_extd = '.'
    ldir = 1
  END IF

  IF( dir_extd(ldir:ldir) /= '/' .AND.  ldir < lenrun ) THEN
    ldir = ldir + 1
    dir_extd(ldir:ldir) = '/'
  END IF

  lenrun = LEN( extdname )
  CALL strlnth( extdname, lenrun )

  IF( extdfmt == 1 ) THEN
    fmtn = 'bin'
  ELSE IF ( extdfmt == 2 ) THEN
    fmtn = 'asc'
  ELSE IF ( extdfmt == 3 ) THEN
    fmtn = 'hdf'
  ELSE IF ( extdfmt == 4 ) THEN
    fmtn = 'pak'
  ELSE IF ( extdfmt == 6 ) THEN
    fmtn = 'bn2'
  ELSE IF ( extdfmt == 7 ) THEN
    fmtn = 'net'
  ELSE IF ( extdfmt == 8 ) THEN
    fmtn = 'npk'
  ELSE IF ( extdfmt == 9 ) THEN
    fmtn = 'gad'
  ELSE IF ( extdfmt == 10 ) THEN
    fmtn = 'grb'
  ELSE
    WRITE(6,'(a,a,a)')                                                  &
        'Unknown format, ', extdfmt, '. Program stopped in GETARPS.'
!   STOP
    CALL arpsstop("GETARPS:  unknown format",1)
  END IF

  READ(extdfcst,'(i3,1x,i2,1x,i2)') ihr,imin,isec

  time_ext = REAL( (ihr*3600)+(imin*60)+isec )
  CALL cvttsnd( time_ext, timsnd, tmstrln )

  grdbasfn_ext = dir_extd(1:ldir)//extdname(1:lenrun)                   &
                                 //'.'//fmtn//'grdbas'
  lengbf_ext = ldir + lenrun + 10

  datafn_ext = dir_extd(1:ldir)//extdname(1:lenrun)                     &
                               //'.'//fmtn//timsnd(1:tmstrln)
  lendtf_ext = ldir + lenrun + 4 + tmstrln

  WRITE(6,*) 'The external grid and base file, grdbasfn = ',            &
             grdbasfn_ext(1:lendtf_ext)

  WRITE(6,*) 'The external time dependent file,  datafn = ',            &
             datafn_ext(1:lengbf_ext)
!
!-----------------------------------------------------------------------
!
!  Since the data reader will change certain parameters stored
!  in common, they need to be saved and restored in common
!  after reading is done.
!
!-----------------------------------------------------------------------
!
  runname_org=runname
  nocmnt_org=nocmnt
  IF( nocmnt > 0 ) THEN
    DO i=1,nocmnt
      cmnt_org(i)=cmnt(i)
    END DO
  END IF
  mapproj_org=mapproj
  month_org=month
  day_org=day
  year_org=year
  hour_org=hour
  minute_org=minute
  second_org=second
!
  umove_org=umove
  vmove_org=vmove
  trulat1_org=trulat1
  trulat2_org=trulat2
  trulon_org=trulon
  sclfct_org=sclfct
  latitud_org=latitud
  ctrlat_org=ctrlat
  ctrlon_org=ctrlon
  dx_org=dx
  dy_org=dy
  xgrdorg_org=xgrdorg
  ygrdorg_org=ygrdorg
  dz_org=dz
  dzmin_org=dzmin
  zrefsfc_org=zrefsfc
  dlayer1_org=dlayer1
  dlayer2_org=dlayer2
  zflat_org=zflat
  strhtune_org=strhtune
  strhopt_org=strhopt
  CALL xytoll(1,1,0.,0.,lat_org,lon_org)
!
!-----------------------------------------------------------------------
!
!  Read ARPS data file
!
!-----------------------------------------------------------------------
!
! uatv_ext & vatu_ext used as temporary variables here
!

!
! We need to disable MPI here, as we DON'T want the data scattered to the
! individual processors!  The current code calls for the entire external
! data to be passed to all processors.  The assumption is that external
! data is much lower resolution, so this can be accomplished.  Eventually,
! this will have to be changed, but not at this release.
!

  mp_opt_save = mp_opt
  mp_opt = 0

  CALL dtaread(nx_ext,ny_ext,nz_ext,nzsoil_ext,nstyps,                  &
               extdfmt, nchanl_ext, grdbasfn_ext,lengbf_ext,            &
               datafn_ext, lendtf_ext, time_ext,                        &
               x_ext,y_ext,z_ext,zp_ext,zpsoil_ext,                     &
               u_ext,v_ext,w_ext,t_ext,p_ext,                           &
               qv_ext, qscalar_ext, vatu_ext,vatu_ext,vatu_ext,         &
               ubar_ext, vbar_ext, wbar_ext,                            &
               ptbar_ext, pbar_ext, vatu_ext, qvbar_ext,                &
               soiltyp_ext,stypfrct_ext,vegtyp_ext,                     &
               lai_ext,roufns_ext,veg_ext,                              &
               tsoil_ext,qsoil_ext,wetcanp_ext,                         &
               snowdpth_ext,                                            &
               tem1_ext(1,1,1),tem1_ext(1,1,2),tem2_ext,                &
               tem2_ext,tem1_ext(1,1,1),tem1_ext(1,1,2),                &
               tem2_ext,tem2_ext,                                       &
               tem1_ext(1,1,1),tem1_ext(1,1,2),                         &
               tem1_ext(1,1,3),tem1_ext(1,1,4),                         &
               ireturn,tem1_ext,tem2_ext,uatv_ext)
!
!-----------------------------------------------------------------------
!
!  Set maprojection parameters
!
!-----------------------------------------------------------------------
!
  iproj_ext=mapproj
  scale_ext=sclfct
  trlon_ext=trulon
  latnot_ext(1)=trulat1
  latnot_ext(2)=trulat2
  CALL setmapr(iproj_ext,scale_ext,latnot_ext,trlon_ext)
  CALL lltoxy(1,1,ctrlat,ctrlon,xctr,yctr)
  dx_ext=x_ext(2)-x_ext(1)
  x0_ext=xctr - 0.5*(nx_ext-3)*dx_ext
  dy_ext=y_ext(2)-y_ext(1)
  y0_ext=yctr - 0.5*(ny_ext-3)*dy_ext
  CALL setorig(1,x0_ext,y0_ext)
!
!-----------------------------------------------------------------------
!
!  Find lat,lon of scalar points
!
!-----------------------------------------------------------------------
!
  DO i=1,nx_ext-1
    xsc(i)=0.5*(x_ext(i)+x_ext(i+1))
  END DO
  xsc(nx_ext)=2.*xsc(nx_ext-1)-xsc(nx_ext-2)
  DO j=1,ny_ext-1
    ysc(j)=0.5*(y_ext(j)+y_ext(j+1))
  END DO
  ysc(ny_ext)=2.*ysc(ny_ext-1)-ysc(ny_ext-2)

  CALL xytoll(nx_ext,ny_ext,xsc,ysc,lat_ext,lon_ext)
  CALL xytoll(nx_ext,ny_ext,x_ext,ysc,latu_ext,lonu_ext)
  CALL xytoll(nx_ext,ny_ext,xsc,y_ext,latv_ext,lonv_ext)
!
!-----------------------------------------------------------------------
!
!  Move z field onto the scalar grid.
!
!-----------------------------------------------------------------------
!
  DO k=1,nz_ext-1
    DO j=1,ny_ext-1
      DO i=1,nx_ext-1
        hgt_ext(i,j,k)=0.5*(zp_ext(i,j,k)+zp_ext(i,j,k+1))
      END DO
    END DO
  END DO
  DO j=1,ny_ext-1
    DO i=1,nx_ext-1
      hgt_ext(i,j,nz_ext)=(2.*hgt_ext(i,j,nz_ext-1))                    &
                             -hgt_ext(i,j,nz_ext-2)
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Combine wind perturbations and mean fields.
!
!-----------------------------------------------------------------------
!
  u_ext = u_ext + ubar_ext
  v_ext = v_ext + vbar_ext
  w_ext = w_ext + wbar_ext
!
!-----------------------------------------------------------------------
!
!  Orient u & v to true north.
!
!-----------------------------------------------------------------------
!
  DO k=1,nz_ext

    ! get u at v grid point locations
    DO j=2,ny_ext
      DO i=1,nx_ext-1
        uatv_ext(i,j,k) = 0.25*(u_ext(i,j-1,k)+u_ext(i+1,j-1,k)  &
                               +u_ext(i,j,k)+u_ext(i+1,j,k))
      END DO
    END DO
    DO i=1,nx_ext-1
      uatv_ext(i,1,k) = 0.5*(u_ext(i,1,k)+u_ext(i+1,1,k))
    END DO
    DO j=2,ny_ext
      uatv_ext(nx_ext,j,k) = 0.5*(u_ext(nx_ext,j-1,k)+u_ext(nx_ext,j,k))
    END DO
    uatv_ext(nx_ext,1,k) = u_ext(nx_ext,1,k)

    ! get v at u grid point locations
    DO j=1,ny_ext-1
      DO i=2,nx_ext
        vatu_ext(i,j,k) = 0.25*(v_ext(i-1,j,k)+v_ext(i,j,k)  &
                               +v_ext(i-1,j+1,k)+v_ext(i,j+1,k))
      END DO
    END DO
    DO j=1,ny_ext-1
      vatu_ext(1,j,k) = 0.5*(v_ext(1,j,k)+v_ext(1,j+1,k))
    END DO
    DO i=2,nx_ext
      vatu_ext(i,ny_ext,k) = 0.5*(v_ext(i-1,ny_ext,k)+v_ext(i,ny_ext,k))
    END DO
    vatu_ext(1,ny_ext,k) = v_ext(1,ny_ext,k)

    CALL uvmptoe(nx_ext,ny_ext,u_ext(1,1,k),vatu_ext(1,1,k),lonu_ext,     &
       tem1_ext(1,1,k),tem2_ext(1,1,k))
    u_ext(1:nx_ext,1:ny_ext,k) = tem1_ext(1:nx_ext,1:ny_ext,k)
    vatu_ext(1:nx_ext,1:ny_ext,k) = tem2_ext(1:nx_ext,1:ny_ext,k)

    CALL uvmptoe(nx_ext,ny_ext,uatv_ext(1,1,k),v_ext(1,1,k),lonv_ext,     &
       tem1_ext(1,1,k),tem2_ext(1,1,k))
    uatv_ext(1:nx_ext,1:ny_ext,k) = tem1_ext(1:nx_ext,1:ny_ext,k)
    v_ext(1:nx_ext,1:ny_ext,k) = tem2_ext(1:nx_ext,1:ny_ext,k)

  END DO
!
!-----------------------------------------------------------------------
!
!  Combine perturbations and mean fields of scalars
!  Convert potential temperature to temperature
!
!-----------------------------------------------------------------------
!
  DO k=2,nz_ext-1
    DO j=2,ny_ext-1
      DO i=2,nx_ext-1
        p_ext(i,j,k)=p_ext(i,j,k)+pbar_ext(i,j,k)
        t_ext(i,j,k)=t_ext(i,j,k)+ptbar_ext(i,j,k)
        t_ext(i,j,k)=t_ext(i,j,k)*((p_ext(i,j,k)/p0)**rddcp)
        qv_ext(i,j,k)=qv_ext(i,j,k)+qvbar_ext(i,j,k)
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Supply data at the edge points (zero gradient, where missing)
!
!-----------------------------------------------------------------------
!
  CALL edgfill(hgt_ext,1,nx_ext,2,nx_ext-1,1,ny_ext,2,ny_ext-1,         &
                       1,nz_ext,1,nz_ext)
  CALL edgfill(p_ext,1,nx_ext,2,nx_ext-1,1,ny_ext,2,ny_ext-1,           &
                     1,nz_ext,1,nz_ext)
  CALL edgfill(t_ext,1,nx_ext,2,nx_ext-1,1,ny_ext,2,ny_ext-1,           &
                     1,nz_ext,1,nz_ext)
  CALL edgfill(qv_ext,1,nx_ext,2,nx_ext-1,1,ny_ext,2,ny_ext-1,          &
                     1,nz_ext,2,nz_ext-1)
! u & v now not on scalar points, so don't edgfill
!  CALL edgfill(u_ext,1,nx_ext,2,nx_ext-1,1,ny_ext,2,ny_ext-1,           &
!                     1,nz_ext,2,nz_ext-1)
!  CALL edgfill(v_ext,1,nx_ext,2,nx_ext-1,1,ny_ext,2,ny_ext-1,           &
!                     1,nz_ext,2,nz_ext-1)
!
!-----------------------------------------------------------------------
!
!  Make top and bottom mass fields via hydrostatic extrapolation.
!
!-----------------------------------------------------------------------
!
  govrd=g/rd
  DO j=1,ny_ext-1
    DO i=1,nx_ext-1

      t_ext(i,j,1)=2.*t_ext(i,j,2)-t_ext(i,j,3)
      tvbot=t_ext(i,j,1) * ( 1.0 + 0.608*qv_ext(i,j,1) )
      tvtop=t_ext(i,j,2) * ( 1.0 + 0.608*qv_ext(i,j,2) )
      tvbar=0.5*(tvtop+tvbot)
      p_ext(i,j,1)=p_ext(i,j,2)                                         &
                   *EXP(govrd*(hgt_ext(i,j,2)-hgt_ext(i,j,1))/tvbar)

      t_ext(i,j,nz_ext)=2.*t_ext(i,j,nz_ext-1)-t_ext(i,j,nz_ext-2)
      tvbot=t_ext(i,j,nz_ext-1)*(1.0 + 0.608*qv_ext(i,j,nz_ext-1))
      tvtop=t_ext(i,j,nz_ext)  *(1.0 + 0.608*qv_ext(i,j,nz_ext))
      tvbar=0.5*(tvtop+tvbot)
      p_ext(i,j,nz_ext)=p_ext(i,j,nz_ext-1)                             &
                   *EXP(govrd*(hgt_ext(i,j,nz_ext-1)-                   &
                        hgt_ext(i,j,nz_ext))/tvbar)
    END DO
  END DO

  CALL edgfill(p_ext,1,nx_ext,1,nx_ext-1,1,ny_ext,1,ny_ext-1,           &
                     1,nz_ext,1,nz_ext)
  CALL edgfill(t_ext,1,nx_ext,1,nx_ext-1,1,ny_ext,1,ny_ext-1,1          &
                      ,nz_ext,1,nz_ext)
!
!-----------------------------------------------------------------------
!
!  Reset info in common to original values
!
!-----------------------------------------------------------------------
!
  runname=runname_org
  nocmnt=nocmnt_org
  IF( nocmnt > 0 ) THEN
    DO i=1,nocmnt
      cmnt(i)=cmnt_org(i)
    END DO
  END IF
  mapproj=mapproj_org
  month=month_org
  day=day_org
  year=year_org
  hour=hour_org
  minute=minute_org
  second=second_org

  umove=umove_org
  vmove=vmove_org
  xgrdorg=xgrdorg_org
  ygrdorg=ygrdorg_org
  trulat1=trulat1_org
  trulat2=trulat2_org
  trulon=trulon_org
  sclfct=sclfct_org
  latitud=latitud_org
  ctrlat=ctrlat_org
  ctrlon=ctrlon_org
  dx=dx_org
  dy=dy_org
  dz=dz_org
  dzmin=dzmin_org
  zrefsfc=zrefsfc_org
  dlayer1=dlayer1_org
  dlayer2=dlayer2_org
  zflat=zflat_org
  strhtune=strhtune_org
  strhopt=strhopt_org
  xgrdorg=xgrdorg_org
  ygrdorg=ygrdorg_org

!
! Restore the correct "mp_opt" for processor 0.
!

  mp_opt = mp_opt_save
  END IF

!
!-----------------------------------------------------------------------
!
!  Reset map projection to previous values
!
!-----------------------------------------------------------------------
!
  CALL mpupdatei(mapproj,1)
  CALL mpupdater(trulat1,1)
  CALL mpupdater(trulat2,1)
  CALL mpupdater(trulon, 1)
  CALL mpupdater(sclfct,1)
  CALL mpupdater(lat_org,1)
  CALL mpupdater(lon_org,1)

  latnot(1)=trulat1
  latnot(2)=trulat2
  CALL setmapr(mapproj,sclfct,latnot,trulon)
  CALL setorig(2,lat_org,lon_org)
  CALL xytoll(1,1,0.,0.,lat_org,lon_org)

  CALL mpupdatei(iproj_ext,1)
  CALL mpupdater(scale_ext,1)
  CALL mpupdater(latnot_ext,2)
  CALL mpupdater(trlon_ext,1)
  CALL mpupdater(x0_ext,1)
  CALL mpupdater(y0_ext,1)
  CALL mpupdater(lat_ext,nx_ext*ny_ext)
  CALL mpupdater(lon_ext,nx_ext*ny_ext)

  CALL mpupdater(lat_ext,nx_ext*ny_ext)
  CALL mpupdater(lon_ext,nx_ext*ny_ext)
  CALL mpupdater(latu_ext,nx_ext*ny_ext)
  CALL mpupdater(lonu_ext,nx_ext*ny_ext)
  CALL mpupdater(latv_ext,nx_ext*ny_ext)
  CALL mpupdater(lonv_ext,nx_ext*ny_ext)
  CALL mpupdater(p_ext,nx_ext*ny_ext*nz_ext)
  CALL mpupdater(hgt_ext,nx_ext*ny_ext*nz_ext)
  CALL mpupdater(zp_ext,nx_ext*ny_ext*nz_ext)
  CALL mpupdater(zpsoil_ext,nx_ext*ny_ext*nzsoil_ext)
  CALL mpupdater(t_ext,nx_ext*ny_ext*nz_ext)
  CALL mpupdater(qv_ext,nx_ext*ny_ext*nz_ext)
  CALL mpupdater(u_ext,nx_ext*ny_ext*nz_ext)
  CALL mpupdater(v_ext,nx_ext*ny_ext*nz_ext)
  CALL mpupdater(uatv_ext,nx_ext*ny_ext*nz_ext)
  CALL mpupdater(vatu_ext,nx_ext*ny_ext*nz_ext)
  CALL mpupdater(w_ext,nx_ext*ny_ext*nz_ext)
  CALL mpupdater(qscalar_ext,nx_ext*ny_ext*nz_ext*nscalar)

  CALL mpupdatei(soiltyp_ext,nx_ext*ny_ext*nstyps)
  CALL mpupdater(stypfrct_ext,nx_ext*ny_ext*nstyps)
  CALL mpupdatei(vegtyp_ext,nx_ext*ny_ext)
  CALL mpupdater(lai_ext,nx_ext*ny_ext)
  CALL mpupdater(roufns_ext,nx_ext*ny_ext)
  CALL mpupdater(veg_ext,nx_ext*ny_ext)

  CALL mpupdater(tsoil_ext,nx_ext*ny_ext*nzsoil_ext*(nstyps+1))
  CALL mpupdater(qsoil_ext,nx_ext*ny_ext*nzsoil_ext*(nstyps+1))
  CALL mpupdater(wetcanp_ext,nx_ext*ny_ext*(nstyps+1))
  CALL mpupdater(snowdpth_ext,nx_ext*ny_ext)

  CALL mpupdater(tem1_ext,nx_ext*ny_ext*nz_ext)
  CALL mpupdater(tem2_ext,nx_ext*ny_ext*nz_ext)

  CALL mpupdater(x_ext,nx_ext)
  CALL mpupdater(y_ext,ny_ext)
  CALL mpupdater(z_ext,nz_ext)
  CALL mpupdater(xsc,nx_ext)
  CALL mpupdater(ysc,ny_ext)
  CALL mpupdater(ubar_ext,nx_ext*ny_ext*nz_ext)
  CALL mpupdater(vbar_ext,nx_ext*ny_ext*nz_ext)
  CALL mpupdater(wbar_ext,nx_ext*ny_ext*nz_ext)
  CALL mpupdater(ptbar_ext,nx_ext*ny_ext*nz_ext)
  CALL mpupdater(pbar_ext,nx_ext*ny_ext*nz_ext)
  CALL mpupdater(qvbar_ext,nx_ext*ny_ext*nz_ext)

  istatus=1
  RETURN
!
!-----------------------------------------------------------------------
!
!  Error destination
!
!-----------------------------------------------------------------------
!
!  598 CONTINUE
!  WRITE(6,'(a)')  ' Error reading last field, returning'
!  RETURN
END SUBROUTINE getarps
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETNMCRUC87                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE getnmcruc87(nx_ext,ny_ext,nz_ext,                            &
           dir_extd,extdname,extdopt,extdfmt,                           &
           extdinit,extdfcst,julfname,                                  &
           iproj_ext,scale_ext,                                         &
           trlon_ext,latnot_ext,x0_ext,y0_ext,                          &
           lat_ext,lon_ext,                                             &
           p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                      &
           tsfc_ext,tdeep_ext,wetsfc_ext,wetdp_ext,wetcanp_ext,         &
           trn_ext,psfc_ext,                                            &
           istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  ARPS version.
!
!  Reads an ARPS file for processing by ext2arps, a program
!  that converts external files to ARPS variables and format.
!  This version is useful when you want to use an ARPS file
!  with a different orientation or terrain than your final
!  ARPS product so arpsr2h does not work.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu and Jinxing Zong
!  01/18/1996
!
!  MODIFICATION HISTORY:
!    6/5/1997 Jinxing Zong and Yuhe Liu
!    Virtualization of temperature is accounted for in variable
!    conversion. Values of constants cp, rd/cp and g used in RUC are
!    adopted in making the variable conversion. Subroutine TV2TQ and
!    function ESW are added to facilitate the conversion.
!
!    01/26/1998 (Yuhe Liu)
!    Removed function esw, instead, used unified function f_es in
!    file thermolib3d.f
!
!    1999/11/30 (Gene Bassett)
!    Changed deep soil temperature and moisture to be an average from
!    0-100cm.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    dir_extd      Directory name for external file
!    extdname      Prefix string of external file name
!    extdopt       Option of external data sources
!    extdfmt       Option of external data format

!    extdinit      Initialized time in mm-dd-yyyy:hh:mm:ss format
!    extdfcst      Forecast hour in HHH:MM:SS format
!    julfname      File name in yyjjjhhmm format
!
!  OUTPUT:
!
!    iproj_ext     Map projection number of external data
!    scale_ext     Scale factor of external data
!    trlon_ext     True longitude of external data (degrees E)
!    latnot_ext(2) True latitude(s) of external data (degrees N)
!    x0_ext        x coordinate of origin of external data
!    y0_ext        y coordinate of origin of external data
!    lat_ext       latitude of external data points (degrees N)
!    lon_ext       longitude of external data points (degrees E)
!    p_ext         pressure (Pascal)
!    hgt_ext       height (m)
!    t_ext         temperature (K)
!    qv_ext        specific humidity (kg/kg)
!    u_ext         u wind component (m/s)
!    v_ext         v wind component (m/s)
!    qc_ext        Cloud water mixing ratio (kg/kg)
!    qr_ext        Rain water mixing ratio (kg/kg)
!    qi_ext        Ice mixing ratio (kg/kg)
!    qs_ext        Snow mixing ratio (kg/kg)
!    qh_ext        Hail mixing ratio (kg/kg)
!
!    tsfc_ext      Surface temperature
!    tdeep_ext     Soil temperature
!    wetsfc_ext    Top layer soil moisture (fraction)
!    wetdp_ext     Deep soil moisture (fraction)
!    wetcanp_ext   Water content on canopy
!
!    trn_ext       External terrain (m)
!    psfc_ext      Surface pressure (Pa)
!
!    istatus       status indicator
!
!  WORK ARRAYS:
!
!    var_grb3d     Arrays to store the GRIB 3-D variables:
!                  var_grb3d(nxgrb,nygrb,nzgrb,1,1) - pressure (Pa)
!                  var_grb3d(nxgrb,nygrb,nzgrb,2,1) - Montgomery stream
!                                                     function (m**2/s**2)
!                  var_grb3d(nxgrb,nygrb,nzgrb,3,1) - Potential
!                                                     temperature (K)
!                  var_grb3d(nxgrb,nygrb,nzgrb,4,1) - Condensation pressure
!                                                     of lifted parcel (Pa)
!                  var_grb3d(nxgrb,nygrb,nzgrb,5,1) - u wind (m/s)
!                  var_grb3d(nxgrb,nygrb,nzgrb,6,1) - v wind (m/s)
!                  var_grb3d(nxgrb,nygrb,nzgrb,1,2) - soil temp. (K)
!                  var_grb3d(nxgrb,nygrb,nzgrb,2,2) - soil moist.
!                                                     (fraction)
!
!    var_grb2d     Arrays to store the GRIB 2-D variables:
!                  var_grb2d(nxgrb,nygrb,1) - Ground surface temperature (K)
!                  var_grb2d(nxgrb,nygrb,2) - Geometric height (m)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'gribcst.inc'

  CHARACTER (LEN=*) :: dir_extd
  CHARACTER (LEN=*) :: extdname

  INTEGER :: extdopt
  INTEGER :: extdfmt

  CHARACTER (LEN=19) :: extdinit
  CHARACTER (LEN=9) :: extdfcst
  CHARACTER (LEN=9) :: julfname
!
!-----------------------------------------------------------------------
!
!  External grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj_ext
  REAL :: scale_ext,trlon_ext
  REAL :: latnot_ext(2)
  REAL :: x0_ext,y0_ext
  REAL :: dx_ext,dy_ext
!
!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx_ext, ny_ext, nz_ext

  REAL :: lat_ext(nx_ext,ny_ext)
  REAL :: lon_ext(nx_ext,ny_ext)
  REAL :: p_ext  (nx_ext,ny_ext,nz_ext)   ! Pressure (Pascals)
  REAL :: hgt_ext(nx_ext,ny_ext,nz_ext)   ! Height (m)
  REAL :: t_ext  (nx_ext,ny_ext,nz_ext)   ! Temperature (K)
  REAL :: qv_ext (nx_ext,ny_ext,nz_ext)   ! Specific humidity (kg/kg)
  REAL :: u_ext  (nx_ext,ny_ext,nz_ext)   ! Eastward wind component
  REAL :: v_ext  (nx_ext,ny_ext,nz_ext)   ! Northward wind component
!  REAL :: qc_ext (nx_ext,ny_ext,nz_ext)   ! Cloud H2O mixing ratio (kg/kg)
!  REAL :: qr_ext (nx_ext,ny_ext,nz_ext)   ! Rain  H2O mixing ratio (kg/kg)
!  REAL :: qi_ext (nx_ext,ny_ext,nz_ext)   ! Ice   mixing ratio (kg/kg)
!  REAL :: qs_ext (nx_ext,ny_ext,nz_ext)   ! Snow  mixing ratio (kg/kg)
!  REAL :: qh_ext (nx_ext,ny_ext,nz_ext)   ! Hail  mixing ratio (kg/kg)

  REAL :: tsfc_ext   (nx_ext,ny_ext)      ! Temperature at surface (K)
  REAL :: tdeep_ext  (nx_ext,ny_ext)      ! Deep soil temperature (K)
  REAL :: wetsfc_ext (nx_ext,ny_ext)      ! Surface soil moisture (fraction)
  REAL :: wetdp_ext  (nx_ext,ny_ext)      ! Deep soil moisture (fraction)
  REAL :: wetcanp_ext(nx_ext,ny_ext)      ! Canopy water amount

  REAL :: trn_ext    (nx_ext,ny_ext)          ! Geometrical heights
  REAL :: psfc_ext   (nx_ext,ny_ext)      ! Surface pressure (Pa)

!
!-----------------------------------------------------------------------
!
!  Work arrays for storing grib data
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: var_grb2d(:,:,:,:)   ! GRIB variables
  REAL, ALLOCATABLE :: var_grb3d(:,:,:,:,:) ! GRIB 3-D variables
  INTEGER, ALLOCATABLE :: var_lev3d(:,:,:)  ! Levels (hybrid) for
                                            ! each 3-D variable
  REAL, ALLOCATABLE :: rcdata(:)            ! temporary data array
!
!-----------------------------------------------------------------------
!
!  Other  external variable arrays
!
!-----------------------------------------------------------------------
!
  REAL :: x_ext(nx_ext)
  REAL :: y_ext(ny_ext)

  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Original grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj
  REAL :: scale,trlon,x0,y0
  REAL :: latnot(2)
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: gribfile
  CHARACTER (LEN=14)  :: gribtime
  INTEGER :: i,j,k
  INTEGER :: grbflen, grbtlen

  INTEGER :: m,n,nz1,max_nr2d,max_nr3d

  REAL :: pilev

  REAL :: tv_ext, tvc_ext
  REAL :: rovcp_p, cpd_p, g0_p, rd_p

  INTEGER :: chklev, lvscan, kk, jj

  REAL :: tema, temb

  INTEGER :: iret             ! Return flag

  REAL, ALLOCATABLE :: utmp(:,:), vtmp(:,:)
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  GRIB grid information
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=42) :: gridesc ! Grid description

  INTEGER :: iproj_grb    ! Map projection indicator
  INTEGER :: gthin        ! Indicator of whether the grid is "thinned"

  INTEGER :: ni_grb       ! Number of points along x-axis
  INTEGER :: nj_grb       ! Number of points along y-axis
  INTEGER :: np_grb       ! Total number of horizontal grid points

  INTEGER :: nk_grb       ! Number of vertical parameters
  REAL :: zk_grb(nz_ext)   ! Vertical coordinate parameters

  INTEGER :: npeq         ! Number of lat circles from pole to equator
  INTEGER :: nit(nz_ext)   ! Number of x-points for thinned grid

  REAL :: pi_grb          ! x-coordinate of pole point
  REAL :: pj_grb          ! y-coordinate of pole point
  INTEGER :: ipole        ! Projection center flag

  REAL :: di_grb          ! x-direction increment or grid length
  REAL :: dj_grb          ! y-direction increment or grid length

  REAL :: latsw           ! Latitude  of South West corner point
  REAL :: lonsw           ! Longitude of South West corner point
  REAL :: latne           ! Latitude  of North East corner point
  REAL :: lonne           ! Longitude of North East corner point

  REAL :: lattru1         ! Latitude (1st) at which projection is true
  REAL :: lattru2         ! Latitude (2nd) at which projection is true
  REAL :: lontrue         ! Longitude      at which projection is true

  REAL :: latrot          ! Latitude  of southern pole of rotation
  REAL :: lonrot          ! Longitude of southern pole of rotation
  REAL :: angrot          ! Angle of rotation

  REAL :: latstr          ! Latitude  of the pole of stretching
  REAL :: lonstr          ! Longitude of the pole of stretching
  REAL :: facstr          ! Stretching factor

  INTEGER :: scanmode     ! Scanning indicator
  INTEGER :: iscan        ! x-direction   scanning indicator
  INTEGER :: jscan        ! y-direction   scanning indicator
  INTEGER :: kscan        ! FORTRAN index scanning indicator

  INTEGER :: ires         ! Resolution direction increments indicator
  INTEGER :: iearth       ! Earth shape indicator: spherical or oblate?
  INTEGER :: icomp        ! (u,v) components decomposition indicator

  INTEGER :: jpenta       ! J-Pentagonal resolution parameter
  INTEGER :: kpenta       ! K-Pentagonal resolution parameter
  INTEGER :: mpenta       ! M-Pentagonal resolution parameter
  INTEGER :: ispect       ! Spectral representation type
  INTEGER :: icoeff       ! Spectral coefficient storage mode

  REAL :: xp_grb          ! X coordinate of sub-satellite point
  REAL :: yp_grb          ! Y coordinate of sub-satellite point
  REAL :: xo_grb          ! X coordinate of image sector origin
  REAL :: yo_grb          ! Y coordinate of image sector origin
  REAL :: zo_grb          ! Camera altitude from center of Earth
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ALLOCATE(var_grb2d(nx_ext,ny_ext,n2dvs,n2dlvt))
  ALLOCATE(var_grb3d(nx_ext,ny_ext,nz_ext,n3dvs,n3dlvt))
  ALLOCATE(rcdata(nx_ext*ny_ext))
  ALLOCATE(var_lev3d(nz_ext,n3dvs,n3dlvt))
  ALLOCATE(utmp(nx_ext,ny_ext))
  ALLOCATE(vtmp(nx_ext,ny_ext))

  CALL getgrbfname(extdopt,extdfmt,dir_extd,extdname,                   &
                   extdinit,extdfcst,'    ',                            &
                   gribfile,grbflen,gribtime,grbtlen,istatus)
  IF (istatus /=0 ) RETURN
!
!-----------------------------------------------------------------------
!
!  RDNMCGRB reads NMC GRIB data
!
!-----------------------------------------------------------------------
!
  gridtyp   = ruc87grid
  mproj_grb = ruc87proj

  n2dvars  = ruc87nvs2d
  n2dlvtps = ruc87nlvt2d

  DO m=1,n2dlvtps
    DO n=1,n2dvars
      var_id2d(n,m) = ruc87var_id2d(n,m)
    END DO
    levtyp2d(m) = ruc87levs2d(m)
  END DO

  n3dvars  = ruc87nvs3d
  n3dlvtps = ruc87nlvt3d

  DO m=1,n3dlvtps
    DO n=1,n3dvars
      var_id3d(n,m) = ruc87var_id3d(n,m)
    END DO
    levtyp3d(m) = ruc87levs3d(m)
  END DO

  CALL rdnmcgrb(nx_ext,ny_ext,nz_ext,gribfile,grbflen, gribtime,        &
                gridesc, iproj_grb, gthin,                              &
                ni_grb,nj_grb,np_grb, nk_grb,zk_grb, npeq,nit,          &
                pi_grb,pj_grb,ipole, di_grb,dj_grb,                     &
                latsw,lonsw, latne,lonne,                               &
                latrot,lonrot,angrot,                                   &
                latstr,lonstr,facstr,                                   &
                lattru1,lattru2,lontrue,                                &
                scanmode, iscan,jscan,kscan,                            &
                ires,iearth,icomp,                                      &
                jpenta,kpenta,mpenta,ispect,icoeff,                     &
                xp_grb,yp_grb, xo_grb,yo_grb,zo_grb,                    &
                rcdata,var_grb2d,var_grb3d,var_lev3d,lvldbg,iret)

  IF (iret /= 0)  THEN
    ISTATUS = -888
    GOTO 999
  END IF

  max_nr2d = 0
  DO n=1,n2dvars
    DO m=1,n2dlvtps
      max_nr2d = MAX( max_nr2d, var_nr2d(n,m) )
    END DO
  END DO

  max_nr3d = 0
  DO n=1,n3dvars
    max_nr3d = MAX( max_nr3d, var_nr3d(n,1) )
  END DO

  IF ( max_nr3d == 0 ) THEN
    WRITE (6,'(a)')                                                     &
        'No 3-D variables was found in the GRIB file',                  &
        'Program stopped in GETNMCRUC.'
!   STOP
    ISTATUS = -888
    GOTO 999
  END IF

  IF ( max_nr2d == 0 ) THEN
    WRITE (6,'(a)')                                                     &
        'No 2-D variables was found in the GRIB file'
  END IF

!  write (6,'(/a7,2x,6(i7))')
!    :      'Lev\\VID',(var_id3d(n,1),n=1,n3dvars)

!  DO 60 k=1,max_nr3d
!     write (6,'(i5,4x,6(i7))')
!    :      k,(var_lev3d(k,n,1),n=1,n3dvars)
  60    CONTINUE

  DO k=1,max_nr3d
    DO n=2,n3dvars
      IF ( var_lev3d(k,1,1) /= var_lev3d(k,n,1) ) THEN
        WRITE (6,'(a)')                                                 &
            'Variables were not at the same level.',                    &
            'Program stopped in GETNMCRUC.'
!       STOP
        ISTATUS = -888
        GOTO 999
      END IF
    END DO
  END DO

  IF ( iproj_grb == 5 .AND. ipole == 0 ) THEN      ! Center in NP
    iproj_ext = 1
  ELSE IF ( iproj_grb == 5 .AND. ipole == 1 ) THEN  ! Center in SP
    iproj_ext = -1
  ELSE IF ( iproj_grb == 3 .AND. ipole == 0 ) THEN  ! Center in NP
    iproj_ext = 2
  ELSE IF ( iproj_grb == 3 .AND. ipole == 1 ) THEN  ! Center in SP
    iproj_ext = -2
  ELSE IF ( iproj_grb == 1 ) THEN
    iproj_ext = 3
  ELSE IF ( iproj_grb == 0 ) THEN
    iproj_ext = 4
  ELSE
    WRITE (6,'(a)')                                                     &
        'Unknown map projection. Set to non-projection.'
    iproj_ext = 0
  END IF

  scale_ext = 1.0
  latnot_ext(1) = lattru1
  latnot_ext(2) = lattru2
  trlon_ext = lontrue

  dx_ext = di_grb
  dy_ext = dj_grb

  CALL getmapr(iproj,scale,latnot,trlon,x0,y0)
  CALL setmapr(iproj_ext,scale_ext,latnot_ext,trlon_ext)
  CALL lltoxy(1,1,latsw,lonsw,x0_ext,y0_ext)

  DO i=1,nx_ext
    x_ext(i)=x0_ext+(i-1)*dx_ext
  END DO

  DO j=1,ny_ext
    y_ext(j)=y0_ext+(j-1)*dy_ext
  END DO

  CALL xytoll(nx_ext,ny_ext,x_ext,y_ext,lat_ext,lon_ext)
!
!-----------------------------------------------------------------------
!
!  Retrieve 2-D variables
!
!-----------------------------------------------------------------------
!

  DO j=1,ny_ext
    DO i=1,nx_ext
      IF ( var_nr2d(1,1) == 0 ) THEN
        tsfc_ext   (i,j) = -999.0
      ELSE
        tsfc_ext   (i,j) = var_grb2d(i,j,1,1)
      END IF

      IF ( var_nr2d(2,1) == 0 ) THEN
        trn_ext    (i,j) = -999.0
      ELSE
        trn_ext    (i,j) = var_grb2d(i,j,2,1)
      END IF

      IF ( var_nr3d(1,2) == 0 ) THEN
        tsfc_ext   (i,j) = var_grb3d(i,j,1,1,2)
             ! note that this is the 5cm value and not the surface value
        tdeep_ext  (i,j) = 0.1 * var_grb3d(i,j,1,1,2)  & !   5cm
                     + 0.2 * var_grb3d(i,j,2,1,2)  & !  20cm
                     + 0.4 * var_grb3d(i,j,3,1,2)  & !  40cm
                     + 0.3 * var_grb3d(i,j,4,1,2)  ! 160cm
        wetsfc_ext (i,j) = var_grb3d(i,j,1,2,2)
             ! note that this is the 5cm value and not the surface value
        wetdp_ext  (i,j) = 0.1 * var_grb3d(i,j,1,2,2)  & !   5cm
                     + 0.2 * var_grb3d(i,j,2,2,2)  & !  20cm
                     + 0.4 * var_grb3d(i,j,3,2,2)  & !  40cm
                     + 0.3 * var_grb3d(i,j,4,2,2)  ! 160cm
        wetcanp_ext(i,j) = var_grb2d(i,j,3,1)*1.e-3  ! in meters
      ELSE
        tsfc_ext   (i,j) = -999.0
        tdeep_ext  (i,j) = -999.0
        wetsfc_ext (i,j) = -999.0
        wetdp_ext  (i,j) = -999.0
        wetcanp_ext(i,j) = -999.0
      END IF

      psfc_ext   (i,j) = -999.0
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Retrieve 3-D variables
!
!-----------------------------------------------------------------------
!

  cpd_p   = 1004.686    ! cp in RUC
  rovcp_p = 0.285714    ! rd/cp used in RUC
  g0_p    = 9.80665     ! gravity in RUC

  nz1 = MIN(var_nr3d(1,1),nz_ext)

  IF ( var_lev3d(1,1,1) < var_lev3d(nz1,1,1) ) THEN  ! 1st level at sfc
    chklev = 1
    lvscan = 0
  ELSE
    chklev = -1
    lvscan = nz1+1
  END IF

  DO k=1,nz1
    kk = chklev * k + lvscan
    DO j=1,ny_ext
      DO i=1,nx_ext

        p_ext(i,j,kk) = var_grb3d(i,j,k,1,1)       ! Pressure (Pa)

        pilev = (p_ext(i,j,kk)/100000.)**rovcp_p
        tv_ext = var_grb3d(i,j,k,3,1)*pilev        ! Virtual Temperature (K)
        hgt_ext(i,j,kk) = (var_grb3d(i,j,k,2,1)-cpd_p*tv_ext)/g0_p
                                ! Height (m) from Mongmery function with
! M = CpTv + gz

        u_ext(i,j,kk) = var_grb3d(i,j,k,5,1)       ! u wind (m/s)
        v_ext(i,j,kk) = var_grb3d(i,j,k,6,1)       ! v wind (m/s)

        tvc_ext = var_grb3d(i,j,k,3,1)                                  &
                     * (var_grb3d(i,j,k,4,1)/100000.)**rovcp_p
                              ! Virtual Temperature (K) at LCL
        CALL tv2tq( tvc_ext, var_grb3d(i,j,k,4,1), tema, temb )
        t_ext(i,j,kk)  = tema                                           &
            * (p_ext(i,j,kk)/var_grb3d(i,j,k,4,1))**rovcp_p  ! Temperature (K)
        qv_ext(i,j,kk) = temb                      ! Specific humidity
!        qc_ext(i,j,kk) = -999.
!        qr_ext(i,j,kk) = -999.
!        qi_ext(i,j,kk) = -999.
!        qs_ext(i,j,kk) = -999.
!        qh_ext(i,j,kk) = -999.
      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Rotate winds to be relative to true north.
!  The RUC data are sent as grid-relative.
!
!-----------------------------------------------------------------------
!
  DO k=1,nz1
    CALL uvmptoe(nx_ext,ny_ext,u_ext(1,1,k),v_ext(1,1,k),               &
                 lon_ext,utmp,vtmp)
    u_ext(:,:,k) = utmp(:,:)
    v_ext(:,:,k) = vtmp(:,:)
  END DO

  istatus = 1

!
!-----------------------------------------------------------------------
!
!  Reset map projection to previous values
!
!-----------------------------------------------------------------------
!

  999  CONTINUE
  CALL setmapr(iproj,scale,latnot,trlon)
  CALL setorig(1,x0,y0)


  DEALLOCATE(var_grb2d,var_grb3d,rcdata,var_lev3d)
  DEALLOCATE(utmp)
  DEALLOCATE(vtmp)

  RETURN
END SUBROUTINE getnmcruc87
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETNMCRUC211               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getnmcruc211(nx_ext,ny_ext,nz_ext,                           &
           dir_extd,extdname,extdopt,extdfmt,                           &
           extdinit,extdfcst,julfname,                                  &
           iproj_ext,scale_ext,                                         &
           trlon_ext,latnot_ext,x0_ext,y0_ext,                          &
           lat_ext,lon_ext,                                             &
           p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                      &
           tsfc_ext,tdeep_ext,wetsfc_ext,wetdp_ext,wetcanp_ext,         &
           trn_ext,psfc_ext, t_2m_ext, rh_2m_ext, u_10m_ext,            &
           v_10m_ext, istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  ARPS version.
!
!  Reads an ARPS file for processing by ext2arps, a program
!  that converts external files to ARPS variables and format.
!  This version is useful when you want to use an ARPS file
!  with a different orientation or terrain than your final
!  ARPS product so arpsr2h does not work.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dennis Moon, SSESCO
!  06/19/1996
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    dir_extd      Directory name for external file
!    extdname      Prefix string of external file name
!    extdopt       Option of external data sources
!    extdfmt       Option of external data format

!    extdinit      Initialized time in mm-dd-yyyy:hh:mm:ss format
!    extdfcst      Forecast hour in HHH:MM:SS format
!    julfname      File name in yyjjjhhmm format
!
!  OUTPUT:
!
!    iproj_ext     Map projection number of external data
!    scale_ext     Scale factor of external data
!    trlon_ext     True longitude of external data (degrees E)
!    latnot_ext(2) True latitude(s) of external data (degrees N)
!    x0_ext        x coordinate of origin of external data
!    y0_ext        y coordinate of origin of external data
!    lat_ext       latitude of external data points (degrees N)
!    lon_ext       longitude of external data points (degrees E)
!    p_ext         pressure (Pascal)
!    hgt_ext       height (m)
!    t_ext         temperature (K)
!    qv_ext        specific humidity (kg/kg)
!    u_ext         u wind component (m/s)
!    v_ext         v wind component (m/s)
!    qc_ext        Cloud water mixing ratio (kg/kg)
!    qr_ext        Rain water mixing ratio (kg/kg)
!    qi_ext        Ice mixing ratio (kg/kg)
!    qs_ext        Snow mixing ratio (kg/kg)
!    qh_ext        Hail mixing ratio (kg/kg)
!
!    tsfc_ext      Surface temperature
!    tdeep_ext     Soil temperature
!    wetsfc_ext    Top layer soil moisture (fraction)
!    wetdp_ext     Deep soil moisture (fraction)
!    wetcanp_ext   Water content on canopy
!
!    trn_ext       External terrain (m)
!    psfc_ext      Surface pressure (Pa)
!
!    T_2m_ext      Temperature at 2m AGL
!    rh_2m_ext     Specific Humidity at 2m AGL
!    U_10m_ext     U at 10m AGL
!    V_10m_ext     V at 10m AGL
!
!    istatus       status indicator
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'gribcst.inc'

  CHARACTER (LEN=*) :: dir_extd
  CHARACTER (LEN=*) :: extdname

  INTEGER :: extdopt
  INTEGER :: extdfmt

  CHARACTER (LEN=19) :: extdinit
  CHARACTER (LEN=9) :: extdfcst
  CHARACTER (LEN=9) :: julfname
!
!-----------------------------------------------------------------------
!
!  External grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj_ext
  REAL :: scale_ext,trlon_ext
  REAL :: latnot_ext(2)
  REAL :: x0_ext,y0_ext
  REAL :: dx_ext,dy_ext
!
!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx_ext,ny_ext,nz_ext

  REAL :: lat_ext(nx_ext,ny_ext)
  REAL :: lon_ext(nx_ext,ny_ext)
  REAL :: p_ext  (nx_ext,ny_ext,nz_ext)   ! Pressure (Pascals)
  REAL :: hgt_ext(nx_ext,ny_ext,nz_ext)   ! Height (m)
  REAL :: t_ext  (nx_ext,ny_ext,nz_ext)   ! Temperature (K)
  REAL :: qv_ext (nx_ext,ny_ext,nz_ext)   ! Specific humidity (kg/kg)
  REAL :: u_ext  (nx_ext,ny_ext,nz_ext)   ! Eastward wind component
  REAL :: v_ext  (nx_ext,ny_ext,nz_ext)   ! Northward wind component
!  REAL :: qc_ext (nx_ext,ny_ext,nz_ext)   ! Cloud H2O mixing ratio (kg/kg)
!  REAL :: qr_ext (nx_ext,ny_ext,nz_ext)   ! Rain  H2O mixing ratio (kg/kg)
!  REAL :: qi_ext (nx_ext,ny_ext,nz_ext)   ! Ice   mixing ratio (kg/kg)
!  REAL :: qs_ext (nx_ext,ny_ext,nz_ext)   ! Snow  mixing ratio (kg/kg)
!  REAL :: qh_ext (nx_ext,ny_ext,nz_ext)   ! Hail  mixing ratio (kg/kg)

  REAL :: tsfc_ext   (nx_ext,ny_ext)      ! Temperature at surface (K)
  REAL :: tdeep_ext  (nx_ext,ny_ext)      ! Deep soil temperature (K)
  REAL :: wetsfc_ext (nx_ext,ny_ext)      ! Surface soil moisture (fraction)
  REAL :: wetdp_ext  (nx_ext,ny_ext)      ! Deep soil moisture (fraction)
  REAL :: wetcanp_ext(nx_ext,ny_ext)      ! Canopy water amount

  REAL :: trn_ext    (nx_ext,ny_ext)          ! Geometrical heights
  REAL :: psfc_ext   (nx_ext,ny_ext)      ! Surface pressure (Pa)

  REAL :: t_2m_ext (nx_ext,ny_ext)
  REAL :: rh_2m_ext(nx_ext,ny_ext)
  REAL :: u_10m_ext(nx_ext,ny_ext)
  REAL :: v_10m_ext(nx_ext,ny_ext)
!
!----------------------------------------------------------------------
!
!  Other  external variable arrays
!
!----------------------------------------------------------------------
!
  REAL :: x_ext(nx_ext)
  REAL :: y_ext(ny_ext)
  REAL :: z_ext(nz_ext)

  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Work arrays for storing grib data
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: var_grb2d(:,:,:,:)   ! GRIB variables
  REAL, ALLOCATABLE :: var_grb3d(:,:,:,:,:) ! GRIB 3-D variables
  INTEGER, ALLOCATABLE :: var_lev3d(:,:,:)  ! Levels (hybrid) for
                                            ! each 3-D variable
  REAL, ALLOCATABLE :: rcdata(:)            ! temporary data array
!
!-----------------------------------------------------------------------
!
!  Original grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj
  REAL :: scale,trlon,x0,y0
  REAL :: latnot(2)
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: gribfile
  CHARACTER (LEN=14 ) :: gribtime
  CHARACTER (LEN=10 ) :: runstr
  CHARACTER (LEN=3  ) :: fmtn
  INTEGER :: ichr,bchar,echar
  INTEGER :: i,j,k,ldir,ireturn
  INTEGER :: grbflen, grbtlen

  INTEGER :: m,n,nz1,max_nr2d,max_nr3d

  REAL :: qvsat, pilev, qvsatice

  REAL :: tema, temb

  INTEGER :: chklev, lvscan, kk, jj

  INTEGER :: iret             ! Return flag

  REAL, ALLOCATABLE :: utmp(:,:), vtmp(:,:)
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  GRIB grid information
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=42) :: gridesc ! Grid description

  INTEGER :: iproj_grb    ! Map projection indicator
  INTEGER :: gthin        ! Indicator of whether the grid is "thinned"

  INTEGER :: ni_grb       ! Number of points along x-axis
  INTEGER :: nj_grb       ! Number of points along y-axis
  INTEGER :: np_grb       ! Total number of horizontal grid points

  INTEGER :: nk_grb       ! Number of vertical parameters
  REAL :: zk_grb(nz_ext)  ! Vertical coordinate parameters

  INTEGER :: npeq         ! Number of lat circles from pole to equator
  INTEGER :: nit(nz_ext)  ! Number of x-points for thinned grid

  REAL :: pi_grb          ! x-coordinate of pole point
  REAL :: pj_grb          ! y-coordinate of pole point
  INTEGER :: ipole        ! Projection center flag

  REAL :: di_grb          ! x-direction increment or grid length
  REAL :: dj_grb          ! y-direction increment or grid length

  REAL :: lrb             ! y-direction increment or grid length

  REAL :: latsw           ! Latitude  of South West corner point
  REAL :: lonsw           ! Longitude of South West corner point
  REAL :: latne           ! Latitude  of North East corner point
  REAL :: lonne           ! Longitude of North East corner point

  REAL :: lattru1         ! Latitude (1st) at which projection is true
  REAL :: lattru2         ! Latitude (2nd) at which projection is true
  REAL :: lontrue         ! Longitude      at which projection is true

  REAL :: latrot          ! Latitude  of southern pole of rotation
  REAL :: lonrot          ! Longitude of southern pole of rotation
  REAL :: angrot          ! Angle of rotation

  REAL :: latstr          ! Latitude  of the pole of stretching
  REAL :: lonstr          ! Longitude of the pole of stretching
  REAL :: facstr          ! Stretching factor

  INTEGER :: scanmode     ! Scanning indicator
  INTEGER :: iscan        ! x-direction   scanning indicator
  INTEGER :: jscan        ! y-direction   scanning indicator
  INTEGER :: kscan        ! FORTRAN index scanning indicator

  INTEGER :: ires         ! Resolution direction increments indicator
  INTEGER :: iearth       ! Earth shape indicator: spherical or oblate?
  INTEGER :: icomp        ! (u,v) components decomposition indicator

  INTEGER :: jpenta       ! J-Pentagonal resolution parameter
  INTEGER :: kpenta       ! K-Pentagonal resolution parameter
  INTEGER :: mpenta       ! M-Pentagonal resolution parameter
  INTEGER :: ispect       ! Spectral representation type
  INTEGER :: icoeff       ! Spectral coefficient storage mode

  REAL :: xp_grb          ! X coordinate of sub-satellite point
  REAL :: yp_grb          ! Y coordinate of sub-satellite point
  REAL :: xo_grb          ! X coordinate of image sector origin
  REAL :: yo_grb          ! Y coordinate of image sector origin
  REAL :: zo_grb          ! Camera altitude from center of Earth
!
!-----------------------------------------------------------------------
!
!  Function f_qvsatl and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_qvsatl
!fpp$ expand (f_desdt)
!fpp$ expand (f_qvsatl)
!!dir$ inline always f_desdt, f_qvsatl
!*$*  inline routine (f_desdt, f_qvsatl)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ALLOCATE(var_grb2d(nx_ext,ny_ext,n2dvs,n2dlvt))
  ALLOCATE(var_grb3d(nx_ext,ny_ext,nz_ext,n3dvs,n3dlvt))
  ALLOCATE(rcdata(nx_ext*ny_ext))
  ALLOCATE(var_lev3d(nz_ext,n3dvs,n3dlvt))
  ALLOCATE(utmp(nx_ext,ny_ext))
  ALLOCATE(vtmp(nx_ext,ny_ext))

  CALL getgrbfname(extdopt,extdfmt,dir_extd,extdname,                   &
                   extdinit,extdfcst,'    ',                            &
                   gribfile,grbflen,gribtime,grbtlen,istatus)
  IF (istatus /=0 ) RETURN
!
!-----------------------------------------------------------------------
!
!  RDNMCGRB reads NMC GRIB data
!
!-----------------------------------------------------------------------
!
  gridtyp   = ruc211grid
  mproj_grb = ruc211proj

  n2dvars  = ruc211nvs2d
  n2dlvtps = ruc211nlvt2d

  DO k=1,n2dlvtps
    DO n=1,n2dvars
      var_id2d(n,k) = ruc211var_id2d(n,k)
    END DO
    levtyp2d(k) = ruc211levs2d(k)
  END DO

  n3dvars  = ruc211nvs3d
  n3dlvtps = ruc211nlvt3d

  DO m=1,n3dlvtps
    DO n=1,n3dvars
      var_id3d(n,m) = ruc211var_id3d(n,m)
    END DO
    levtyp3d(m) = ruc211levs3d(m)
  END DO

  CALL rdnmcgrb(nx_ext,ny_ext,nz_ext,gribfile,grbflen, gribtime,        &
                gridesc, iproj_grb, gthin,                              &
                ni_grb,nj_grb,np_grb, nk_grb,zk_grb, npeq,nit,          &
                pi_grb,pj_grb,ipole, di_grb,dj_grb,                     &
                latsw,lonsw, latne,lonne,                               &
                latrot,lonrot,angrot,                                   &
                latstr,lonstr,facstr,                                   &
                lattru1,lattru2,lontrue,                                &
                scanmode, iscan,jscan,kscan,                            &
                ires,iearth,icomp,                                      &
                jpenta,kpenta,mpenta,ispect,icoeff,                     &
                xp_grb,yp_grb, xo_grb,yo_grb,zo_grb,                    &
                rcdata,var_grb2d,var_grb3d,var_lev3d,lvldbg,iret)

  IF (iret /= 0)  THEN
    ISTATUS = -888
    GOTO 999
  END IF

  max_nr2d = 0
  DO n=1,n2dvars
    DO m=1,n2dlvtps
      max_nr2d = MAX( max_nr2d, var_nr2d(n,m) )
    END DO
  END DO

  max_nr3d = 0
  DO n=1,n3dvars
    max_nr3d = MAX( max_nr3d, var_nr3d(n,1))
  END DO

  IF ( max_nr3d == 0 ) THEN
    WRITE (6,'(a)')                                                     &
        'No 3-D variables was found in the GRIB file',                  &
        'Program stopped in GETNMCRUC211.'
!   STOP
    ISTATUS = -888
    GOTO 999
  END IF

  IF ( max_nr2d == 0 ) THEN
    WRITE (6,'(a)')                                                     &
        'No 2-D variables was found in the GRIB file'
  END IF

!  write (6,'(/a7,2x,6(i7))')
!    :  'Lev\\VID',(var_id3d(n,1),n=1,n3dvars)

!  DO 60 k=1,max_nr3d
!    write (6,'(i5,4x,6(i7))')
!    :  k,(var_lev3d(k,n,1),n=1,n3dvars)
  60    CONTINUE

  DO k=1,max_nr3d
    DO n=2,n3dvars
      IF ( var_lev3d(k,1,1) /= var_lev3d(k,n,1) ) THEN
        WRITE (6,'(a)')                                                 &
            'Variables were not at the same level.',                    &
            'Program stopped in GETNMCRUC211.'
!       STOP
        ISTATUS = -888
        GOTO 999
      END IF
    END DO
  END DO

  IF ( iproj_grb == 5 .AND. ipole == 0 ) THEN      ! Center in NP
    iproj_ext = 1
  ELSE IF ( iproj_grb == 5 .AND. ipole == 1 ) THEN  ! Center in SP
    iproj_ext = -1
  ELSE IF ( iproj_grb == 3 .AND. ipole == 0 ) THEN  ! Center in NP
    iproj_ext = 2
  ELSE IF ( iproj_grb == 3 .AND. ipole == 1 ) THEN  ! Center in SP
    iproj_ext = -2
  ELSE IF ( iproj_grb == 1 ) THEN
    iproj_ext = 3
  ELSE IF ( iproj_grb == 0 ) THEN
    iproj_ext = 4
  ELSE
    WRITE (6,'(a)')                                                     &
        'Unknown map projection. Set to non-projection.'
    iproj_ext = 0
  END IF

  scale_ext = 1.0
  latnot_ext(1) = lattru1
  latnot_ext(2) = lattru2
  trlon_ext = lontrue

  dx_ext = di_grb
  dy_ext = dj_grb

  CALL getmapr(iproj,scale,latnot,trlon,x0,y0)
  CALL setmapr(iproj_ext,scale_ext,latnot_ext,trlon_ext)
  CALL lltoxy(1,1,latsw,lonsw,x0_ext,y0_ext)

  DO i=1,nx_ext
    x_ext(i)=x0_ext+(i-1)*dx_ext
  END DO

  DO j=1,ny_ext
    y_ext(j)=y0_ext+(j-1)*dy_ext
  END DO

  CALL xytoll(nx_ext,ny_ext,x_ext,y_ext,lat_ext,lon_ext)
!
!-----------------------------------------------------------------------
!
!  Retrieve 2-D variables
!
!-----------------------------------------------------------------------
!
  DO j=1,ny_ext
    DO i=1,nx_ext
      IF ( var_nr2d(1,1) == 0 ) THEN
        psfc_ext   (i,j) = -999.0
      ELSE
        psfc_ext   (i,j) = var_grb2d(i,j,1,1)
      END IF

      IF ( var_nr2d(2,1) == 0 ) THEN
        trn_ext    (i,j) = -999.0
      ELSE
        trn_ext    (i,j) = var_grb2d(i,j,2,1)
      END IF

      IF( var_nr2d(1,2) == 0 ) THEN
        t_2m_ext(i,j)= -999.0
      ELSE
        t_2m_ext(i,j)= var_grb2d(i,j,1,2)
      END IF

!   at this point we are reading in the relative humidity
!   later we'll convert to specific humidity

      IF( var_nr2d(2,2) == 0 ) THEN
        rh_2m_ext(i,j)= -999.0
      ELSE
        rh_2m_ext(i,j)= var_grb2d(i,j,2,2)
      END IF

      IF( var_nr2d(3,2) == 0 ) THEN
        u_10m_ext(i,j)= -999.0
      ELSE
        u_10m_ext(i,j)= var_grb2d(i,j,3,2)
      END IF

      IF( var_nr2d(4,2) == 0 ) THEN
        v_10m_ext(i,j)= -999.0
      ELSE
        v_10m_ext(i,j)= var_grb2d(i,j,4,2)
      END IF

      tsfc_ext   (i,j) = -999.0
      tdeep_ext  (i,j) = -999.0
      wetsfc_ext (i,j) = -999.0
      wetdp_ext  (i,j) = -999.0
      wetcanp_ext(i,j) = -999.0

    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Retrieve 3-D variables
!
!-----------------------------------------------------------------------
!
  nz1 = MIN(var_nr3d(1,1),nz_ext)

  IF ( var_lev3d(1,1,1) > var_lev3d(nz1,1,1) ) THEN  ! 1st level at sfc
    chklev = 1
    lvscan = 0
  ELSE
    chklev = -1
    lvscan = nz1+1
  END IF

  DO k=1,nz1
    kk = chklev * k + lvscan
    DO j=1,ny_ext
      DO i=1,nx_ext
        p_ext(i,j,kk) = 100.0 * REAL(var_lev3d(k,1,1)) ! Pressure
        t_ext(i,j,kk) = var_grb3d(i,j,k,2,1)    ! Temperature (K)
        u_ext(i,j,kk) = var_grb3d(i,j,k,4,1)    ! u wind (m/s)
        v_ext(i,j,kk) = var_grb3d(i,j,k,5,1)    ! v wind (m/s)
        hgt_ext(i,j,kk) = var_grb3d(i,j,k,1,1)    ! height (m)

!   check for portions of constant pressure grids that are below the surface

!       IF(((kk == 1) .OR. (p_ext(i,j,kk) > psfc_ext(i,j))) .AND.       &
!             (psfc_ext(i,j) > 0.) ) THEN
!         p_ext(i,j,kk)= psfc_ext(i,j) - 1. * (kk - 1)
!         t_ext(i,j,kk)= t_2m_ext(i,j)
!         u_ext(i,j,kk)= u_10m_ext(i,j)
!         v_ext(i,j,kk)= v_10m_ext(i,j)
!         hgt_ext(i,j,kk)= trn_ext(i,j) + 0.1 * (kk - 1)
!       END IF

        IF( (p_ext(i,j,kk) > 0.0) .AND. (t_ext(i,j,kk) > 0.0) ) THEN
          qvsat = f_qvsatl( p_ext(i,j,kk), t_ext(i,j,kk) )
          qv_ext(i,j,kk)= var_grb3d(i,j,k,3,1) * qvsat * 0.01

          IF((kk == 1) .OR. (p_ext(i,j,kk) > psfc_ext(i,j)))THEN
            qv_ext(i,j,kk)= rh_2m_ext(i,j) * qvsat * 0.01
          END IF
        ELSE
          qv_ext(i,j,kk)= 0.0
        END IF

!        qc_ext(i,j,kk)= 0.0
!        qi_ext(i,j,kk)= 0.0

!        qr_ext (i,j,kk) = -999.
!        qs_ext (i,j,kk) = -999.
!        qh_ext (i,j,kk) = -999.
      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Rotate winds to be relative to true north.
!  The RUCawips data are sent as grid-relative.
!
!-----------------------------------------------------------------------
!
  DO k=1,nz1
    CALL uvmptoe(nx_ext,ny_ext,u_ext(1,1,k),v_ext(1,1,k),               &
                 lon_ext,utmp,vtmp)
    u_ext(:,:,k) = utmp(:,:)
    v_ext(:,:,k) = vtmp(:,:)
  END DO

  istatus = 1

!
!-----------------------------------------------------------------------
!
!  Reset map projection to previous values
!
!-----------------------------------------------------------------------
!
  999  CONTINUE
  CALL setmapr(iproj,scale,latnot,trlon)
  CALL setorig(1,x0,y0)

  DEALLOCATE(var_grb2d,var_grb3d,rcdata,var_lev3d)
  DEALLOCATE(utmp)
  DEALLOCATE(vtmp)

  RETURN
END SUBROUTINE getnmcruc211
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETREANALT62               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getreanalt62(nx_ext,ny_ext,nz_ext,                           &
           dir_extd,extdname,extdopt,extdfmt,                           &
           extdinit,extdfcst,julfname,                                  &
           iproj_ext,scale_ext,                                         &
           trlon_ext,latnot_ext,x0_ext,y0_ext,                          &
           lat_ext,lon_ext,                                             &
           p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                      &
           tsfc_ext,tdeep_ext,wetsfc_ext,wetdp_ext,wetcanp_ext,         &
           trn_ext,psfc_ext, istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!
!  Reads a Global ReAnalaysis file for processing by ext2arps,
!  a program that converts external files to ARPS variables and format.
!  This version is useful when you want to use an ARPS file
!  with a different orientation or terrain than your final
!  ARPS product so arpsr2h does not work.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dennis Moon, SSESCO
!  06/19/1998
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    dir_extd      Directory name for external file
!    extdname      Prefix string of external file name
!    extdopt       Option of external data sources
!    extdfmt       Option of external data format
!
!    extdinit      Initialized time in mm-dd-yyyy:hh:mm:ss format
!    extdfcst      Forecast hour in HHH:MM:SS format
!    julfname      File name in yyjjjhhmm format
!
!  OUTPUT:
!
!    iproj_ext     Map projection number of external data
!    scale_ext     Scale factor of external data
!    trlon_ext     True longitude of external data (degrees E)
!    latnot_ext(2) True latitude(s) of external data (degrees N)
!    x0_ext        x coordinate of origin of external data
!    y0_ext        y coordinate of origin of external data
!    lat_ext       latitude of external data points (degrees N)
!    lon_ext       longitude of external data points (degrees E)
!    p_ext         pressure (Pascal)
!    hgt_ext       height (m)
!    t_ext         temperature (K)
!    qv_ext        specific humidity (kg/kg)
!    u_ext         u wind component (m/s)
!    v_ext         v wind component (m/s)
!    qc_ext        Cloud water mixing ratio (kg/kg)
!    qr_ext        Rain water mixing ratio (kg/kg)
!    qi_ext        Ice mixing ratio (kg/kg)
!    qs_ext        Snow mixing ratio (kg/kg)
!    qh_ext        Hail mixing ratio (kg/kg)
!
!    tsfc_ext      Surface temperature
!    tdeep_ext     Soil temperature
!    wetsfc_ext    Top layer soil moisture (fraction)
!    wetdp_ext     Deep soil moisture (fraction)
!    wetcanp_ext   Water content on canopy
!
!    trn_ext       External terrain (m)
!    psfc_ext      Surface pressure (Pa)
!
!    istatus       status indicator
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'gribcst.inc'

  CHARACTER (LEN=*) :: dir_extd
  CHARACTER (LEN=*) :: extdname

  INTEGER :: extdopt
  INTEGER :: extdfmt

  CHARACTER (LEN=19) :: extdinit
  CHARACTER (LEN=9) :: extdfcst
  CHARACTER (LEN=9) :: julfname
!
!-----------------------------------------------------------------------
!
!  External grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj_ext
  REAL :: scale_ext,trlon_ext
  REAL :: latnot_ext(2)
  REAL :: x0_ext,y0_ext
  REAL :: dx_ext,dy_ext
!
!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx_ext,ny_ext,nz_ext

  REAL :: lat_ext(nx_ext,ny_ext)
  REAL :: lon_ext(nx_ext,ny_ext)
  REAL :: p_ext  (nx_ext,ny_ext,nz_ext)   ! Pressure (Pascals)
  REAL :: hgt_ext(nx_ext,ny_ext,nz_ext)   ! Height (m)
  REAL :: t_ext  (nx_ext,ny_ext,nz_ext)   ! Temperature (K)
  REAL :: qv_ext (nx_ext,ny_ext,nz_ext)   ! Specific humidity (kg/kg)
  REAL :: u_ext  (nx_ext,ny_ext,nz_ext)   ! Eastward wind component
  REAL :: v_ext  (nx_ext,ny_ext,nz_ext)   ! Northward wind component
!  REAL :: qc_ext (nx_ext,ny_ext,nz_ext)   ! Cloud H2O mixing ratio (kg/kg)
!  REAL :: qr_ext (nx_ext,ny_ext,nz_ext)   ! Rain  H2O mixing ratio (kg/kg)
!  REAL :: qi_ext (nx_ext,ny_ext,nz_ext)   ! Ice   mixing ratio (kg/kg)
!  REAL :: qs_ext (nx_ext,ny_ext,nz_ext)   ! Snow  mixing ratio (kg/kg)
!  REAL :: qh_ext (nx_ext,ny_ext,nz_ext)   ! Hail  mixing ratio (kg/kg)
!
  REAL :: tsfc_ext   (nx_ext,ny_ext)      ! Temperature at surface (K)
  REAL :: tdeep_ext  (nx_ext,ny_ext)      ! Deep soil temperature (K)
  REAL :: wetsfc_ext (nx_ext,ny_ext)      ! Surface soil moisture (fraction)
  REAL :: wetdp_ext  (nx_ext,ny_ext)      ! Deep soil moisture (fraction)
  REAL :: wetcanp_ext(nx_ext,ny_ext)      ! Canopy water amount

  REAL :: trn_ext    (nx_ext,ny_ext)          ! Geometrical heights
  REAL :: psfc_ext   (nx_ext,ny_ext)      ! Surface pressure (Pa)

!
!----------------------------------------------------------------------
!
!  Other  external variable arrays
!
!----------------------------------------------------------------------
!
  REAL :: x_ext(nx_ext)
  REAL :: y_ext(ny_ext)
  REAL :: z_ext(nz_ext)

  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Work arrays for storing grib data
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: var_grb2d(:,:,:,:)   ! GRIB variables
  REAL, ALLOCATABLE :: var_grb3d(:,:,:,:,:) ! GRIB 3-D variables
  INTEGER, ALLOCATABLE :: var_lev3d(:,:,:)  ! Levels (hybrid) for
                                            ! each 3-D variable
  REAL, ALLOCATABLE :: rcdata(:)            ! temporary data array
!
!-----------------------------------------------------------------------
!
!  Original grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj
  REAL :: scale,trlon,x0,y0
  REAL :: latnot(2)
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: gribfile
  CHARACTER (LEN=14 ) :: gribtime
  CHARACTER (LEN=10 ) :: runstr
  CHARACTER (LEN=3  ) :: fmtn
  INTEGER :: ichr,bchar,echar
  INTEGER :: i,j,k,ldir,ireturn
  INTEGER :: grbflen, grbtlen

  INTEGER :: m,n,nz1,max_nr2d,max_nr3d

  REAL :: qvsat, pilev, qvsatice

  REAL :: tema, temb, qvavg, mixrat, tvirtbar, tmpval
  REAL :: gausslat(200)

  INTEGER :: chklev, lvscan, kk, jj

  INTEGER :: iret             ! Return flag
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  GRIB grid information
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=42) :: gridesc ! Grid description

  INTEGER :: iproj_grb    ! Map projection indicator
  INTEGER :: gthin        ! Indicator of whether the grid is "thinned"

  INTEGER :: ni_grb       ! Number of points along x-axis
  INTEGER :: nj_grb       ! Number of points along y-axis
  INTEGER :: np_grb       ! Total number of horizontal grid points

  INTEGER :: nk_grb       ! Number of vertical parameters
  REAL :: zk_grb(nz_ext)  ! Vertical coordinate parameters

  INTEGER :: npeq         ! Number of lat circles from pole to equator
  INTEGER :: nit(nz_ext)  ! Number of x-points for thinned grid

  REAL :: pi_grb          ! x-coordinate of pole point
  REAL :: pj_grb          ! y-coordinate of pole point
  INTEGER :: ipole        ! Projection center flag

  REAL :: di_grb          ! x-direction increment or grid length
  REAL :: dj_grb          ! y-direction increment or grid length

  REAL :: lrb             ! y-direction increment or grid length

  REAL :: latsw           ! Latitude  of South West corner point
  REAL :: lonsw           ! Longitude of South West corner point
  REAL :: latne           ! Latitude  of North East corner point
  REAL :: lonne           ! Longitude of North East corner point

  REAL :: lattru1         ! Latitude (1st) at which projection is true
  REAL :: lattru2         ! Latitude (2nd) at which projection is true
  REAL :: lontrue         ! Longitude      at which projection is true

  REAL :: latrot          ! Latitude  of southern pole of rotation
  REAL :: lonrot          ! Longitude of southern pole of rotation
  REAL :: angrot          ! Angle of rotation

  REAL :: latstr          ! Latitude  of the pole of stretching
  REAL :: lonstr          ! Longitude of the pole of stretching
  REAL :: facstr          ! Stretching factor

  INTEGER :: scanmode     ! Scanning indicator
  INTEGER :: iscan        ! x-direction   scanning indicator
  INTEGER :: jscan        ! y-direction   scanning indicator
  INTEGER :: kscan        ! FORTRAN index scanning indicator

  INTEGER :: ires         ! Resolution direction increments indicator
  INTEGER :: iearth       ! Earth shape indicator: spherical or oblate?
  INTEGER :: icomp        ! (u,v) components decomposition indicator

  INTEGER :: jpenta       ! J-Pentagonal resolution parameter
  INTEGER :: kpenta       ! K-Pentagonal resolution parameter
  INTEGER :: mpenta       ! M-Pentagonal resolution parameter
  INTEGER :: ispect       ! Spectral representation type
  INTEGER :: icoeff       ! Spectral coefficient storage mode

  REAL :: xp_grb          ! X coordinate of sub-satellite point
  REAL :: yp_grb          ! Y coordinate of sub-satellite point
  REAL :: xo_grb          ! X coordinate of image sector origin
  REAL :: yo_grb          ! Y coordinate of image sector origin
  REAL :: zo_grb          ! Camera altitude from center of Earth
!
!-----------------------------------------------------------------------
!
!  T62 Gaussian grid latitudes
!
!-----------------------------------------------------------------------
!
  DATA gausslat /  88.5420,  86.6532,  84.7532,  82.8508,  80.9474,     &
                   79.0435,  77.1393,  75.2351,  73.3307,  71.4262,     &
                   69.5217,  67.6171,  65.7125,  63.8079,  61.9033,     &
                   59.9986,  58.0940,  56.1893,  54.2846,  52.3799,     &
                   50.4752,  48.5705,  46.6658,  44.7611,  42.8564,     &
                   40.9517,  39.0470,  37.1422,  35.2375,  33.3328,     &
                   31.4281,  29.5234,  27.6186,  25.7139,  23.8092,     &
                   21.9044,  19.9997,  18.0950,  16.1902,  14.2855,     &
                   12.3808,  10.4760,   8.5713,   6.6666,   4.7618,     &
                    2.8571,   0.9524,  -0.9524,  -2.8571,  -4.7618,     &
                   -6.6666,  -8.5713, -10.4760, -12.3808, -14.2855,     &
                  -16.1902, -18.0950, -19.9997, -21.9044, -23.8092,     &
                  -25.7139, -27.6186, -29.5234, -31.4281, -33.3328,     &
                  -35.2375, -37.1422, -39.0470, -40.9517, -42.8564,     &
                  -44.7611, -46.6658, -48.5705, -50.4752, -52.3799,     &
                  -54.2846, -56.1893, -58.0940, -59.9986, -61.9033,     &
                  -63.8079, -65.7125, -67.6171, -69.5217, -71.4262,     &
                  -73.3307, -75.2351, -77.1393, -79.0435, -80.9474,     &
                  -82.8508, -84.7532, -86.6532, -88.5420, 106*0.0/
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ALLOCATE(var_grb2d(nx_ext,ny_ext,n2dvs,n2dlvt))
  ALLOCATE(var_grb3d(nx_ext,ny_ext,nz_ext,n3dvs,n3dlvt))
  ALLOCATE(rcdata(nx_ext*ny_ext))
  ALLOCATE(var_lev3d(nz_ext,n3dvs,n3dlvt))

  CALL getgrbfname(extdopt,extdfmt,dir_extd,extdname,                   &
                   extdinit,extdfcst,'    ',                            &
                   gribfile,grbflen,gribtime,grbtlen,istatus)
  IF (istatus /=0 ) RETURN

!-----------------------------------------------------------------------
!
!  RDNMCGRB reads NMC GRIB data
!
!-----------------------------------------------------------------------
!
  gridtyp   = reanalt62grid
  mproj_grb = reanalt62proj

  n2dvars  = reanalt62nvs2d
  n2dlvtps = reanalt62nlvt2d

  DO k=1,n2dlvtps
    DO n=1,n2dvars
      var_id2d(n,k) = reanalt62var_id2d(n,k)
    END DO
    levtyp2d(k) = reanalt62levs2d(k)
  END DO

  n3dvars  = reanalt62nvs3d
  n3dlvtps = reanalt62nlvt3d

  DO m=1,n3dlvtps
    DO n=1,n3dvars
      var_id3d(n,m) = reanalt62var_id3d(n,m)
    END DO
    levtyp3d(m) = reanalt62levs3d(m)
  END DO

  CALL rdnmcgrb(nx_ext,ny_ext,nz_ext,gribfile,grbflen, gribtime,        &
                gridesc, iproj_grb, gthin,                              &
                ni_grb,nj_grb,np_grb, nk_grb,zk_grb, npeq,nit,          &
                pi_grb,pj_grb,ipole, di_grb,dj_grb,                     &
                latsw,lonsw, latne,lonne,                               &
                latrot,lonrot,angrot,                                   &
                latstr,lonstr,facstr,                                   &
                lattru1,lattru2,lontrue,                                &
                scanmode, iscan,jscan,kscan,                            &
                ires,iearth,icomp,                                      &
                jpenta,kpenta,mpenta,ispect,icoeff,                     &
                xp_grb,yp_grb, xo_grb,yo_grb,zo_grb,                    &
                rcdata,var_grb2d,var_grb3d,var_lev3d,lvldbg,iret)

  IF (iret /= 0)  THEN
    ISTATUS = -888
    GOTO 999
  END IF

  max_nr2d = 0
  DO n=1,n2dvars
    DO m=1,n2dlvtps
      max_nr2d = MAX( max_nr2d, var_nr2d(n,m) )
    END DO
  END DO

  max_nr3d = 0
  DO n=1,n3dvars
    max_nr3d = MAX( max_nr3d, var_nr3d(n,1))
  END DO

  IF ( max_nr3d == 0 ) THEN
    WRITE (6,'(a)')                                                     &
        'No 3-D variables was found in the GRIB file',                  &
        'Program stopped in GETREANALT62.'
!   STOP
    ISTATUS = -888
    GOTO 999
  END IF

  IF ( max_nr2d == 0 ) THEN
    WRITE (6,'(a)')                                                     &
        'No 2-D variables was found in the GRIB file'
  END IF

!  write (6,'(/a7,2x,6(i7))')
!    :  'Lev\\VID',(var_id3d(n,1),n=1,n3dvars)

!  DO 60 k=1,max_nr3d
!    write (6,'(i5,4x,6(i7))')
!    :  k,(var_lev3d(k,n,1),n=1,n3dvars)
  60    CONTINUE

  DO k=1,max_nr3d
    DO n=2,n3dvars
      IF ( var_lev3d(k,1,1) /= var_lev3d(k,n,1) ) THEN
        WRITE (6,'(a)')                                                 &
            'Variables were not at the same level.',                    &
            'Program stopped in GETREANALT62.'
!       STOP
        ISTATUS = -888
        GOTO 999
      END IF
    END DO
  END DO

  IF ( iproj_grb == 5 .AND. ipole == 0 ) THEN      ! Center in NP
    iproj_ext = 1
  ELSE IF ( iproj_grb == 5 .AND. ipole == 1 ) THEN  ! Center in SP
    iproj_ext = -1
  ELSE IF ( iproj_grb == 3 .AND. ipole == 0 ) THEN  ! Center in NP
    iproj_ext = 2
  ELSE IF ( iproj_grb == 3 .AND. ipole == 1 ) THEN  ! Center in SP
    iproj_ext = -2
  ELSE IF ( iproj_grb == 1 ) THEN
    iproj_ext = 3
  ELSE IF ( (iproj_grb == 0) .OR. (iproj_grb == 4) ) THEN
    iproj_ext = 4
  ELSE
    WRITE (6,'(a)')                                                     &
        'Unknown map projection. Set to non-projection.'
    iproj_ext = 0
  END IF

  scale_ext = 1.0
  latnot_ext(1) = lattru1
  latnot_ext(2) = lattru2
  trlon_ext = lontrue

  latnot_ext(1)= 0.
  trlon_ext= 0.

  dx_ext = di_grb
  dy_ext = dj_grb

  DO i=1, nx_ext
    DO j=1, ny_ext
      lat_ext(i,j)= gausslat(j)
      lon_ext(i,j)= lonsw + (i-1) * dx_ext
    END DO
  END DO

!   swap the latitude values so that +j moves to the north

  DO i=1, nx_ext
    DO j=1, ny_ext/2
      tmpval= lat_ext(i,ny_ext-j+1)
      lat_ext(i,ny_ext-j+1)= lat_ext(i,j)
      lat_ext(i,j)= tmpval
    END DO
  END DO

!   call getmapr(iproj,scale,latnot,trlon,x0,y0)
!   call setmapr(iproj_ext,scale_ext,latnot_ext,trlon_ext)
!   call lltoxy(1,1,LatSW,LonSW,x0_ext,y0_ext)

!   DO 100 i=1,nx_ext
!     x_ext(i)=x0_ext+(i-1)*dx_ext
!100   CONTINUE

!   DO 110 j=1,ny_ext
!     y_ext(j)=y0_ext+(j-1)*dy_ext
!110   CONTINUE

!   CALL xytoll(nx_ext,ny_ext,x_ext,y_ext,lat_ext,lon_ext)
!
!-----------------------------------------------------------------------
!
!  Retrieve 2-D variables
!
!-----------------------------------------------------------------------
!
  DO j=1,ny_ext
    DO i=1,nx_ext
      IF ( var_nr2d(1,1) == 0 ) THEN
        psfc_ext   (i,j) = -999.0
      ELSE
        psfc_ext   (i,j) = var_grb2d(i,j,1,1)
      END IF

      IF ( var_nr2d(2,1) == 0 ) THEN
        trn_ext    (i,j) = -999.0
      ELSE
        trn_ext    (i,j) = var_grb2d(i,j,2,1)
      END IF

      tsfc_ext   (i,j) = -999.0
      tdeep_ext  (i,j) = -999.0
      wetsfc_ext (i,j) = -999.0
      wetdp_ext  (i,j) = -999.0
      wetcanp_ext(i,j) = -999.0

    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Retrieve 3-D variables
!
!-----------------------------------------------------------------------
!
  nz1 = MIN(var_nr3d(1,1),nz_ext)

  IF ( var_lev3d(1,1,1) > var_lev3d(nz1,1,1) ) THEN  ! 1st level at sfc
    chklev = 1
    lvscan = 0
  ELSE
    chklev = -1
    lvscan = nz1+1
  END IF

  DO k=1,nz1
    kk = chklev * k + lvscan
    DO j=1,ny_ext
      DO i=1,nx_ext
        p_ext(i,j,kk) = 1.e-4 * REAL(var_lev3d(k,1,1))                 &
                         * psfc_ext(i,j)        ! Pressure
        t_ext(i,j,kk) = var_grb3d(i,j,k,1,1)    ! Temperature (K)
        qv_ext(i,j,kk)= var_grb3d(i,j,k,2,1)    ! Specific Humidity
        u_ext(i,j,kk) = var_grb3d(i,j,k,3,1)    ! u wind (m/s)
        v_ext(i,j,kk) = var_grb3d(i,j,k,4,1)    ! v wind (m/s)

!   calculate height from sigma-P values

        IF( k == 1) THEN
!          mixrat= 1.0 / ( 1./qv_ext(i,j,kk) - 1.)
          mixrat= qv_ext(i,j,kk) / ( 1. - qv_ext(i,j,kk) )
          tvirtbar= t_ext(i,j,kk) * ( 1.0 + .61* mixrat)
          hgt_ext(i,j,kk)= trn_ext(i,j) + 29.3 * tvirtbar *             &
                          LOG(psfc_ext(i,j)/p_ext(i,j,kk))
        ELSE
          qvavg= 0.5 * (qv_ext(i,j,kk) + qv_ext(i,j,kk-1))
!          mixrat= 1.0 / (1./qvavg - .1)
          mixrat= qvavg / (1. - qvavg)
          tvirtbar= 0.5 * (t_ext(i,j,kk) + t_ext(i,j,kk-1)) *           &
                          ( 1.0 + .61 * mixrat)
          hgt_ext(i,j,kk)= hgt_ext(i,j,kk-1) + 29.3 * tvirtbar *        &
                          LOG(p_ext(i,j,kk-1)/p_ext(i,j,kk))
        END IF

!        qc_ext(i,j,kk)= 0.0
!        qi_ext(i,j,kk)= 0.0
!
!        qr_ext (i,j,kk) = -999.
!        qs_ext (i,j,kk) = -999.
!        qh_ext (i,j,kk) = -999.
      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  No need to rotate winds to be relative to true north.
!  The Re-analysis grid winds are true N-S, E-W
!
!-----------------------------------------------------------------------
!
!   DO 300 k=1,nz1
!     CALL uvmptoe(nx_ext,ny_ext,u_ext(1,1,k),v_ext(1,1,k),
!  :               lon_ext,u_ext(1,1,k),v_ext(1,1,k))
!300   CONTINUE

  istatus = 1

!
!-----------------------------------------------------------------------
!
!  Reset map projection to previous values
!
!-----------------------------------------------------------------------
!
  999  CONTINUE

!   call setmapr(iproj,scale,latnot,trlon)
!   call setorig(1,x0,y0)

  DEALLOCATE(var_grb2d,var_grb3d,rcdata,var_lev3d)

  RETURN
END SUBROUTINE getreanalt62

!
!##################################################################
!##################################################################
!######                                                      ######
!######                   SUBROUTINE TV2TQ                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE tv2tq(tv,p,t,q)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  REAL :: tv       ! virtual temperature (K)
  REAL :: p        ! pressure (pa)
  REAL :: t        ! temperature (K)
  REAL :: q        ! specific humidity (kg/kg)

  REAL :: tg,qg,fac,tvg,dtvdt,tgnew
  REAL :: es, es1

  INTEGER :: it
!
!-----------------------------------------------------------------------
!
!  Function and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_es

!fpp$ expand (f_es)
!!dir$ inline always f_es
!*$*  inline routine (f_es)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  it = 0
!
! -- TG - guess at non-virtual temperature
!
  es = f_es(p,tv)
  tg = tv/( 1.0 + 0.6078 * 0.62197*es                                   &
                / ( p - es ) )
  10    CONTINUE
  it = it + 1
!
! -- QG - guess at specific humidity
!
  es = f_es(p,tg)
  qg = 0.62197*es/( p - es )
  fac = 1.+0.6078*qg
!
! -- TVG - guess at virtual temperature
!
  tvg = tg*fac
!
! -- DTVDT - d(Tv)/dT estimated over 0.1 K interval
!    from 1st term Taylor series expansion
!
  es1 = f_es(p,tg+0.1)
  dtvdt = fac + tg * 0.6078 * 0.62197*p/(p-es)**2                       &
       * ( es1 - es ) / 0.1
  tgnew = tg + (tv-tvg)/dtvdt
!
!    write(12,*) ' '
!    write(12,*) 'IT=',IT,'  TV=',TV,'  TVG=',TVG,
!    :              '  TV-TVG=',TV-TVG,'  DTVDT=',DTVDT
!
  IF (ABS(tv-tvg) < 1.0E-4) GO TO 20
  tg = tgnew
  GO TO 10

  20    CONTINUE

  t = tgnew
  es = f_es(p,t)
  q = 0.62197*es / ( p - es )

  RETURN
END SUBROUTINE tv2tq
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETNMCETA212               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getnmceta212(nx_ext,ny_ext,nz_ext,nzsoil_ext,                &
           nstyps_ext,dir_extd,extdname,extdopt,extdfmt,                &
           extdinit,extdfcst,julfname,                                  &
           iproj_ext,scale_ext,                                         &
           trlon_ext,latnot_ext,x0_ext,y0_ext,                          &
           lat_ext,lon_ext,zpsoil_ext,                                  &
           p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                      &
           qc_ext,                                                      &
           tsoil_ext,qsoil_ext,wetcanp_ext,                             &
           snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext,                   &
           t_2m_ext, qv_2m_ext, u_10m_ext, v_10m_ext, rain_ext,         &
           istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads a NMC GRIB (Grid #212, 40km) ETA file for processing by
!  ext2arps, a program that converts external files to ARPS variables
!  and format.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  02/22/1996
!
!  MODIFICATION HISTORY:
!
!  1999/08/03 (Gene Bassett)
!  Modified the deep soil moisture and temperature to be an average of
!  the three soil layers covering 0 to 100 cm in depth.
!
!  2003/10/23 (F. KONG)
!  Modified to read in from grid 212 four more near surface variables:
!  2-m temperature and humidity, and 10-m wind (u, v)
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    dir_extd      Directory name for external file
!    extdname      Prefix string of external file name
!    extdopt       Option of external data sources
!    extdfmt       Option of external data format

!    extdinit      Initialized time in mm-dd-yyyy:hh:mm:ss format
!    extdfcst      Forecast hour in HHH:MM:SS format
!    julfname      File name in yyjjjhhmm format
!
!  OUTPUT:
!
!    iproj_ext     Map projection number of external data
!    scale_ext     Scale factor of external data
!    trlon_ext     True longitude of external data (degrees E)
!    latnot_ext(2) True latitude(s) of external data (degrees N)
!    x0_ext        x coordinate of origin of external data
!    y0_ext        y coordinate of origin of external data
!    lat_ext       latitude of external data points (degrees N)
!    lon_ext       longitude of external data points (degrees E)
!    p_ext         pressure (Pascal)
!    hgt_ext       height (m)
!    t_ext         temperature (K)
!    qv_ext        specific humidity (kg/kg)
!    u_ext         u wind component (m/s)
!    v_ext         v wind component (m/s)
!    qc_ext        Cloud water mixing ratio (kg/kg)
!    qr_ext        Rain water mixing ratio (kg/kg)
!    qi_ext        Ice mixing ratio (kg/kg)
!    qs_ext        Snow mixing ratio (kg/kg)
!    qh_ext        Hail mixing ratio (kg/kg)
!
!    tsoil_ext     Soil temperature (K)
!    qsoil_ext     Soil moisture (m**3/m**3)
!    wetcanp_ext   Water content on canopy
!
!    trn_ext       External terrain (m)
!    psfc_ext      Surface pressure (Pa)
!
! F.KONG add four near surface variables
!    T_2m_ext      Temperature at 2m AGL
!    qv_2m_ext     Specific Humidity at 2m AGL
!    U_10m_ext     U at 10m AGL
!    V_10m_ext     V at 10m AGL
!    rain_ext
! end F.KONG mod
!
!    istatus       status indicator
!
!  WORK ARRAYS:
!
!    var_grb3d     Arrays to store the GRIB 3-D variables:
!                  var_grb3d(nxgrb,nygrb,nzgrb,1,1) - Temperature (K)
!                  var_grb3d(nxgrb,nygrb,nzgrb,2,1) - Specific humidity
!                                                     (kg/kg)
!                  var_grb3d(nxgrb,nygrb,nzgrb,3,1) - u wind (m/s)
!                  var_grb3d(nxgrb,nygrb,nzgrb,4,1) - v wind (m/s)
!                  var_grb3d(nxgrb,nygrb,nzgrb,5,1) - Geopotential
!                                                     height (gpm)
!                  var_grb3d(nxgrb,nygrb,nzgrb,6,1) - Pressure vertical
!                                                     velocity (Pa/s)
!                                                     (if applied)
!                  var_grb3d(nxgrb,nygrb,nzgrb,1,2) - soil temp. (K)
!                  var_grb3d(nxgrb,nygrb,nzgrb,2,2) - vol. soil moist.
!                                                     (m**3/m**3)
!
!    var_grb2d     Arrays to store the GRIB 2-D variables:
!                  var_grb2d(nxgrb,nygrb,1) - Surface pressure (Pa)
!                  var_grb2d(nxgrb,nygrb,2) - Geopotential height (gpm)
!                  var_grb2d(nxgrb,nygrb,3) - Surface temperature (K)
!                  var_grb2d(nxgrb,nygrb,4) - Plant canopy surface
!                                             water (kg/m**2)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'gribcst.inc'

  CHARACTER (LEN=*) :: dir_extd
  CHARACTER (LEN=*) :: extdname

  INTEGER :: extdopt
  INTEGER :: extdfmt

  CHARACTER (LEN=19) :: extdinit
  CHARACTER (LEN=9) :: extdfcst
  CHARACTER (LEN=9) :: julfname
!
!-----------------------------------------------------------------------
!
!  External grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj_ext
  REAL :: scale_ext,trlon_ext
  REAL :: latnot_ext(2)
  REAL :: x0_ext,y0_ext
  REAL :: dx_ext,dy_ext
!
!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx_ext,ny_ext,nz_ext
  INTEGER :: nzsoil_ext, nstyps_ext

  REAL :: lat_ext(nx_ext,ny_ext)
  REAL :: lon_ext(nx_ext,ny_ext)
  REAL :: p_ext  (nx_ext,ny_ext,nz_ext)   ! Pressure (Pascals)
  REAL :: hgt_ext(nx_ext,ny_ext,nz_ext)   ! Height (m)
  REAL :: t_ext  (nx_ext,ny_ext,nz_ext)   ! Temperature (K)
  REAL :: qv_ext (nx_ext,ny_ext,nz_ext)   ! Specific humidity (kg/kg)
  REAL :: u_ext  (nx_ext,ny_ext,nz_ext)   ! Eastward wind component
  REAL :: v_ext  (nx_ext,ny_ext,nz_ext)   ! Northward wind component
  REAL :: qc_ext (nx_ext,ny_ext,nz_ext)   ! Cloud H2O mixing ratio (kg/kg)
!  REAL :: qr_ext (nx_ext,ny_ext,nz_ext)   ! Rain  H2O mixing ratio (kg/kg)
!  REAL :: qi_ext (nx_ext,ny_ext,nz_ext)   ! Ice   mixing ratio (kg/kg)
!  REAL :: qs_ext (nx_ext,ny_ext,nz_ext)   ! Snow  mixing ratio (kg/kg)
!  REAL :: qh_ext (nx_ext,ny_ext,nz_ext)   ! Hail  mixing ratio (kg/kg)
  REAL :: zpsoil_ext(nx_ext,ny_ext,nzsoil_ext) !Soil depths (m)

  REAL :: tsoil_ext (nx_ext,ny_ext,nzsoil_ext) ! Soil temperature (K)
  REAL :: qsoil_ext (nx_ext,ny_ext,nzsoil_ext) ! Soil moisture (m**3/m**3)
  REAL :: wetcanp_ext(nx_ext,ny_ext)      ! Canopy water amount
  REAL :: snowdpth_ext(nx_ext,ny_ext)     ! Snow depth (m)

  REAL :: trn_ext    (nx_ext,ny_ext)      ! External terrain (m)
  REAL :: psfc_ext   (nx_ext,ny_ext)      ! Surface pressure (Pa)

  INTEGER soiltyp_ext (nx_ext,ny_ext)     ! Soil type

  REAL :: t_2m_ext (nx_ext,ny_ext)
  REAL :: qv_2m_ext(nx_ext,ny_ext)
  REAL :: u_10m_ext(nx_ext,ny_ext)
  REAL :: v_10m_ext(nx_ext,ny_ext)
  REAL :: rain_ext (nx_ext,ny_ext)
!
!-----------------------------------------------------------------------
!
!  Other  external variable arrays
!
!-----------------------------------------------------------------------
!
  REAL :: x_ext(nx_ext)
  REAL :: y_ext(ny_ext)

  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Work arrays for storing grib data
!
!-----------------------------------------------------------------------
!
  REAL,    ALLOCATABLE :: var_grb2d(:,:,:,:)   ! GRIB variables
  REAL,    ALLOCATABLE :: var_grb3d(:,:,:,:,:) ! GRIB 3-D variables
  INTEGER, ALLOCATABLE :: var_lev3d(:,:,:)  ! Levels (hybrid) for
                                            ! each 3-D variable
  REAL, ALLOCATABLE :: rcdata(:)            ! temporary data array

  REAL, ALLOCATABLE :: qvs_ext(:,:,:)
!
!-----------------------------------------------------------------------
!
!  Original grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj
  REAL    :: scale,trlon,x0,y0
  REAL    :: latnot(2)
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: gribfile
  CHARACTER (LEN=14)  :: gribtime
  INTEGER :: i,j,k,kk
  INTEGER :: grbflen, grbtlen

  INTEGER :: m,n,nz1,max_nr2d,max_nr3d

  REAL :: govrd
  INTEGER, PARAMETER :: nzsoilin_ext = 4
  REAL,    PARAMETER :: soildepth_ext(nzsoilin_ext+1) = (/0.0, 0.1, 0.4, 1.0, 2.0/)
           ! grid #212 contains 4 soil layers 0-1cm, 1-4cm, 4-100cm, 100-200cm.

  INTEGER :: chklev, lvscan

  LOGICAL :: fexist
  INTEGER :: igrbfmt

  REAL, ALLOCATABLE :: utmp(:,:), vtmp(:,:)
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'grid.inc'
!
!-----------------------------------------------------------------------
!
!  GRIB grid information
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=42) :: gridesc ! Grid description

  INTEGER :: iproj_grb    ! Map projection indicator
  INTEGER :: gthin        ! Indicator of whether the grid is "thinned"

  INTEGER :: ni_grb       ! Number of points along x-axis
  INTEGER :: nj_grb       ! Number of points along y-axis
  INTEGER :: np_grb       ! Total number of horizontal grid points

  INTEGER :: nk_grb       ! Number of vertical parameters
  REAL :: zk_grb(nz_ext)  ! Vertical coordinate parameters

  INTEGER :: npeq         ! Number of lat circles from pole to equator
  INTEGER :: nit(nz_ext)  ! Number of x-points for thinned grid

  REAL :: pi_grb          ! x-coordinate of pole point
  REAL :: pj_grb          ! y-coordinate of pole point
  INTEGER :: ipole        ! Projection center flag

  REAL :: di_grb          ! x-direction increment or grid length
  REAL :: dj_grb          ! y-direction increment or grid length

  REAL :: latsw           ! Latitude  of South West corner point
  REAL :: lonsw           ! Longitude of South West corner point
  REAL :: latne           ! Latitude  of North East corner point
  REAL :: lonne           ! Longitude of North East corner point

  REAL :: lattru1         ! Latitude (1st) at which projection is true
  REAL :: lattru2         ! Latitude (2nd) at which projection is true
  REAL :: lontrue         ! Longitude      at which projection is true

  REAL :: latrot          ! Latitude  of southern pole of rotation
  REAL :: lonrot          ! Longitude of southern pole of rotation
  REAL :: angrot          ! Angle of rotation

  REAL :: latstr          ! Latitude  of the pole of stretching
  REAL :: lonstr          ! Longitude of the pole of stretching
  REAL :: facstr          ! Stretching factor

  INTEGER :: scanmode     ! Scanning indicator
  INTEGER :: iscan        ! x-direction   scanning indicator
  INTEGER :: jscan        ! y-direction   scanning indicator
  INTEGER :: kscan        ! FORTRAN index scanning indicator

  INTEGER :: ires         ! Resolution direction increments indicator
  INTEGER :: iearth       ! Earth shape indicator: spherical or oblate?
  INTEGER :: icomp        ! (u,v) components decomposition indicator

  INTEGER :: jpenta       ! J-Pentagonal resolution parameter
  INTEGER :: kpenta       ! K-Pentagonal resolution parameter
  INTEGER :: mpenta       ! M-Pentagonal resolution parameter
  INTEGER :: ispect       ! Spectral representation type
  INTEGER :: icoeff       ! Spectral coefficient storage mode

  REAL :: xp_grb          ! X coordinate of sub-satellite point
  REAL :: yp_grb          ! Y coordinate of sub-satellite point
  REAL :: xo_grb          ! X coordinate of image sector origin
  REAL :: yo_grb          ! Y coordinate of image sector origin
  REAL :: zo_grb          ! Camera altitude from center of Earth
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  CALL getmapr(iproj,scale,latnot,trlon,x0,y0)

  ALLOCATE(utmp(nx_ext,ny_ext))
  ALLOCATE(vtmp(nx_ext,ny_ext))

  CALL getgrbfname(extdopt,extdfmt,dir_extd,extdname,                   &
                   extdinit,extdfcst,'    ',                            &
                   gribfile,grbflen,gribtime,grbtlen,istatus)
  IF (istatus /=0 ) RETURN

  CALL chkgrb(gribfile(1:grbflen),grbflen,igrbfmt,istatus)
  IF (istatus /=0 ) RETURN

  IF (igrbfmt == 2) THEN

    !WRITE(6,'(1x,2a,I4)') 'GRIB file to be read is ',gribfile(1:grbflen),grbflen

    CALL  getnam212_grb2(nx_ext,ny_ext,nz_ext,nzsoil_ext,nzsoilin_ext,  &
             gribfile,gribtime,grbflen,soildepth_ext,                   &
             dx_ext,dy_ext,                                             &
             iproj_ext,scale_ext,trlon_ext,latnot_ext,latsw,lonsw,      &
             p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                    &
             qc_ext,                                                    &
             tsoil_ext,qsoil_ext,wetcanp_ext,                           &
             snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext,                 &
             t_2m_ext, qv_2m_ext, u_10m_ext, v_10m_ext, rain_ext,       &
             istatus)

    IF (istatus /= 0) RETURN

  ELSE

    ALLOCATE(var_grb2d(nx_ext,ny_ext,n2dvs,n2dlvt))
    ALLOCATE(var_grb3d(nx_ext,ny_ext,nz_ext,n3dvs,n3dlvt))
    ALLOCATE(rcdata(nx_ext*ny_ext))
    ALLOCATE(var_lev3d(nz_ext,n3dvs,n3dlvt))
!
!-----------------------------------------------------------------------
!
!  RDNMCGRB reads NMC GRIB data
!
!-----------------------------------------------------------------------
!
    govrd = g/rd

    gridtyp  = eta212grid
    mproj_grb = eta212proj

    n2dvars  = eta212nvs2d
    n2dlvtps = eta212nlvt2d

    DO k=1,n2dlvtps
      DO n=1,n2dvars
        var_id2d(n,k) = eta212var_id2d(n,k)
      END DO
      levtyp2d(k) = eta212levs2d(k)
    END DO

    n3dvars  = eta212nvs3d
    n3dlvtps = eta212nlvt3d

    DO m=1,n3dlvtps
      DO n=1,n3dvars
        IF (extdfmt == 101) THEN
          var_id3d(n,m) = eta212var_id3d_sref(n,m)
        ELSE
          var_id3d(n,m) = eta212var_id3d(n,m)
        END IF
      END DO
      levtyp3d(m) = eta212levs3d(m)
    END DO

    CALL rdnmcgrb(nx_ext,ny_ext,nz_ext,gribfile,grbflen, gribtime,        &
                  gridesc, iproj_grb, gthin,                              &
                  ni_grb,nj_grb,np_grb, nk_grb,zk_grb, npeq,nit,          &
                  pi_grb,pj_grb,ipole, di_grb,dj_grb,                     &
                  latsw,lonsw, latne,lonne,                               &
                  latrot,lonrot,angrot,                                   &
                  latstr,lonstr,facstr,                                   &
                  lattru1,lattru2,lontrue,                                &
                  scanmode, iscan,jscan,kscan,                            &
                  ires,iearth,icomp,                                      &
                  jpenta,kpenta,mpenta,ispect,icoeff,                     &
                  xp_grb,yp_grb, xo_grb,yo_grb,zo_grb,                    &
                  rcdata,var_grb2d,var_grb3d,var_lev3d,lvldbg,istatus)

    IF (istatus /= 0)  THEN
      ISTATUS = -888
      GOTO 999
    END IF

    max_nr2d = 0
    DO n=1,n2dvars
      DO m=1,n2dlvtps
        max_nr2d = MAX( max_nr2d, var_nr2d(n,m) )
      END DO
    END DO

    max_nr3d = 0
    DO n=1,n3dvars
      max_nr3d = MAX( max_nr3d, var_nr3d(n,1) )
    END DO

    IF ( max_nr3d == 0 ) THEN
      WRITE (6,'(a)')                                                     &
          'No 3-D variable was found in the GRIB file',                   &
          'Dataset problem in GETNMCETA212.',                             &
          ' '
      !STOP
      ISTATUS = -888
      GOTO 999
    END IF

    IF ( max_nr2d == 0 ) THEN
      WRITE (6,'(a)')                                                     &
          'No 2-D variables was found in the GRIB file'
    END IF

  !  WRITE (6,'(/a7,2x,6(i7))') 'Lev\VID',(var_id3d(n,1),n=1,n3dvars)
  !  DO k=1,max_nr3d
  !    write (6,'(/i5,4x,6(i7))') k,(var_lev3d(k,n,1),n=1,n3dvars)
  !   END DO

    DO k=1,max_nr3d
      DO n=2,n3dvars
        IF ( var_lev3d(k,1,1) /= var_lev3d(k,n,1) ) THEN
          WRITE (6,'(a,/,a,I3,a,I7,/,a,I7,a,I3,/,a,/)')                   &
          'Variables were not at the same level.',                        &
          'The K (k= ',k,')th level for variable 1 is ',var_lev3d(k,1,1), &
          'but it is ',var_lev3d(k,n,1),' for variable n = ',n,           &
          'Dataset problem in GETNMCETA212.'
          !STOP
          ISTATUS = -888
          GOTO 999
        END IF
      END DO
    END DO

    IF ( iproj_grb == 5 .AND. ipole == 0 ) THEN      ! Center in NP
      iproj_ext = 1
    ELSE IF ( iproj_grb == 5 .AND. ipole == 1 ) THEN  ! Center in SP
      iproj_ext = -1
    ELSE IF ( iproj_grb == 3 .AND. ipole == 0 ) THEN  ! Center in NP
      iproj_ext = 2
    ELSE IF ( iproj_grb == 3 .AND. ipole == 1 ) THEN  ! Center in SP
      iproj_ext = -2
    ELSE IF ( iproj_grb == 1 ) THEN
      iproj_ext = 3
    ELSE IF ( iproj_grb == 0 ) THEN
      iproj_ext = 4
    ELSE
      WRITE (6,'(a)')                                                     &
          'Unknown map projection. Set to non-projection.'
      iproj_ext = 0
    END IF

    scale_ext = 1.0
    latnot_ext(1) = lattru1
    latnot_ext(2) = lattru2
    trlon_ext = lontrue

    dx_ext = di_grb
    dy_ext = dj_grb
!
!-----------------------------------------------------------------------
!
!  Retrieve 2-D variables
!
!-----------------------------------------------------------------------
!

    DO j=1,ny_ext
      DO i=1,nx_ext

        IF ( var_nr2d(1,1) == 0 ) THEN
          psfc_ext   (i,j) = -999.0
        ELSE
  !        psfc_ext   (i,j) = var_grb2d(i,j,1,1) * 100.0
          psfc_ext   (i,j) = var_grb2d(i,j,1,1) ! already with Pa unit (F.KONG)
        END IF
        IF ( var_nr2d(2,1) == 0 ) THEN
          trn_ext    (i,j) = -999.0
        ELSE
          trn_ext    (i,j) = var_grb2d(i,j,2,1)
        END IF

        IF ( var_nr3d(1,2) == 0 ) THEN
          DO k=1,nzsoil_ext
             tsoil_ext (i,j,k) = -999.0
             qsoil_ext (i,j,k) = -999.0
          END DO

        ELSE

          IF (soilmodel_option == 1) THEN ! Old ARPS Force-Restore Soil Model
            tsoil_ext(i,j,1) = var_grb2d(i,j,3,1)

            IF ( nint(var_grb2d(i,j,5,1)) == 1 ) THEN  ! soil temp over land
              tsoil_ext(i,j,2) = 0.1 * var_grb3d(i,j,1,1,2) & !   0-10cm
                               + 0.3 * var_grb3d(i,j,2,1,2) & !  10-40cm
                               + 0.6 * var_grb3d(i,j,3,1,2)   ! 40-100cm
              soiltyp_ext(i,j) = 0
            ELSE                                       ! sfc temp over sea
              tsoil_ext(i,j,2) = var_grb2d(i,j,3,1)
              soiltyp_ext(i,j) = 13 ! Set soil type to water
            END IF

            qsoil_ext(i,j,1) = var_grb3d(i,j,1,2,2)
            qsoil_ext(i,j,2) = 0.1 * var_grb3d(i,j,1,2,2)  & !   0-10cm
                              + 0.3 * var_grb3d(i,j,2,2,2) & !  10-40cm
                              + 0.6 * var_grb3d(i,j,3,2,2)   ! 40-100cm

          ELSE ! OU Soil model

  !EMK 15 June 2002...Reorganized, added another soil level, and added
  !water temperature.
            IF ( nint(var_grb2d(i,j,5,1)) == 1 ) THEN  ! Land
              tsoil_ext (i,j,1) = var_grb2d(i,j,3,1)   ! Ground temperature
              qsoil_ext (i,j,1) = var_grb3d(i,j,1,2,2) ! Assumed to be same
                                                       ! as first below ground
                                                       ! level.
              DO k=2,nzsoil_ext
                ! "TSOIL" in GRIB is below ground, treated as separate
                ! variable from ground temperature.
                tsoil_ext (i,j,k) = var_grb3d(i,j,k-1,1,2) ! Not a bug;
                qsoil_ext (i,j,k) = var_grb3d(i,j,k-1,2,2)
              END DO


            ELSE  ! Water
              DO k=1,nzsoil_ext
                tsoil_ext (i,j,k) = var_grb2d(i,j,7,1) ! Water temperature
                qsoil_ext (i,j,k) = 1.                 ! 100% water
              END DO
            END IF ! Land or water?

  !        qsoil_ext (i,j,1) = var_grb3d(i,j,1,2,2)  ! 0 cm
  !        DO k=2,nzsoil_ext
  !          qsoil_ext (i,j,k) = var_grb3d(i,j,k-1,2,2)
  !        END DO

  !EMK END 15 June 2002

          END IF  ! soilmodel_option

        END IF

        IF ( var_nr2d(4,1) == 0 ) THEN
          wetcanp_ext(i,j) = -999.0
        ELSE
          wetcanp_ext(i,j) = var_grb2d(i,j,4,1)*1.e-3     ! in meter
        END IF

        IF ( var_nr2d(6,1) == 0 ) THEN
          snowdpth_ext(i,j) = -999
        ELSE

  !      Convert water equiv. of accum. snow depth (kg/m**2) to meters
  !      (where 1 meter liquid water is set equivqlent to 10 meters snow).
  !          0.01 = 10. (m snow/m liquid) / (1000 kg/m**3)

          snowdpth_ext(i,j) = 0.01 * var_grb2d(i,j,6,1)  ! in meters

        END IF
        ! F.KONG add
        IF( var_nr2d(1,2) == 0 ) THEN
          t_2m_ext(i,j)= -999.0
        ELSE
          t_2m_ext(i,j)= var_grb2d(i,j,1,2)
        END IF
        IF( var_nr2d(2,2) == 0 ) THEN
          qv_2m_ext(i,j)= -999.0
        ELSE
          qv_2m_ext(i,j)= var_grb2d(i,j,2,2)
        END IF

        IF( var_nr2d(3,2) == 0 ) THEN
          u_10m_ext(i,j)= -999.0
        ELSE
          u_10m_ext(i,j)= var_grb2d(i,j,3,2)
        END IF

        IF( var_nr2d(4,2) == 0 ) THEN
          v_10m_ext(i,j)= -999.0
        ELSE
          v_10m_ext(i,j)= var_grb2d(i,j,4,2)
        END IF

        IF ( var_nr2d(8,1) == 0 ) THEN
          rain_ext(i,j)= -999.0
        ELSE
          rain_ext(i,j)= var_grb2d(i,j,8,1)
        END IF

        ! end F.KONG mod

      END DO
    END DO

!
!-----------------------------------------------------------------------
!
!  Retrieve 3-D variables
!
!-----------------------------------------------------------------------
!
    nz1 = MIN(var_nr3d(1,1),nz_ext)

    IF ( var_lev3d(1,1,1) > var_lev3d(nz1,1,1) ) THEN  ! 1st level at sfc
      chklev = 1
      lvscan = 0
    ELSE
      chklev = -1
      lvscan = nz1+1
    END IF

    DO k=1,nz1
      kk = chklev * k + lvscan
      DO j=1,ny_ext
        DO i=1,nx_ext
  !      p_ext  (i,j,kk) = psfc_ext(i,j)
  !    :                         * float(var_lev3d(k,1))/10000.
                ! Pressure derived from sigma ccordinates

          p_ext  (i,j,kk) = 100.0 * REAL(var_lev3d(k,1,1)) ! Pressure
          t_ext  (i,j,kk) = var_grb3d(i,j,k,1,1)    ! Temperature (K)
          qv_ext (i,j,kk) = var_grb3d(i,j,k,2,1)    ! Spec. humidity
                                                    ! (kg/kg)
          u_ext  (i,j,kk) = var_grb3d(i,j,k,3,1)    ! u wind (m/s)
          v_ext  (i,j,kk) = var_grb3d(i,j,k,4,1)    ! v wind (m/s)
          hgt_ext(i,j,kk) = var_grb3d(i,j,k,5,1)    ! height (m)

          qc_ext (i,j,kk) = -999.
!          qr_ext (i,j,kk) = -999.
!          qi_ext (i,j,kk) = -999.
!          qs_ext (i,j,kk) = -999.
!          qh_ext (i,j,kk) = -999.
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Calculate the height by using hytrostatical equation
!
!-----------------------------------------------------------------------
!
!   DO 220 j=1,ny_ext
!   DO 220 i=1,nx_ext
!     tema = 0.5 * ( tsfc_ext (i,j) + t_ext (i,j,1) )
!     temb = 0.5 * ( qvsfc_ext(i,j) + qv_ext(i,j,1) )
!     temc = tema * (1.0 + 0.608*temb)         ! Virtual temperature
!     tema = 0.5 * ( psfc_ext (i,j) + p_ext (i,j,1) )
!     temb = temc / tema
!
!     hgt_ext(i,j,1) = trn_ext(i,j) + govrd * temb      ! Height (m)
!  :                 * ( psfc_ext(i,j) - p_ext(i,j,1) )
!220   CONTINUE
!
!   DO 230 k=2,nz1
!   DO 230 j=1,ny_ext
!   DO 230 i=1,nx_ext
!     tema = 0.5 * ( t_ext (i,j,k-1) + t_ext (i,j,k) )
!     temb = 0.5 * ( qv_ext(i,j,k-1) + qv_ext(i,j,k) )
!     temc = tema * (1.0 + 0.608*temb)         ! Virtual temperature
!     tema = 0.5 * ( p_ext (i,j,k-1) + p_ext (i,j,k) )
!     temb = temc / tema
!
!     hgt_ext(i,j,k) = hgt_ext(i,j,k-1) + govrd * temb   ! Height (m)
!  :                 * ( p_ext (i,j,k-1) - p_ext (i,j,k) )
!230   CONTINUE
!
!-----------------------------------------------------------------------
!
! Convert relative humidity to specific humidity
!
!-----------------------------------------------------------------------

    IF (var_id3d(2,1) == 52) THEN    ! 52 is relative humidity
                                     ! 51 is specific humidity
      ALLOCATE(qvs_ext(nx_ext,ny_ext,nz_ext), STAT = istatus)
      CALL check_alloc_status(istatus, "getnmceta212:qvs_ext")

      CALL getqvs(nx_ext,ny_ext,nz_ext, 1,nx_ext,1,ny_ext,1,nz_ext,       &
                  p_ext, t_ext, qvs_ext )

      DO k=1,nz_ext
        DO j=1,ny_ext
          DO i=1,nx_ext
            qv_ext(i,j,k) = 0.01*qv_ext(i,j,k)*qvs_ext(i,j,k)
          END DO
        END DO
      END DO
      DEALLOCATE(qvs_ext)
    END IF

    DEALLOCATE(var_grb2d,var_grb3d,rcdata,var_lev3d)
  END IF

  IF (lonsw > 180) lonsw = lonsw-360.
  CALL setmapr(iproj_ext,scale_ext,latnot_ext,trlon_ext)
  CALL lltoxy(1,1,latsw,lonsw,x0_ext,y0_ext)

  DO i=1,nx_ext
    x_ext(i)=x0_ext+(i-1)*dx_ext
  END DO

  DO j=1,ny_ext
    y_ext(j)=y0_ext+(j-1)*dy_ext
  END DO

  CALL xytoll(nx_ext,ny_ext,x_ext,y_ext,lat_ext,lon_ext)

!---------------------------------------------------------------------
!
!  Define soil depths
!
!---------------------------------------------------------------------

  IF(soilmodel_option == 1 .AND. nzsoil_ext /= 2) THEN
    WRITE(6,'(2a,I3,a/2a,I2,a/,a,a)')                                   &
               'NCEP GFS grid #212 only provides 4 soil layers, However,',&
               ' You are trying to extract ',nzsoil_ext, ' layers.',    &
               ' for ARPS two-layer force-restore model ',              &
               '(soilmodel_option = ',soilmodel_option,')',             &
               ' Please check the code ext2arps.f90 for NAM212 grid.',  &
               ' Terminating ...'
    CALL arpsstop('Wrong number soil layers.',1)
  ELSE IF( soilmodel_option == 2 .AND. nzsoil_ext /= nzsoilin_ext+1) THEN
    WRITE(6,'(2a,I3,a/2a,I2,a/,a,a)')                                   &
               'NCEP GFS grid #212 only provides 4 soil layers, However,',&
               ' You are trying to extract ',nzsoil_ext, ' layers.',    &
               ' for ARPS multi-layer OUSoil model ',                   &
               '(soilmodel_option = ',soilmodel_option,')',             &
               ' Please check the code ext2arps.f90 for NAM212 grid.',  &
               ' Terminating ...'
    CALL arpsstop('Wrong number of soil layers.',1)
  END IF

  IF(soilmodel_option == 1) THEN
    DO j=1,ny_ext
      DO i=1,nx_ext
        zpsoil_ext(i,j,1) = 0.05  ! any values should work because
        zpsoil_ext(i,j,2) = 1.0   ! the ARPS system will ignore these values
      END DO
    END DO
  ELSE
    DO j=1,ny_ext
      DO i=1,nx_ext
        zpsoil_ext(i,j,1) = 0.0
        DO k=2,nzsoil_ext
          zpsoil_ext(i,j,k) = - (soildepth_ext(k)-soildepth_ext(k-1))/2 - soildepth_ext(k-1)
                                  ! The middle point for each layer
        END DO
      END DO
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Rotate winds to be relative to true north.
!  The RUC data are sent as grid-relative.
!
!-----------------------------------------------------------------------
!
  DO k=1,nz_ext
    CALL uvmptoe(nx_ext,ny_ext,u_ext(1,1,k),v_ext(1,1,k),               &
                 lon_ext,utmp,vtmp)
    u_ext(:,:,k) = utmp(:,:)
    v_ext(:,:,k) = vtmp(:,:)
  END DO

  istatus = 1

!
!-----------------------------------------------------------------------
!
!  Reset map projection to previous values
!
!-----------------------------------------------------------------------
!

  999  CONTINUE
  CALL setmapr(iproj,scale,latnot,trlon)
  CALL setorig(1,x0,y0)

  DEALLOCATE(utmp)
  DEALLOCATE(vtmp)

  RETURN
END SUBROUTINE getnmceta212

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETNMCETA212               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getnmceta132(nx_ext,ny_ext,nz_ext,nzsoil_ext,                &
           nstyps_ext,dir_extd,extdname,extdopt,extdfmt,                &
           extdinit,extdfcst,julfname,                                  &
           iboxs,iboxe,jboxs,jboxe,                                     &
           iproj_ext,scale_ext,                                         &
           trlon_ext,latnot_ext,x0_ext,y0_ext,                          &
           lat_ext,lon_ext,zpsoil_ext,                                  &
           p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                      &
           qc_ext,                                                      &
           tsoil_ext,qsoil_ext,wetcanp_ext,                             &
           snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext,                   &
           t_2m_ext, qv_2m_ext, u_10m_ext, v_10m_ext, rain_ext,         &
           istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads a NMC GRIB (Grid #132, 16km) NAM file for processing by
!  ext2arps, a program that converts external files to ARPS variables
!  and format.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  01/24/2013
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  !INCLUDE 'gribcst.inc'

  CHARACTER (LEN=*), INTENT(INOUT) :: dir_extd
  CHARACTER (LEN=*), INTENT(IN)    :: extdname

  INTEGER,           INTENT(IN)    :: extdopt
  INTEGER,           INTENT(IN)    :: extdfmt

  CHARACTER (LEN=19),INTENT(IN)    :: extdinit
  CHARACTER (LEN=9) ,INTENT(IN)    :: extdfcst
  CHARACTER (LEN=9) ,INTENT(IN)    :: julfname

  INTEGER,           INTENT(IN)    :: iboxs, iboxe, jboxs, jboxe
!
!-----------------------------------------------------------------------
!
!  External grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER, INTENT(OUT) :: iproj_ext
  REAL,    INTENT(OUT) :: scale_ext,trlon_ext
  REAL,    INTENT(OUT) :: latnot_ext(2)
  REAL,    INTENT(OUT) :: x0_ext,y0_ext
!
!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------
!
  INTEGER, INTENT(IN)  :: nx_ext,ny_ext,nz_ext
  INTEGER, INTENT(IN)  :: nzsoil_ext, nstyps_ext

  REAL,    INTENT(OUT) :: lat_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: lon_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: p_ext  (nx_ext,ny_ext,nz_ext)   ! Pressure (Pascals)
  REAL,    INTENT(OUT) :: hgt_ext(nx_ext,ny_ext,nz_ext)   ! Height (m)
  REAL,    INTENT(OUT) :: t_ext  (nx_ext,ny_ext,nz_ext)   ! Temperature (K)
  REAL,    INTENT(OUT) :: qv_ext (nx_ext,ny_ext,nz_ext)   ! Specific humidity (kg/kg)
  REAL,    INTENT(OUT) :: u_ext  (nx_ext,ny_ext,nz_ext)   ! Eastward wind component
  REAL,    INTENT(OUT) :: v_ext  (nx_ext,ny_ext,nz_ext)   ! Northward wind component
  REAL,    INTENT(OUT) :: qc_ext (nx_ext,ny_ext,nz_ext)   ! Cloud H2O mixing ratio (kg/kg)
!  REAL,    INTENT(OUT) :: qr_ext (nx_ext,ny_ext,nz_ext)   ! Rain  H2O mixing ratio (kg/kg)
!  REAL,    INTENT(OUT) :: qi_ext (nx_ext,ny_ext,nz_ext)   ! Ice   mixing ratio (kg/kg)
!  REAL,    INTENT(OUT) :: qs_ext (nx_ext,ny_ext,nz_ext)   ! Snow  mixing ratio (kg/kg)
!  REAL,    INTENT(OUT) :: qh_ext (nx_ext,ny_ext,nz_ext)   ! Hail  mixing ratio (kg/kg)
  REAL,    INTENT(OUT) :: zpsoil_ext(nx_ext,ny_ext,nzsoil_ext) !Soil depths (m)

  REAL,    INTENT(OUT) :: tsoil_ext (nx_ext,ny_ext,nzsoil_ext) ! Soil temperature (K)
  REAL,    INTENT(OUT) :: qsoil_ext (nx_ext,ny_ext,nzsoil_ext) ! Soil moisture (m**3/m**3)
  REAL,    INTENT(OUT) :: wetcanp_ext(nx_ext,ny_ext)      ! Canopy water amount
  REAL,    INTENT(OUT) :: snowdpth_ext(nx_ext,ny_ext)     ! Snow depth (m)

  REAL,    INTENT(OUT) :: trn_ext    (nx_ext,ny_ext)      ! External terrain (m)
  REAL,    INTENT(OUT) :: psfc_ext   (nx_ext,ny_ext)      ! Surface pressure (Pa)

  INTEGER, INTENT(OUT) :: soiltyp_ext (nx_ext,ny_ext)     ! Soil type

  REAL,    INTENT(OUT) :: t_2m_ext (nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: qv_2m_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: u_10m_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: v_10m_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: rain_ext (nx_ext,ny_ext)

  INTEGER, INTENT(OUT) :: istatus
!
!-----------------------------------------------------------------------
!
!  Other  external variable arrays
!
!-----------------------------------------------------------------------
!
  REAL :: dx_ext,dy_ext

  REAL :: x_ext(nx_ext)
  REAL :: y_ext(ny_ext)
!
!-----------------------------------------------------------------------
!
!  Work arrays for storing grib data
!
!-----------------------------------------------------------------------
!
  !REAL,    ALLOCATABLE :: var_grb2d(:,:,:,:)   ! GRIB variables
  !REAL,    ALLOCATABLE :: var_grb3d(:,:,:,:,:) ! GRIB 3-D variables
  !INTEGER, ALLOCATABLE :: var_lev3d(:,:,:)     ! Levels (hybrid) for
  !                                             ! each 3-D variable
  !REAL,    ALLOCATABLE :: rcdata(:)            ! temporary data array
  !
  !REAL,    ALLOCATABLE :: qvs_ext(:,:,:)
  REAL,    ALLOCATABLE :: utmp(:,:), vtmp(:,:)

!
!-----------------------------------------------------------------------
!
!  Original grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj
  REAL    :: scale,trlon,x0,y0
  REAL    :: latnot(2)
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: gribfile
  CHARACTER (LEN=14)  :: gribtime
  INTEGER :: i,j,k,kk
  INTEGER :: grbflen, grbtlen

  !INTEGER :: m,n,nz1,max_nr2d,max_nr3d
  !
  !REAL :: govrd
  INTEGER, PARAMETER :: nzsoilin_ext = 4
  REAL,    PARAMETER :: soildepth_ext(nzsoilin_ext+1) = (/0.0, 0.1, 0.4, 1.0, 3.0/)
           ! grid #212 contains 4 soil layers 0-1cm, 1-4cm, 4-100cm, 100-200cm, 300cm

  !INTEGER :: chklev, lvscan

  INTEGER :: igrbfmt
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'grid.inc'
!
!-----------------------------------------------------------------------
!
!  GRIB grid information
!
!-----------------------------------------------------------------------
!
  !CHARACTER (LEN=42) :: gridesc ! Grid description
  !
  !INTEGER :: iproj_grb    ! Map projection indicator
  !INTEGER :: gthin        ! Indicator of whether the grid is "thinned"
  !
  !INTEGER :: ni_grb       ! Number of points along x-axis
  !INTEGER :: nj_grb       ! Number of points along y-axis
  !INTEGER :: np_grb       ! Total number of horizontal grid points
  !
  !INTEGER :: nk_grb       ! Number of vertical parameters
  !REAL :: zk_grb(nz_ext)  ! Vertical coordinate parameters
  !
  !INTEGER :: npeq         ! Number of lat circles from pole to equator
  !INTEGER :: nit(nz_ext)  ! Number of x-points for thinned grid
  !
  !REAL :: pi_grb          ! x-coordinate of pole point
  !REAL :: pj_grb          ! y-coordinate of pole point
  !INTEGER :: ipole        ! Projection center flag
  !
  !REAL :: di_grb          ! x-direction increment or grid length
  !REAL :: dj_grb          ! y-direction increment or grid length
  !
  REAL :: latsw           ! Latitude  of South West corner point
  REAL :: lonsw           ! Longitude of South West corner point
  !REAL :: latne           ! Latitude  of North East corner point
  !REAL :: lonne           ! Longitude of North East corner point
  !
  !REAL :: lattru1         ! Latitude (1st) at which projection is true
  !REAL :: lattru2         ! Latitude (2nd) at which projection is true
  !REAL :: lontrue         ! Longitude      at which projection is true
  !
  !REAL :: latrot          ! Latitude  of southern pole of rotation
  !REAL :: lonrot          ! Longitude of southern pole of rotation
  !REAL :: angrot          ! Angle of rotation
  !
  !REAL :: latstr          ! Latitude  of the pole of stretching
  !REAL :: lonstr          ! Longitude of the pole of stretching
  !REAL :: facstr          ! Stretching factor
  !
  !INTEGER :: scanmode     ! Scanning indicator
  !INTEGER :: iscan        ! x-direction   scanning indicator
  !INTEGER :: jscan        ! y-direction   scanning indicator
  !INTEGER :: kscan        ! FORTRAN index scanning indicator
  !
  !INTEGER :: ires         ! Resolution direction increments indicator
  !INTEGER :: iearth       ! Earth shape indicator: spherical or oblate?
  !INTEGER :: icomp        ! (u,v) components decomposition indicator
  !
  !INTEGER :: jpenta       ! J-Pentagonal resolution parameter
  !INTEGER :: kpenta       ! K-Pentagonal resolution parameter
  !INTEGER :: mpenta       ! M-Pentagonal resolution parameter
  !INTEGER :: ispect       ! Spectral representation type
  !INTEGER :: icoeff       ! Spectral coefficient storage mode
  !
  !REAL :: xp_grb          ! X coordinate of sub-satellite point
  !REAL :: yp_grb          ! Y coordinate of sub-satellite point
  !REAL :: xo_grb          ! X coordinate of image sector origin
  !REAL :: yo_grb          ! Y coordinate of image sector origin
  !REAL :: zo_grb          ! Camera altitude from center of Earth
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  CALL getmapr(iproj,scale,latnot,trlon,x0,y0)

  CALL getgrbfname(extdopt,extdfmt,dir_extd,extdname,                   &
                   extdinit,extdfcst,'    ',                            &
                   gribfile,grbflen,gribtime,grbtlen,istatus)
  IF (istatus /=0 ) RETURN

  CALL chkgrb(gribfile(1:grbflen),grbflen,igrbfmt,istatus)
  IF (istatus /=0 ) RETURN

  IF (igrbfmt == 2) THEN

    !WRITE(6,'(1x,2a,I4)') 'GRIB file to be read is ',gribfile(1:grbflen),grbflen

    CALL  getnam132_grb2(nx_ext,ny_ext,nz_ext,nzsoil_ext,nzsoilin_ext,  &
             gribfile,gribtime,grbflen,soildepth_ext,                   &
             dx_ext,dy_ext,iboxs,iboxe,jboxs,jboxe,                     &
             iproj_ext,scale_ext,trlon_ext,latnot_ext,latsw,lonsw,      &
             lat_ext,lon_ext,p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,    &
             qc_ext,                                                    &
             tsoil_ext,qsoil_ext,wetcanp_ext,                           &
             snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext,                 &
             t_2m_ext, qv_2m_ext, u_10m_ext, v_10m_ext, rain_ext,       &
             istatus)

    IF (istatus /= 0) RETURN

  ELSE

    WRITE(6,'(1x,2a,I4)') 'GRIB file to be read is ',gribfile(1:grbflen),grbflen
    WRITE(*,'(1x,a)') 'ERROR: No support for GRIB 1 edition.'
    CALL arpsstop('ERROR: Not supported at present.',1)

  END IF

  !WHERE(lon_ext > 180) lon_ext = lon_ext-360.

  CALL setmapr(iproj_ext,scale_ext,latnot_ext,trlon_ext)
  IF (lonsw > 180) lonsw = lonsw-360.
  CALL lltoxy(1,1,latsw,lonsw,x0_ext,y0_ext)

  x0_ext = x0_ext + (iboxs-1)*dx_ext
  y0_ext = y0_ext + (jboxs-1)*dy_ext

  DO i=1,nx_ext
    x_ext(i)=x0_ext+(i-1)*dx_ext
  END DO

  DO j=1,ny_ext
    y_ext(j)=y0_ext+(j-1)*dy_ext
  END DO

  CALL xytoll(nx_ext,ny_ext,x_ext,y_ext,lat_ext,lon_ext)

  !CALL lltoxy(1,1,lat_ext(1,1),lon_ext(1,1),x0_ext,y0_ext)

!---------------------------------------------------------------------
!
!  Define soil depths
!
!---------------------------------------------------------------------

  IF(soilmodel_option == 1 .AND. nzsoil_ext /= 2) THEN
    WRITE(6,'(2a,I3,a/2a,I2,a/,a,a)')                                   &
               'NCEP GFS grid #212 only provides 4 soil layers, However,',&
               ' You are trying to extract ',nzsoil_ext, ' layers.',    &
               ' for ARPS two-layer force-restore model ',              &
               '(soilmodel_option = ',soilmodel_option,')',             &
               ' Please check the code ext2arps.f90 for NAM212 grid.',  &
               ' Terminating ...'
    CALL arpsstop('Wrong number soil layers.',1)
  ELSE IF( soilmodel_option == 2 .AND. nzsoil_ext /= nzsoilin_ext+1) THEN
    WRITE(6,'(2a,I3,a/2a,I2,a/,a,a)')                                   &
               'NCEP GFS grid #212 only provides 4 soil layers, However,',&
               ' You are trying to extract ',nzsoil_ext, ' layers.',    &
               ' for ARPS multi-layer OUSoil model ',                   &
               '(soilmodel_option = ',soilmodel_option,')',             &
               ' Please check the code ext2arps.f90 for NAM212 grid.',  &
               ' Terminating ...'
    CALL arpsstop('Wrong number of soil layers.',1)
  END IF

  IF(soilmodel_option == 1) THEN
    DO j=1,ny_ext
      DO i=1,nx_ext
        zpsoil_ext(i,j,1) = 0.05  ! any values should work because
        zpsoil_ext(i,j,2) = 1.0   ! the ARPS system will ignore these values
      END DO
    END DO
  ELSE
    DO j=1,ny_ext
      DO i=1,nx_ext
        zpsoil_ext(i,j,1) = 0.0
        DO k=2,nzsoil_ext
          zpsoil_ext(i,j,k) = - (soildepth_ext(k)-soildepth_ext(k-1))/2 - soildepth_ext(k-1)
                                  ! The middle point for each layer
        END DO
      END DO
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Rotate winds to be relative to true north.
!  The RUC data are sent as grid-relative.
!
!-----------------------------------------------------------------------
!
  ALLOCATE(utmp(nx_ext,ny_ext), STAT = istatus)
  ALLOCATE(vtmp(nx_ext,ny_ext), STAT = istatus)

  DO k=1,nz_ext
    CALL uvmptoe(nx_ext,ny_ext,u_ext(1,1,k),v_ext(1,1,k),               &
                 lon_ext,utmp,vtmp)
    u_ext(:,:,k) = utmp(:,:)
    v_ext(:,:,k) = vtmp(:,:)
  END DO

  DEALLOCATE(utmp)
  DEALLOCATE(vtmp)
!
!-----------------------------------------------------------------------
!
!  Reset map projection to previous values
!
!-----------------------------------------------------------------------
!
  istatus = 1

  999  CONTINUE
  CALL setmapr(iproj,scale,latnot,trlon)
  CALL setorig(1,x0,y0)

  RETURN
END SUBROUTINE getnmceta132
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETNMCETA218               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getnmceta218(notiles,nx_ext,ny_ext,nz_ext,nzsoil_ext,        &
           nstyps_ext,dir_extd,extdname,extdopt,extdfmt,                &
           extdinit,extdfcst,domain_tile,npr,npc,                       &
           iproj_ext,scale_ext,trlon_ext,latnot_ext,x0_ext,y0_ext,      &
           lat_ext,lon_ext,zpsoil_ext,                                  &
           p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                      &
           qc_ext,                                                      &
           tsoil_ext,qsoil_ext,wetcanp_ext,                             &
           snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext,                   &
           t_2m_ext, qv_2m_ext, u_10m_ext, v_10m_ext, rain_ext,         &
           istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads a NMC GRIB (Grid #218, 12km) ETA file for processing by
!  ext2arps, a program that converts external files to ARPS variables
!  and format. This subroutine calls getnmceta218_data for each tile.
!
!  NOTE:
!
!   The description about Eta 12km tiled output is available at
!   http://www.emc.ncep.noaa.gov/mmb/research/tiles.218.html.
!
!   The data can be downloaded from
!   ftpprd.ncep.noaa.gov/pub/data/nccf/com/nam/prod/nam.YYYYMMDD/tiles.tHHz
!
!   This subroutine assumes that data has been downloaded to local
!   disk and stored in "dir_extd"/"extdname".tCCz.awip218HH.YY
!
!   where:
!
!     o "dir_extd" is specified inside arps.input;
!     o "extdname" supposed to be "eta" or others;
!     o CC is the model cycle, extracted from "extdinit" (00,06,12,18);
!     o HH is the forecast hour, extracted from "extdfcst";
!     o YY is the tile number, determined from "domain_tile".
!
!   This subroutine also requires ASCII files of the grid lat/lon
!   for each tile, which can be got from
!   ftp://ftpprd.ncep.noaa.gov/pub/emc/mmb/mmbpll/g218tiles.latlon
!   and they should be stored to local disk inside directory
!   "dir_extd"/.
!
!   If the data is downloaded from NOMADS, it is a large file contains
!   data over all tiles. The data file name should be:
!   "dir_extd"/"extdname".YYYYMMDDHHfhh
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  07/07/2004
!
!  MODIFICATION HISTORY:
!
!  Yunheng Wang (03/14/2006)
!  Added support for one large file contains data for all tiles, which
!  is usually downloaded from NOMADS.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx_ext
!    ny_ext
!    nz_ext
!    nzsoil_ext
!    nstyp_ext
!    dir_extd      Directory name for external file
!    extdname      Prefix string of external file name
!    extdinit      Initialized time in mm-dd-yyyy:hh:mm:ss format
!    extdfcst      Forecast hour in HHH:MM:SS format
!
!    domain_tile   A 2D integer array which specifies the tiles to be read
!                  and their order.
!    npr           tile number per row of the above 2D array
!    npc           tile number per column of the above 2D array
!
!  OUTPUT:
!
!    iproj_ext     Map projection number of external data
!    scale_ext     Scale factor of external data
!    trlon_ext     True longitude of external data (degrees E)
!    latnot_ext(2) True latitude(s) of external data (degrees N)
!    x0_ext        x coordinate of origin of external data
!    y0_ext        y coordinate of origin of external data
!    lat_ext       latitude of external data points (degrees N)
!    lon_ext       longitude of external data points (degrees E)
!    p_ext         pressure (Pascal)
!    hgt_ext       height (m)
!    t_ext         temperature (K)
!    qv_ext        specific humidity (kg/kg)
!    u_ext         u wind component (m/s)
!    v_ext         v wind component (m/s)
!    qc_ext        Cloud water mixing ratio (kg/kg)
!    qr_ext        Rain water mixing ratio (kg/kg)
!    qi_ext        Ice mixing ratio (kg/kg)
!    qs_ext        Snow mixing ratio (kg/kg)
!    qh_ext        Hail mixing ratio (kg/kg)
!
!    tsoil_ext     Soil temperature (K)
!    qsoil_ext     Soil moisture (m**3/m**3)
!    wetcanp_ext   Water content on canopy
!
!    trn_ext       External terrain (m)
!    psfc_ext      Surface pressure (Pa)
!
!    T_2m_ext      Temperature at 2m AGL
!    qv_2m_ext     Specific Humidity at 2m AGL
!    U_10m_ext     U at 10m AGL
!    V_10m_ext     V at 10m AGL
!
!    istatus       status indicator
!
!  WORK ARRAYS:
!
!    var_grb3d     Arrays to store the GRIB 3-D variables:
!                  var_grb3d(nxgrb,nygrb,nzgrb,1,1) - Temperature (K)
!                  var_grb3d(nxgrb,nygrb,nzgrb,2,1) - Specific humidity
!                                                     (kg/kg)
!                  var_grb3d(nxgrb,nygrb,nzgrb,3,1) - u wind (m/s)
!                  var_grb3d(nxgrb,nygrb,nzgrb,4,1) - v wind (m/s)
!                  var_grb3d(nxgrb,nygrb,nzgrb,5,1) - Geopotential
!                                                     height (gpm)
!                  var_grb3d(nxgrb,nygrb,nzgrb,6,1) - Pressure vertical
!                                                     velocity (Pa/s)
!                                                     (if applied)
!                  var_grb3d(nxgrb,nygrb,nzgrb,1,2) - soil temp. (K)
!                  var_grb3d(nxgrb,nygrb,nzgrb,2,2) - vol. soil moist.
!                                                     (m**3/m**3)
!
!    var_grb2d     Arrays to store the GRIB 2-D variables:
!                  var_grb2d(nxgrb,nygrb,1) - Surface pressure (Pa)
!                  var_grb2d(nxgrb,nygrb,2) - Geopotential height (gpm)
!                  var_grb2d(nxgrb,nygrb,3) - Surface temperature (K)
!                  var_grb2d(nxgrb,nygrb,4) - Plant canopy surface
!                                             water (kg/m**2)
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'gribcst.inc'

!-----------------------------------------------------------------------
!
! Variable declarations
!
!-----------------------------------------------------------------------

  CHARACTER (LEN=*),  INTENT(INOUT)  :: dir_extd
  CHARACTER (LEN=*),  INTENT(IN)     :: extdname

  CHARACTER (LEN=19), INTENT(IN)  :: extdinit
  CHARACTER (LEN=9),  INTENT(IN)  :: extdfcst

  INTEGER, INTENT(IN)  :: extdopt, extdfmt

  INTEGER, INTENT(IN)  :: notiles   ! = 0, Tiled files
                                    ! = 1, not in tiles
  INTEGER, INTENT(IN)  :: nx_ext,ny_ext,nz_ext
  INTEGER, INTENT(IN)  :: nzsoil_ext, nstyps_ext

  INTEGER, INTENT(IN)  :: npr,npc
  INTEGER, INTENT(IN)  :: domain_tile(npr,npc)
!
!-----------------------------------------------------------------------
!
!  External grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER, INTENT(OUT) :: iproj_ext
  REAL,    INTENT(OUT) :: scale_ext,trlon_ext
  REAL,    INTENT(OUT) :: latnot_ext(2)
  REAL,    INTENT(OUT) :: x0_ext,y0_ext
!
!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------
!
  REAL, INTENT(OUT)  :: lat_ext(nx_ext,ny_ext)
  REAL, INTENT(OUT)  :: lon_ext(nx_ext,ny_ext)
  REAL, INTENT(OUT)  :: p_ext  (nx_ext,ny_ext,nz_ext)   ! Pressure (Pascals)
  REAL, INTENT(OUT)  :: hgt_ext(nx_ext,ny_ext,nz_ext)   ! Height (m)
  REAL, INTENT(OUT)  :: t_ext  (nx_ext,ny_ext,nz_ext)   ! Temperature (K)
  REAL, INTENT(OUT)  :: qv_ext (nx_ext,ny_ext,nz_ext)   ! Specific humidity (kg/kg)
  REAL, INTENT(OUT)  :: u_ext  (nx_ext,ny_ext,nz_ext)   ! Eastward wind component
  REAL, INTENT(OUT)  :: v_ext  (nx_ext,ny_ext,nz_ext)   ! Northward wind component
  REAL, INTENT(OUT)  :: qc_ext (nx_ext,ny_ext,nz_ext)   ! Cloud H2O mixing ratio (kg/kg)
!  REAL, INTENT(OUT)  :: qr_ext (nx_ext,ny_ext,nz_ext)   ! Rain  H2O mixing ratio (kg/kg)
!  REAL, INTENT(OUT)  :: qi_ext (nx_ext,ny_ext,nz_ext)   ! Ice   mixing ratio (kg/kg)
!  REAL, INTENT(OUT)  :: qs_ext (nx_ext,ny_ext,nz_ext)   ! Snow  mixing ratio (kg/kg)
!  REAL, INTENT(OUT)  :: qh_ext (nx_ext,ny_ext,nz_ext)   ! Hail  mixing ratio (kg/kg)

  REAL, INTENT(OUT)  :: zpsoil_ext(nx_ext,ny_ext,nzsoil_ext) !Soil depths (m)
  REAL, INTENT(OUT)  :: tsoil_ext (nx_ext,ny_ext,nzsoil_ext) ! Soil temperature (K)
  REAL, INTENT(OUT)  :: qsoil_ext (nx_ext,ny_ext,nzsoil_ext) ! Soil moisture (m**3/m**3)

  REAL, INTENT(OUT)  :: wetcanp_ext (nx_ext,ny_ext)      ! Canopy water amount
  REAL, INTENT(OUT)  :: snowdpth_ext(nx_ext,ny_ext)      ! Snow depth (m)
  REAL, INTENT(OUT)  :: trn_ext     (nx_ext,ny_ext)      ! External terrain (m)
  REAL, INTENT(OUT)  :: psfc_ext    (nx_ext,ny_ext)      ! Surface pressure (Pa)

  INTEGER, INTENT(OUT) ::  soiltyp_ext (nx_ext,ny_ext)     ! Soil type

  REAL, INTENT(OUT)  :: t_2m_ext (nx_ext,ny_ext)
  REAL, INTENT(OUT)  :: qv_2m_ext(nx_ext,ny_ext)
  REAL, INTENT(OUT)  :: u_10m_ext(nx_ext,ny_ext)
  REAL, INTENT(OUT)  :: v_10m_ext(nx_ext,ny_ext)
  REAL, INTENT(OUT)  :: rain_ext(nx_ext,ny_ext)
!
!-----------------------------------------------------------------------
!
!  Other  external variable arrays
!
!-----------------------------------------------------------------------
!
  REAL    :: dx_ext,dy_ext
  REAL    :: x_ext(nx_ext)
  REAL    :: y_ext(ny_ext)

  REAL    :: latsw, lonsw
  REAL    :: latsw_t, lonsw_t
  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: gribfile
  CHARACTER (LEN=256) :: tmpstr
  CHARACTER (LEN=256) :: latlonfile
  CHARACTER (LEN=14)  :: gribtime
  INTEGER             :: grbflen, grbtlen, lenstr

  INTEGER :: it,jt
  INTEGER :: i,j,k,kk

  INTEGER :: nxgrb,nygrb
  INTEGER :: offset_x, offset_y
  INTEGER :: inunit, dir_len

  REAL    :: alat, alon
!
!-----------------------------------------------------------------------
!
!  Work arrays for storing grib data
!
!-----------------------------------------------------------------------
!
  REAL,    ALLOCATABLE :: var_grb2d(:,:,:,:)   ! GRIB variables
  REAL,    ALLOCATABLE :: var_grb3d(:,:,:,:,:) ! GRIB 3-D variables
  INTEGER, ALLOCATABLE :: var_lev3d(:,:,:)     ! Levels (hybrid) for
                                               ! each 3-D variable
  REAL,    ALLOCATABLE :: rcdata(:)            ! temporary data array

  REAL,    ALLOCATABLE :: utmp(:,:), vtmp(:,:)

  INTEGER, PARAMETER :: nxmax = 69, nymax = 72 ! max size of the tiles
  INTEGER, PARAMETER :: nxbdy = 62, nybdy = 68 ! size of tile at boundary

  REAL,    ALLOCATABLE :: qvs_ext(:,:,:)
  REAL,    ALLOCATABLE :: lat(:,:), lon(:,:)

  INTEGER   :: igrbfmt
  LOGICAL   :: fexist

  INTEGER, PARAMETER :: nzsoilin_ext = 4
  REAL,    PARAMETER :: soildepth_ext(nzsoilin_ext+1) = (/0.0, 0.1, 0.4, 1.0, 2.0/)
           ! grid #218 contains 4 soil layers 0-1cm, 1-4cm, 4-100cm, 100-200cm.
!
!-----------------------------------------------------------------------
!
!  Original grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj
  REAL    :: scale,trlon,x0,y0
  REAL    :: latnot(2)

  INTEGER :: nxextin,nyextin
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begining of executable code .....
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (notiles == 1) THEN
    nxextin = nx_ext
    nyextin = ny_ext
  ELSE
    nxextin = nxmax
    nyextin = nymax
  END IF

  ALLOCATE(utmp(nx_ext,ny_ext), STAT = istatus)
  CALL check_alloc_status(istatus, "getnmceta218:utmp")
  ALLOCATE(vtmp(nx_ext,ny_ext), STAT = istatus)
  CALL check_alloc_status(istatus, "getnmceta218:vtmp")

!-----------------------------------------------------------------------
!
!  Save orignal grid map projection information
!
!-----------------------------------------------------------------------

  CALL getmapr(iproj,scale,latnot,trlon,x0,y0)

!---------------------------------------------------------------------
!
!  Define soil depths
!
!---------------------------------------------------------------------

  IF(soilmodel_option == 1 .AND. nzsoil_ext /= 2) THEN
    WRITE(6,'(2a,I3,a/2a,I2,a/,a,a)')                                   &
               'NCEP NAM grid #218 only provides 4 soil layers, However,',&
               ' You are trying to extract ',nzsoil_ext, ' layers.',    &
               ' for ARPS two-layer force-restore model ',              &
               '(soilmodel_option = ',soilmodel_option,')',             &
               ' Please check the code ext2arps.f90 for NAM218 grid.',  &
               ' Terminating ...'
    CALL arpsstop('Wrong number soil layers.',1)
  ELSE IF( soilmodel_option == 2 .AND. nzsoil_ext /= nzsoilin_ext+1) THEN
    WRITE(6,'(2a,I3,a/2a,I2,a/,a,a)')                                   &
               'NCEP NAM grid #218 only provides 4 soil layers, However,',&
               ' You are trying to extract ',nzsoil_ext, ' layers.',    &
               ' for ARPS multi-layer OUSoil model ',                   &
               '(soilmodel_option = ',soilmodel_option,')',             &
               ' Please check the code ext2arps.f90 for NAM218 grid.',  &
               ' Terminating ...'
    CALL arpsstop('Wrong number of soil layers.',1)
  END IF

  IF(soilmodel_option == 1) THEN
    DO j=1,ny_ext
      DO i=1,nx_ext
        zpsoil_ext(i,j,1) = 0.05  ! any values should work because
        zpsoil_ext(i,j,2) = 1.0   ! the ARPS system will ignore these values
      END DO
    END DO
  ELSE
    DO j=1,ny_ext
      DO i=1,nx_ext
        DO k=2,nzsoil_ext
          zpsoil_ext(i,j,k) = - (soildepth_ext(k)-soildepth_ext(k-1))/2 - soildepth_ext(k-1)
                                  ! The middle point for each layer
        END DO
        zpsoil_ext(i,j,1) = 0.0
      END DO
    END DO
  END IF

!-----------------------------------------------------------------------
!
! Find GRIB or GRIB2
!
!-----------------------------------------------------------------------

  tmpstr = ' '
  IF (notiles /= 1) THEN
     WRITE(tmpstr,'(a,i2.2)') ".",domain_tile(1,1)
  END IF

  CALL getgrbfname(extdopt,extdfmt,dir_extd,extdname,                   &
                   extdinit,extdfcst, tmpstr,                           &
                   gribfile,grbflen,gribtime,grbtlen,istatus)
  IF (istatus /=0 ) RETURN

  tmpstr = gribfile
  CALL chkgrb(tmpstr,grbflen,igrbfmt,istatus)
  IF (istatus /=0 ) RETURN

  IF (notiles /= 1) THEN     ! Get rid of the tile no from file name
                             ! Since we will loop over tiles later
     gribfile = tmpstr(1:grbflen-3)
     grbflen = grbflen-3
  END IF

  IF (igrbfmt == 1) THEN
    ALLOCATE(var_grb2d(nxextin,nyextin,n2dvs,n2dlvt),        STAT = istatus)
    CALL check_alloc_status(istatus, "getnmceta218:var_grb2d")
    ALLOCATE(var_grb3d(nxextin,nyextin,nz_ext,n3dvs,n3dlvt), STAT = istatus)
    CALL check_alloc_status(istatus, "getnmceta218:var_grb3d")
    ALLOCATE(rcdata   (nxextin*nyextin),                     STAT = istatus)
    CALL check_alloc_status(istatus, "getnmceta218:rcdata")
    ALLOCATE(var_lev3d(nz_ext,n3dvs,n3dlvt),                 STAT = istatus)
    CALL check_alloc_status(istatus, "getnmceta218:var_lev3d")
  END IF

!-----------------------------------------------------------------------
!
! Begining loop over each tile
! Note for non-tiled data, npc=npr=1
!
!-----------------------------------------------------------------------

  DO jt = 1,npc
    DO it = 1,npr

      IF (notiles == 1) THEN
        nxgrb = nx_ext
        nygrb = ny_ext

        tmpstr = gribfile
        lenstr = grbflen
      ELSE

        WRITE(tmpstr,'(2a,I2.2)') TRIM(gribfile),'.',domain_tile(it,jt)
        lenstr = grbflen + 3

        nxgrb = nxextin
        nygrb = nyextin

        IF ( MOD(domain_tile(it,jt), 9) == 0) nxgrb = nxbdy ! East boudnary
        IF ( domain_tile(it,jt) >= 46 )       nygrb = nybdy ! North boudnary
      END IF

      IF ( igrbfmt == 2 ) THEN
!-----------------------------------------------------------------------
!
!     Call getnmceta218_grb2
!
!-----------------------------------------------------------------------

        CALL  getnam218_grb2(notiles,nx_ext,ny_ext,nz_ext,nzsoil_ext,   &
                 nzsoilin_ext,nxgrb,nygrb,it,jt,                        &
                 tmpstr,lenstr,gribtime,                                &
                 soildepth_ext,dx_ext,dy_ext,                           &
                 iproj_ext,scale_ext,                                   &
                 trlon_ext,latnot_ext,latsw_t,lonsw_t,                  &
                 p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                &
                 qc_ext,                                                &
                 tsoil_ext,qsoil_ext,wetcanp_ext,                       &
                 snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext,             &
                 t_2m_ext, qv_2m_ext, u_10m_ext, v_10m_ext, rain_ext,   &
                 istatus)

      ELSE
!-----------------------------------------------------------------------
!
!     Call getnmceta218_data for GRIB file
!
!-----------------------------------------------------------------------

        CALL  getnmceta218_data(notiles,nx_ext,ny_ext,nz_ext,nzsoil_ext,&
                  nstyps_ext,nxgrb,nygrb,it,jt,                         &
                  tmpstr,lenstr,gribtime,dx_ext,dy_ext,                 &
                  iproj_ext,scale_ext,                                  &
                  trlon_ext,latnot_ext,latsw_t,lonsw_t,                 &
                  p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,               &
                  qc_ext,                                               &
                  tsoil_ext,qsoil_ext,wetcanp_ext,                      &
                  snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext,            &
                  t_2m_ext, qv_2m_ext, u_10m_ext, v_10m_ext, rain_ext,  &
                  var_grb2d,var_grb3d,rcdata,var_lev3d,istatus)
      END IF

      IF (istatus < 0) GOTO 999

      IF (it == 1 .AND. jt == 1) THEN
        latsw = latsw_t
        lonsw = lonsw_t
      END IF

    END DO
  END DO

  !
  ! Convert relative humidity to specific humidity
  !

  ALLOCATE(qvs_ext(nx_ext,ny_ext,nz_ext), STAT = istatus)
  CALL check_alloc_status(istatus, "getnmceta218:qvs_ext")

  CALL getqvs(nx_ext,ny_ext,nz_ext, 1,nx_ext,1,ny_ext,1,nz_ext,         &
              p_ext, t_ext, qvs_ext )

  DO k=1,nz_ext
    DO j=1,ny_ext
      DO i=1,nx_ext
        qv_ext(i,j,k) = 0.01*qv_ext(i,j,k)*qvs_ext(i,j,k)
      END DO
    END DO
  END DO

  !
  ! Convert surface 2-m RH to QV
  ! Assume 2-m pressure is the same as psfc
  !

  CALL getqvs(nx_ext,ny_ext,1, 1,nx_ext,1,ny_ext,1,1,                   &
              psfc_ext, t_2m_ext, qvs_ext )
  DO j=1,ny_ext
    DO i=1,nx_ext
      qv_2m_ext(i,j) = 0.01*qv_2m_ext(i,j)*qvs_ext(i,j,1)
    END DO
  END DO

  DEALLOCATE(qvs_ext)

  IF (igrbfmt == 1) THEN   ! Deallocate working arrays for GRIB1 only

    DEALLOCATE(var_grb2d,var_grb3d,rcdata,var_lev3d)

  END IF

!-----------------------------------------------------------------------
!
!  Set map projection and check for lat_ext, lon_ext etc.
!
!-----------------------------------------------------------------------

  CALL setmapr(iproj_ext,scale_ext,latnot_ext,trlon_ext)
  CALL lltoxy(1,1,latsw,lonsw,x0_ext,y0_ext)
  CALL setorig( 1, x0_ext, y0_ext)

  DO i=1,nx_ext
    x_ext(i)=(i-1)*dx_ext
  END DO

  DO j=1,ny_ext
    y_ext(j)=(j-1)*dy_ext
  END DO

  CALL xytoll(nx_ext,ny_ext,x_ext,y_ext,lat_ext,lon_ext)

!-----------------------------------------------------------------------
!
!  Check for data consistence.
!
!-----------------------------------------------------------------------

!    ALLOCATE(lat(nx_ext,ny_ext), STAT = iret)
!    ALLOCATE(lon(nx_ext,ny_ext), STAT = iret)

!    lat(:,:) = lat_ext(:,:)
!    lon(:,:) = lon_ext(:,:)


!-----------------------------------------------------------------------
!
!  Get grid lat/lon from external ASCII file (For GRIB1 tiled data only)
!
!  NOTE: This block is not necessary after we make sure everything
!        is working correctly above.
!
!-----------------------------------------------------------------------
!
  !  lat_ext(:,:) = 0.0
  !  lon_ext(:,:) = 0.0
  !
  !  inunit = 123
  !  DO jt = 1, npc
  !    DO it = 1, npr
  !      offset_x = (it-1)*nxmax
  !      offset_y = (jt-1)*nymax
  !
  !      WRITE(latlonfile,'(a,a,I2.2)') dir_extd(1:dir_len),               &
  !                                 'latlon.grid218.',domain_tile(it,jt)
  !      OPEN(inunit, FILE=TRIM(latlonfile),FORM='FORMATTED',STATUS='OLD')
  !
  !      DO WHILE(.TRUE.)
  !        READ(inunit,*, END=995) i,j, alat, alon
  !        lat_ext(offset_x+i,offset_y+j) = alat
  !        lon_ext(offset_x+i,offset_y+j) = -1*alon
  !      END DO
  !      995 CONTINUE
  !
  !      CLOSE(inunit)
  !    END DO
  !  END DO
  !
  !  IF( ANY(lat_ext <= 12.0)  .OR. ANY(lat_ext >= 62.0) .OR.              &
  !      ANY(lon_ext >= -49.0) .OR. ANY(lon_ext <= -160.0) ) THEN
  !    WRITE(6,*) 'Some lat_ext or lon_ext have wrong value. ',            &
  !               'Please check the code.'
  !    STOP
  !  END IF

!-----------------------------------------------------------------------
!
! Check grid latitude and longitude
!
!-----------------------------------------------------------------------
!
  !  DO j = 1, ny_ext
  !    DO i = 1, nx_ext
  !
  !      IF ( ABS(lat(i,j)-lat_ext(i,j)) > 0.01 .OR.                       &
  !           ABS(lon(i,j)-lon_ext(i,j)) > 0.01 )  THEN
  !        WRITE(6,'(2(a,I3),/,2(a,2F7.2,/))')                             &
  !           'Unmatched grid at i = ',i,' j = ',j,                        &
  !           'Expecting lat,lon = ',lat_ext(i,j), lon_ext(i,j),           &
  !           'Found     lat,lon = ',lat(i,j),lon(i,j)
  !        STOP
  !      END IF
  !
  !    END DO
  !  END DO

!  DEALLOCATE(lat,lon)

!-----------------------------------------------------------------------
!
!  Rotate winds to be relative to true north.
!  The Eta data are sent as grid-relative.
!
!-----------------------------------------------------------------------
!
  DO k=1,nz_ext
    CALL uvmptoe(nx_ext,ny_ext,u_ext(1,1,k),v_ext(1,1,k),               &
                 lon_ext,utmp,vtmp)
    u_ext(:,:,k) = utmp(:,:)
    v_ext(:,:,k) = vtmp(:,:)
  END DO

  istatus = 1

!
!-----------------------------------------------------------------------
!
!  Reset map projection to previous values
!
!-----------------------------------------------------------------------
!

  999  CONTINUE

  CALL setmapr(iproj,scale,latnot,trlon)
  CALL setorig(1,x0,y0)

  DEALLOCATE(utmp, vtmp)

  RETURN
END SUBROUTINE getnmceta218


SUBROUTINE getnmceta218_data(notiles,nx_ext,ny_ext,nz_ext,nzsoil_ext,   &
           nstyps_ext,nxgrb,nygrb,tile_x,tile_y,                        &
           gribfile,grbflen,gribtime,dx_ext,dy_ext,                     &
           iproj_ext,scale_ext,                                         &
           trlon_ext,latnot_ext,latsw_t,lonsw_t,                        &
           p_ext,hgt_ext,t_ext,rh_ext,u_ext,v_ext,                      &
           qc_ext,                                                      &
           tsoil_ext,qsoil_ext,wetcanp_ext,                             &
           snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext,                   &
           t_2m_ext, rh_2m_ext, u_10m_ext, v_10m_ext, rain_ext,         &
           var_grb2d,var_grb3d,rcdata,var_lev3d,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads a NMC GRIB (Grid #218, 12km) ETA file for processing by
!  EXT2ARPS. This subroutine reads one tile data and insert the data
!  into the global domain (maybe multiple tiles or one tile) arrays.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  07/06/2004
!  Based on subroutine getnmceta212.
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx_ext        X dimension size of the global domain (maybe mupltiple tiles)
!    ny_ext        Y dimension size of the global domain (maybe mupltiple tiles)
!    nz_ext        Z dimension size (39 pressure levels)
!    nzsoil_ext    5 layers (Surface, 0-10cm, 10-40cm, 40-100cm, 100-200cm)
!    nstyps_ext    =1
!    nxgrb         Tile size in X direction (69/62)
!    nygrb         Tile size in Y direction (72/68)
!    tile_x        This tile location in the multiple tile domain in X direction
!    tile_y        This tile location in the multiple tile domain in Y direction
!    dir_extd      Directory name for external file
!    extdname      External file name
!    gribtime      Data valid time, YYYYMMDDHHfHH
!
!
!  OUTPUT:
!
!    iproj_ext     Map projection number of external data
!    scale_ext     Scale factor of external data
!    trlon_ext     True longitude of external data (degrees E)
!    latnot_ext(2) True latitude(s) of external data (degrees N)
!    latsw_t       latitude of the southwest corner of this tile (degrees N)
!    lonsw_t       longitude of the southwest corner of this tile (degrees E)
!    p_ext         pressure (Pascal)
!    hgt_ext       height (m)
!    t_ext         temperature (K)
!    rh_ext        relative humidity (%)
!    u_ext         u wind component (m/s)
!    v_ext         v wind component (m/s)
!    qc_ext        Cloud water mixing ratio (kg/kg)
!    qr_ext        Rain water mixing ratio (kg/kg)
!    qi_ext        Ice mixing ratio (kg/kg)
!    qs_ext        Snow mixing ratio (kg/kg)
!    qh_ext        Hail mixing ratio (kg/kg)
!
!    tsoil_ext     Soil temperature (K)
!    qsoil_ext     Soil moisture (m**3/m**3)
!    wetcanp_ext   Water content on canopy
!
!    trn_ext       External terrain (m)
!    psfc_ext      Surface pressure (Pa)
!
!    T_2m_ext      Temperature at 2m AGL
!    Rh_2m_ext     Relative Humidity at 2m AGL
!    U_10m_ext     U at 10m AGL
!    V_10m_ext     V at 10m AGL
!
!    istatus       status indicator
!
!  WORK ARRAYS:
!
!    var_grb3d     Arrays to store the GRIB 3-D variables:
!                  var_grb3d(nxgrb,nygrb,nzgrb,1,1) - Temperature (K)
!                  var_grb3d(nxgrb,nygrb,nzgrb,2,1) - Specific humidity
!                                                     (kg/kg)
!                  var_grb3d(nxgrb,nygrb,nzgrb,3,1) - u wind (m/s)
!                  var_grb3d(nxgrb,nygrb,nzgrb,4,1) - v wind (m/s)
!                  var_grb3d(nxgrb,nygrb,nzgrb,5,1) - Geopotential
!                                                     height (gpm)
!                  var_grb3d(nxgrb,nygrb,nzgrb,6,1) - Pressure vertical
!                                                     velocity (Pa/s)
!                                                     (if applied)
!                  var_grb3d(nxgrb,nygrb,nzgrb,1,2) - soil temp. (K)
!                  var_grb3d(nxgrb,nygrb,nzgrb,2,2) - vol. soil moist.
!                                                     (m**3/m**3)
!
!    var_grb2d     Arrays to store the GRIB 2-D variables:
!                  var_grb2d(nxgrb,nygrb,1,1) - Surface pressure (Pa)
!                  var_grb2d(nxgrb,nygrb,2,1) - Geopotential height (gpm)
!                  var_grb2d(nxgrb,nygrb,3,1) - Surface temperature (K)
!                  var_grb2d(nxgrb,nygrb,4,1) - Plant canopy surface
!                                               water (kg/m**2)
!                  var_grb2d(nxgrb,nygrb,1,2) - 2-m temperature (Pa)
!                  var_grb2d(nxgrb,nygrb,2,2) - 2-m specific humidity (kg/kg)
!                  var_grb2d(nxgrb,nygrb,3,2) - 10-m U velocity (m/s)
!                  var_grb2d(nxgrb,nygrb,4,2) - 10-m V velocity (m/s)
!
!   rcdata
!
!   var_lev3d
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'gribcst.inc'

  CHARACTER(LEN=256), INTENT(IN) :: gribfile
  INTEGER,            INTENT(IN) :: grbflen
  CHARACTER(LEN=*),   INTENT(IN) :: gribtime

!-----------------------------------------------------------------------
!
!  External grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER, INTENT(IN) :: notiles
  INTEGER, INTENT(IN) :: nx_ext,ny_ext     ! multiple tiles size
  INTEGER, INTENT(IN) :: nxgrb, nygrb      ! one tile size
  INTEGER, INTENT(IN) :: nz_ext
  INTEGER, INTENT(IN) :: nzsoil_ext, nstyps_ext

  INTEGER, INTENT(IN) :: tile_x, tile_y

  INTEGER, INTENT(OUT) :: iproj_ext
  REAL,    INTENT(OUT) :: scale_ext,trlon_ext
  REAL,    INTENT(OUT) :: latnot_ext(2)
  REAL,    INTENT(OUT) :: latsw_t
  REAL,    INTENT(OUT) :: lonsw_t
  REAL,    INTENT(OUT) :: dx_ext,dy_ext
!
!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------
!
  REAL, INTENT(OUT) :: p_ext  (nx_ext,ny_ext,nz_ext)   ! Pressure (Pascals)
  REAL, INTENT(OUT) :: hgt_ext(nx_ext,ny_ext,nz_ext)   ! Height (m)
  REAL, INTENT(OUT) :: t_ext  (nx_ext,ny_ext,nz_ext)   ! Temperature (K)
  REAL, INTENT(OUT) :: rh_ext (nx_ext,ny_ext,nz_ext)   ! Relative humidity (%)
  REAL, INTENT(OUT) :: u_ext  (nx_ext,ny_ext,nz_ext)   ! Eastward wind component
  REAL, INTENT(OUT) :: v_ext  (nx_ext,ny_ext,nz_ext)   ! Northward wind component
  REAL, INTENT(OUT) :: qc_ext (nx_ext,ny_ext,nz_ext)   ! Cloud H2O mixing ratio (kg/kg)
!  REAL, INTENT(OUT) :: qr_ext (nx_ext,ny_ext,nz_ext)   ! Rain  H2O mixing ratio (kg/kg)
!  REAL, INTENT(OUT) :: qi_ext (nx_ext,ny_ext,nz_ext)   ! Ice   mixing ratio (kg/kg)
!  REAL, INTENT(OUT) :: qs_ext (nx_ext,ny_ext,nz_ext)   ! Snow  mixing ratio (kg/kg)
!  REAL, INTENT(OUT) :: qh_ext (nx_ext,ny_ext,nz_ext)   ! Hail  mixing ratio (kg/kg)

  REAL, INTENT(OUT) :: tsoil_ext (nx_ext,ny_ext,nzsoil_ext) ! Soil temperature (K)
  REAL, INTENT(OUT) :: qsoil_ext (nx_ext,ny_ext,nzsoil_ext) ! Soil moisture (m**3/m**3)

  REAL, INTENT(OUT) :: wetcanp_ext (nx_ext,ny_ext)      ! Canopy water amount
  REAL, INTENT(OUT) :: snowdpth_ext(nx_ext,ny_ext)      ! Snow depth (m)

  REAL, INTENT(OUT) :: trn_ext     (nx_ext,ny_ext)      ! External terrain (m)
  REAL, INTENT(OUT) :: psfc_ext    (nx_ext,ny_ext)      ! Surface pressure (Pa)

  INTEGER, INTENT(OUT) :: soiltyp_ext (nx_ext,ny_ext)   ! Soil type

  REAL, INTENT(OUT) :: t_2m_ext  (nx_ext,ny_ext)
  REAL, INTENT(OUT) :: rh_2m_ext (nx_ext,ny_ext)
  REAL, INTENT(OUT) :: u_10m_ext (nx_ext,ny_ext)
  REAL, INTENT(OUT) :: v_10m_ext (nx_ext,ny_ext)
  REAL, INTENT(OUT) :: rain_ext (nx_ext,ny_ext)
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'grid.inc'
!
!-----------------------------------------------------------------------
!
!  Other  external variable arrays
!
!-----------------------------------------------------------------------
!
  INTEGER, INTENT(OUT) :: istatus
!
!-----------------------------------------------------------------------
!
!  Work arrays for storing grib data
!
!-----------------------------------------------------------------------
!
  REAL,    INTENT(OUT) :: var_grb2d(nxgrb,nygrb,n2dvs,n2dlvt)   ! GRIB variables
  REAL,    INTENT(OUT) :: var_grb3d(nxgrb,nygrb,nz_ext,n3dvs,n3dlvt)
                                                             ! GRIB 3-D variables
  INTEGER, INTENT(OUT) :: var_lev3d(nz_ext,n3dvs,n3dlvt)     ! Levels (hybrid) for
                                                             ! each 3-D variable
  REAL,    INTENT(OUT) :: rcdata(nxgrb*nygrb)                ! temporary data array
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!

  INTEGER :: offset_x, offset_y

  INTEGER :: i,j,k,kk

  INTEGER :: m,n,nz1,max_nr2d,max_nr3d

  INTEGER :: chklev, lvscan
!
!-----------------------------------------------------------------------
!
!  GRIB grid information
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=42) :: gridesc ! Grid description

  INTEGER :: iproj_grb    ! Map projection indicator
  INTEGER :: gthin        ! Indicator of whether the grid is "thinned"

  INTEGER :: ni_grb       ! Number of points along x-axis
  INTEGER :: nj_grb       ! Number of points along y-axis
  INTEGER :: np_grb       ! Total number of horizontal grid points

  INTEGER :: nk_grb       ! Number of vertical parameters
  REAL :: zk_grb(nz_ext)  ! Vertical coordinate parameters

  INTEGER :: npeq         ! Number of lat circles from pole to equator
  INTEGER :: nit(nz_ext)  ! Number of x-points for thinned grid

  REAL :: pi_grb          ! x-coordinate of pole point
  REAL :: pj_grb          ! y-coordinate of pole point
  INTEGER :: ipole        ! Projection center flag

  REAL :: di_grb          ! x-direction increment or grid length
  REAL :: dj_grb          ! y-direction increment or grid length

  REAL :: latsw           ! Latitude  of South West corner point
  REAL :: lonsw           ! Longitude of South West corner point
  REAL :: latne           ! Latitude  of North East corner point
  REAL :: lonne           ! Longitude of North East corner point

  REAL :: lattru1         ! Latitude (1st) at which projection is true
  REAL :: lattru2         ! Latitude (2nd) at which projection is true
  REAL :: lontrue         ! Longitude      at which projection is true

  REAL :: latrot          ! Latitude  of southern pole of rotation
  REAL :: lonrot          ! Longitude of southern pole of rotation
  REAL :: angrot          ! Angle of rotation

  REAL :: latstr          ! Latitude  of the pole of stretching
  REAL :: lonstr          ! Longitude of the pole of stretching
  REAL :: facstr          ! Stretching factor

  INTEGER :: scanmode     ! Scanning indicator
  INTEGER :: iscan        ! x-direction   scanning indicator
  INTEGER :: jscan        ! y-direction   scanning indicator
  INTEGER :: kscan        ! FORTRAN index scanning indicator

  INTEGER :: ires         ! Resolution direction increments indicator
  INTEGER :: iearth       ! Earth shape indicator: spherical or oblate?
  INTEGER :: icomp        ! (u,v) components decomposition indicator

  INTEGER :: jpenta       ! J-Pentagonal resolution parameter
  INTEGER :: kpenta       ! K-Pentagonal resolution parameter
  INTEGER :: mpenta       ! M-Pentagonal resolution parameter
  INTEGER :: ispect       ! Spectral representation type
  INTEGER :: icoeff       ! Spectral coefficient storage mode

  REAL :: xp_grb          ! X coordinate of sub-satellite point
  REAL :: yp_grb          ! Y coordinate of sub-satellite point
  REAL :: xo_grb          ! X coordinate of image sector origin
  REAL :: yo_grb          ! Y coordinate of image sector origin
  REAL :: zo_grb          ! Camera altitude from center of Earth


  LOGICAL, PARAMETER :: verbose = .FALSE.
  INTEGER, PARAMETER :: nxmax = 69, nymax = 72
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
!  RDNMCGRB reads NMC GRIB data
!
!-----------------------------------------------------------------------
!
  IF (notiles == 1) THEN
    gridtyp  = eta218grid
  ELSE
    gridtyp = 255
  END IF
  mproj_grb = eta218proj

  n2dvars  = eta218nvs2d
  n2dlvtps = eta218nlvt2d

  DO k=1,n2dlvtps
    DO n=1,n2dvars
      var_id2d(n,k) = eta218var_id2d(n,k)
    END DO
    levtyp2d(k) = eta218levs2d(k)
  END DO

  n3dvars  = eta218nvs3d
  n3dlvtps = eta218nlvt3d

  DO m=1,n3dlvtps
    DO n=1,n3dvars
      var_id3d(n,m) = eta218var_id3d(n,m)
    END DO
    levtyp3d(m) = eta218levs3d(m)
  END DO

  CALL rdnmcgrb(nxgrb,nygrb,nz_ext,gribfile,grbflen, gribtime,          &
                gridesc, iproj_grb, gthin,                              &
                ni_grb,nj_grb,np_grb, nk_grb,zk_grb, npeq,nit,          &
                pi_grb,pj_grb,ipole, di_grb,dj_grb,                     &
                latsw,lonsw, latne,lonne,                               &
                latrot,lonrot,angrot,                                   &
                latstr,lonstr,facstr,                                   &
                lattru1,lattru2,lontrue,                                &
                scanmode, iscan,jscan,kscan,                            &
                ires,iearth,icomp,                                      &
                jpenta,kpenta,mpenta,ispect,icoeff,                     &
                xp_grb,yp_grb, xo_grb,yo_grb,zo_grb,                    &
                rcdata,var_grb2d,var_grb3d,var_lev3d,lvldbg,istatus)

  IF (istatus /= 0)  RETURN

  max_nr2d = 0
  DO n=1,n2dvars
    DO m=1,n2dlvtps
      max_nr2d = MAX( max_nr2d, var_nr2d(n,m) )
    END DO
  END DO

  max_nr3d = 0
  DO n=1,n3dvars
    max_nr3d = MAX( max_nr3d, var_nr3d(n,1) )
  END DO

  IF ( max_nr3d == 0 ) THEN
    WRITE (6,'(a)')                                                     &
        'No 3-D variable was found in the GRIB file ',TRIM(gribfile),   &
        ' Dataset problem in GETNMCETA218.',                            &
        ' '
    istatus = -888
    RETURN
  END IF

  IF ( max_nr2d == 0 ) THEN
    WRITE (6,'(a)')                                                     &
        'No 2-D variables was found in the GRIB file'
  END IF

  IF (verbose) THEN
    WRITE (6,'(/a7,2x,6(i7))')  'Lev\\VID',(var_id3d(n,1),n=1,n3dvars)
    DO k=1,max_nr3d
      WRITE (6,'(/i5,4x,6(i7))') k,(var_lev3d(k,n,1),n=1,n3dvars)
    END DO
  END IF

  DO k=1,max_nr3d
    DO n=2,n3dvars
      IF ( var_lev3d(k,1,1) /= var_lev3d(k,n,1) ) THEN
        WRITE (6,'(a)')                                                 &
            'Variables were not at the same level.',                    &
            'Dataset problem in GETNMCETA218.',                         &
            ' '
        istatus = -888
        RETURN
      END IF
    END DO
  END DO

  IF ( iproj_grb == 5 .AND. ipole == 0 ) THEN       ! Center in NP
    iproj_ext = 1
  ELSE IF ( iproj_grb == 5 .AND. ipole == 1 ) THEN  ! Center in SP
    iproj_ext = -1
  ELSE IF ( iproj_grb == 3 .AND. ipole == 0 ) THEN  ! Center in NP
    iproj_ext = 2
  ELSE IF ( iproj_grb == 3 .AND. ipole == 1 ) THEN  ! Center in SP
    iproj_ext = -2
  ELSE IF ( iproj_grb == 1 ) THEN
    iproj_ext = 3
  ELSE IF ( iproj_grb == 0 ) THEN
    iproj_ext = 4
  ELSE
    WRITE (6,'(a)') 'Unknown map projection. Set to non-projection.'
    iproj_ext = 0
  END IF

  scale_ext = 1.0
  latnot_ext(1) = lattru1
  latnot_ext(2) = lattru2
  trlon_ext     = lontrue

!
!  Do some check here, such as nxgrd, nygrd.
!
  IF (iproj_grb /= mproj_grb) THEN
    WRITE(6,*) 'Map projection unmatch, expecting ',mproj_grb,' Found ',iproj_ext
!   STOP
    ISTATUS = -888
    RETURN
  END IF

  IF(ni_grb /= nxgrb .OR. nj_grb /= nygrb) THEN
    WRITE(6,'(2(a,I4),/,2(a,I4))')                                      &
      'Data array size inconsistent, expecting  nxgrb = ',nxgrb,        &
      '  nygrb = ',nygrb,                                               &
      '                              Found:    ni_grb = ',ni_grb,       &
      ' nj_grb = ',nj_grb
    istatus = -888
    RETURN
  END IF

  dx_ext = di_grb
  dy_ext = dj_grb

  latsw_t  = latsw
  lonsw_t  = lonsw

  offset_x = (tile_x-1)*nxmax       ! Offset of this tile in the multiple
  offset_y = (tile_y-1)*nymax       ! tile domain.
!
!-----------------------------------------------------------------------
!
!  Retrieve 2-D variables
!
!-----------------------------------------------------------------------
!

  DO j=1,nygrb
    DO i=1,nxgrb

      IF ( var_nr2d(1,1) == 0 ) THEN
        psfc_ext   (offset_x+i,offset_y+j) = -999.0
      ELSE
        psfc_ext   (offset_x+i,offset_y+j) = var_grb2d(i,j,1,1)
      END IF
      IF ( var_nr2d(2,1) == 0 ) THEN
        trn_ext    (offset_x+i,offset_y+j) = -999.0
      ELSE
        trn_ext    (offset_x+i,offset_y+j) = var_grb2d(i,j,2,1)
      END IF

      IF ( var_nr2d(4,1) == 0 ) THEN
        wetcanp_ext(offset_x+i,offset_y+j) = -999.0
      ELSE
        wetcanp_ext(offset_x+i,offset_y+j) = var_grb2d(i,j,4,1)*1.e-3     ! in meter
      END IF

      IF ( var_nr2d(6,1) == 0 ) THEN
        snowdpth_ext(offset_x+i,offset_y+j) = -999
      ELSE
!      Convert water equiv. of accum. snow depth (kg/m**2) to meters
!      (where 1 meter liquid water is set equivqlent to 10 meters snow).
!          0.01 = 10. (m snow/m liquid) / (1000 kg/m**3)

        snowdpth_ext(offset_x+i,offset_y+j) = 0.01 * var_grb2d(i,j,6,1)  ! in meters
      END IF

      IF ( var_nr2d(7,1) == 0 ) THEN
        rain_ext(offset_x+i,offset_y+j)= -999.0
      ELSE
        rain_ext(offset_x+i,offset_y+j)= var_grb2d(i,j,7,1)
      END IF

      IF( var_nr2d(1,2) == 0 ) THEN
        t_2m_ext(offset_x+i,offset_y+j)= -999.0
      ELSE
        t_2m_ext(offset_x+i,offset_y+j)= var_grb2d(i,j,1,2)     ! 2-m temperature
      END IF
      IF( var_nr2d(2,2) == 0 ) THEN
        rh_2m_ext(offset_x+i,offset_y+j)= -999.0
      ELSE
        rh_2m_ext(offset_x+i,offset_y+j)= var_grb2d(i,j,2,2)    ! 2-m relative humidity
      END IF

      IF( var_nr2d(3,2) == 0 ) THEN
        u_10m_ext(offset_x+i,offset_y+j)= -999.0
      ELSE
        u_10m_ext(offset_x+i,offset_y+j)= var_grb2d(i,j,3,2)    ! 10-m U
      END IF

      IF( var_nr2d(4,2) == 0 ) THEN
        v_10m_ext(offset_x+i,offset_y+j)= -999.0
      ELSE
        v_10m_ext(offset_x+i,offset_y+j)= var_grb2d(i,j,4,2)    ! 10-m V
      END IF

    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Retrieve 3-D variables
!
!-----------------------------------------------------------------------
!
  nz1 = MIN(var_nr3d(1,1),nz_ext)

  IF ( var_lev3d(1,1,1) > var_lev3d(nz1,1,1) ) THEN  ! 1st level at sfc
    chklev = 1
    lvscan = 0
  ELSE
    chklev = -1
    lvscan = nz1+1
  END IF

  DO k=1,nz1
    kk = chklev * k + lvscan
    DO j=1,nygrb
      DO i=1,nxgrb
        p_ext  (offset_x+i,offset_y+j,kk) = 100.0 * REAL(var_lev3d(k,1,1))
                                                  ! Pressure
        t_ext  (offset_x+i,offset_y+j,kk) = var_grb3d(i,j,k,1,1)
                                                  ! Temperature (K)
        rh_ext (offset_x+i,offset_y+j,kk) = var_grb3d(i,j,k,2,1)
                                                  ! Relative humidity (%)
        u_ext  (offset_x+i,offset_y+j,kk) = var_grb3d(i,j,k,3,1)    ! u wind (m/s)
        v_ext  (offset_x+i,offset_y+j,kk) = var_grb3d(i,j,k,4,1)    ! v wind (m/s)
        hgt_ext(offset_x+i,offset_y+j,kk) = var_grb3d(i,j,k,5,1)    ! height (m)

        qc_ext (offset_x+i,offset_y+j,kk) = -999.
!        qr_ext (offset_x+i,offset_y+j,kk) = -999.
!        qi_ext (offset_x+i,offset_y+j,kk) = -999.
!        qs_ext (offset_x+i,offset_y+j,kk) = -999.
!        qh_ext (offset_x+i,offset_y+j,kk) = -999.
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Retrive land surface variables
!
!-----------------------------------------------------------------------

  DO j = 1,nygrb
    DO i = 1, nxgrb

      IF ( var_nr3d(1,2) == 0 ) THEN
        DO k=1,nzsoil_ext
           tsoil_ext (offset_x+i,offset_y+j,k) = -999.0
           qsoil_ext (offset_x+i,offset_y+j,k) = -999.0
        END DO

      ELSE

        IF (soilmodel_option == 1) THEN   ! Old ARPS Force-Restore Soil Model
          tsoil_ext(offset_x+i,offset_y+j,1) = var_grb2d(i,j,3,1)

          IF ( nint(var_grb2d(i,j,5,1)) == 1 ) THEN  ! Over land
            tsoil_ext(offset_x+i,offset_y+j,2) =          &
                               0.1 * var_grb3d(i,j,1,1,2) & !   0-10cm
                             + 0.3 * var_grb3d(i,j,2,1,2) & !  10-40cm
                             + 0.6 * var_grb3d(i,j,3,1,2)   ! 40-100cm
            soiltyp_ext(offset_x+i,offset_y+j) = 0
          ELSE                                       ! Over sea
            tsoil_ext(offset_x+i,offset_y+j,2) = var_grb2d(i,j,3,1)
            soiltyp_ext(offset_x+i,offset_y+j) = 13  ! Set soil type to water
          END IF

          qsoil_ext(offset_x+i,offset_y+j,1) = var_grb3d(i,j,1,2,2)
          qsoil_ext(offset_x+i,offset_y+j,2) = 0.1 * var_grb3d(i,j,1,2,2)  & !   0-10cm
                           + 0.3 * var_grb3d(i,j,2,2,2)  & !  10-40cm
                           + 0.6 * var_grb3d(i,j,3,2,2)    ! 40-100cm

        ELSE    ! OU Soil model

          IF ( nint(var_grb2d(i,j,5,1)) == 1 ) THEN  ! Over land
            tsoil_ext (offset_x+i,offset_y+j,1) = var_grb2d(i,j,3,1)   ! Ground temperature
            qsoil_ext (offset_x+i,offset_y+j,1) = var_grb3d(i,j,1,2,2)
                                    ! Assumed to be same as first below ground level.
            DO k=2,nzsoil_ext
              ! "TSOIL" in GRIB is below ground, treated as separate
              ! variable from ground temperature.
              tsoil_ext (offset_x+i,offset_y+j,k) = var_grb3d(i,j,k-1,1,2) ! Not a bug;
              qsoil_ext (offset_x+i,offset_y+j,k) = var_grb3d(i,j,k-1,2,2)
            END DO
          ELSE                                       ! Over water
            DO k=1,nzsoil_ext
              tsoil_ext (offset_x+i,offset_y+j,k) = var_grb2d(i,j,3,1) ! Ground temperature
              qsoil_ext (offset_x+i,offset_y+j,k) = 1.                 ! 100% water
            END DO
          END IF

       END IF  ! soilmodel_option

      END IF
    END DO
  END DO

  istatus = 1

  RETURN
END SUBROUTINE getnmceta218_data
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE GETNMCRUCN236               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getnmcrucn236(nx_ext,ny_ext,nz_ext,                          &
           dir_extd,extdname,extdopt,extdfmt,                           &
           extdinit,extdfcst,julfname,                                  &
           iproj_ext,scale_ext,                                         &
           trlon_ext,latnot_ext,x0_ext,y0_ext,                          &
           lat_ext,lon_ext,                                             &
           p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                      &
           qc_ext,qr_ext,qi_ext,qs_ext,qh_ext,                          &
           tsfc_ext,tdeep_ext,wetsfc_ext,wetdp_ext,wetcanp_ext,         &
           trn_ext,psfc_ext,snowdpth_ext,t_2m_ext,qv_2m_ext,            &
           u_10m_ext,v_10m_ext,rain_ext,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads in a GRIB file containing native coordinate RUC2 data
!  (Grid 236 or Grid 252) and extracts/converts selected variables
!  for use by EXT2ARPS.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Eric Kemp
!  09/17/1999
!  Based on subroutine GETNMCRUC87.
!
!  MODIFICATION HISTORY:
!
!  11/01/1999 Eric Kemp
!  Corrected several bugs in variable retrievals, and added
!  retrieval of RUC2 hydrometeor fields.
!
!  1999/11/30 Gene Bassett
!  Corrected RUC2 soil moisture and changed deep soil moisture &
!  temperature to be an average from 0 to 100cm.
!
!  2009/06/17 Y. Wang
!  Added GRIB2 format and near surface fields.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    dir_extd      Directory name for external file
!    extdname      Prefix string of external file name
!    extdopt       Option of external data sources
!    extdfmt       Option of external data format

!    extdinit      Initialized time in mm-dd-yyyy:hh:mm:ss format
!    extdfcst      Forecast hour in HHH:MM:SS format
!    julfname      File name in yyjjjhhmm format
!
!  OUTPUT:
!
!    iproj_ext     Map projection number of external data
!    scale_ext     Scale factor of external data
!    trlon_ext     True longitude of external data (degrees E)
!    latnot_ext(2) True latitude(s) of external data (degrees N)
!    x0_ext        x coordinate of origin of external data
!    y0_ext        y coordinate of origin of external data
!    lat_ext       latitude of external data points (degrees N)
!    lon_ext       longitude of external data points (degrees E)
!    p_ext         pressure (Pascal)
!    hgt_ext       height (m)
!    t_ext         temperature (K)
!    qv_ext        specific humidity (kg/kg)
!    u_ext         u wind component (m/s)
!    v_ext         v wind component (m/s)
!    qc_ext        Cloud water mixing ratio (kg/kg)
!    qr_ext        Rain water mixing ratio (kg/kg)
!    qi_ext        Ice mixing ratio (kg/kg)
!    qs_ext        Snow mixing ratio (kg/kg)
!    qh_ext        Hail mixing ratio (kg/kg)
!
!    tsfc_ext      Surface temperature
!    tdeep_ext     Soil temperature
!    wetsfc_ext    Top layer soil moisture (fraction)
!    wetdp_ext     Deep soil moisture (fraction)
!    wetcanp_ext   Water content on canopy
!
!    trn_ext       External terrain (m)
!    psfc_ext      Surface pressure (Pa)
!    snowdpth_ext  Snow depth (m)
!
!    istatus       status indicator
!
!  WORK ARRAYS:
!
!    var_grb3d     Arrays to store the GRIB 3-D variables:
!                  var_grb3d(nxgrb,nygrb,nzgrb,1,1) - pressure (Pa)
!                  var_grb3d(nxgrb,nygrb,nzgrb,2,1) - height (m)
!                  var_grb3d(nxgrb,nygrb,nzgrb,3,1) - Virtual Potential
!                                                     temperature (K)
!                  var_grb3d(nxgrb,nygrb,nzgrb,4,1) - Water vapor
!                                                     mixing ratio
!                                                     (kg/kg)
!                  var_grb3d(nxgrb,nygrb,nzgrb,5,1) - u wind (m/s)
!                  var_grb3d(nxgrb,nygrb,nzgrb,6,1) - v wind (m/s)
!                  var_grb3d(nxgrb,nygrb,nzgrb,7,1)
!                        - cloud water mixing ratio (kg/kg)
!                  var_grb3d(nxgrb,nygrb,nzgrb,8,1)
!                        - rain water mixing ratio (kg/kg)
!                  var_grb3d(nxgrb,nygrb,nzgrb,9,1)
!                        - ice mixing ratio (kg/kg)
!                  var_grb3d(nxgrb,nygrb,nzgrb,10,1)
!                        - snow mixing ratio (kg/kg)
!                  var_grb3d(nxgrb,nygrb,nzgrb,11,1)
!                        - graupel mixing ratio (kg/kg)
!                  var_grb3d(nxgrb,nygrb,nzgrb,1,2) - soil temp. (K)
!                  var_grb3d(nxgrb,nygrb,nzgrb,2,2) - soil moist.
!                                                    (fraction)
!
!    var_grb2d     Arrays to store the GRIB 2-D variables:
!                  var_grb2d(nxgrb,nygrb,1) - Canopy Water
!                                             (kg/m**2)
!                  var_grb2d(nxgrb,nygrb,2) - Snow depth (m)
!                  var_grb2d(nxgrb,nygrb,3) - Surface soil temperature (K)
!                  var_grb2d(nxgrb,nygrb,4) - Surface soil moisture
!                                             (fraction)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'gribcst.inc'

  CHARACTER (LEN=*) :: dir_extd
  CHARACTER (LEN=*) :: extdname

  INTEGER :: extdopt
  INTEGER :: extdfmt

  CHARACTER (LEN=19) :: extdinit
  CHARACTER (LEN=9) :: extdfcst
  CHARACTER (LEN=9) :: julfname
!
!
!-----------------------------------------------------------------------
!
!  External grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj_ext
  REAL :: scale_ext,trlon_ext
  REAL :: latnot_ext(2)
  REAL :: x0_ext,y0_ext
  REAL :: dx_ext,dy_ext
!
!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx_ext, ny_ext, nz_ext

  REAL :: lat_ext(nx_ext,ny_ext)
  REAL :: lon_ext(nx_ext,ny_ext)
  REAL :: p_ext  (nx_ext,ny_ext,nz_ext)   ! Pressure (Pascals)
  REAL :: hgt_ext(nx_ext,ny_ext,nz_ext)   ! Height (m)
  REAL :: t_ext  (nx_ext,ny_ext,nz_ext)   ! Temperature (K)
  REAL :: qv_ext (nx_ext,ny_ext,nz_ext)   ! Specific humidity (kg/kg)
  REAL :: u_ext  (nx_ext,ny_ext,nz_ext)   ! Eastward wind component
  REAL :: v_ext  (nx_ext,ny_ext,nz_ext)   ! Northward wind component
  REAL :: qc_ext (nx_ext,ny_ext,nz_ext)   ! Cloud H2O mixing ratio (kg/kg)
  REAL :: qr_ext (nx_ext,ny_ext,nz_ext)   ! Rain  H2O mixing ratio (kg/kg)
  REAL :: qi_ext (nx_ext,ny_ext,nz_ext)   ! Ice   mixing ratio (kg/kg)
  REAL :: qs_ext (nx_ext,ny_ext,nz_ext)   ! Snow  mixing ratio (kg/kg)
  REAL :: qh_ext (nx_ext,ny_ext,nz_ext)   ! Hail  mixing ratio (kg/kg)

  REAL :: tsfc_ext   (nx_ext,ny_ext)      ! Temperature at surface (K)
  REAL :: tdeep_ext  (nx_ext,ny_ext)      ! Deep soil temperature (K)
  REAL :: wetsfc_ext (nx_ext,ny_ext)      ! Surface soil moisture (fraction)
  REAL :: wetdp_ext  (nx_ext,ny_ext)      ! Deep soil moisture (fraction)
  REAL :: wetcanp_ext(nx_ext,ny_ext)      ! Canopy water amount

  REAL :: trn_ext    (nx_ext,ny_ext)          ! Geometrical heights
  REAL :: psfc_ext   (nx_ext,ny_ext)      ! Surface pressure (Pa)
  REAL :: snowdpth_ext (nx_ext,ny_ext)    ! Snow depth (m)

  REAL, INTENT(OUT) ::  t_2m_ext (nx_ext,ny_ext)  ! Temperature(K) at 2m
  REAL, INTENT(OUT) :: qv_2m_ext (nx_ext,ny_ext)  ! Humidity Mixing Ratio (kg kg-1) at 2m
  REAL, INTENT(OUT) :: u_10m_ext (nx_ext,ny_ext)  ! U-Component of Wind (m s-1) at 10m
  REAL, INTENT(OUT) :: v_10m_ext (nx_ext,ny_ext)  ! V-Component of Wind (m s-1) at 2m
  REAL, INTENT(OUT) :: rain_ext (nx_ext,ny_ext)
!
!-----------------------------------------------------------------------
!
!  Other  external variable arrays
!
!-----------------------------------------------------------------------
!
  REAL :: x_ext(nx_ext)
  REAL :: y_ext(ny_ext)

  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Work arrays for storing grib data
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: var_grb2d(:,:,:,:)   ! GRIB variables
  REAL, ALLOCATABLE :: var_grb3d(:,:,:,:,:) ! GRIB 3-D variables
  INTEGER, ALLOCATABLE :: var_lev3d(:,:,:)  ! Levels (hybrid) for
                                            ! each 3-D variable
  REAL, ALLOCATABLE :: rcdata(:)            ! temporary data array
!
!-----------------------------------------------------------------------
!
!  Original grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj
  REAL :: scale,trlon,x0,y0
  REAL :: latnot(2)
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: gribfile
  CHARACTER (LEN=14) :: gribtime
  INTEGER :: i,j,k
  INTEGER :: grbflen, grbtlen

  INTEGER :: m,n,nz1,max_nr2d,max_nr3d

  INTEGER, PARAMETER :: nzsoilin_ext = 6
  REAL,    PARAMETER :: soildepth_ext(nzsoilin_ext) = (/0.0, 0.05, 0.2, 0.4, 1.6, 3.0/)

  REAL    :: pilev
  REAL    :: tv_ext, tvc_ext
  REAL    :: rovcp_p, cpd_p, g0_p, rd_p

  INTEGER :: chklev, lvscan, kk, jj

  REAL    :: tema, temb
  REAL    :: a,b
  INTEGER :: iret             ! Return flag

  REAL, ALLOCATABLE :: utmp(:,:), vtmp(:,:)
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  GRIB grid information
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=42) :: gridesc ! Grid description

  INTEGER :: iproj_grb    ! Map projection indicator
  INTEGER :: gthin        ! Indicator of whether the grid is "thinned"

  INTEGER :: ni_grb       ! Number of points along x-axis
  INTEGER :: nj_grb       ! Number of points along y-axis
  INTEGER :: np_grb       ! Total number of horizontal grid points

  INTEGER :: nk_grb       ! Number of vertical parameters
  REAL :: zk_grb(nz_ext)  ! Vertical coordinate parameters

  INTEGER :: npeq         ! Number of lat circles from pole to equator
  INTEGER :: nit(nz_ext)  ! Number of x-points for thinned grid

  REAL :: pi_grb          ! x-coordinate of pole point
  REAL :: pj_grb          ! y-coordinate of pole point
  INTEGER :: ipole        ! Projection center flag

  REAL :: di_grb          ! x-direction increment or grid length
  REAL :: dj_grb          ! y-direction increment or grid length

  REAL :: latsw           ! Latitude  of South West corner point
  REAL :: lonsw           ! Longitude of South West corner point
  REAL :: latne           ! Latitude  of North East corner point
  REAL :: lonne           ! Longitude of North East corner point

  REAL :: lattru1         ! Latitude (1st) at which projection is true
  REAL :: lattru2         ! Latitude (2nd) at which projection is true
  REAL :: lontrue         ! Longitude      at which projection is true

  REAL :: latrot          ! Latitude  of southern pole of rotation
  REAL :: lonrot          ! Longitude of southern pole of rotation
  REAL :: angrot          ! Angle of rotation

  REAL :: latstr          ! Latitude  of the pole of stretching
  REAL :: lonstr          ! Longitude of the pole of stretching
  REAL :: facstr          ! Stretching factor

  INTEGER :: scanmode     ! Scanning indicator
  INTEGER :: iscan        ! x-direction   scanning indicator
  INTEGER :: jscan        ! y-direction   scanning indicator
  INTEGER :: kscan        ! FORTRAN index scanning indicator

  INTEGER :: ires         ! Resolution direction increments indicator
  INTEGER :: iearth       ! Earth shape indicator: spherical or oblate?
  INTEGER :: icomp        ! (u,v) components decomposition indicator

  INTEGER :: jpenta       ! J-Pentagonal resolution parameter
  INTEGER :: kpenta       ! K-Pentagonal resolution parameter
  INTEGER :: mpenta       ! M-Pentagonal resolution parameter
  INTEGER :: ispect       ! Spectral representation type
  INTEGER :: icoeff       ! Spectral coefficient storage mode

  REAL :: xp_grb          ! X coordinate of sub-satellite point
  REAL :: yp_grb          ! Y coordinate of sub-satellite point
  REAL :: xo_grb          ! X coordinate of image sector origin
  REAL :: yo_grb          ! Y coordinate of image sector origin
  REAL :: zo_grb          ! Camera altitude from center of Earth

  INTEGER :: igrbfmt
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ALLOCATE(utmp(nx_ext,ny_ext))
  ALLOCATE(vtmp(nx_ext,ny_ext))

  CALL getmapr(iproj,scale,latnot,trlon,x0,y0)

  CALL getgrbfname(extdopt,extdfmt,dir_extd,extdname,                   &
                   extdinit,extdfcst,'    ',                            &
                   gribfile,grbflen,gribtime,grbtlen,istatus)
  IF (istatus /=0 ) RETURN

  CALL chkgrb(gribfile(1:grbflen),grbflen,igrbfmt,istatus)

  IF (igrbfmt == 2) THEN

    WRITE(6,'(1x,2a,I4)') 'GRIB2 file to be read is ',gribfile(1:grbflen),grbflen

    CALL  getruc236n_grb2(nx_ext,ny_ext,nz_ext,nzsoilin_ext,            &
             gribfile,gribtime,grbflen,soildepth_ext,                   &
             dx_ext,dy_ext,                                             &
             iproj_ext,scale_ext,trlon_ext,latnot_ext,latsw,lonsw,      &
             p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                    &
             qc_ext,qr_ext,qi_ext,qs_ext,qh_ext,                        &
             tsfc_ext,tdeep_ext,wetsfc_ext,wetdp_ext,wetcanp_ext,       &
             snowdpth_ext,trn_ext,psfc_ext,t_2m_ext,qv_2m_ext,          &
             u_10m_ext,v_10m_ext,rain_ext,istatus)

    IF (istatus /= 0) RETURN

  ELSE

    ALLOCATE(var_grb2d(nx_ext,ny_ext,n2dvs,n2dlvt))
    ALLOCATE(var_grb3d(nx_ext,ny_ext,nz_ext,n3dvs,n3dlvt))
    ALLOCATE(rcdata(nx_ext*ny_ext))
    ALLOCATE(var_lev3d(nz_ext,n3dvs,n3dlvt))

!
!-----------------------------------------------------------------------
!
!  RDNMCGRB reads NMC GRIB data
!
!-----------------------------------------------------------------------
!
    IF (extdopt == 18) THEN
      gridtyp   = rucn252grid
    ELSE
      gridtyp   = rucn236grid
    END IF
    mproj_grb = rucn236proj

    n2dvars  = rucn236nvs2d
    n2dlvtps = rucn236nlvt2d

    DO m=1,n2dlvtps
      DO n=1,n2dvars
        var_id2d(n,m) = rucn236var_id2d(n,m)
      END DO
      levtyp2d(m) = rucn236levs2d(m)
    END DO

    n3dvars  = rucn236nvs3d
    n3dlvtps = rucn236nlvt3d

    DO m=1,n3dlvtps
      DO n=1,n3dvars
        var_id3d(n,m) = rucn236var_id3d(n,m)
      END DO
      levtyp3d(m) = rucn236levs3d(m)
    END DO

    CALL rdnmcgrb(nx_ext,ny_ext,nz_ext,gribfile,grbflen, gribtime,        &
                  gridesc, iproj_grb, gthin,                              &
                  ni_grb,nj_grb,np_grb, nk_grb,zk_grb, npeq,nit,          &
                  pi_grb,pj_grb,ipole, di_grb,dj_grb,                     &
                  latsw,lonsw, latne,lonne,                               &
                  latrot,lonrot,angrot,                                   &
                  latstr,lonstr,facstr,                                   &
                  lattru1,lattru2,lontrue,                                &
                  scanmode, iscan,jscan,kscan,                            &
                  ires,iearth,icomp,                                      &
                  jpenta,kpenta,mpenta,ispect,icoeff,                     &
                  xp_grb,yp_grb, xo_grb,yo_grb,zo_grb,                    &
                  rcdata,var_grb2d,var_grb3d,var_lev3d,lvldbg,iret)

    IF (iret /= 0)  THEN
      ISTATUS = -888
      GOTO 999
    END IF

    max_nr2d = 0
    DO n=1,n2dvars
      DO m=1,n2dlvtps
        max_nr2d = MAX( max_nr2d, var_nr2d(n,m) )
      END DO
    END DO

    max_nr3d = 0
    DO n=1,n3dvars
      max_nr3d = MAX( max_nr3d, var_nr3d(n,1) )
    END DO

    IF ( max_nr3d == 0 ) THEN
      WRITE (6,'(a)')                                                     &
          'No 3-D variables was found in the GRIB file',                  &
          'Program stopped in GETNMCRUC.'
  !   STOP
      ISTATUS = -888
      GOTO 999
    END IF

    IF ( max_nr2d == 0 ) THEN
      WRITE (6,'(a)')                                                     &
          'No 2-D variables was found in the GRIB file'
    END IF

  !  write (6,'(/a7,2x,6(i7))')
  !    :      'Lev\\VID',(var_id3d(n,1),n=1,n3dvars)

  !  DO 60 k=1,max_nr3d
  !     write (6,'(i5,4x,6(i7))')
  !    :      k,(var_lev3d(k,n,1),n=1,n3dvars)
    60    CONTINUE

    DO k=1,max_nr3d
      DO n=2,n3dvars
        IF ( var_lev3d(k,1,1) /= var_lev3d(k,n,1) ) THEN
          WRITE (6,'(a)')                                                 &
              'Variables were not at the same level.',                    &
              'Program stopped in GETNMCRUC.'
          WRITE(6,*)                                                      &
              'var_lev3d(k,1,1) = ',var_lev3d(k,1,1),                     &
              'var_lev3d(k,n,1) = ',var_lev3d(k,n,1)
  !       STOP
          ISTATUS = -888
          GOTO 999
        END IF
      END DO
    END DO

    IF ( iproj_grb == 5 .AND. ipole == 0 ) THEN      ! Center in NP
      iproj_ext = 1
    ELSE IF ( iproj_grb == 5 .AND. ipole == 1 ) THEN  ! Center in SP
      iproj_ext = -1
    ELSE IF ( iproj_grb == 3 .AND. ipole == 0 ) THEN  ! Center in NP
      iproj_ext = 2
    ELSE IF ( iproj_grb == 3 .AND. ipole == 1 ) THEN  ! Center in SP
      iproj_ext = -2
    ELSE IF ( iproj_grb == 1 ) THEN
      iproj_ext = 3
    ELSE IF ( iproj_grb == 0 ) THEN
      iproj_ext = 4
    ELSE
      WRITE (6,'(a)')                                                     &
          'Unknown map projection. Set to non-projection.'
      iproj_ext = 0
    END IF

    scale_ext = 1.0
    latnot_ext(1) = lattru1
    latnot_ext(2) = lattru2
    trlon_ext = lontrue

    dx_ext = di_grb
    dy_ext = dj_grb

!
!-----------------------------------------------------------------------
!
!  Retrieve 2-D variables
!
!-----------------------------------------------------------------------
!
    DO j=1,ny_ext
      DO i=1,nx_ext

        tsfc_ext(i,j) = -999.0

        trn_ext(i,j) = -999.0

        IF ( var_nr3d(1,2) == 0 ) THEN
          tsfc_ext   (i,j) = -999.0
          tdeep_ext  (i,j) = -999.0
          wetsfc_ext (i,j) = -999.0
          wetdp_ext  (i,j) = -999.0
          wetcanp_ext(i,j) = -999.0
        ELSE
  !      tsfc_ext   (i,j) = var_grb3d(i,j,1,1,2)
          tsfc_ext   (i,j) = var_grb2d(i,j,3,1)
          tdeep_ext  (i,j) = 0.1 * var_grb3d(i,j,1,1,2)  & !   5cm
                       + 0.2 * var_grb3d(i,j,2,1,2)      & !  20cm
                       + 0.4 * var_grb3d(i,j,3,1,2)      & !  40cm
                       + 0.3 * var_grb3d(i,j,4,1,2)        ! 160cm
  !      wetsfc_ext (i,j) = var_grb3d(i,j,1,2,2)
          wetsfc_ext (i,j) = var_grb2d(i,j,4,1)
          wetdp_ext  (i,j) = 0.1 * var_grb3d(i,j,1,2,2)  & !   5cm
                       + 0.2 * var_grb3d(i,j,2,2,2)      & !  20cm
                       + 0.4 * var_grb3d(i,j,3,2,2)      & !  40cm
                       + 0.3 * var_grb3d(i,j,4,2,2)        ! 160cm
          wetcanp_ext(i,j) = var_grb2d(i,j,1,1)*1.e-3    ! in meters
        END IF

        psfc_ext   (i,j) = -999.0

        IF ( var_nr2d(2,1) == 0 ) THEN
          snowdpth_ext(i,j) = -999
        ELSE
          snowdpth_ext(i,j) = var_grb2d(i,j,2,1)  ! in meters
        END IF

        t_2m_ext (i,j) = var_grb2d(i,j,1,2)    ! T
        qv_2m_ext(i,j) = var_grb2d(i,j,2,2)    ! QV
        u_10m_ext(i,j) = var_grb2d(i,j,3,2)    ! U
        v_10m_ext(i,j) = var_grb2d(i,j,4,2)    ! V

      END DO
    END DO

!
!-----------------------------------------------------------------------
!
!  Retrieve 3-D variables
!
!-----------------------------------------------------------------------
!

    cpd_p   = 1004.686    ! cp in RUC
    rovcp_p = 0.285714    ! rd/cp used in RUC
    g0_p    = 9.80665     ! gravity in RUC

    nz1 = MIN(var_nr3d(1,1),nz_ext)

    IF ( var_lev3d(1,1,1) < var_lev3d(nz1,1,1) ) THEN  ! 1st level at sfc
      chklev = 1
      lvscan = 0
    ELSE
      chklev = -1
      lvscan = nz1+1
    END IF

    DO k=1,nz1
      kk = chklev * k + lvscan
      DO j=1,ny_ext
        DO i=1,nx_ext

          p_ext(i,j,kk) = var_grb3d(i,j,k,1,1)       ! Pressure (Pa)

          hgt_ext(i,j,kk) = var_grb3d(i,j,k,2,1)     ! Height (m)

          u_ext(i,j,kk) = var_grb3d(i,j,k,5,1)       ! u wind (m/s)
          v_ext(i,j,kk) = var_grb3d(i,j,k,6,1)       ! v wind (m/s)

          a = REAL(100000)/var_grb3d(i,j,k,1,1)
          a = a**rovcp_p
          tvc_ext = var_grb3d(i,j,k,3,1)/a ! Virtual Temperature
          b = 0.61*var_grb3d(i,j,k,4,1)
          b = REAL(1) + b
          t_ext(i,j,kk) = tvc_ext/b ! Temperature (K)

          a = var_grb3d(i,j,k,4,1)*var_grb3d(i,j,k,1,1)
          a = a/(0.622 - var_grb3d(i,j,k,4,1))
          qv_ext(i,j,kk) = a*0.622/var_grb3d(i,j,k,1,1) ! SpecificHumidity
  !
  !-----------------------------------------------------------------------
  !
  !      Retrieve hydrometeor data.
  !
  !-----------------------------------------------------------------------
  !
          IF (var_nr3d(7,1) == 0) THEN
            qc_ext(i,j,kk) = -999.
          ELSE
            qc_ext(i,j,kk) = var_grb3d(i,j,k,7,1)
          END IF

          IF (var_nr3d(8,1) == 0) THEN
            qr_ext(i,j,kk) = -999.
          ELSE
            qr_ext(i,j,kk) = var_grb3d(i,j,k,8,1)
          END IF

          IF (var_nr3d(9,1) == 0) THEN
            qi_ext(i,j,kk) = -999.
          ELSE
            qi_ext(i,j,kk) = var_grb3d(i,j,k,9,1)
          END IF

          IF (var_nr3d(10,1) == 0) THEN
            qs_ext(i,j,kk) = -999.
          ELSE
            qs_ext(i,j,kk) = var_grb3d(i,j,k,10,1)
          END IF

          IF (var_nr3d(11,1) == 0) THEN
            qh_ext(i,j,kk) = -999.
          ELSE
            qh_ext(i,j,kk) = var_grb3d(i,j,k,11,1)
          END IF

        END DO
      END DO
    END DO

    DEALLOCATE(var_grb2d,var_grb3d,rcdata,var_lev3d)
  END IF

  CALL setmapr(iproj_ext,scale_ext,latnot_ext,trlon_ext)
  CALL lltoxy(1,1,latsw,lonsw,x0_ext,y0_ext)

  DO i=1,nx_ext
    x_ext(i)=x0_ext+(i-1)*dx_ext
  END DO

  DO j=1,ny_ext
    y_ext(j)=y0_ext+(j-1)*dy_ext
  END DO

  CALL xytoll(nx_ext,ny_ext,x_ext,y_ext,lat_ext,lon_ext)
!
!-----------------------------------------------------------------------
!
!  Rotate winds to be relative to true north.
!  The RUC data are sent as grid-relative.
!
!-----------------------------------------------------------------------
!
  DO k=1,nz_ext
    CALL uvmptoe(nx_ext,ny_ext,u_ext(1,1,k),v_ext(1,1,k),               &
                 lon_ext,utmp,vtmp)
    u_ext(:,:,k) = utmp(:,:)
    v_ext(:,:,k) = vtmp(:,:)
  END DO

  istatus = 1

!
!-----------------------------------------------------------------------
!
!  Reset map projection to previous values
!
!-----------------------------------------------------------------------
!
  999  CONTINUE
  CALL setmapr(iproj,scale,latnot,trlon)
  CALL setorig(1,x0,y0)

  DEALLOCATE(utmp)
  DEALLOCATE(vtmp)

  RETURN
END SUBROUTINE getnmcrucn236
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE GETNMCRUCP236               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getnmcrucp236(nx_ext,ny_ext,nz_ext,                          &
           dir_extd,extdname,extdopt,extdfmt,                           &
           extdinit,extdfcst,julfname,                                  &
           iproj_ext,scale_ext,                                         &
           trlon_ext,latnot_ext,x0_ext,y0_ext,                          &
           lat_ext,lon_ext,                                             &
           p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                      &
           tsfc_ext,tdeep_ext,wetsfc_ext,wetdp_ext,wetcanp_ext,         &
           trn_ext,psfc_ext, t_2m_ext, rh_2m_ext,                       &
           u_10m_ext, v_10m_ext, snowdpth_ext,                          &
           istatus)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads in a GRIB file containing isobaric RUC2 data (Grid 236)
!  and extracts/converts selected variables for use by EXT2ARPS.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Eric Kemp
!  09/17/1999
!  Based on subroutine GETNMCRUC211.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    dir_extd      Directory name for external file
!    extdname      Prefix string of external file name
!    extdopt       Option of external data sources
!    extdfmt       Option of external data format

!    extdinit      Initialized time in mm-dd-yyyy:hh:mm:ss format
!    extdfcst      Forecast hour in HHH:MM:SS format
!    julfname      File name in yyjjjhhmm format
!
!  OUTPUT:
!
!    iproj_ext     Map projection number of external data
!    scale_ext     Scale factor of external data
!    trlon_ext     True longitude of external data (degrees E)
!    latnot_ext(2) True latitude(s) of external data (degrees N)
!    x0_ext        x coordinate of origin of external data
!    y0_ext        y coordinate of origin of external data
!    lat_ext       latitude of external data points (degrees N)
!    lon_ext       longitude of external data points (degrees E)
!    p_ext         pressure (Pascal)
!    hgt_ext       height (m)
!    t_ext         temperature (K)
!    qv_ext        specific humidity (kg/kg)
!    u_ext         u wind component (m/s)
!    v_ext         v wind component (m/s)
!    qc_ext        Cloud water mixing ratio (kg/kg)
!    qr_ext        Rain water mixing ratio (kg/kg)
!    qi_ext        Ice mixing ratio (kg/kg)
!    qs_ext        Snow mixing ratio (kg/kg)
!    qh_ext        Hail mixing ratio (kg/kg)
!
!    tsfc_ext      Surface temperature
!    tdeep_ext     Soil temperature
!    wetsfc_ext    Top layer soil moisture (fraction)
!    wetdp_ext     Deep soil moisture (fraction)
!    wetcanp_ext   Water content on canopy
!
!    trn_ext       External terrain (m)
!    psfc_ext      Surface pressure (Pa)
!
!    T_2m_ext      Temperature at 2m AGL
!    rh_2m_ext     Specific Humidity at 2m AGL
!    U_10m_ext     U at 10m AGL
!    V_10m_ext     V at 10m AGL
!
!    snowdpth_ext  Snow depth (m)
!
!    istatus       status indicator
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'gribcst.inc'

  CHARACTER (LEN=*) :: dir_extd
  CHARACTER (LEN=*) :: extdname

  INTEGER :: extdopt
  INTEGER :: extdfmt

  CHARACTER (LEN=19) :: extdinit
  CHARACTER (LEN=9) :: extdfcst
  CHARACTER (LEN=9) :: julfname
!
!-----------------------------------------------------------------------
!
!  External grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj_ext
  REAL :: scale_ext,trlon_ext
  REAL :: latnot_ext(2)
  REAL :: x0_ext,y0_ext
  REAL :: dx_ext,dy_ext
!
!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx_ext,ny_ext,nz_ext

  REAL :: lat_ext(nx_ext,ny_ext)
  REAL :: lon_ext(nx_ext,ny_ext)
  REAL :: p_ext  (nx_ext,ny_ext,nz_ext)   ! Pressure (Pascals)
  REAL :: hgt_ext(nx_ext,ny_ext,nz_ext)   ! Height (m)
  REAL :: t_ext  (nx_ext,ny_ext,nz_ext)   ! Temperature (K)
  REAL :: qv_ext (nx_ext,ny_ext,nz_ext)   ! Specific humidity (kg/kg)
  REAL :: u_ext  (nx_ext,ny_ext,nz_ext)   ! Eastward wind component
  REAL :: v_ext  (nx_ext,ny_ext,nz_ext)   ! Northward wind component
!  REAL :: qc_ext (nx_ext,ny_ext,nz_ext)   ! Cloud H2O mixing ratio
!                                          ! (kg/kg)
!  REAL :: qr_ext (nx_ext,ny_ext,nz_ext)   ! Rain  H2O mixing ratio
!                                          ! (kg/kg)
!  REAL :: qi_ext (nx_ext,ny_ext,nz_ext)   ! Ice   mixing ratio (kg/kg)
!  REAL :: qs_ext (nx_ext,ny_ext,nz_ext)   ! Snow  mixing ratio (kg/kg)
!  REAL :: qh_ext (nx_ext,ny_ext,nz_ext)   ! Hail  mixing ratio (kg/kg)
!
  REAL :: tsfc_ext   (nx_ext,ny_ext)      ! Temperature at surface (K)
  REAL :: tdeep_ext  (nx_ext,ny_ext)      ! Deep soil temperature (K)
  REAL :: wetsfc_ext (nx_ext,ny_ext)      ! Surface soil moisture (fraction)
  REAL :: wetdp_ext  (nx_ext,ny_ext)      ! Deep soil moisture (fraction)
  REAL :: wetcanp_ext(nx_ext,ny_ext)      ! Canopy water amount

  REAL :: trn_ext    (nx_ext,ny_ext)          ! Geometrical heights
  REAL :: psfc_ext   (nx_ext,ny_ext)      ! Surface pressure (Pa)

  REAL :: t_2m_ext (nx_ext,ny_ext)
  REAL :: rh_2m_ext(nx_ext,ny_ext)
  REAL :: u_10m_ext(nx_ext,ny_ext)
  REAL :: v_10m_ext(nx_ext,ny_ext)

  REAL :: snowdpth_ext(nx_ext,ny_ext)     ! Snow depth (m)
!
!----------------------------------------------------------------------
!
!  Other  external variable arrays
!
!----------------------------------------------------------------------
!
  REAL :: x_ext(nx_ext)
  REAL :: y_ext(ny_ext)
  REAL :: z_ext(nz_ext)

  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Work arrays for storing grib data
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: var_grb2d(:,:,:,:)   ! GRIB variables
  REAL, ALLOCATABLE :: var_grb3d(:,:,:,:,:) ! GRIB 3-D variables
  INTEGER, ALLOCATABLE :: var_lev3d(:,:,:)  ! Levels (hybrid) for
                                            ! each 3-D variable
  REAL, ALLOCATABLE :: rcdata(:)            ! temporary data array
!
!-----------------------------------------------------------------------
!
!  Original grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj
  REAL :: scale,trlon,x0,y0
  REAL :: latnot(2)
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: gribfile
  CHARACTER (LEN=14) :: gribtime
  CHARACTER (LEN=10) :: runstr
  CHARACTER (LEN=3) :: fmtn
  INTEGER :: ichr,bchar,echar
  INTEGER :: i,j,k,ldir,ireturn
  INTEGER :: grbflen, grbtlen

  INTEGER :: m,n,nz1,max_nr2d,max_nr3d

  REAL :: qvsat, pilev, qvsatice

  REAL :: tema, temb

  INTEGER :: chklev, lvscan, kk, jj

  INTEGER :: iret             ! Return flag

  REAL, ALLOCATABLE :: utmp(:,:), vtmp(:,:)
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
!
!
!-----------------------------------------------------------------------
!
!  GRIB grid information
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=42) :: gridesc ! Grid description

  INTEGER :: iproj_grb    ! Map projection indicator
  INTEGER :: gthin        ! Indicator of whether the grid is "thinned"

  INTEGER :: ni_grb       ! Number of points along x-axis
  INTEGER :: nj_grb       ! Number of points along y-axis
  INTEGER :: np_grb       ! Total number of horizontal grid points

  INTEGER :: nk_grb       ! Number of vertical parameters
  REAL :: zk_grb(nz_ext)  ! Vertical coordinate parameters

  INTEGER :: npeq         ! Number of lat circles from pole to equator
  INTEGER :: nit(nz_ext)  ! Number of x-points for thinned grid

  REAL :: pi_grb          ! x-coordinate of pole point
  REAL :: pj_grb          ! y-coordinate of pole point
  INTEGER :: ipole        ! Projection center flag

  REAL :: di_grb          ! x-direction increment or grid length
  REAL :: dj_grb          ! y-direction increment or grid length

  REAL :: lrb             ! y-direction increment or grid length

  REAL :: latsw           ! Latitude  of South West corner point
  REAL :: lonsw           ! Longitude of South West corner point
  REAL :: latne           ! Latitude  of North East corner point
  REAL :: lonne           ! Longitude of North East corner point

  REAL :: lattru1         ! Latitude (1st) at which projection is true
  REAL :: lattru2         ! Latitude (2nd) at which projection is true
  REAL :: lontrue         ! Longitude      at which projection is true

  REAL :: latrot          ! Latitude  of southern pole of rotation
  REAL :: lonrot          ! Longitude of southern pole of rotation
  REAL :: angrot          ! Angle of rotation

  REAL :: latstr          ! Latitude  of the pole of stretching
  REAL :: lonstr          ! Longitude of the pole of stretching
  REAL :: facstr          ! Stretching factor

  INTEGER :: scanmode     ! Scanning indicator
  INTEGER :: iscan        ! x-direction   scanning indicator
  INTEGER :: jscan        ! y-direction   scanning indicator
  INTEGER :: kscan        ! FORTRAN index scanning indicator

  INTEGER :: ires         ! Resolution direction increments indicator
  INTEGER :: iearth       ! Earth shape indicator: spherical or oblate?
  INTEGER :: icomp        ! (u,v) components decomposition indicator

  INTEGER :: jpenta       ! J-Pentagonal resolution parameter
  INTEGER :: kpenta       ! K-Pentagonal resolution parameter
  INTEGER :: mpenta       ! M-Pentagonal resolution parameter
  INTEGER :: ispect       ! Spectral representation type
  INTEGER :: icoeff       ! Spectral coefficient storage mode

  REAL :: xp_grb          ! X coordinate of sub-satellite point
  REAL :: yp_grb          ! Y coordinate of sub-satellite point
  REAL :: xo_grb          ! X coordinate of image sector origin
  REAL :: yo_grb          ! Y coordinate of image sector origin
  REAL :: zo_grb          ! Camera altitude from center of Earth

  INTEGER :: igrbfmt

  INTEGER, PARAMETER :: nzsoilin_ext = 1
  REAL,    PARAMETER :: soildepth_ext(nzsoilin_ext) = (/0.0/)
!
!
!-----------------------------------------------------------------------
!
!  Function f_qvsatl and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_qvsatl
!fpp$ expand (f_desdt)
!fpp$ expand (f_qvsatl)
!!dir$ inline always f_desdt, f_qvsatl
!*$*  inline routine (f_desdt, f_qvsatl)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ALLOCATE(utmp(nx_ext,ny_ext))
  ALLOCATE(vtmp(nx_ext,ny_ext))

  CALL getmapr(iproj,scale,latnot,trlon,x0,y0)

  CALL getgrbfname(extdopt,extdfmt,dir_extd,extdname,                   &
                   extdinit,extdfcst,'    ',                            &
                   gribfile,grbflen,gribtime,grbtlen,istatus)
  IF (istatus /=0 ) RETURN

  CALL chkgrb(gribfile(1:grbflen),grbflen,igrbfmt,istatus)

  IF (igrbfmt == 2) THEN

    WRITE(6,'(1x,2a,I4)') 'GRIB2 file to be read is ',gribfile(1:grbflen),grbflen

    CALL  getruc236p_grb2(nx_ext,ny_ext,nz_ext,nzsoilin_ext,            &
             gribfile,gribtime,grbflen,soildepth_ext,                   &
             dx_ext,dy_ext,                                             &
             iproj_ext,scale_ext,trlon_ext,latnot_ext,latsw,lonsw,      &
             p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                    &
             tsfc_ext,tdeep_ext,wetsfc_ext,wetdp_ext,wetcanp_ext,       &
             snowdpth_ext,trn_ext,psfc_ext,t_2m_ext,rh_2m_ext,          &
             u_10m_ext,v_10m_ext,istatus)

    IF (istatus /= 0) RETURN

  ELSE

    ALLOCATE(var_grb2d(nx_ext,ny_ext,n2dvs,n2dlvt))
    ALLOCATE(var_grb3d(nx_ext,ny_ext,nz_ext,n3dvs,n3dlvt))
    ALLOCATE(rcdata(nx_ext*ny_ext))
    ALLOCATE(var_lev3d(nz_ext,n3dvs,n3dlvt))

!
!-----------------------------------------------------------------------
!
!  RDNMCGRB reads NMC GRIB data
!
!-----------------------------------------------------------------------
!
    IF (extdopt == 19) THEN
      gridtyp   = rucp252grid
    ELSE
      gridtyp   = rucp236grid
    END IF
    mproj_grb = rucp236proj

    n2dvars  = rucp236nvs2d
    n2dlvtps = rucp236nlvt2d

    DO k=1,n2dlvtps
      DO n=1,n2dvars
        var_id2d(n,k) = rucp236var_id2d(n,k)
      END DO
      levtyp2d(k) = rucp236levs2d(k)
    END DO

    n3dvars  = rucp236nvs3d
    n3dlvtps = rucp236nlvt3d

    DO m=1,n3dlvtps
      DO n=1,n3dvars
        var_id3d(n,m) = rucp236var_id3d(n,m)
      END DO
      levtyp3d(m) = rucp236levs3d(m)
    END DO

    CALL rdnmcgrb(nx_ext,ny_ext,nz_ext,gribfile,grbflen, gribtime,      &
                  gridesc, iproj_grb, gthin,                            &
                  ni_grb,nj_grb,np_grb, nk_grb,zk_grb, npeq,nit,        &
                  pi_grb,pj_grb,ipole, di_grb,dj_grb,                   &
                  latsw,lonsw, latne,lonne,                             &
                  latrot,lonrot,angrot,                                 &
                  latstr,lonstr,facstr,                                 &
                  lattru1,lattru2,lontrue,                              &
                  scanmode, iscan,jscan,kscan,                          &
                  ires,iearth,icomp,                                    &
                  jpenta,kpenta,mpenta,ispect,icoeff,                   &
                  xp_grb,yp_grb, xo_grb,yo_grb,zo_grb,                  &
                  rcdata,var_grb2d,var_grb3d,var_lev3d,lvldbg,iret)

    IF (iret /= 0)  THEN
      ISTATUS = -888
      GOTO 999
    END IF

    max_nr2d = 0
    DO n=1,n2dvars
      DO m=1,n2dlvtps
        max_nr2d = MAX( max_nr2d, var_nr2d(n,m) )
      END DO
    END DO

    max_nr3d = 0
    DO n=1,n3dvars
      max_nr3d = MAX( max_nr3d, var_nr3d(n,1))
    END DO

    IF ( max_nr3d == 0 ) THEN
      WRITE (6,'(a)')                                                   &
          'No 3-D variables was found in the GRIB file',                &
          'Program returned from GETNMCRUCP236.'
      ISTATUS = -888
      GOTO 999
    END IF

    IF ( max_nr2d == 0 ) THEN
      WRITE (6,'(a)')                                                   &
          'No 2-D variables was found in the GRIB file'
    END IF

!  WRITE (6,'(/a7,2x,6(i7))') 'Lev\VID',(var_id3d(n,1),n=1,n3dvars)
!  DO  k=1,max_nr3d
!    WRITE (6,'(i5,4x,6(i7))') k,(var_lev3d(k,n,1),n=1,n3dvars)
!  END DO

    DO k=1,max_nr3d
      DO n=2,n3dvars
        IF ( var_lev3d(k,1,1) /= var_lev3d(k,n,1) ) THEN
          WRITE (6,'(a)')                                               &
              'Variables were not at the same level.',                  &
              'Program stopped in GETNMCRUCP236.'
!         STOP
          ISTATUS = -888
          GOTO 999
        END IF
      END DO
    END DO

    IF ( iproj_grb == 5 .AND. ipole == 0 ) THEN      ! Center in NP
      iproj_ext = 1
    ELSE IF ( iproj_grb == 5 .AND. ipole == 1 ) THEN  ! Center in SP
      iproj_ext = -1
    ELSE IF ( iproj_grb == 3 .AND. ipole == 0 ) THEN  ! Center in NP
      iproj_ext = 2
    ELSE IF ( iproj_grb == 3 .AND. ipole == 1 ) THEN  ! Center in SP
      iproj_ext = -2
    ELSE IF ( iproj_grb == 1 ) THEN
      iproj_ext = 3
    ELSE IF ( iproj_grb == 0 ) THEN
      iproj_ext = 4
    ELSE
      WRITE (6,'(a)')                                                     &
          'Unknown map projection. Set to non-projection.'
      iproj_ext = 0
    END IF

    scale_ext = 1.0
    latnot_ext(1) = lattru1
    latnot_ext(2) = lattru2
    trlon_ext = lontrue

    dx_ext = di_grb
    dy_ext = dj_grb
!
!-----------------------------------------------------------------------
!
!  Retrieve 2-D variables
!
!-----------------------------------------------------------------------
!
    DO j=1,ny_ext
      DO i=1,nx_ext
        IF ( var_nr2d(1,1) == 0 ) THEN
          psfc_ext   (i,j) = -999.0
        ELSE
          psfc_ext   (i,j) = var_grb2d(i,j,1,1)
        END IF

        IF ( var_nr2d(2,1) == 0 ) THEN
          trn_ext    (i,j) = -999.0
        ELSE
          trn_ext    (i,j) = var_grb2d(i,j,2,1)
        END IF

        IF( var_nr2d(1,2) == 0 ) THEN
          t_2m_ext(i,j)= -999.0
        ELSE
          t_2m_ext(i,j)= var_grb2d(i,j,1,2)
        END IF

!     at this point we are reading in the relative humidity
!     later we'll convert to specific humidity

        IF( var_nr2d(2,2) == 0 ) THEN
          rh_2m_ext(i,j)= -999.0
        ELSE
          rh_2m_ext(i,j)= var_grb2d(i,j,2,2)
        END IF

        IF( var_nr2d(3,2) == 0 ) THEN
          u_10m_ext(i,j)= -999.0
        ELSE
          u_10m_ext(i,j)= var_grb2d(i,j,3,2)
        END IF

        IF( var_nr2d(4,2) == 0 ) THEN
          v_10m_ext(i,j)= -999.0
        ELSE
          v_10m_ext(i,j)= var_grb2d(i,j,4,2)
        END IF


        tsfc_ext   (i,j) = -999.0
        tdeep_ext  (i,j) = -999.0
        wetsfc_ext (i,j) = -999.0
        wetdp_ext  (i,j) = -999.0
        wetcanp_ext(i,j) = -999.0


        IF ( var_nr2d(3,1) == 0 ) THEN
          snowdpth_ext(i,j) = -999
        ELSE
          snowdpth_ext(i,j) = var_grb2d(i,j,3,1)   ! in meters
        END IF

      END DO
    END DO

!
!-----------------------------------------------------------------------
!
!  Retrieve 3-D variables
!
!-----------------------------------------------------------------------
!
    nz1 = MIN(var_nr3d(1,1),nz_ext)

    IF ( var_lev3d(1,1,1) > var_lev3d(nz1,1,1) ) THEN  ! 1st level at
                                                       !sfc
      chklev = 1
      lvscan = 0
    ELSE
      chklev = -1
      lvscan = nz1+1
    END IF

    DO k=1,nz1
      kk = chklev * k + lvscan
      DO j=1,ny_ext
        DO i=1,nx_ext
          p_ext(i,j,kk) = 100.0 * REAL(var_lev3d(k,1,1)) ! Pressure
          t_ext(i,j,kk) = var_grb3d(i,j,k,2,1)    ! Temperature (K)
          u_ext(i,j,kk) = var_grb3d(i,j,k,4,1)    ! u wind (m/s)
          v_ext(i,j,kk) = var_grb3d(i,j,k,5,1)    ! v wind (m/s)
          hgt_ext(i,j,kk) = var_grb3d(i,j,k,1,1)    ! height (m)

  !   check for portions of constant pressure grids that are below the
  !   surface
  !
  !       IF(((kk == 1) .OR. (p_ext(i,j,kk) > psfc_ext(i,j))) .AND.       &
  !             (psfc_ext(i,j) > 0.) ) THEN
  !         p_ext(i,j,kk)= psfc_ext(i,j) - 1. * (kk - 1)
  !         t_ext(i,j,kk)= t_2m_ext(i,j)
  !         u_ext(i,j,kk)= u_10m_ext(i,j)
  !         v_ext(i,j,kk)= v_10m_ext(i,j)
  !         hgt_ext(i,j,kk)= trn_ext(i,j) + 0.1 * (kk - 1)
  !       END IF

          IF( (p_ext(i,j,kk) > 0.0) .AND. (t_ext(i,j,kk) > 0.0) ) THEN
            qvsat = f_qvsatl( p_ext(i,j,kk), t_ext(i,j,kk) )
            qv_ext(i,j,kk)= var_grb3d(i,j,k,3,1) * qvsat * 0.01

          IF((kk == 1) .OR. (p_ext(i,j,kk) > psfc_ext(i,j)))THEN
            qv_ext(i,j,kk)= rh_2m_ext(i,j) * qvsat * 0.01
          END IF
!          qc_ext(i,j,kk)= 0.0
!          qi_ext(i,j,kk)= 0.0
        ELSE
          qv_ext(i,j,kk)= 0.0
!          qi_ext(i,j,kk)= 0.0
!          qc_ext(i,j,kk)= 0.0
        END IF

!        qr_ext (i,j,kk) = -999.
!        qs_ext (i,j,kk) = -999.
!        qh_ext (i,j,kk) = -999.
        END DO
      END DO
    END DO

    DEALLOCATE(var_grb2d,var_grb3d,rcdata,var_lev3d)
  END IF

  CALL setmapr(iproj_ext,scale_ext,latnot_ext,trlon_ext)
  CALL lltoxy(1,1,latsw,lonsw,x0_ext,y0_ext)

  DO i=1,nx_ext
    x_ext(i)=x0_ext+(i-1)*dx_ext
  END DO

  DO j=1,ny_ext
    y_ext(j)=y0_ext+(j-1)*dy_ext
  END DO

  CALL xytoll(nx_ext,ny_ext,x_ext,y_ext,lat_ext,lon_ext)

!
!-----------------------------------------------------------------------
!
!  Rotate winds to be relative to true north.
!  The RUCawips data are sent as grid-relative.
!
!-----------------------------------------------------------------------
!
  DO k=1,nz_ext
    CALL uvmptoe(nx_ext,ny_ext,u_ext(1,1,k),v_ext(1,1,k),               &
                 lon_ext,utmp,vtmp)
    u_ext(:,:,k) = utmp(:,:)
    v_ext(:,:,k) = vtmp(:,:)
  END DO

  istatus = 1

!
!-----------------------------------------------------------------------
!
!  Reset map projection to previous values
!
!-----------------------------------------------------------------------
!
  999  CONTINUE
  CALL setmapr(iproj,scale,latnot,trlon)
  CALL setorig(1,x0,y0)

  DEALLOCATE(utmp)
  DEALLOCATE(vtmp)

  RETURN
END SUBROUTINE getnmcrucp236
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE GETNMCRUCN130               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getnmcrucn130(nx_ext,ny_ext,nz_ext,nzsoil_ext,               &
           dir_extd,extdname,extdopt,extdfmt,                           &
           extdinit,extdfcst,julfname,                                  &
           iproj_ext,scale_ext,                                         &
           trlon_ext,latnot_ext,x0_ext,y0_ext,                          &
           lat_ext,lon_ext,zpsoil_ext,                                  &
           p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                      &
           qc_ext,qr_ext,qi_ext,qs_ext,qh_ext,                          &
           tsoil_ext,qsoil_ext,wetcanp_ext,                             &
           trn_ext,psfc_ext,snowdpth_ext,soiltyp_ext,                   &
           t_2m_ext,qv_2m_ext,u_10m_ext,v_10m_ext,rain_ext,             &
           istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads in a GRIB file containing native coordinate RUC2 data
!  (Grid 130) and extracts/converts selected variables
!  for use by EXT2ARPS.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  05/06/2009
!  Based on subroutine GETNMCRUCN236.
!
!  MODIFICATION HISTORY:
!  Keith Brewster 06/29/2011  Added GRIB2 processing for RUC native #130.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    dir_extd      Directory name for external file
!    extdname      Prefix string of external file name
!    extdopt       Option of external data sources
!    extdfmt       Option of external data format

!    extdinit      Initialized time in mm-dd-yyyy:hh:mm:ss format
!    extdfcst      Forecast hour in HHH:MM:SS format
!    julfname      File name in yyjjjhhmm format
!
!  OUTPUT:
!
!    iproj_ext     Map projection number of external data
!    scale_ext     Scale factor of external data
!    trlon_ext     True longitude of external data (degrees E)
!    latnot_ext(2) True latitude(s) of external data (degrees N)
!    x0_ext        x coordinate of origin of external data
!    y0_ext        y coordinate of origin of external data
!    lat_ext       latitude of external data points (degrees N)
!    lon_ext       longitude of external data points (degrees E)
!    p_ext         pressure (Pascal)
!    hgt_ext       height (m)
!    t_ext         temperature (K)
!    qv_ext        specific humidity (kg/kg)
!    u_ext         u wind component (m/s)
!    v_ext         v wind component (m/s)
!    qc_ext        Cloud water mixing ratio (kg/kg)
!    qr_ext        Rain water mixing ratio (kg/kg)
!    qi_ext        Ice mixing ratio (kg/kg)
!    qs_ext        Snow mixing ratio (kg/kg)
!    qh_ext        Hail mixing ratio (kg/kg)
!
!    tsfc_ext      Surface temperature
!    tdeep_ext     Soil temperature
!    wetsfc_ext    Top layer soil moisture (fraction)
!    wetdp_ext     Deep soil moisture (fraction)
!    wetcanp_ext   Water content on canopy
!
!    trn_ext       External terrain (m)
!    psfc_ext      Surface pressure (Pa)
!    snowdpth_ext  Snow depth (m)
!
!    istatus       status indicator
!
!  WORK ARRAYS:
!
!    var_grb3d     Arrays to store the GRIB 3-D variables:
!                  var_grb3d(nxgrb,nygrb,nzgrb,1,1) - pressure (Pa)
!                  var_grb3d(nxgrb,nygrb,nzgrb,2,1) - height (m)
!                  var_grb3d(nxgrb,nygrb,nzgrb,3,1) - Virtual Potential
!                                                     temperature (K)
!                  var_grb3d(nxgrb,nygrb,nzgrb,4,1) - Water vapor
!                                                     mixing ratio
!                                                     (kg/kg)
!                  var_grb3d(nxgrb,nygrb,nzgrb,5,1) - u wind (m/s)
!                  var_grb3d(nxgrb,nygrb,nzgrb,6,1) - v wind (m/s)
!                  var_grb3d(nxgrb,nygrb,nzgrb,7,1)
!                        - cloud water mixing ratio (kg/kg)
!                  var_grb3d(nxgrb,nygrb,nzgrb,8,1)
!                        - rain water mixing ratio (kg/kg)
!                  var_grb3d(nxgrb,nygrb,nzgrb,9,1)
!                        - ice mixing ratio (kg/kg)
!                  var_grb3d(nxgrb,nygrb,nzgrb,10,1)
!                        - snow mixing ratio (kg/kg)
!                  var_grb3d(nxgrb,nygrb,nzgrb,11,1)
!                        - graupel mixing ratio (kg/kg)
!                  var_grb3d(nxgrb,nygrb,nzgrb,1,2) - soil temp. (K)
!                  var_grb3d(nxgrb,nygrb,nzgrb,2,2) - soil moist.
!                                                    (fraction)
!
!    var_grb2d     Arrays to store the GRIB 2-D variables:
!                  var_grb2d(nxgrb,nygrb,1) - Canopy Water
!                                             (kg/m**2)
!                  var_grb2d(nxgrb,nygrb,2) - Snow depth (m)
!                  var_grb2d(nxgrb,nygrb,3) - Surface soil temperature (K)
!                  var_grb2d(nxgrb,nygrb,4) - Surface soil moisture
!                                             (fraction)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'gribcst.inc'

  CHARACTER (LEN=*) :: dir_extd
  CHARACTER (LEN=*) :: extdname

  INTEGER :: extdopt
  INTEGER :: extdfmt

  CHARACTER (LEN=19) :: extdinit
  CHARACTER (LEN=9) :: extdfcst
  CHARACTER (LEN=9) :: julfname
!
!
!-----------------------------------------------------------------------
!
!  External grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj_ext
  REAL :: scale_ext,trlon_ext
  REAL :: latnot_ext(2)
  REAL :: x0_ext,y0_ext
  REAL :: dx_ext,dy_ext
!
!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------
!
  INTEGER, INTENT(IN) :: nx_ext, ny_ext, nz_ext, nzsoil_ext

  REAL :: lat_ext(nx_ext,ny_ext)
  REAL :: lon_ext(nx_ext,ny_ext)
  REAL :: p_ext  (nx_ext,ny_ext,nz_ext)   ! Pressure (Pascals)
  REAL :: hgt_ext(nx_ext,ny_ext,nz_ext)   ! Height (m)
  REAL :: t_ext  (nx_ext,ny_ext,nz_ext)   ! Temperature (K)
  REAL :: qv_ext (nx_ext,ny_ext,nz_ext)   ! Specific humidity (kg/kg)
  REAL :: u_ext  (nx_ext,ny_ext,nz_ext)   ! Eastward wind component
  REAL :: v_ext  (nx_ext,ny_ext,nz_ext)   ! Northward wind component
  REAL :: qc_ext (nx_ext,ny_ext,nz_ext)   ! Cloud H2O mixing ratio (kg/kg)
  REAL :: qr_ext (nx_ext,ny_ext,nz_ext)   ! Rain  H2O mixing ratio (kg/kg)
  REAL :: qi_ext (nx_ext,ny_ext,nz_ext)   ! Ice   mixing ratio (kg/kg)
  REAL :: qs_ext (nx_ext,ny_ext,nz_ext)   ! Snow  mixing ratio (kg/kg)
  REAL :: qh_ext (nx_ext,ny_ext,nz_ext)   ! Hail  mixing ratio (kg/kg)

  REAL :: tsoil_ext  (nx_ext,ny_ext,nzsoil_ext)      ! Temperature at surface (K)
  REAL :: qsoil_ext  (nx_ext,ny_ext,nzsoil_ext)      ! soil moisture
  REAL :: wetcanp_ext(nx_ext,ny_ext)      ! Canopy water amount

  REAL :: trn_ext    (nx_ext,ny_ext)      ! Geometrical heights
  REAL :: psfc_ext   (nx_ext,ny_ext)      ! Surface pressure (Pa)
  REAL :: snowdpth_ext (nx_ext,ny_ext)    ! Snow depth (m)

  REAL,    INTENT(OUT) :: zpsoil_ext(nx_ext,ny_ext,nzsoil_ext)
  INTEGER, INTENT(OUT) :: soiltyp_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: t_2m_ext (nx_ext,ny_ext), qv_2m_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: u_10m_ext(nx_ext,ny_ext), v_10m_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: rain_ext(nx_ext,ny_ext)

  INTEGER, INTENT(OUT) :: istatus
!
!-----------------------------------------------------------------------
!
!  Other  external variable arrays
!
!-----------------------------------------------------------------------
!
  REAL :: x_ext(nx_ext)
  REAL :: y_ext(ny_ext)
!
!-----------------------------------------------------------------------
!
!  Work arrays for storing grib data
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: var_grb2d(:,:,:,:)   ! GRIB variables
  REAL, ALLOCATABLE :: var_grb3d(:,:,:,:,:) ! GRIB 3-D variables
  INTEGER, ALLOCATABLE :: var_lev3d(:,:,:)  ! Levels (hybrid) for
                                            ! each 3-D variable
  REAL, ALLOCATABLE :: rcdata(:)            ! temporary data array
!
!-----------------------------------------------------------------------
!
!  Original grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj
  REAL    :: scale,trlon,x0,y0
  REAL    :: latnot(2)
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: gribfile
  CHARACTER (LEN=14) :: gribtime
  INTEGER :: i,j,k
  INTEGER :: grbflen, grbtlen

  INTEGER :: m,n,nz1,max_nr2d,max_nr3d

  INTEGER, PARAMETER :: nzsoilin_ext = 6
  REAL,    PARAMETER :: soildepth_ext(nzsoilin_ext) = (/0.0, 0.05, 0.2, 0.4, 1.6, 3.0/)

  REAL :: pilev

  REAL :: tv_ext, tvc_ext
  REAL :: rovcp_p, cpd_p, g0_p, rd_p

  INTEGER :: chklev, lvscan, kk, jj

  REAL :: tema, temb

  REAL :: a,b

  INTEGER :: iret             ! Return flag

  REAL, ALLOCATABLE :: utmp(:,:), vtmp(:,:)
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  GRIB grid information
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=42) :: gridesc ! Grid description

  INTEGER :: iproj_grb    ! Map projection indicator
  INTEGER :: gthin        ! Indicator of whether the grid is "thinned"

  INTEGER :: ni_grb       ! Number of points along x-axis
  INTEGER :: nj_grb       ! Number of points along y-axis
  INTEGER :: np_grb       ! Total number of horizontal grid points

  INTEGER :: nk_grb       ! Number of vertical parameters
  REAL :: zk_grb(nz_ext)  ! Vertical coordinate parameters

  INTEGER :: npeq         ! Number of lat circles from pole to equator
  INTEGER :: nit(nz_ext)  ! Number of x-points for thinned grid

  REAL :: pi_grb          ! x-coordinate of pole point
  REAL :: pj_grb          ! y-coordinate of pole point
  INTEGER :: ipole        ! Projection center flag

  REAL :: di_grb          ! x-direction increment or grid length
  REAL :: dj_grb          ! y-direction increment or grid length

  REAL :: latsw           ! Latitude  of South West corner point
  REAL :: lonsw           ! Longitude of South West corner point
  REAL :: latne           ! Latitude  of North East corner point
  REAL :: lonne           ! Longitude of North East corner point

  REAL :: lattru1         ! Latitude (1st) at which projection is true
  REAL :: lattru2         ! Latitude (2nd) at which projection is true
  REAL :: lontrue         ! Longitude      at which projection is true

  REAL :: latrot          ! Latitude  of southern pole of rotation
  REAL :: lonrot          ! Longitude of southern pole of rotation
  REAL :: angrot          ! Angle of rotation

  REAL :: latstr          ! Latitude  of the pole of stretching
  REAL :: lonstr          ! Longitude of the pole of stretching
  REAL :: facstr          ! Stretching factor

  INTEGER :: scanmode     ! Scanning indicator
  INTEGER :: iscan        ! x-direction   scanning indicator
  INTEGER :: jscan        ! y-direction   scanning indicator
  INTEGER :: kscan        ! FORTRAN index scanning indicator

  INTEGER :: ires         ! Resolution direction increments indicator
  INTEGER :: iearth       ! Earth shape indicator: spherical or oblate?
  INTEGER :: icomp        ! (u,v) components decomposition indicator

  INTEGER :: jpenta       ! J-Pentagonal resolution parameter
  INTEGER :: kpenta       ! K-Pentagonal resolution parameter
  INTEGER :: mpenta       ! M-Pentagonal resolution parameter
  INTEGER :: ispect       ! Spectral representation type
  INTEGER :: icoeff       ! Spectral coefficient storage mode

  REAL :: xp_grb          ! X coordinate of sub-satellite point
  REAL :: yp_grb          ! Y coordinate of sub-satellite point
  REAL :: xo_grb          ! X coordinate of image sector origin
  REAL :: yo_grb          ! Y coordinate of image sector origin
  REAL :: zo_grb          ! Camera altitude from center of Earth

  INTEGER :: igrbfmt
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ALLOCATE(utmp(nx_ext,ny_ext))
  ALLOCATE(vtmp(nx_ext,ny_ext))

!---------------------------------------------------------------------
!
!  Define soil depths
!
!---------------------------------------------------------------------

  IF(soilmodel_option == 1 .AND. nzsoil_ext /= 2) THEN
    WRITE(6,'(2a,I3,a/2a,I2,a/,a,a)')                                   &
               'NCEP RUC grid #130 only provides 5 soil layers, However,',&
               ' You are trying to extract ',nzsoil_ext, ' layers.',    &
               ' for ARPS two-layer force-restore model ',              &
               '(soilmodel_option = ',soilmodel_option,')',             &
               ' Please check the code ext2arps.f90 for RUCN130 grid.',  &
               ' Terminating ...'
    CALL arpsstop('Wrong number soil layers.',1)
  ELSE IF( soilmodel_option == 2 .AND. nzsoil_ext /= nzsoilin_ext+1) THEN
    WRITE(6,'(2a,I3,a/2a,I2,a/,a,a)')                                   &
               'NCEP RUC grid #130 only provides 5 soil layers, However,',&
               ' You are trying to extract ',nzsoil_ext, ' layers.',    &
               ' for ARPS multi-layer OUSoil model ',                   &
               '(soilmodel_option = ',soilmodel_option,')',             &
               ' Please check the code ext2arps.f90 for NAM218 grid.',  &
               ' Terminating ...'
    CALL arpsstop('Wrong number of soil layers.',1)
  END IF

  IF(soilmodel_option == 1) THEN
    DO j=1,ny_ext
      DO i=1,nx_ext
        zpsoil_ext(i,j,1) = 0.05  ! any values should work because
        zpsoil_ext(i,j,2) = 1.0   ! the ARPS system will ignore these values
      END DO
    END DO
  ELSE
    DO j=1,ny_ext
      DO i=1,nx_ext
        DO k=2,nzsoil_ext
          zpsoil_ext(i,j,k) = - (soildepth_ext(k)-soildepth_ext(k-1))/2 - soildepth_ext(k-1)
                                  ! The middle point for each layer
        END DO
        zpsoil_ext(i,j,1) = 0.0
      END DO
    END DO
  END IF

!-----------------------------------------------------------------------
!
! Start processing and dispatch job based on whether it is GRIB or GRIB2
!
!-----------------------------------------------------------------------

  CALL getmapr(iproj,scale,latnot,trlon,x0,y0)

  CALL getgrbfname(extdopt,extdfmt,dir_extd,extdname,                   &
                   extdinit,extdfcst,'    ',                            &
                   gribfile,grbflen,gribtime,grbtlen,istatus)
  IF (istatus /=0 ) RETURN

  CALL chkgrb(gribfile(1:grbflen),grbflen,igrbfmt,istatus)

  IF (igrbfmt == 2) THEN

    WRITE(6,'(1x,2a,I4)') 'GRIB2 file to be read is ',gribfile(1:grbflen),grbflen

    CALL getruc130n_grb2(nx_ext,ny_ext,nz_ext,nzsoil_ext,nzsoilin_ext, &
           gribfile,gribtime,grbflen,soildepth_ext,                    &
           dx_ext,dy_ext,                                              &
           iproj_ext,scale_ext,trlon_ext,latnot_ext,latsw,lonsw,       &
           p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                     &
           qc_ext,qr_ext,qi_ext,qs_ext,qh_ext,                         &
           tsoil_ext,qsoil_ext,wetcanp_ext,                            &
           snowdpth_ext,trn_ext,psfc_ext,t_2m_ext,qv_2m_ext,           &
           u_10m_ext,v_10m_ext,rain_ext,istatus)
!
    IF (istatus /= 0) RETURN

  ELSE

    ALLOCATE(var_grb2d(nx_ext,ny_ext,n2dvs,n2dlvt))
    ALLOCATE(var_grb3d(nx_ext,ny_ext,nz_ext,n3dvs,n3dlvt))
    ALLOCATE(rcdata(nx_ext*ny_ext))
    ALLOCATE(var_lev3d(nz_ext,n3dvs,n3dlvt))

!
!-----------------------------------------------------------------------
!
!  RDNMCGRB reads NMC GRIB data
!
!-----------------------------------------------------------------------
!
    gridtyp   = rucn130grid
    mproj_grb = rucn130proj

    n2dvars  = rucn130nvs2d
    n2dlvtps = rucn130nlvt2d

    DO m=1,n2dlvtps
      DO n=1,n2dvars
        var_id2d(n,m) = rucn130var_id2d(n,m)
      END DO
      levtyp2d(m) = rucn130levs2d(m)
    END DO

    n3dvars  = rucn130nvs3d
    n3dlvtps = rucn130nlvt3d

    DO m=1,n3dlvtps
      DO n=1,n3dvars
        var_id3d(n,m) = rucn130var_id3d(n,m)
      END DO
      levtyp3d(m) = rucn130levs3d(m)
    END DO

    CALL rdnmcgrb(nx_ext,ny_ext,nz_ext,gribfile,grbflen, gribtime,        &
                  gridesc, iproj_grb, gthin,                              &
                  ni_grb,nj_grb,np_grb, nk_grb,zk_grb, npeq,nit,          &
                  pi_grb,pj_grb,ipole, di_grb,dj_grb,                     &
                  latsw,lonsw, latne,lonne,                               &
                  latrot,lonrot,angrot,                                   &
                  latstr,lonstr,facstr,                                   &
                  lattru1,lattru2,lontrue,                                &
                  scanmode, iscan,jscan,kscan,                            &
                  ires,iearth,icomp,                                      &
                  jpenta,kpenta,mpenta,ispect,icoeff,                     &
                  xp_grb,yp_grb, xo_grb,yo_grb,zo_grb,                    &
                  rcdata,var_grb2d,var_grb3d,var_lev3d,lvldbg,iret)

    IF (iret /= 0)  THEN
      ISTATUS = -888
      GOTO 999
    END IF

    max_nr2d = 0
    DO n=1,n2dvars
      DO m=1,n2dlvtps
        max_nr2d = MAX( max_nr2d, var_nr2d(n,m) )
      END DO
    END DO

    max_nr3d = 0
    DO n=1,n3dvars
      max_nr3d = MAX( max_nr3d, var_nr3d(n,1) )
    END DO

    IF ( max_nr3d == 0 ) THEN
      WRITE (6,'(a)')                                                     &
          'No 3-D variables was found in the GRIB file',                  &
          'Program stopped in GETNMCRUC.'
  !   STOP
      ISTATUS = -888
      GOTO 999
    END IF

    IF ( max_nr2d == 0 ) THEN
      WRITE (6,'(a)')                                                     &
          'No 2-D variables was found in the GRIB file'
    END IF

!   write (6,'(/a7,2x,6(i7))')
!     :      'Lev\\VID',(var_id3d(n,1),n=1,n3dvars)

!   DO 60 k=1,max_nr3d
!      write (6,'(i5,4x,6(i7))')
!     :      k,(var_lev3d(k,n,1),n=1,n3dvars)
!   60    CONTINUE

    DO k=1,max_nr3d
      DO n=2,n3dvars
        IF ( var_lev3d(k,1,1) /= var_lev3d(k,n,1) ) THEN
          WRITE (6,'(a)')                                                 &
              'Variables were not at the same level.',                    &
              'Program stopped in GETNMCRUC.'
          WRITE(6,*)                                                      &
              'var_lev3d(k,1,1) = ',var_lev3d(k,1,1),                     &
              'var_lev3d(k,n,1) = ',var_lev3d(k,n,1)
  !       STOP
          ISTATUS = -888
          GOTO 999
        END IF
      END DO
    END DO

    IF ( iproj_grb == 5 .AND. ipole == 0 ) THEN      ! Center in NP
      iproj_ext = 1
    ELSE IF ( iproj_grb == 5 .AND. ipole == 1 ) THEN  ! Center in SP
      iproj_ext = -1
    ELSE IF ( iproj_grb == 3 .AND. ipole == 0 ) THEN  ! Center in NP
      iproj_ext = 2
    ELSE IF ( iproj_grb == 3 .AND. ipole == 1 ) THEN  ! Center in SP
      iproj_ext = -2
    ELSE IF ( iproj_grb == 1 ) THEN
      iproj_ext = 3
    ELSE IF ( iproj_grb == 0 ) THEN
      iproj_ext = 4
    ELSE
      WRITE (6,'(a)')                                                     &
          'Unknown map projection. Set to non-projection.'
      iproj_ext = 0
    END IF

    scale_ext = 1.0
    latnot_ext(1) = lattru1
    latnot_ext(2) = lattru2
    trlon_ext = lontrue

    dx_ext = di_grb
    dy_ext = dj_grb
!
!-----------------------------------------------------------------------
!
!  Retrieve 2-D variables
!
!-----------------------------------------------------------------------
!
    DO j=1,ny_ext
      DO i=1,nx_ext

!        trn_ext(i,j) = -999.0
!        psfc_ext   (i,j) = -999.0

        IF ( var_nr2d(1,1) == 0 ) THEN
          rain_ext(i,j) = -999.0
        ELSE
          rain_ext(i,j) = var_grb2d(i,j,1,1) + var_grb2d(i,j,2,1)  ! in kg/m**2
        END IF

        IF ( var_nr2d(3,1) == 0 ) THEN
          snowdpth_ext(i,j) = -999
        ELSE
          snowdpth_ext(i,j) = var_grb2d(i,j,3,1)  ! in meters
        END IF

        wetcanp_ext(i,j) = var_grb2d(i,j,7,1)*1.e-3  ! in meters

        IF ( var_nr2d(1,2) == 0 ) THEN
          t_2m_ext(i,j) = -999.0
        ELSE
          t_2m_ext(i,j) = var_grb2d(i,j,1,2)
        END IF

        IF ( var_nr2d(2,2) == 0 ) THEN
          qv_2m_ext(i,j) = -999.0
        ELSE
          qv_2m_ext(i,j) = var_grb2d(i,j,2,2)
        END IF

        IF ( var_nr2d(3,2) == 0 ) THEN
          u_10m_ext(i,j) = -999.0
        ELSE
          u_10m_ext(i,j) = var_grb2d(i,j,3,2)
        END IF

        IF ( var_nr2d(4,2) == 0 ) THEN
          v_10m_ext(i,j) = -999.0
        ELSE
          v_10m_ext(i,j) = var_grb2d(i,j,4,2)
        END IF

      END DO
    END DO

!
!-----------------------------------------------------------------------
!
!  Retrieve 3-D variables
!
!-----------------------------------------------------------------------
!

    cpd_p   = 1004.686    ! cp in RUC
    rovcp_p = 0.285714    ! rd/cp used in RUC
    g0_p    = 9.80665     ! gravity in RUC

    nz1 = MIN(var_nr3d(1,1),nz_ext)

    IF ( var_lev3d(1,1,1) < var_lev3d(nz1,1,1) ) THEN  ! 1st level at sfc
      chklev = 1
      lvscan = 0
    ELSE
      chklev = -1
      lvscan = nz1+1
    END IF

    DO k=1,nz1
      kk = chklev * k + lvscan
      DO j=1,ny_ext
        DO i=1,nx_ext

          p_ext(i,j,kk) = var_grb3d(i,j,k,1,1)       ! Pressure (Pa)

          hgt_ext(i,j,kk) = var_grb3d(i,j,k,2,1)     ! Height (m)
          u_ext(i,j,kk) = var_grb3d(i,j,k,5,1)       ! u wind (m/s)
          v_ext(i,j,kk) = var_grb3d(i,j,k,6,1)       ! v wind (m/s)

          a = REAL(100000)/var_grb3d(i,j,k,1,1)
          a = a**rovcp_p
          tvc_ext = var_grb3d(i,j,k,3,1)/a           ! Virtual Temperature
          b = 0.61*var_grb3d(i,j,k,4,1)
          b = REAL(1) + b
          t_ext(i,j,kk) = tvc_ext/b                  ! Temperature (K)

          a = var_grb3d(i,j,k,4,1)*var_grb3d(i,j,k,1,1)
          a = a/(0.622 - var_grb3d(i,j,k,4,1))
          qv_ext(i,j,kk) = a*0.622/var_grb3d(i,j,k,1,1) ! Specific Humidity
!
!-----------------------------------------------------------------------
!
!      Retrieve hydrometeor data.
!
!-----------------------------------------------------------------------
!
          IF (var_nr3d(7,1) == 0) THEN
            qc_ext(i,j,kk) = -999.
          ELSE
            qc_ext(i,j,kk) = var_grb3d(i,j,k,7,1)
          END IF

          IF (var_nr3d(8,1) == 0) THEN
            qr_ext(i,j,kk) = -999.
          ELSE
            qr_ext(i,j,kk) = var_grb3d(i,j,k,8,1)
          END IF

          IF (var_nr3d(9,1) == 0) THEN
            qi_ext(i,j,kk) = -999.
          ELSE
            qi_ext(i,j,kk) = var_grb3d(i,j,k,9,1)
          END IF

          IF (var_nr3d(10,1) == 0) THEN
            qs_ext(i,j,kk) = -999.
          ELSE
            qs_ext(i,j,kk) = var_grb3d(i,j,k,10,1)
          END IF

          IF (var_nr3d(11,1) == 0) THEN
            qh_ext(i,j,kk) = -999.
          ELSE
            qh_ext(i,j,kk) = var_grb3d(i,j,k,11,1)
          END IF

        END DO
      END DO
    END DO

    trn_ext (:,:) = hgt_ext(:,:,1)
    psfc_ext(:,:) = p_ext(:,:,1)
!
!-----------------------------------------------------------------------
!
!  Retrieve soil variables
!
!-----------------------------------------------------------------------
!
    IF (soilmodel_option == 1) THEN ! Old ARPS Force-Restore Soil Model
      DO j = 1, ny_ext
        DO i = 1, nx_ext
          tsoil_ext(i,j,1) = var_grb2d(i,j,5,1)
          IF ( nint(var_grb2d(i,j,4,1)) == 1 ) THEN  ! soil temp over land
            tsoil_ext(i,j,2) = 0.1 * var_grb3d(i,j,1,1,2)  & !   5cm
                             + 0.2 * var_grb3d(i,j,2,1,2)  & !  20cm
                             + 0.4 * var_grb3d(i,j,3,1,2)  & !  40cm
                             + 0.3 * var_grb3d(i,j,4,1,2)  ! 160cm

            soiltyp_ext(i,j) = 0
          ELSE
            tsoil_ext(i,j,2) = var_grb2d(i,j,5,1)
            soiltyp_ext(i,j) = 13    ! Set soil type to water
          END IF

          qsoil_ext(i,j,1) = var_grb2d(i,j,6,1)
          qsoil_ext(i,j,2) = 0.1 * var_grb3d(i,j,1,2,2)  & !   5cm
                           + 0.2 * var_grb3d(i,j,2,2,2)  & !  20cm
                           + 0.4 * var_grb3d(i,j,3,2,2)  & !  40cm
                           + 0.3 * var_grb3d(i,j,4,2,2)  ! 160cm

        END DO
      END DO
    ELSE    ! OU Soil Model

      DO j = 1, ny_ext
        DO i = 1, nx_ext

          IF ( nint(var_grb2d(i,j,4,1)) == 1 ) THEN  ! Land

            tsoil_ext (i,j,1) = var_grb2d(i,j,5,1)   ! Ground temperature
            qsoil_ext (i,j,1) = var_grb2d(i,j,6,1)

            DO k=2,nzsoil_ext
              ! "TSOIL" in GRIB is below ground, treated as separate
              ! variable from ground temperature.
              tsoil_ext (i,j,k) = var_grb3d(i,j,k-1,1,2)
              qsoil_ext (i,j,k) = var_grb3d(i,j,k-1,2,2)
            END DO

            soiltyp_ext(i,j) = 0

          ELSE  ! Water

            DO k=1,nzsoil_ext
              tsoil_ext (i,j,k) = var_grb2d(i,j,5,1) ! Water temperature?
              qsoil_ext (i,j,k) = 1.                 ! 100% water
            END DO
            soiltyp_ext(i,j) = 13

          END IF ! Land or water?

          END DO
        END DO

    END IF

    DEALLOCATE(var_grb2d,var_grb3d,rcdata,var_lev3d)
  END IF

  CALL setmapr(iproj_ext,scale_ext,latnot_ext,trlon_ext)
  CALL lltoxy(1,1,latsw,lonsw,x0_ext,y0_ext)

  DO i=1,nx_ext
    x_ext(i)=x0_ext+(i-1)*dx_ext
  END DO

  DO j=1,ny_ext
    y_ext(j)=y0_ext+(j-1)*dy_ext
  END DO

  CALL xytoll(nx_ext,ny_ext,x_ext,y_ext,lat_ext,lon_ext)
!
!-----------------------------------------------------------------------
!
!  Rotate winds to be relative to true north.
!  The RUC data are sent as grid-relative.
!
!-----------------------------------------------------------------------
!
  DO k=1,nz_ext
    CALL uvmptoe(nx_ext,ny_ext,u_ext(1,1,k),v_ext(1,1,k),               &
                 lon_ext,utmp,vtmp)
    u_ext(:,:,k) = utmp(:,:)
    v_ext(:,:,k) = vtmp(:,:)
  END DO

  istatus = 1

!
!-----------------------------------------------------------------------
!
!  Reset map projection to previous values
!
!-----------------------------------------------------------------------
!
  999  CONTINUE
  CALL setmapr(iproj,scale,latnot,trlon)
  CALL setorig(1,x0,y0)

  DEALLOCATE(utmp)
  DEALLOCATE(vtmp)

  RETURN
END SUBROUTINE getnmcrucn130
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETNCEPAVN3                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE getncepavn3(nx_ext,ny_ext,nz_ext,nzsoil_ext,                 &
           dir_extd,extdname,extdopt,extdfmt,nofixdim,lon_0_360,        &
           extdinit,extdfcst,julfname,                                  &
           iproj_ext,scale_ext,                                         &
           trlon_ext,latnot_ext,x0_ext,y0_ext,                          &
           lat_ext,lon_ext,zpsoil_ext,                                  &
           p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                      &
           qc_ext,                                                      &
           tsoil_ext,qsoil_ext,wetcanp_ext,                             &
           snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext,                   &
           t_2m_ext, qv_2m_ext, u_10m_ext, v_10m_ext,                   &
           rain_ext,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads a NCEP AVN GRIB (Grid #3, 1x1 degree) data file for
!  processing by ext2arps.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Donghai Wang and Yuhe Liu
!  05/19/1999
!
!  MODIFICATION HISTORY:
!
!  03/23/2000 (Donghai Wang)
!  Fixed some bugs and modified the code for ARPS4.5.1 version.
!
!  2004/01/15 (F. KONG)
!  Modified to read in from grid 3 four more near surface variables:
!  2-m temperature and humidity, and 10-m wind (u, v)
!
!  08/11/2005 (Y. Wang)
!  Modified to handle cases with Greenwich Meridian in the domain. Note
!  that a new parameter, lon_0_360 was passed in and an integer array
!  ii(nx_ext) was decalared.
!
!  Modified to work with data over subregion instead of global data. A
!  new parameter nofixdim was passed in.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    dir_extd      Directory name for external file
!    extdname      Prefix string of external file name
!    extdopt       Option of external data sources
!    extdfmt       Option of external data format
!
!    extdinit      Initialized time in mm-dd-yyyy:hh:mm:ss format
!    extdfcst      Forecast hour in HHH:MM:SS format
!    julfname      File name in yyjjjhhmm format
!
!  OUTPUT:
!
!    iproj_ext     Map projection number of external data
!    scale_ext     Scale factor of external data
!    trlon_ext     True longitude of external data (degrees E)
!    latnot_ext(2) True latitude(s) of external data (degrees N)
!    x0_ext        x coordinate of origin of external data
!    y0_ext        y coordinate of origin of external data
!    lat_ext       latitude of external data points (degrees N)
!    lon_ext       longitude of external data points (degrees E)
!    p_ext         pressure (Pascal)
!    hgt_ext       height (m)
!    t_ext         temperature (K)
!    qv_ext        specific humidity (kg/kg)
!    u_ext         u wind component (m/s)
!    v_ext         v wind component (m/s)
!    qc_ext        Cloud water mixing ratio (kg/kg)
!    qr_ext        Rain water mixing ratio (kg/kg)
!    qi_ext        Ice mixing ratio (kg/kg)
!    qs_ext        Snow mixing ratio (kg/kg)
!    qh_ext        Hail mixing ratio (kg/kg)
!
!    tsoil_ext     Soil temperature
!    qsoil_ext     Soil moisture
!    wetcanp_ext   Water content on canopy
!
!    trn_ext       External terrain (m)
!    psfc_ext      Surface pressure (Pa)
!
!    T_2m_ext      Temperature at 2m AGL
!    qv_2m_ext     Specific Humidity at 2m AGL
!    U_10m_ext     U at 10m AGL
!    V_10m_ext     V at 10m AGL
!
!    istatus       status indicator
!
!  WORK ARRAYS:
!
!    var_grb3d     Arrays to store the GRIB 3-D variables:
!                  var_grb3d(nxgrb,nygrb,nzgrb,1,1) - Temperature (K)
!                  var_grb3d(nxgrb,nygrb,nzgrb,2,1) - Specific humidity
!                                                     (kg/kg)
!                  var_grb3d(nxgrb,nygrb,nzgrb,3,1) - u wind (m/s)
!                  var_grb3d(nxgrb,nygrb,nzgrb,4,1) - v wind (m/s)
!                  var_grb3d(nxgrb,nygrb,nzgrb,5,1) - Geopotential
!                                                     height (gpm)
!                  var_grb3d(nxgrb,nygrb,nzgrb,6,1) - Pressure vertical
!                                                     velocity (Pa/s)
!                                                     (if applied)
!                  var_grb3d(nxgrb,nygrb,nzgrb,1,2) - soil temp. (K)
!                  var_grb3d(nxgrb,nygrb,nzgrb,2,2) - vol. soil moist.
!                                                     (m**3/m**3)
!
!    var_grb2d     Arrays to store the GRIB 2-D variables:
!                  var_grb2d(nxgrb,nygrb,1) - Surface pressure (Pa)
!                  var_grb2d(nxgrb,nygrb,2) - Geopotential height (gpm)
!                  var_grb2d(nxgrb,nygrb,3) - Surface temperature (K)
!                  var_grb2d(nxgrb,nygrb,4) - Plant canopy surface
!                                             water (kg/m**2)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'gribcst.inc'

  CHARACTER (LEN=*) :: dir_extd
  CHARACTER (LEN=*) :: extdname

  INTEGER, INTENT(IN) :: extdopt, extdfmt
  INTEGER, INTENT(IN) :: nofixdim
  LOGICAL, INTENT(IN) :: lon_0_360

  CHARACTER (LEN=19) :: extdinit
  CHARACTER (LEN=9)  :: extdfcst
  CHARACTER (LEN=9)  :: julfname
!
!-----------------------------------------------------------------------
!
!  External grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj_ext
  REAL    :: scale_ext,trlon_ext
  REAL    :: latnot_ext(2)
  REAL    :: x0_ext,y0_ext
  REAL    :: dx_ext,dy_ext
!
!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx_ext,ny_ext,nz_ext,nzsoil_ext

  REAL :: lat_ext(nx_ext,ny_ext)
  REAL :: lon_ext(nx_ext,ny_ext)
  REAL :: zpsoil_ext(nx_ext,ny_ext,nzsoil_ext)
  REAL :: p_ext  (nx_ext,ny_ext,nz_ext)   ! Pressure (Pascals)
  REAL :: hgt_ext(nx_ext,ny_ext,nz_ext)   ! Height (m)
  REAL :: t_ext  (nx_ext,ny_ext,nz_ext)   ! Temperature (K)
  REAL :: qv_ext (nx_ext,ny_ext,nz_ext)   ! Specific humidity (kg/kg)
  REAL :: u_ext  (nx_ext,ny_ext,nz_ext)   ! Eastward wind component
  REAL :: v_ext  (nx_ext,ny_ext,nz_ext)   ! Northward wind component
  REAL :: qc_ext (nx_ext,ny_ext,nz_ext)   ! Cloud H2O mixing ratio (kg/kg)
!  REAL :: qr_ext (nx_ext,ny_ext,nz_ext)   ! Rain  H2O mixing ratio (kg/kg)
!  REAL :: qi_ext (nx_ext,ny_ext,nz_ext)   ! Ice   mixing ratio (kg/kg)
!  REAL :: qs_ext (nx_ext,ny_ext,nz_ext)   ! Snow  mixing ratio (kg/kg)
!  REAL :: qh_ext (nx_ext,ny_ext,nz_ext)   ! Hail  mixing ratio (kg/kg)

  REAL :: tsoil_ext  (nx_ext,ny_ext,nzsoil_ext)      ! Temperature at surface (K)
  REAL :: qsoil_ext  (nx_ext,ny_ext,nzsoil_ext)      ! Deep soil temperature (K)
  REAL :: wetcanp_ext(nx_ext,ny_ext)                 ! Canopy water amount
  REAL :: snowdpth_ext(nx_ext,ny_ext)                ! Snow depth (m)

  REAL :: trn_ext    (nx_ext,ny_ext)      ! External terrain (m)
  REAL :: psfc_ext   (nx_ext,ny_ext)      ! Surface pressure (Pa)
  REAL :: t_2m_ext (nx_ext,ny_ext)
  REAL :: qv_2m_ext(nx_ext,ny_ext)
  REAL :: u_10m_ext(nx_ext,ny_ext)
  REAL :: v_10m_ext(nx_ext,ny_ext)
  REAL :: rain_ext (nx_ext,ny_ext)

  INTEGER :: soiltyp_ext(nx_ext,ny_ext)
!
!-----------------------------------------------------------------------
!
!  Other  external variable arrays
!
!-----------------------------------------------------------------------
!
!  REAL :: x_ext(nx_ext)
!  REAL :: y_ext(ny_ext)

  INTEGER :: ii(nx_ext)

  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Work arrays for storing grib data
!
!-----------------------------------------------------------------------
!
  REAL,    ALLOCATABLE :: var_grb2d(:,:,:,:)   ! GRIB variables
  REAL,    ALLOCATABLE :: var_grb3d(:,:,:,:,:) ! GRIB 3-D variables
  INTEGER, ALLOCATABLE :: var_lev3d(:,:,:)     ! Levels (hybrid) for
                                               ! each 3-D variable
  REAL,    ALLOCATABLE :: rcdata(:)            ! temporary data array
!
!-----------------------------------------------------------------------
!
!  Original grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj
  REAL    :: scale,trlon,x0,y0
  REAL    :: latnot(2)
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: gribfile
  CHARACTER (LEN=14)  :: gribtime
  INTEGER :: i,j,k,kk
  INTEGER :: grbflen, grbtlen

  INTEGER :: igrbfmt

  INTEGER :: m,n,nz1,max_nr2d,max_nr3d,min_nr3d,nz2

  REAL    :: govrd

  INTEGER :: chklev, lvscan

  INTEGER :: iret             ! Return flag

!  REAL, PARAMETER :: soildepth_ext(4) = (/0.0, 0.05, 1.05/)
!                             ! Surface and two soil layers  0-10cm & 10-200cm
  REAL, PARAMETER :: soildepth_ext(4) = (/0.0, 0.1, 0.4, 1.0/)

  REAL    :: rtmp
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  GRIB grid information
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=42) :: gridesc ! Grid description

  INTEGER :: iproj_grb    ! Map projection indicator
  INTEGER :: gthin        ! Indicator of whether the grid is "thinned"

  INTEGER :: ni_grb       ! Number of points along x-axis
  INTEGER :: nj_grb       ! Number of points along y-axis
  INTEGER :: np_grb       ! Total number of horizontal grid points

  INTEGER :: nk_grb       ! Number of vertical parameters
  REAL :: zk_grb(nz_ext)  ! Vertical coordinate parameters

  INTEGER :: npeq         ! Number of lat circles from pole to equator
  INTEGER :: nit(nz_ext)  ! Number of x-points for thinned grid

  REAL :: pi_grb          ! x-coordinate of pole point
  REAL :: pj_grb          ! y-coordinate of pole point
  INTEGER :: ipole        ! Projection center flag

  REAL :: di_grb          ! x-direction increment or grid length
  REAL :: dj_grb          ! y-direction increment or grid length

  REAL :: latsw           ! Latitude  of South West corner point
  REAL :: lonsw           ! Longitude of South West corner point
  REAL :: latne           ! Latitude  of North East corner point
  REAL :: lonne           ! Longitude of North East corner point

  REAL :: lattru1         ! Latitude (1st) at which projection is true
  REAL :: lattru2         ! Latitude (2nd) at which projection is true
  REAL :: lontrue         ! Longitude      at which projection is true

  REAL :: latrot          ! Latitude  of southern pole of rotation
  REAL :: lonrot          ! Longitude of southern pole of rotation
  REAL :: angrot          ! Angle of rotation

  REAL :: latstr          ! Latitude  of the pole of stretching
  REAL :: lonstr          ! Longitude of the pole of stretching
  REAL :: facstr          ! Stretching factor

  INTEGER :: scanmode     ! Scanning indicator
  INTEGER :: iscan        ! x-direction   scanning indicator
  INTEGER :: jscan        ! y-direction   scanning indicator
  INTEGER :: kscan        ! FORTRAN index scanning indicator

  INTEGER :: ires         ! Resolution direction increments indicator
  INTEGER :: iearth       ! Earth shape indicator: spherical or oblate?
  INTEGER :: icomp        ! (u,v) components decomposition indicator

  INTEGER :: jpenta       ! J-Pentagonal resolution parameter
  INTEGER :: kpenta       ! K-Pentagonal resolution parameter
  INTEGER :: mpenta       ! M-Pentagonal resolution parameter
  INTEGER :: ispect       ! Spectral representation type
  INTEGER :: icoeff       ! Spectral coefficient storage mode

  REAL :: xp_grb          ! X coordinate of sub-satellite point
  REAL :: yp_grb          ! Y coordinate of sub-satellite point
  REAL :: xo_grb          ! X coordinate of image sector origin
  REAL :: yo_grb          ! Y coordinate of image sector origin
  REAL :: zo_grb          ! Camera altitude from center of Earth
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  x0_ext = 0.0
  y0_ext = 0.0

!-----------------------------------------------------------------------
!
!  Define soil depths
!
!-----------------------------------------------------------------------
  IF(soilmodel_option == 1 .AND. nzsoil_ext /= 2) THEN
    WRITE(6,'(2a,I3,a/2a,I2,a/,a,a)')                                   &
               'NCEP GFS grid #3 only provides 2 soil layers, However,',&
               ' You are trying to extract ',nzsoil_ext, ' layers.',    &
               ' for ARPS two-layer force-restore model ',              &
               '(soilmodel_option = ',soilmodel_option,')',             &
               ' Please check the code ext2arps.f90 for NCEPAVN3 grid.',&
               ' Terminating ...'
!   STOP
    CALL arpsstop("NCEPAVN3 problem",1)
  ELSE IF( soilmodel_option == 2 .AND. nzsoil_ext /= 4) THEN
    WRITE(6,'(2a,I3,a/2a,I2,a/,a,a)')                                   &
               'NCEP GFS grid #3 only provides 4 soil layers, However,',&
               ' You are trying to extract ',nzsoil_ext, ' layers.',    &
               ' for ARPS multi-layer OUSoil model ',                   &
               '(soilmodel_option = ',soilmodel_option,')',             &
               ' Please check the code ext2arps.f90 for NCEPAVN3 grid.',&
               ' Terminating ...'
!   STOP
    CALL arpsstop("NCEPAVN3 problem",1)
  END IF

  IF(soilmodel_option == 1) THEN
    DO j=1,ny_ext
      DO i=1,nx_ext
        zpsoil_ext(i,j,1) = soildepth_ext(2)/2.0
        zpsoil_ext(i,j,2) = soildepth_ext(4)
      END DO
    END DO
  ELSE
    DO k=1,nzsoil_ext
      DO j=1,ny_ext
        DO i=1,nx_ext
          zpsoil_ext(i,j,k) = soildepth_ext(k)
        END DO
      END DO
    END DO
  END IF

!-----------------------------------------------------------------------
!
! Check GRIB file format
!
!-----------------------------------------------------------------------

  CALL getgrbfname(extdopt,extdfmt,dir_extd,extdname,                   &
                   extdinit,extdfcst,'    ',                            &
                   gribfile,grbflen,gribtime,grbtlen,istatus)
  IF (istatus /=0 ) RETURN

  CALL chkgrb(gribfile,grbflen,igrbfmt,istatus)
  IF (istatus /=0 ) RETURN

  IF (igrbfmt == 2) THEN  ! Decode GRIB2 file

    CALL getncepgfs_grb2(nx_ext,ny_ext,nz_ext,nzsoil_ext,4,             &
           gribfile,grbflen,gribtime,soildepth_ext,                     &
           lon_0_360,iproj_ext,scale_ext,                               &
           trlon_ext,latnot_ext,lat_ext,lon_ext,                        &
           p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                      &
           qc_ext,                                                      &
           tsoil_ext,qsoil_ext,wetcanp_ext,                             &
           snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext,                   &
           t_2m_ext, qv_2m_ext, u_10m_ext, v_10m_ext,                   &
           rain_ext,istatus)

    IF (istatus /= 0) RETURN
    istatus = 1

  ELSE                    ! Decode GRIB file

    ALLOCATE(var_grb2d(nx_ext,ny_ext,n2dvs,n2dlvt))
    ALLOCATE(var_grb3d(nx_ext,ny_ext,nz_ext,n3dvs,n3dlvt))
    ALLOCATE(rcdata(nx_ext*ny_ext))
    ALLOCATE(var_lev3d(nz_ext,n3dvs,n3dlvt))
    var_grb2d = 0.0
    var_grb3d = 0.0
    rcdata    = 0.0
    var_lev3d = 0.0
!
!-----------------------------------------------------------------------
!
!  RDNMCGRB reads NMC GRIB data
!
!-----------------------------------------------------------------------
!
    gridtyp   = avn3grid
    mproj_grb = avn3proj

    IF (extdopt == 51 .OR. nofixdim == 1) THEN
      WRITE(6,*) 'Special handling for extdopt ', extdopt
      gridtyp  = 255
    END IF

    n2dvars  = avn3nvs2d
    n2dlvtps = avn3nlvt2d

    DO k=1,n2dlvtps
      DO n=1,n2dvars
        var_id2d(n,k) = avn3var_id2d(n,k)
      END DO
      levtyp2d(k) = avn3levs2d(k)
    END DO

    n3dvars  = avn3nvs3d
    n3dlvtps = avn3nlvt3d

    DO m=1,n3dlvtps
      DO n=1,n3dvars
        var_id3d(n,m) = avn3var_id3d(n,m)
      END DO
      levtyp3d(m) = avn3levs3d(m)
    END DO

    CALL rdnmcgrb(nx_ext,ny_ext,nz_ext,gribfile,grbflen, gribtime,        &
                  gridesc, iproj_grb, gthin,                              &
                  ni_grb,nj_grb,np_grb, nk_grb,zk_grb, npeq,nit,          &
                  pi_grb,pj_grb,ipole, di_grb,dj_grb,                     &
                  latsw,lonsw, latne,lonne,                               &
                  latrot,lonrot,angrot,                                   &
                  latstr,lonstr,facstr,                                   &
                  lattru1,lattru2,lontrue,                                &
                  scanmode, iscan,jscan,kscan,                            &
                  ires,iearth,icomp,                                      &
                  jpenta,kpenta,mpenta,ispect,icoeff,                     &
                  xp_grb,yp_grb, xo_grb,yo_grb,zo_grb,                    &
                  rcdata,var_grb2d,var_grb3d,var_lev3d,lvldbg,istatus)

    IF (istatus /= 0)  THEN
      ISTATUS = -888
      GOTO 999
    END IF

    max_nr2d = 0
    DO n=1,n2dvars
      DO m=1,n2dlvtps
        max_nr2d = MAX( max_nr2d, var_nr2d(n,m) )
      END DO
    END DO

    max_nr3d = 0
    min_nr3d = nz_ext
    DO n=1,n3dvars
      max_nr3d = MAX( max_nr3d, var_nr3d(n,1) )
      min_nr3d = MIN( min_nr3d, var_nr3d(n,1) )
    END DO

    IF ( max_nr3d == 0 ) THEN
      WRITE (6,'(a)')                                                     &
          'No 3-D variable was found in the GRIB file',                   &
           'Program stopped in GETNCEPAVN3.'
!     STOP
      ISTATUS = -888
      GOTO 999
    END IF

    IF ( max_nr2d == 0 ) THEN
      WRITE (6,'(a)')                                                     &
          'No 2-D variables was found in the GRIB file'
    END IF

!    WRITE (6,'(/a7,2x,6(i7))')                                          &
!           'Lev\\VID',(var_id3d(n,1),n=1,n3dvars)
!
!    DO k=1,max_nr3d
!      var_lev3d(k,5,1) = var_lev3d(k,1,1)
!      var_lev3d(k,6,1) = var_lev3d(k,1,1)
!      WRITE (6,'(/i5,4x,6(i7))') k,(var_lev3d(k,n,1),n=1,n3dvars)
!    END DO

    DO k=1,min_nr3d
      DO n=2,n3dvars
        IF ( var_lev3d(k,1,1) /= var_lev3d(k,n,1) ) THEN
          WRITE (6,'(a)')                                                 &
              'Variables were not at the same level.',                    &
              'Program stopped in GETNCEPAVN3.'
!         STOP
          ISTATUS = -888
          GOTO 999
        END IF
      END DO
    END DO

    IF ( iproj_grb == 5 .AND. ipole == 0 ) THEN      ! Center in NP
      iproj_ext = 1
    ELSE IF ( iproj_grb == 5 .AND. ipole == 1 ) THEN  ! Center in SP
      iproj_ext = -1
    ELSE IF ( iproj_grb == 3 .AND. ipole == 0 ) THEN  ! Center in NP
      iproj_ext = 2
    ELSE IF ( iproj_grb == 3 .AND. ipole == 1 ) THEN  ! Center in SP
      iproj_ext = -2
    ELSE IF ( iproj_grb == 1 ) THEN
      iproj_ext = 3
    ELSE IF ( iproj_grb == 0 ) THEN
      iproj_ext = 4
    ELSE
      WRITE (6,'(a)')                                                     &
          'Unknown map projection. Set to non-projection.'
      iproj_ext = 0
    END IF

    scale_ext = 1.0
    latnot_ext(1) = lattru1
    latnot_ext(2) = lattru2
    trlon_ext = lontrue

    dx_ext = di_grb
    dy_ext = dj_grb

    DO j=1, ny_ext
      DO i=1, nx_ext
        lon_ext(i,j)= lonsw + (i-1) * dx_ext
        lat_ext(i,j)= latsw + (j-1) * dy_ext
      END DO
    END DO

    IF ( lon_0_360 ) THEN

      DO i = 1,nx_ext
        ii(i) = i
      END DO

    ELSE
      IF (mod(nx_ext,2) /= 0) THEN
        WRITE(6,'(a/)') 'Wrong size of nx_ext in GETNCEPAVN3 for lon_0_360.'
        CALL arpsstop('Wrong nx_ext size.',1)
      END IF

      DO i = 1,nx_ext/2              ! map 1-180 to 181-360
        ii(i) = i + nx_ext/2
      END DO

      DO i = nx_ext/2+1, nx_ext      ! map 181-360 to 1-180
        ii(i) = i - nx_ext/2
      END DO

      WHERE (lon_ext >= 180) lon_ext = lon_ext - 360

      DO j = 1,ny_ext                ! swap lat_ext & lon_ext
        DO i = 1,nx_ext/2
          rtmp = lat_ext(ii(i),j)
          lat_ext(ii(i),j) = lat_ext(i,j)
          lat_ext(i,    j) = rtmp

          rtmp = lon_ext(ii(i),j)
          lon_ext(ii(i),j) = lon_ext(i,j)
          lon_ext(i,    j) = rtmp
        END DO
      END DO
    END IF

    IF (myproc == 0) PRINT *,'LatSW = ',latsw,' LonSW = ',lonsw
!
!-----------------------------------------------------------------------
!
!  Retrieve 2-D variables
!
!-----------------------------------------------------------------------
!

    DO j=1,ny_ext
      DO i=1,nx_ext
        IF ( var_nr2d(1,1) == 0 ) THEN
          psfc_ext   (i,j) = -999.0
        ELSE
          !psfc_ext   (i,j) = var_grb2d(i,j,1,1) * 100.0
          psfc_ext   (i,j) = var_grb2d(ii(i),j,1,1) ! already Pa (F.KONG)
        END IF
        IF ( var_nr2d(2,1) == 0 ) THEN
          trn_ext    (i,j) = -999.0
        ELSE
          !trn_ext    (i,j) = var_grb2d(i,j,2,1)/g
          trn_ext    (i,j) = var_grb2d(ii(i),j,2,1)   ! why divided by g? (F.KONG)
        END IF

        IF ( var_nr3d(1,2) == 0 ) THEN
          tsoil_ext (i,j,:) = -999.0
          qsoil_ext (i,j,:) = -999.0
        ELSE

!          IF(soilmodel_option == 1) THEN       ! ARPS: two-layer force-restore model)
!
!            tsoil_ext(i,j,1) = var_grb2d(ii(i),j,3,1)      ! sfc temp.
!
!            IF ( nint(var_grb2d(ii(i),j,5,1)) == 1 ) THEN  ! soil temp over land
!              tsoil_ext(i,j,2) = var_grb3d(ii(i),j,1,1,2)
!              IF ( tsoil_ext (i,j,2) <= 200. ) THEN
!                tsoil_ext(i,j,2) = tsoil_ext(i,j,1)
!              END IF
!              qsoil_ext(i,j,1) = var_grb3d(ii(i),j,2,2,2)
!              qsoil_ext(i,j,2) = var_grb3d(ii(i),j,1,2,2)
!              soiltyp_ext(i,j)= 0          ! We do not know the soil type
!
!            ELSE                                       ! sfc temp over sea
!              tsoil_ext(i,j,2) = tsoil_ext(i,j,1)
!              qsoil_ext(i,j,:) = 1.0
!              soiltyp_ext(i,j) = 13                    ! Water
!            END IF
!
!          ELSE                                ! ARPS: multi-layer OUSoil model
!
!            tsoil_ext(i,j,1) = var_grb2d(ii(i),j,3,1)      ! sfc temp.
!
!            IF ( nint(var_grb2d(ii(i),j,5,1)) == 1 ) THEN  ! soil temp over land
!              tsoil_ext(i,j,2) = var_grb3d(ii(i),j,2,1,2)  !  0 - 10  cm
!              tsoil_ext(i,j,3) = var_grb3d(ii(i),j,1,1,2)  ! 10 - 200 cm
!
!              qsoil_ext(i,j,1) = var_grb3d(ii(i),j,2,2,2)  !  0 - 10 cm
!                                                       ! Value at surface assume equals
!                                                       ! the first layer below ground
!              qsoil_ext(i,j,2) = var_grb3d(ii(i),j,2,2,2)  !  0 - 10 cm
!              qsoil_ext(i,j,3) = var_grb3d(ii(i),j,1,2,2)  ! 10 - 200 cm
!
!              soiltyp_ext(i,j)= 0          ! We do not know the soil type
!
!            ELSE                                       ! sfc temp over sea
!              tsoil_ext(i,j,2:3) = tsoil_ext(i,j,1)    ! over water, All soil layers = SST
!              qsoil_ext(i,j,1:3) = 1.0                 !             All soil layers = 1.0
!              soiltyp_ext(i,j)   = 13        ! Water
!            END IF
!
!          END IF    ! soilmodel_option

          IF ( soilmodel_option == 1) THEN       ! ARPS: two-layer force-restore model)

            tsoil_ext(i,j,1) = var_grb2d(ii(i),j,3,1)      ! sfc temp.

            IF ( nint(var_grb2d(ii(i),j,5,1)) == 1 ) THEN  ! soil temp over land
              tsoil_ext(i,j,2) = 0.1 * var_grb3d(ii(i),j,1,1,2) & !   0-10cm
                               + 0.3 * var_grb3d(ii(i),j,2,1,2) & !  10-40cm
                               + 0.6 * var_grb3d(ii(i),j,3,1,2)   ! 40-100cm
              IF ( tsoil_ext (i,j,2) <= 200. ) tsoil_ext(i,j,2) = tsoil_ext(i,j,1)

              qsoil_ext(i,j,1) = var_grb3d(ii(i),j,1,2,2)
              qsoil_ext(i,j,2) = 0.1 * var_grb3d(ii(i),j,1,2,2) & !   0-10cm
                               + 0.3 * var_grb3d(ii(i),j,2,2,2) & !  10-40cm
                               + 0.6 * var_grb3d(ii(i),j,3,2,2)   ! 40-100cm
              soiltyp_ext(i,j)= 0          ! We do not know the soil type

            ELSE                                       ! sfc temp over sea
              tsoil_ext(i,j,2) = tsoil_ext(i,j,1)
              qsoil_ext(i,j,:) = 1.0
              soiltyp_ext(i,j) = 13                    ! Water
            END IF

          ELSE                                ! ARPS: multi-layer OUSoil model

            tsoil_ext(i,j,1) = 0.5*(  var_grb2d(ii(i),j,3,1)         &     ! sfc temp.
                                    + var_grb3d(ii(i),j,1,1,2) )

            IF ( nint(var_grb2d(ii(i),j,5,1)) == 1 ) THEN  ! soil temp over land
              qsoil_ext(i,j,1) = var_grb3d(ii(i),j,1,2,2)

              DO k = 2, nzsoil_ext
                tsoil_ext(i,j,k) = var_grb3d(ii(i),j,k,1,2)
                qsoil_ext(i,j,k) = var_grb3d(ii(i),j,k,2,2)
              END DO

              soiltyp_ext(i,j)= 0          ! We do not know the soil type

            ELSE                                       ! sfc temp over sea
              tsoil_ext(i,j,2:nzsoil_ext) = tsoil_ext(i,j,1)    ! over water, All soil layers = SST
              qsoil_ext(i,j,1:nzsoil_ext) = 1.0                 !             All soil layers = 1.0
              soiltyp_ext(i,j)   = 13        ! Water
            END IF

          END IF    ! soilmodel_option

        END IF

        IF ( var_nr2d(4,1) == 0 ) THEN
          wetcanp_ext(i,j) = -999.0
        ELSE
          wetcanp_ext(i,j) = var_grb2d(ii(i),j,4,1)*1.e-3     ! in meter
        END IF

        IF ( var_nr2d(6,1) == 0 ) THEN
          snowdpth_ext(i,j) = -999.
        ELSE

!        Convert water equiv. of accum. snow depth (kg/m**2) to meters
!        (where 1 meter liquid water is set equivqlent to 10 meters snow).
!            0.01 = 10. (m snow/m liquid) / (1000 kg/m**3)

          snowdpth_ext(i,j) = 0.01 * var_grb2d(ii(i),j,6,1)  ! in meters

        END IF

!   F.KONG add
        IF( var_nr2d(1,4) == 0 ) THEN
          t_2m_ext(i,j)= -999.0
        ELSE
          t_2m_ext(i,j)= var_grb2d(ii(i),j,1,4)
        END IF
        IF( var_nr2d(2,4) == 0 ) THEN
          qv_2m_ext(i,j)= -999.0
        ELSE
          qv_2m_ext(i,j)= var_grb2d(ii(i),j,2,4)
        END IF

        IF( var_nr2d(3,4) == 0 ) THEN
          u_10m_ext(i,j)= -999.0
        ELSE
          u_10m_ext(i,j)= var_grb2d(ii(i),j,3,4)
        END IF

        IF( var_nr2d(4,4) == 0 ) THEN
          v_10m_ext(i,j)= -999.0
        ELSE
          v_10m_ext(i,j)= var_grb2d(ii(i),j,4,4)
        END IF
! end F.KONG mod

      END DO
    END DO

!
!-----------------------------------------------------------------------
!
!  Retrieve 3-D variables
!
!-----------------------------------------------------------------------
!
    nz1 = MIN(var_nr3d(1,1),nz_ext)

    IF ( var_lev3d(1,1,1) > var_lev3d(nz1,1,1) ) THEN  ! 1st level at sfc
      chklev = 1
      lvscan = 0
    ELSE
      chklev = -1
      lvscan = nz1+1
    END IF

    DO k=1,nz1
      kk = chklev * k + lvscan
      DO j=1,ny_ext
        DO i=1,nx_ext
          p_ext  (i,j,kk) = 100.0 * REAL(var_lev3d(k,1,1)) ! Pressure
          hgt_ext(i,j,kk) = var_grb3d(ii(i),j,k,1,1)
          u_ext  (i,j,kk) = var_grb3d(ii(i),j,k,2,1)    ! u wind (m/s)
          v_ext  (i,j,kk) = var_grb3d(ii(i),j,k,3,1)    ! v wind (m/s)

          t_ext  (i,j,kk) = var_grb3d(ii(i),j,k,4,1)    ! Temperature (K)

          qc_ext (i,j,kk) = var_grb3d(ii(i),j,k,7,1)
          !qr_ext (i,j,kk) = -999.
          !qi_ext (i,j,kk) = -999.
          !qs_ext (i,j,kk) = -999.
          !qh_ext (i,j,kk) = -999.
        END DO
      END DO
    END DO

    CALL getqvs(nx_ext,ny_ext,nz1, 1,nx_ext,1,ny_ext,1,nz1,               &
                p_ext, t_ext, qv_ext )

    nz2 = MIN( nz1, min_nr3d )

    DO k=1,nz2
      kk = chklev * k + lvscan
      DO j=1,ny_ext
        DO i=1,nx_ext
         qv_ext(i,j,kk) = 0.01*var_grb3d(ii(i),j,k,5,1)*qv_ext(i,j,kk)
        END DO
      END DO
    END DO

    IF ( nz2 < nz1 ) THEN
      DO k=nz2+1,nz1
        kk = chklev * k + lvscan
        DO j=1,ny_ext
          DO i=1,nx_ext
            qv_ext(i,j,kk) = 0.0
            var_lev3d(k,5,1) = var_lev3d(k,1,1)
          END DO
        END DO
      END DO
    END IF

    istatus = 1

    999  CONTINUE

    DEALLOCATE(var_grb2d,var_grb3d,rcdata,var_lev3d)

  END IF

  RETURN
END SUBROUTINE getncepavn3

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETNCEPAVN2                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE getncepavn2(nx_ext,ny_ext,nz_ext,nzsoil_ext,                 &
           dir_extd,extdname,extdopt,extdfmt,lon_0_360,                 &
           extdinit,extdfcst,julfname,                                  &
           iproj_ext,scale_ext,                                         &
           trlon_ext,latnot_ext,x0_ext,y0_ext,                          &
           lat_ext,lon_ext,zpsoil_ext,                                  &
           p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                      &
           qc_ext,                                                      &
           tsoil_ext,qsoil_ext,wetcanp_ext,                             &
           snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext,                   &
           istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads a NCEP GFS or NCAR/NCEP Reanalysis GRIB (Grid #2, 2.5 x 2.5 degree)
!  data file for processing by ext2arps.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Richard Carpenter
!  2003-07-03
!
!  MODIFICATION HISTORY:
!
!  08/11/2005 (Y. Wang)
!  Modified to handle cases with Greenwich Meridian in the domain. Note
!  that a new parameter, lon_0_360 was passed in and an integer array
!  ii(nx_ext) was decalared. Still not tested.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    dir_extd      Directory name for external file
!    extdname      Prefix string of external file name
!    extdopt       Option of external data sources
!    extdfmt       Option of external data format
!
!    extdinit      Initialized time in mm-dd-yyyy:hh:mm:ss format
!    extdfcst      Forecast hour in HHH:MM:SS format
!    julfname      File name in yyjjjhhmm format
!
!  OUTPUT:
!
!    iproj_ext     Map projection number of external data
!    scale_ext     Scale factor of external data
!    trlon_ext     True longitude of external data (degrees E)
!    latnot_ext(2) True latitude(s) of external data (degrees N)
!    x0_ext        x coordinate of origin of external data
!    y0_ext        y coordinate of origin of external data
!    lat_ext       latitude of external data points (degrees N)
!    lon_ext       longitude of external data points (degrees E)
!    p_ext         pressure (Pascal)
!    hgt_ext       height (m)
!    t_ext         temperature (K)
!    qv_ext        specific humidity (kg/kg)
!    u_ext         u wind component (m/s)
!    v_ext         v wind component (m/s)
!    qc_ext        Cloud water mixing ratio (kg/kg)
!    qr_ext        Rain water mixing ratio (kg/kg)
!    qi_ext        Ice mixing ratio (kg/kg)
!    qs_ext        Snow mixing ratio (kg/kg)
!    qh_ext        Hail mixing ratio (kg/kg)
!
!    tsoil_ext     Soil temperature
!    qsoil_ext     Soil moisture
!    wetcanp_ext   Water content on canopy
!
!    trn_ext       External terrain (m)
!    psfc_ext      Surface pressure (Pa)
!
!    istatus       status indicator
!
!  WORK ARRAYS:
!
!    var_grb3d     Arrays to store the GRIB 3-D variables:
!                  var_grb3d(nxgrb,nygrb,nzgrb,1,1) - Temperature (K)
!                  var_grb3d(nxgrb,nygrb,nzgrb,2,1) - Specific humidity
!                                                     (kg/kg)
!                  var_grb3d(nxgrb,nygrb,nzgrb,3,1) - u wind (m/s)
!                  var_grb3d(nxgrb,nygrb,nzgrb,4,1) - v wind (m/s)
!                  var_grb3d(nxgrb,nygrb,nzgrb,5,1) - Geopotential
!                                                     height (gpm)
!                  var_grb3d(nxgrb,nygrb,nzgrb,6,1) - Pressure vertical
!                                                     velocity (Pa/s)
!                                                     (if applied)
!                  var_grb3d(nxgrb,nygrb,nzgrb,1,2) - soil temp. (K)
!                  var_grb3d(nxgrb,nygrb,nzgrb,2,2) - vol. soil moist.
!                                                     (m**3/m**3)
!
!    var_grb2d     Arrays to store the GRIB 2-D variables:
!                  var_grb2d(nxgrb,nygrb,1) - Surface pressure (Pa)
!                  var_grb2d(nxgrb,nygrb,2) - Geopotential height (gpm)
!                  var_grb2d(nxgrb,nygrb,3) - Surface temperature (K)
!                  var_grb2d(nxgrb,nygrb,4) - Plant canopy surface
!                                             water (kg/m**2)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'gribcst.inc'
  INCLUDE 'mp.inc'

  CHARACTER (LEN=*) :: dir_extd
  CHARACTER (LEN=*) :: extdname

  INTEGER :: extdopt
  INTEGER :: extdfmt
  LOGICAL, INTENT(IN) :: lon_0_360

  CHARACTER (LEN=19) :: extdinit
  CHARACTER (LEN=9) :: extdfcst
  CHARACTER (LEN=9) :: julfname
!
!-----------------------------------------------------------------------
!
!  External grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj_ext
  REAL :: scale_ext,trlon_ext
  REAL :: latnot_ext(2)
  REAL :: x0_ext,y0_ext
  REAL :: dx_ext,dy_ext
!
!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx_ext,ny_ext,nz_ext,nzsoil_ext

  REAL :: lat_ext(nx_ext,ny_ext)
  REAL :: lon_ext(nx_ext,ny_ext)
  REAL :: zpsoil_ext(nx_ext,ny_ext,nzsoil_ext)
  REAL :: p_ext  (nx_ext,ny_ext,nz_ext)   ! Pressure (Pascals)
  REAL :: hgt_ext(nx_ext,ny_ext,nz_ext)   ! Height (m)
  REAL :: t_ext  (nx_ext,ny_ext,nz_ext)   ! Temperature (K)
  REAL :: qv_ext (nx_ext,ny_ext,nz_ext)   ! Specific humidity (kg/kg)
  REAL :: u_ext  (nx_ext,ny_ext,nz_ext)   ! Eastward wind component
  REAL :: v_ext  (nx_ext,ny_ext,nz_ext)   ! Northward wind component
  REAL :: qc_ext (nx_ext,ny_ext,nz_ext)   ! Cloud H2O mixing ratio (kg/kg)
!  REAL :: qr_ext (nx_ext,ny_ext,nz_ext)   ! Rain  H2O mixing ratio (kg/kg)
!  REAL :: qi_ext (nx_ext,ny_ext,nz_ext)   ! Ice   mixing ratio (kg/kg)
!  REAL :: qs_ext (nx_ext,ny_ext,nz_ext)   ! Snow  mixing ratio (kg/kg)
!  REAL :: qh_ext (nx_ext,ny_ext,nz_ext)   ! Hail  mixing ratio (kg/kg)
!
  REAL :: tsoil_ext  (nx_ext,ny_ext,nzsoil_ext)      ! Temperature at surface (K)
  REAL :: qsoil_ext  (nx_ext,ny_ext,nzsoil_ext)      ! Deep soil temperature (K)
  REAL :: wetcanp_ext(nx_ext,ny_ext)                 ! Canopy water amount
  REAL :: snowdpth_ext(nx_ext,ny_ext)                ! Snow depth (m)

  REAL :: trn_ext    (nx_ext,ny_ext)      ! External terrain (m)
  REAL :: psfc_ext   (nx_ext,ny_ext)      ! Surface pressure (Pa)

  INTEGER :: soiltyp_ext(nx_ext,ny_ext)
!
!-----------------------------------------------------------------------
!
!  Other  external variable arrays
!
!-----------------------------------------------------------------------
!
!  REAL :: x_ext(nx_ext)
!  REAL :: y_ext(ny_ext)

  INTEGER :: ii(nx_ext)
  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Work arrays for storing grib data
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: var_grb2d(:,:,:,:)   ! GRIB variables
  REAL, ALLOCATABLE :: var_grb3d(:,:,:,:,:) ! GRIB 3-D variables
  INTEGER, ALLOCATABLE :: var_lev3d(:,:,:)  ! Levels (hybrid) for
                                            ! each 3-D variable
  REAL, ALLOCATABLE :: rcdata(:)            ! temporary data array
!
!-----------------------------------------------------------------------
!
!  Original grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj
  REAL :: scale,trlon,x0,y0
  REAL :: latnot(2)
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: gribfile
  CHARACTER (LEN=14) :: gribtime
  INTEGER :: i,j,k,kk
  INTEGER :: grbflen, grbtlen

  INTEGER :: m,n,nz1,max_nr2d,max_nr3d,min_nr3d,nz2

  REAL :: govrd

  INTEGER :: chklev, lvscan

  INTEGER :: iret             ! Return flag

  REAL    :: rtmp

  REAL, PARAMETER :: soildepth_ext(3) = (/0.0, 0.05, 1.05/)
                             ! Surface and two soil layers  0-10cm & 10-200cm

!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  GRIB grid information
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=42) :: gridesc ! Grid description

  INTEGER :: iproj_grb    ! Map projection indicator
  INTEGER :: gthin        ! Indicator of whether the grid is "thinned"

  INTEGER :: ni_grb       ! Number of points along x-axis
  INTEGER :: nj_grb       ! Number of points along y-axis
  INTEGER :: np_grb       ! Total number of horizontal grid points

  INTEGER :: nk_grb       ! Number of vertical parameters
  REAL :: zk_grb(nz_ext)  ! Vertical coordinate parameters

  INTEGER :: npeq         ! Number of lat circles from pole to equator
  INTEGER :: nit(nz_ext)  ! Number of x-points for thinned grid

  REAL :: pi_grb          ! x-coordinate of pole point
  REAL :: pj_grb          ! y-coordinate of pole point
  INTEGER :: ipole        ! Projection center flag

  REAL :: di_grb          ! x-direction increment or grid length
  REAL :: dj_grb          ! y-direction increment or grid length

  REAL :: latsw           ! Latitude  of South West corner point
  REAL :: lonsw           ! Longitude of South West corner point
  REAL :: latne           ! Latitude  of North East corner point
  REAL :: lonne           ! Longitude of North East corner point

  REAL :: lattru1         ! Latitude (1st) at which projection is true
  REAL :: lattru2         ! Latitude (2nd) at which projection is true
  REAL :: lontrue         ! Longitude      at which projection is true

  REAL :: latrot          ! Latitude  of southern pole of rotation
  REAL :: lonrot          ! Longitude of southern pole of rotation
  REAL :: angrot          ! Angle of rotation

  REAL :: latstr          ! Latitude  of the pole of stretching
  REAL :: lonstr          ! Longitude of the pole of stretching
  REAL :: facstr          ! Stretching factor

  INTEGER :: scanmode     ! Scanning indicator
  INTEGER :: iscan        ! x-direction   scanning indicator
  INTEGER :: jscan        ! y-direction   scanning indicator
  INTEGER :: kscan        ! FORTRAN index scanning indicator

  INTEGER :: ires         ! Resolution direction increments indicator
  INTEGER :: iearth       ! Earth shape indicator: spherical or oblate?
  INTEGER :: icomp        ! (u,v) components decomposition indicator

  INTEGER :: jpenta       ! J-Pentagonal resolution parameter
  INTEGER :: kpenta       ! K-Pentagonal resolution parameter
  INTEGER :: mpenta       ! M-Pentagonal resolution parameter
  INTEGER :: ispect       ! Spectral representation type
  INTEGER :: icoeff       ! Spectral coefficient storage mode

  REAL :: xp_grb          ! X coordinate of sub-satellite point
  REAL :: yp_grb          ! Y coordinate of sub-satellite point
  REAL :: xo_grb          ! X coordinate of image sector origin
  REAL :: yo_grb          ! Y coordinate of image sector origin
  REAL :: zo_grb          ! Camera altitude from center of Earth
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  x0_ext = 0.0
  y0_ext = 0.0

!---------------------------------------------------------------------
!
!  Define soil depths
!
!---------------------------------------------------------------------

  IF(soilmodel_option == 1 .AND. nzsoil_ext /= 2) THEN
    WRITE(6,'(2a,I3,a/2a,I2,a/,a,a)')                                   &
               'NCEP GFS grid #2 only provides 2 soil layers, However,',&
               ' You are trying to extract ',nzsoil_ext, ' layers.',    &
               ' for ARPS two-layer force-restore model ',              &
               '(soilmodel_option = ',soilmodel_option,')',             &
               ' Please check the code ext2arps.f90 for NCEPAVN2 grid.',&
               ' Terminating ...'
!   STOP
    CALL arpsstop("NCEPAVN2 problem",1)
  ELSE IF( soilmodel_option == 2 .AND. nzsoil_ext /= 3) THEN
    WRITE(6,'(2a,I3,a/2a,I2,a/,a,a)')                                   &
               'NCEP GFS grid #2 only provides 2 soil layers, However,',&
               ' You are trying to extract ',nzsoil_ext, ' layers.',    &
               ' for ARPS multi-layer OUSoil model ',                   &
               '(soilmodel_option = ',soilmodel_option,')',             &
               ' Please check the code ext2arps.f90 for NCEPAVN2 grid.',&
               ' Terminating ...'
!   STOP
    CALL arpsstop("NCEPAVN2 problem",1)
  END IF

  IF(soilmodel_option == 1) THEN
    DO k=1,nzsoil_ext
      DO j=1,ny_ext
        DO i=1,nx_ext
          zpsoil_ext(i,j,k) = soildepth_ext(k+1)
        END DO
      END DO
    END DO
  ELSE
    DO k=1,nzsoil_ext
      DO j=1,ny_ext
        DO i=1,nx_ext
          zpsoil_ext(i,j,k) = soildepth_ext(k)
        END DO
      END DO
    END DO
  END IF

  ALLOCATE(var_grb2d(nx_ext,ny_ext,n2dvs,n2dlvt))
  ALLOCATE(var_grb3d(nx_ext,ny_ext,nz_ext,n3dvs,n3dlvt))
  ALLOCATE(rcdata(nx_ext*ny_ext))
  ALLOCATE(var_lev3d(nz_ext,n3dvs,n3dlvt))
  var_grb2d = 0.0
  var_grb3d = 0.0
  rcdata    = 0.0
  var_lev3d = 0.0

  CALL getgrbfname(extdopt,extdfmt,dir_extd,extdname,                   &
                   extdinit,extdfcst,'    ',                            &
                   gribfile,grbflen,gribtime,grbtlen,istatus)
  IF (istatus /=0 ) RETURN
!
!-----------------------------------------------------------------------
!
!  RDNMCGRB reads NMC GRIB data
!
!-----------------------------------------------------------------------
!
  gridtyp  = avn2grid
  mproj_grb = avn2proj

  IF (extdopt == 51) THEN
    WRITE(6,*) 'Special handling for extdopt ', extdopt
    gridtyp  = 255
  END IF

  n2dvars  = avn2nvs2d
  n2dlvtps = avn2nlvt2d

  DO k=1,n2dlvtps
    DO n=1,n2dvars
      var_id2d(n,k) = avn2var_id2d(n,k)
    END DO
    levtyp2d(k) = avn2levs2d(k)
  END DO

  n3dvars  = avn2nvs3d
  n3dlvtps = avn2nlvt3d

  DO m=1,n3dlvtps
    DO n=1,n3dvars
      var_id3d(n,m) = avn2var_id3d(n,m)
    END DO
    levtyp3d(m) = avn2levs3d(m)
  END DO

  CALL rdnmcgrb(nx_ext,ny_ext,nz_ext,gribfile,grbflen, gribtime,        &
                gridesc, iproj_grb, gthin,                              &
                ni_grb,nj_grb,np_grb, nk_grb,zk_grb, npeq,nit,          &
                pi_grb,pj_grb,ipole, di_grb,dj_grb,                     &
                latsw,lonsw, latne,lonne,                               &
                latrot,lonrot,angrot,                                   &
                latstr,lonstr,facstr,                                   &
                lattru1,lattru2,lontrue,                                &
                scanmode, iscan,jscan,kscan,                            &
                ires,iearth,icomp,                                      &
                jpenta,kpenta,mpenta,ispect,icoeff,                     &
                xp_grb,yp_grb, xo_grb,yo_grb,zo_grb,                    &
                rcdata,var_grb2d,var_grb3d,var_lev3d,lvldbg,iret)

  IF (iret /= 0)  THEN
    ISTATUS = -888
    GOTO 999
  END IF

  max_nr2d = 0
  DO n=1,n2dvars
    DO m=1,n2dlvtps
      max_nr2d = MAX( max_nr2d, var_nr2d(n,m) )
    END DO
  END DO

  max_nr3d = 0
  min_nr3d = nz_ext
  DO n=1,n3dvars
    max_nr3d = MAX( max_nr3d, var_nr3d(n,1) )
    min_nr3d = MIN( min_nr3d, var_nr3d(n,1) )
  END DO

  IF ( max_nr3d == 0 ) THEN
    WRITE (6,'(a)')                                                     &
        'No 3-D variable was found in the GRIB file',                   &
         'Program stopped in GETNCEPAVN2.'
!   STOP
    ISTATUS = -888
    GOTO 999
  END IF

  IF ( max_nr2d == 0 ) THEN
    WRITE (6,'(a)')                                                     &
        'No 2-D variables was found in the GRIB file'
  END IF

!  write (6,'(/a7,2x,6(i7))')
!    :    'Lev\\VID',(var_id3d(n,1),n=1,n3dvars)

!  DO 60 k=1,max_nr3d
!    var_lev3d(k,5,1) = var_lev3d(k,1,1)
!    var_lev3d(k,6,1) = var_lev3d(k,1,1)
!    write (6,'(/i5,4x,6(i7))')
!    :    k,(var_lev3d(k,n,1),n=1,n3dvars)
!  60    CONTINUE

  DO k=1,min_nr3d
    DO n=2,n3dvars
      IF ( var_lev3d(k,1,1) /= var_lev3d(k,n,1) ) THEN
        WRITE (6,'(a)')                                                 &
            'Variables were not at the same level.',                    &
            'Program stopped in GETNCEPAVN2.'
!       STOP
        ISTATUS = -888
        GOTO 999
      END IF
    END DO
  END DO

  IF ( iproj_grb == 5 .AND. ipole == 0 ) THEN      ! Center in NP
    iproj_ext = 1
  ELSE IF ( iproj_grb == 5 .AND. ipole == 1 ) THEN  ! Center in SP
    iproj_ext = -1
  ELSE IF ( iproj_grb == 3 .AND. ipole == 0 ) THEN  ! Center in NP
    iproj_ext = 2
  ELSE IF ( iproj_grb == 3 .AND. ipole == 1 ) THEN  ! Center in SP
    iproj_ext = -2
  ELSE IF ( iproj_grb == 1 ) THEN
    iproj_ext = 3
  ELSE IF ( iproj_grb == 0 ) THEN
    iproj_ext = 4
  ELSE
    WRITE (6,'(a)')                                                     &
        'Unknown map projection. Set to non-projection.'
    iproj_ext = 0
  END IF

  scale_ext = 1.0
  latnot_ext(1) = lattru1
  latnot_ext(2) = lattru2
  trlon_ext = lontrue

  dx_ext = di_grb
  dy_ext = dj_grb

  DO i=1, nx_ext
    DO j=1, ny_ext
      lon_ext(i,j)= lonsw + (i-1) * dx_ext
      lat_ext(i,j)= latsw + (j-1) * dy_ext
    END DO
  END DO

  IF ( lon_0_360 ) THEN

    DO i = 1,nx_ext
      ii(i) = i
    END DO

  ELSE
    IF (mod(nx_ext,2) /= 0) THEN
      WRITE(6,'(a/)') 'Wrong size of nx_ext in GETNCEPAVN2 for lon_0_360.'
      CALL arpsstop('Wrong nx_ext size.',1)
    END IF

    DO i = 1,nx_ext/2              ! map 1-180 to 181-360
      ii(i) = i + nx_ext/2
    END DO

    DO i = nx_ext/2+1, nx_ext      ! map 181-360 to 1-180
      ii(i) = i - nx_ext/2
    END DO

    WHERE (lon_ext >= 180.) lon_ext = lon_ext - 360

    DO j = 1,ny_ext                ! swap lat_ext & lon_ext
      DO i = 1,nx_ext/2
        rtmp = lat_ext(ii(i),j)
        lat_ext(ii(i),j) = lat_ext(i,j)
        lat_ext(i,    j) = rtmp

        rtmp = lon_ext(ii(i),j)
        lon_ext(ii(i),j) = lon_ext(i,j)
        lon_ext(i,    j) = rtmp
      END DO
    END DO
  END IF

  IF (myproc == 0) PRINT *,'LatSW = ',latsw,' LonSW = ',lonsw
!
!-----------------------------------------------------------------------
!
!  Retrieve 2-D variables
!
!-----------------------------------------------------------------------
!

  DO j=1,ny_ext
    DO i=1,nx_ext
      IF ( var_nr2d(1,1) == 0 ) THEN
        psfc_ext   (i,j) = -999.0
      ELSE
        psfc_ext   (i,j) = var_grb2d(ii(i),j,1,1) * 100.0
      END IF
      IF ( var_nr2d(2,1) == 0 ) THEN
        trn_ext    (i,j) = -999.0
      ELSE
        trn_ext    (i,j) = var_grb2d(ii(i),j,2,1)/g
      END IF

      IF ( var_nr3d(1,2) == 0 ) THEN
        tsoil_ext (i,j,:) = -999.0
        qsoil_ext (i,j,:) = -999.0
      ELSE

        IF(soilmodel_option == 1) THEN       ! ARPS: two-layer force-restore model)

          tsoil_ext(i,j,1) = var_grb2d(ii(i),j,3,1)      ! sfc temp.

          IF ( nint(var_grb2d(ii(i),j,5,1)) == 1 ) THEN  ! soil temp over land
            tsoil_ext(i,j,2) = var_grb3d(ii(i),j,1,1,2)
            IF ( tsoil_ext (i,j,2) <= 200. ) THEN
              tsoil_ext(i,j,2) = tsoil_ext(i,j,1)
            END IF
            qsoil_ext(i,j,1) = var_grb3d(ii(i),j,2,2,2)
            qsoil_ext(i,j,2) = var_grb3d(ii(i),j,1,2,2)
            soiltyp_ext(i,j)= 0          ! We do not know the soil type

          ELSE                                       ! sfc temp over sea
            tsoil_ext(i,j,2) = tsoil_ext(i,j,1)
            qsoil_ext(i,j,:) = 1.0
            soiltyp_ext(i,j) = 13                    ! Water
          END IF

        ELSE                                ! ARPS: multi-layer OUSoil model

          tsoil_ext(i,j,1) = var_grb2d(ii(i),j,3,1)      ! sfc temp.

          IF ( nint(var_grb2d(ii(i),j,5,1)) == 1 ) THEN  ! soil temp over land
            tsoil_ext(i,j,2) = var_grb3d(ii(i),j,2,1,2)  !  0 - 10  cm
            tsoil_ext(i,j,3) = var_grb3d(ii(i),j,1,1,2)  ! 10 - 200 cm

            qsoil_ext(i,j,1) = var_grb3d(ii(i),j,2,2,2)  !  0 - 10 cm
                                                     ! Value at surface assume equals
                                                     ! the first layer below ground
            qsoil_ext(i,j,2) = var_grb3d(ii(i),j,2,2,2)  !  0 - 10 cm
            qsoil_ext(i,j,3) = var_grb3d(ii(i),j,1,2,2)  ! 10 - 200 cm

            soiltyp_ext(i,j)= 0          ! We do not know the soil type

          ELSE                                       ! sfc temp over sea
            tsoil_ext(i,j,2:3) = tsoil_ext(i,j,1)    ! over water, All soil layers = SST
            qsoil_ext(i,j,1:3) = 1.0                 !             All soil layers = 1.0
            soiltyp_ext(i,j)   = 13        ! Water
          END IF

        END IF    ! soilmodel_option

      END IF

      IF ( var_nr2d(4,1) == 0 ) THEN
        wetcanp_ext(i,j) = -999.0
      ELSE
        wetcanp_ext(i,j) = var_grb2d(ii(i),j,4,1)*1.e-3     ! in meter
      END IF

      IF ( var_nr2d(6,1) == 0 ) THEN
        snowdpth_ext(i,j) = -999.
      ELSE

!      Convert water equiv. of accum. snow depth (kg/m**2) to meters
!      (where 1 meter liquid water is set equivqlent to 10 meters snow).
!          0.01 = 10. (m snow/m liquid) / (1000 kg/m**3)

        snowdpth_ext(i,j) = 0.01 * var_grb2d(ii(i),j,6,1)  ! in meters

      END IF

    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Retrieve 3-D variables
!
!-----------------------------------------------------------------------
!
  nz1 = MIN(var_nr3d(1,1),nz_ext)

  IF ( var_lev3d(1,1,1) > var_lev3d(nz1,1,1) ) THEN  ! 1st level at sfc
    chklev = 1
    lvscan = 0
  ELSE
    chklev = -1
    lvscan = nz1+1
  END IF

  DO k=1,nz1
    kk = chklev * k + lvscan
    DO j=1,ny_ext
      DO i=1,nx_ext
        p_ext  (i,j,kk) = 100.0 * REAL(var_lev3d(k,1,1)) ! Pressure
        hgt_ext(i,j,kk) = var_grb3d(ii(i),j,k,1,1)
        u_ext  (i,j,kk) = var_grb3d(ii(i),j,k,2,1)    ! u wind (m/s)
        v_ext  (i,j,kk) = var_grb3d(ii(i),j,k,3,1)    ! v wind (m/s)

        t_ext  (i,j,kk) = var_grb3d(ii(i),j,k,4,1)    ! Temperature (K)

        qc_ext (i,j,kk) = var_grb3d(ii(i),j,k,7,1)
!        qr_ext (i,j,kk) = -999.
!        qi_ext (i,j,kk) = -999.
!        qs_ext (i,j,kk) = -999.
!        qh_ext (i,j,kk) = -999.
      END DO
    END DO
  END DO

  CALL getqvs(nx_ext,ny_ext,nz1, 1,nx_ext,1,ny_ext,1,nz1,               &
           p_ext, t_ext, qv_ext )

  nz2 = MIN( nz1, min_nr3d )

  DO k=1,nz2
    kk = chklev * k + lvscan
    DO j=1,ny_ext
      DO i=1,nx_ext
       qv_ext(i,j,kk) = 0.01*var_grb3d(ii(i),j,k,5,1)*qv_ext(i,j,kk)
      END DO
    END DO
  END DO

  IF ( nz2 < nz1 ) THEN
    DO k=nz2+1,nz1
      kk = chklev * k + lvscan
      DO j=1,ny_ext
        DO i=1,nx_ext
          qv_ext(i,j,kk) = 0.0
          var_lev3d(k,5,1) = var_lev3d(k,1,1)
        END DO
      END DO
    END DO
  END IF

  istatus = 1

  999  CONTINUE
  DEALLOCATE(var_grb2d,var_grb3d,rcdata,var_lev3d)

  RETURN
END SUBROUTINE getncepavn2

!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GET_AVN_BIN                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE get_avn_bin(nx_ext,ny_ext,nz_ext,extdopt,extdfmt,            &
           dir_extd,extdname,extdinit,extdfcst,julfname,                &
           iproj_ext,scale_ext,trlon_ext,latnot_ext,x0_ext,y0_ext,      &
           lat_ext,lon_ext,p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,      &
           tsfc_ext,tdeep_ext,wetsfc_ext,wetdp_ext,wetcanp_ext,         &
           snowdpth_ext,trn_ext,psfc_ext,                               &
           istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads and pass out a section of NCEP AVN GRIB
!  (Grid #3, 1x1 degree) data file (created by program EXTRACT_AVN).
!
!-----------------------------------------------------------------------
!
!  AUTHOR: M. Xue
!  07/25/2000
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    dir_extd      Directory name for external file
!    extdname      Prefix string of external file name
!    extdinit      Initialized time in mm-dd-yyyy:hh:mm:ss format
!    extdfcst      Forecast hour in HHH:MM:SS format
!    julfname      File name in yyjjjhhmm format
!
!  OUTPUT:
!
!    iproj_ext     Map projection number of external data
!    scale_ext     Scale factor of external data
!    trlon_ext     True longitude of external data (degrees E)
!    latnot_ext(2) True latitude(s) of external data (degrees N)
!    x0_ext        x coordinate of origin of external data
!    y0_ext        y coordinate of origin of external data
!    lat_ext       latitude of external data points (degrees N)
!    lon_ext       longitude of external data points (degrees E)
!    p_ext         pressure (Pascal)
!    hgt_ext       height (m)
!    t_ext         temperature (K)
!    qv_ext        specific humidity (kg/kg)
!    u_ext         u wind component (m/s)
!    v_ext         v wind component (m/s)
!    qc_ext        Cloud water mixing ratio (kg/kg)
!    qr_ext        Rain water mixing ratio (kg/kg)
!    qi_ext        Ice mixing ratio (kg/kg)
!    qs_ext        Snow mixing ratio (kg/kg)
!    qh_ext        Hail mixing ratio (kg/kg)
!
!    tsfc_ext      Surface temperature
!    tdeep_ext     Soil temperature
!    wetsfc_ext    Top layer soil moisture
!    wetdp_ext     Deep soil moisture
!    wetcanp_ext   Water content on canopy
!
!    trn_ext       External terrain (m)
!    psfc_ext      Surface pressure (Pa)
!
!    istatus       status indicator
!
!  WORK ARRAYS:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx_ext,ny_ext,nz_ext

  INTEGER, INTENT(IN) :: extdopt, extdfmt

  CHARACTER (LEN=*) :: dir_extd
  CHARACTER (LEN=*) :: extdname
  CHARACTER (LEN=19) :: extdinit
  CHARACTER (LEN=9) :: extdfcst
  CHARACTER (LEN=9) :: julfname
!
!-----------------------------------------------------------------------
!
!  External grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj_ext
  REAL :: scale_ext,trlon_ext
  REAL :: latnot_ext(2)
  REAL :: x0_ext,y0_ext
  REAL :: dx_ext,dy_ext
!
!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------
!
  REAL :: lat_ext(nx_ext,ny_ext)
  REAL :: lon_ext(nx_ext,ny_ext)
  REAL :: p_ext  (nx_ext,ny_ext,nz_ext)   ! Pressure (Pascals)
  REAL :: hgt_ext(nx_ext,ny_ext,nz_ext)   ! Height (m)
  REAL :: t_ext  (nx_ext,ny_ext,nz_ext)   ! Temperature (K)
  REAL :: qv_ext (nx_ext,ny_ext,nz_ext)   ! Specific humidity (kg/kg)
  REAL :: u_ext  (nx_ext,ny_ext,nz_ext)   ! Eastward wind component
  REAL :: v_ext  (nx_ext,ny_ext,nz_ext)   ! Northward wind component
!  REAL :: qc_ext (nx_ext,ny_ext,nz_ext)   ! Cloud H2O mixing ratio (kg/kg)
!  REAL :: qr_ext (nx_ext,ny_ext,nz_ext)   ! Rain  H2O mixing ratio (kg/kg)
!  REAL :: qi_ext (nx_ext,ny_ext,nz_ext)   ! Ice   mixing ratio (kg/kg)
!  REAL :: qs_ext (nx_ext,ny_ext,nz_ext)   ! Snow  mixing ratio (kg/kg)
!  REAL :: qh_ext (nx_ext,ny_ext,nz_ext)   ! Hail  mixing ratio (kg/kg)

  REAL :: tsfc_ext   (nx_ext,ny_ext)      ! Temperature at surface (K)
  REAL :: tdeep_ext  (nx_ext,ny_ext)      ! Deep soil temperature (K)
  REAL :: wetsfc_ext (nx_ext,ny_ext)      ! Surface soil moisture
  REAL :: wetdp_ext  (nx_ext,ny_ext)      ! Deep soil moisture
  REAL :: wetcanp_ext(nx_ext,ny_ext)      ! Canopy water amount
  REAL :: snowdpth_ext(nx_ext,ny_ext)     ! Snow depth (m)

  REAL :: trn_ext    (nx_ext,ny_ext)      ! External terrain (m)
  REAL :: psfc_ext   (nx_ext,ny_ext)      ! Surface pressure (Pa)

  REAL :: temp_ext   (nx_ext,ny_ext)      ! Automatic work array

!-----------------------------------------------------------------------
!
!  Local variables
!
!-----------------------------------------------------------------------

  INTEGER :: istatus,ierr
  INTEGER :: nunit, idummy
  REAL :: rdummy
  CHARACTER (LEN=10) :: var_label

  CHARACTER (LEN=256) :: gribfile
  CHARACTER (LEN=14) :: gribtime
  INTEGER :: i,j,k
  INTEGER :: grbflen, grbtlen

  INTEGER :: iuout, ivout, ipout, ihout,itout,                          &
          iqvout, itsfcout,itsoilout,iwsfcout,iwdpout,                  &
          iwcnpout,isnowout,itrnout,ipsfcout,                           &
          iugrdout,ivgrdout,itgrdout,iptgrdout,irhgrdout,ipmslout

  INTEGER :: nx_ext_in, ny_ext_in, nz_ext_in
  REAL :: latbgn,latend,lonbgn,lonend, del_lat, del_lon
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  CALL getgrbfname(extdopt,extdfmt,dir_extd,extdname,                   &
                   extdinit,extdfcst,'    ',                            &
                   gribfile,grbflen,gribtime,grbtlen,istatus)
  IF (istatus /=0 ) RETURN

  psfc_ext   = -999.0
  trn_ext    = -999.0
  tsfc_ext   = -999.0
  tdeep_ext  = -999.0
  wetsfc_ext = -999.0
  wetdp_ext  = -999.0
  wetcanp_ext= -999.0
  snowdpth_ext=-999.0

  u_ext   = -999.0
  v_ext   = -999.0
  p_ext   = -999.0
  hgt_ext = -999.0
  t_ext   = -999.0
  qv_ext  = -999.0

  CALL getunit( nunit)
  CALL asnctl ('NEWLOCAL', 1, ierr)
  CALL asnfile(trim(gribfile), '-F f77 -N ieee', ierr)

  OPEN(UNIT=nunit,FILE=trim(gribfile),                                  &
       STATUS='old',FORM='unformatted',IOSTAT=istatus)

  IF( istatus /= 0 ) THEN
    WRITE(6,'(1x,a,a,/1x,i3,a)')                                        &
        'Error occured when opening file ',trim(gribfile),              &
        'using FORTRAN unit ',nunit,' Program stopped.'
!   STOP
    ISTATUS = -888
    RETURN
  END IF

  PRINT*,'To read file ',trim(gribfile)

  READ (nunit,ERR=910,END=920) nx_ext_in,ny_ext_in, nz_ext_in

  IF( nx_ext_in /= nx_ext.OR.ny_ext_in /= ny_ext .OR. nz_ext_in /= nz_ext ) THEN
    WRITE(6,'(1x,a)')                                                   &
         ' Dimensions in GET_AVN_BIN inconsistent with data.'
    WRITE(6,'(1x,a,3I15)') ' Read were: ',                              &
         nx_ext_in, ny_ext_in, nz_ext_in
    WRITE(6,'(1x,a,3I15)') ' Expected:  ',                              &
         nx_ext, ny_ext, nz_ext
    WRITE(6,'(1x,a)')                                                   &
         ' Program aborted in GET_AVN_BIN.'
!   STOP
    CALL arpsstop("get_avn_bin",1)
  END IF

  READ (nunit,ERR=910,END=920) lonbgn,lonend,latbgn,latend
  READ (nunit,ERR=910,END=920) del_lon, del_lat

  READ (nunit,ERR=910,END=920)                                          &
        iuout, ivout, ipout, ihout,itout,                               &
        iqvout, itsfcout,itsoilout,iwsfcout,iwdpout,                    &
        iwcnpout,isnowout,itrnout,ipsfcout,iugrdout,                    &
        ivgrdout,itgrdout,iptgrdout,irhgrdout,ipmslout,                 &
        iproj_ext,idummy,idummy, idummy, idummy,                        &
        idummy, idummy,  idummy, idummy, idummy

  READ (nunit,ERR=910,END=920)                                          &
        scale_ext,trlon_ext,latnot_ext(1),latnot_ext(2),                &
        x0_ext, y0_ext, rdummy, rdummy, rdummy, rdummy,                 &
        rdummy, rdummy, rdummy, rdummy, rdummy,                         &
        rdummy, rdummy, rdummy, rdummy, rdummy,                         &
        rdummy, rdummy, rdummy, rdummy, rdummy,                         &
        rdummy, rdummy, rdummy, rdummy, rdummy

  IF( iuout == 1 ) THEN
    READ (nunit,ERR=910,END=920) var_label
    READ (nunit,ERR=910,END=920) u_ext
    WRITE(6,'(a,a,a)')'Read ',var_label,' into array u_ext.'
  END IF

  IF( ivout == 1  ) THEN
    READ (nunit,ERR=910,END=920) var_label
    READ (nunit,ERR=910,END=920) v_ext
    WRITE(6,'(a,a,a)')'Read ',var_label,' into array v_ext.'
  END IF

  IF( ipout == 1  ) THEN
    READ (nunit,ERR=910,END=920) var_label
    READ (nunit,ERR=910,END=920) p_ext
    WRITE(6,'(a,a,a)')'Read ',var_label,' into array p_ext.'
  END IF

  IF( ihout == 1  ) THEN
    READ (nunit,ERR=910,END=920) var_label
    READ (nunit,ERR=910,END=920) hgt_ext
    WRITE(6,'(a,a,a)')'Read ',var_label,' into array hgt_ext.'
  END IF

  IF( itout == 1  ) THEN
    READ (nunit,ERR=910,END=920) var_label
    READ (nunit,ERR=910,END=920) t_ext
    WRITE(6,'(a,a,a)')'Read ',var_label,' into array t_ext.'
  END IF

  IF( iqvout== 1  ) THEN
    READ (nunit,ERR=910,END=920) var_label
    READ (nunit,ERR=910,END=920) qv_ext
    WRITE(6,'(a,a,a)')'Read ',var_label,' into array qv_ext.'
  END IF

  IF( itsfcout==1 ) THEN
    READ (nunit,ERR=910,END=920) var_label
    READ (nunit,ERR=910,END=920) tsfc_ext
    WRITE(6,'(a,a,a)')'Read ',var_label,' into array tsfc_ext.'
  END IF

  IF( itsoilout==1) THEN
    READ (nunit,ERR=910,END=920) var_label
    READ (nunit,ERR=910,END=920) tdeep_ext
    WRITE(6,'(a,a,a)')'Read ',var_label,' into array tdeep_ext.'
  END IF

  IF( iwsfcout ==1) THEN
    READ (nunit,ERR=910,END=920) var_label
    READ (nunit,ERR=910,END=920) wetsfc_ext
    WRITE(6,'(a,a,a)')'Read ',var_label,' into array wetsfc_ext.'
  END IF

  IF( iwdpout ==1 ) THEN
    READ (nunit,ERR=910,END=920) var_label
    READ (nunit,ERR=910,END=920) wetdp_ext
    WRITE(6,'(a,a,a)')'Read ',var_label,' into array wetdp_ext.'
  END IF


  IF( iwcnpout ==1) THEN
    READ (nunit,ERR=910,END=920) var_label
    READ (nunit,ERR=910,END=920) wetcanp_ext
    WRITE(6,'(a,a,a)')'Read ',var_label,' into array wetcanp_ext.'
  END IF

  IF( isnowout ==1) THEN
    READ (nunit,ERR=910,END=920) var_label
    READ (nunit,ERR=910,END=920) snowdpth_ext
    WRITE(6,'(a,a,a)')'Read ',var_label,' into array snowdpth_ext.'
  END IF

  IF( itrnout ==1 ) THEN
    READ (nunit,ERR=910,END=920) var_label
    READ (nunit,ERR=910,END=920) trn_ext
    WRITE(6,'(a,a,a)')'Read ',var_label,' into array trn_ext.'
  END IF

  IF( ipsfcout ==1) THEN
    READ (nunit,ERR=910,END=920) var_label
    READ (nunit,ERR=910,END=920) psfc_ext
    WRITE(6,'(a,a,a)')'Read ',var_label,' into array psfc_ext.'
  END IF

  IF( iugrdout ==1) THEN
    READ (nunit,ERR=910,END=920) var_label
    READ (nunit,ERR=910,END=920) temp_ext ! discarded
    WRITE(6,'(a,a,a)')'Read ',var_label,' into array temp_ext.'
  END IF

  IF( ivgrdout ==1) THEN
    READ (nunit,ERR=910,END=920) var_label
    READ (nunit,ERR=910,END=920) temp_ext ! discarded
    WRITE(6,'(a,a,a)')'Read ',var_label,' into array temp_ext.'
  END IF

  IF( itgrdout ==1) THEN
    READ (nunit,ERR=910,END=920) var_label
    READ (nunit,ERR=910,END=920) temp_ext ! discarded
    WRITE(6,'(a,a,a)')'Read ',var_label,' into array temp_ext.'
  END IF

  IF( irhgrdout ==1) THEN
    READ (nunit,ERR=910,END=920) var_label
    READ (nunit,ERR=910,END=920) temp_ext ! discarded
    WRITE(6,'(a,a,a)')'Read ',var_label,' into array temp_ext.'
  END IF

  IF( iptgrdout ==1) THEN
    READ (nunit,ERR=910,END=920) var_label
    READ (nunit,ERR=910,END=920) temp_ext ! discarded
    WRITE(6,'(a,a,a)')'Read ',var_label,' into array temp_ext.'
  END IF

  IF( ipmslout ==1) THEN
    READ (nunit,ERR=910,END=920) var_label
    READ (nunit,ERR=910,END=920) temp_ext ! discarded
    WRITE(6,'(a,a,a)')'Read ',var_label,' into array temp_ext.'
  END IF

  CLOSE (UNIT=nunit)
  CALL retunit( nunit)

  PRINT*,'Finished reading file ',trim(gribfile)
  PRINT*,' '

  900 CONTINUE

  dx_ext = del_lon
  dy_ext = del_lat

  IF( lonend > 180.0 ) THEN
    lonend = lonend - 360.0
  END IF

  IF( lonbgn > 180.0 ) THEN
    lonbgn = lonbgn - 360.0
  END IF

  IF( lonbgn > lonend ) lonbgn = lonbgn - 360.0

  DO i=1,nx_ext
    DO j=1,ny_ext
      lon_ext(i,j)= lonbgn + (i-1) * dx_ext
      lat_ext(i,j)= latbgn + (j-1) * dy_ext - 90.0
    END DO
  END DO

  istatus = 1

  RETURN
!
!-----------------------------------------------------------------------
!
!  Error during read
!
!----------------------------------------------------------------------
!

  910   CONTINUE
  WRITE(6,'(/a/)') ' Error reading data in GET_AVN_BIN.'
  istatus =2
  RETURN
!
!-----------------------------------------------------------------------
!
!  End-of-file during read
!
!----------------------------------------------------------------------
!

  920   CONTINUE
  WRITE(6,'(/a/)') ' End of file reached in GET_AVN_BIN.'
  istatus =3
  RETURN
END SUBROUTINE get_avn_bin

!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE GETNCEPGFS                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE getncepgfs(nx_ext,ny_ext,nz_ext,nzsoil_ext,dir_extd,extdname,&
                      extdopt,extdfmt,nofixdim,lon_0_360,               &
                      extdinit,extdfcst,julfname,                       &
                      iproj_ext,scale_ext,                              &
                      trlon_ext,latnot_ext,x0_ext,y0_ext,               &
                      lat_ext,lon_ext,zpsoil_ext,                       &
                      p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,           &
                      qc_ext,                                           &
                      tsoil_ext,qsoil_ext,wetcanp_ext,                  &
                      snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext,        &
                      t_2m_ext, qv_2m_ext, u_10m_ext, v_10m_ext,        &
                      rain_ext,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads a NCEP GFS 0.5 degree global data in either GRIB format or
!  GRIB2 format for processing by ext2arps.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  08/30/2010
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    dir_extd      Directory name for external file
!    extdname      Prefix string of external file name
!    extdopt       Option of external data sources
!    extdfmt       Option of external data format
!
!    extdinit      Initialized time in mm-dd-yyyy:hh:mm:ss format
!    extdfcst      Forecast hour in HHH:MM:SS format
!    julfname      File name in yyjjjhhmm format
!
!  OUTPUT:
!
!    iproj_ext     Map projection number of external data
!    scale_ext     Scale factor of external data
!    trlon_ext     True longitude of external data (degrees E)
!    latnot_ext(2) True latitude(s) of external data (degrees N)
!    x0_ext        x coordinate of origin of external data
!    y0_ext        y coordinate of origin of external data
!    lat_ext       latitude of external data points (degrees N)
!    lon_ext       longitude of external data points (degrees E)
!    p_ext         pressure (Pascal)
!    hgt_ext       height (m)
!    t_ext         temperature (K)
!    qv_ext        specific humidity (kg/kg)
!    u_ext         u wind component (m/s)
!    v_ext         v wind component (m/s)
!    qc_ext        Cloud water mixing ratio (kg/kg)
!    qr_ext        Rain water mixing ratio (kg/kg)
!    qi_ext        Ice mixing ratio (kg/kg)
!    qs_ext        Snow mixing ratio (kg/kg)
!    qh_ext        Hail mixing ratio (kg/kg)
!
!    tsoil_ext     Soil temperature
!    qsoil_ext     Soil moisture
!    wetcanp_ext   Water content on canopy
!
!    trn_ext       External terrain (m)
!    psfc_ext      Surface pressure (Pa)
!
!    T_2m_ext      Temperature at 2m AGL
!    qv_2m_ext     Specific Humidity at 2m AGL
!    U_10m_ext     U at 10m AGL
!    V_10m_ext     V at 10m AGL
!
!    istatus       status indicator
!
!  WORK ARRAYS:
!
!    var_grb3d     Arrays to store the GRIB 3-D variables:
!
!    var_grb2d     Arrays to store the GRIB 2-D variables:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx_ext,ny_ext,nz_ext,nzsoil_ext

  CHARACTER (LEN=*)   :: dir_extd
  CHARACTER (LEN=*)   :: extdname

  INTEGER, INTENT(IN) :: extdopt, extdfmt
  INTEGER, INTENT(IN) :: nofixdim
  LOGICAL, INTENT(IN) :: lon_0_360

  CHARACTER (LEN=19) :: extdinit
  CHARACTER (LEN=9)  :: extdfcst
  CHARACTER (LEN=9)  :: julfname
!
!-----------------------------------------------------------------------
!
!  External grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj_ext
  REAL    :: scale_ext,trlon_ext
  REAL    :: latnot_ext(2)
  REAL    :: x0_ext,y0_ext
  REAL    :: dx_ext,dy_ext
!
!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------
!
  REAL :: lat_ext(nx_ext,ny_ext)
  REAL :: lon_ext(nx_ext,ny_ext)
  REAL :: zpsoil_ext(nx_ext,ny_ext,nzsoil_ext)
  REAL :: p_ext  (nx_ext,ny_ext,nz_ext)   ! Pressure (Pascals)
  REAL :: hgt_ext(nx_ext,ny_ext,nz_ext)   ! Height (m)
  REAL :: t_ext  (nx_ext,ny_ext,nz_ext)   ! Temperature (K)
  REAL :: qv_ext (nx_ext,ny_ext,nz_ext)   ! Specific humidity (kg/kg)
  REAL :: u_ext  (nx_ext,ny_ext,nz_ext)   ! Eastward wind component
  REAL :: v_ext  (nx_ext,ny_ext,nz_ext)   ! Northward wind component
  REAL :: qc_ext (nx_ext,ny_ext,nz_ext)   ! Cloud H2O mixing ratio (kg/kg)
  !REAL :: qr_ext (nx_ext,ny_ext,nz_ext)   ! Rain  H2O mixing ratio (kg/kg)
  !REAL :: qi_ext (nx_ext,ny_ext,nz_ext)   ! Ice   mixing ratio (kg/kg)
  !REAL :: qs_ext (nx_ext,ny_ext,nz_ext)   ! Snow  mixing ratio (kg/kg)
  !REAL :: qh_ext (nx_ext,ny_ext,nz_ext)   ! Hail  mixing ratio (kg/kg)

  REAL :: tsoil_ext  (nx_ext,ny_ext,nzsoil_ext)      ! Temperature at surface (K)
  REAL :: qsoil_ext  (nx_ext,ny_ext,nzsoil_ext)      ! Deep soil temperature (K)
  REAL :: wetcanp_ext(nx_ext,ny_ext)                 ! Canopy water amount
  REAL :: snowdpth_ext(nx_ext,ny_ext)                ! Snow depth (m)

  REAL :: trn_ext    (nx_ext,ny_ext)      ! External terrain (m)
  REAL :: psfc_ext   (nx_ext,ny_ext)      ! Surface pressure (Pa)
  REAL :: t_2m_ext (nx_ext,ny_ext)
  REAL :: qv_2m_ext(nx_ext,ny_ext)
  REAL :: u_10m_ext(nx_ext,ny_ext)
  REAL :: v_10m_ext(nx_ext,ny_ext)
  REAL :: rain_ext (nx_ext,ny_ext)

  INTEGER :: soiltyp_ext(nx_ext,ny_ext)

  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Work arrays for storing grib data
!
!-----------------------------------------------------------------------
!
  REAL,    ALLOCATABLE :: var_grb2d(:,:,:,:)   ! GRIB 2D variables
  REAL,    ALLOCATABLE :: var_grb3d(:,:,:,:,:) ! GRIB 3-D variables
  INTEGER, ALLOCATABLE :: var_lev3d(:,:,:)   ! Levels (hybrid) for
                                             ! each 3-D variable
  REAL,    ALLOCATABLE :: rcdata(:)            ! temporary data array
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: gribfile
  CHARACTER (LEN=14)  :: gribtime
  INTEGER :: grbflen
  INTEGER :: grbtlen

  INTEGER :: i,j,k,kk
  INTEGER :: m,n,max_nr2d,max_nr3d
  INTEGER :: ii(nx_ext)
  REAL    :: rtmp
  INTEGER :: igrbfmt

!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'gribcst.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  GRIB grid information
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=42) :: gridesc ! Grid description

  INTEGER :: iproj_grb    ! Map projection indicator
  INTEGER :: gthin        ! Indicator of whether the grid is "thinned"

  INTEGER :: ni_grb       ! Number of points along x-axis
  INTEGER :: nj_grb       ! Number of points along y-axis
  INTEGER :: np_grb       ! Total number of horizontal grid points

  INTEGER :: nk_grb       ! Number of vertical parameters
  REAL :: zk_grb(nz_ext)  ! Vertical coordinate parameters

  INTEGER :: npeq         ! Number of lat circles from pole to equator
  INTEGER :: nit(nz_ext)  ! Number of x-points for thinned grid

  REAL :: pi_grb          ! x-coordinate of pole point
  REAL :: pj_grb          ! y-coordinate of pole point
  INTEGER :: ipole        ! Projection center flag

  REAL :: di_grb          ! x-direction increment or grid length
  REAL :: dj_grb          ! y-direction increment or grid length

  REAL :: latsw           ! Latitude  of South West corner point
  REAL :: lonsw           ! Longitude of South West corner point
  REAL :: latne           ! Latitude  of North East corner point
  REAL :: lonne           ! Longitude of North East corner point

  REAL :: lattru1         ! Latitude (1st) at which projection is true
  REAL :: lattru2         ! Latitude (2nd) at which projection is true
  REAL :: lontrue         ! Longitude      at which projection is true

  REAL :: latrot          ! Latitude  of southern pole of rotation
  REAL :: lonrot          ! Longitude of southern pole of rotation
  REAL :: angrot          ! Angle of rotation

  REAL :: latstr          ! Latitude  of the pole of stretching
  REAL :: lonstr          ! Longitude of the pole of stretching
  REAL :: facstr          ! Stretching factor

  INTEGER :: scanmode     ! Scanning indicator
  INTEGER :: iscan        ! x-direction   scanning indicator
  INTEGER :: jscan        ! y-direction   scanning indicator
  INTEGER :: kscan        ! FORTRAN index scanning indicator

  INTEGER :: ires         ! Resolution direction increments indicator
  INTEGER :: iearth       ! Earth shape indicator: spherical or oblate?
  INTEGER :: icomp        ! (u,v) components decomposition indicator

  INTEGER :: jpenta       ! J-Pentagonal resolution parameter
  INTEGER :: kpenta       ! K-Pentagonal resolution parameter
  INTEGER :: mpenta       ! M-Pentagonal resolution parameter
  INTEGER :: ispect       ! Spectral representation type
  INTEGER :: icoeff       ! Spectral coefficient storage mode

  REAL :: xp_grb          ! X coordinate of sub-satellite point
  REAL :: yp_grb          ! Y coordinate of sub-satellite point
  REAL :: xo_grb          ! X coordinate of image sector origin
  REAL :: yo_grb          ! Y coordinate of image sector origin
  REAL :: zo_grb          ! Camera altitude from center of Earth

!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: nzsoilin_ext = 4
  REAL,    PARAMETER :: soildepth_ext(nzsoilin_ext) = (/0.0, 0.1, 0.4, 1.0/)
                             ! Surface and two soil layers  0-10cm & 10-200cm
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  x0_ext = 0.0
  y0_ext = 0.0

!---------------------------------------------------------------------
!
!  Define soil depths
!
!---------------------------------------------------------------------

  IF(soilmodel_option == 1 .AND. nzsoil_ext /= 2) THEN
    WRITE(6,'(2a,I3,a/2a,I2,a/,a,a)')                                   &
               'NCEP GFS grid #4 only provides 4 soil layers, However,',&
               ' You are trying to extract ',nzsoil_ext, ' layers.',    &
               ' for ARPS two-layer force-restore model ',              &
               '(soilmodel_option = ',soilmodel_option,')',             &
               ' Please check the code ext2arps.f90 for NCEPGFS_grb2 grid.',&
               ' Terminating ...'
    CALL arpsstop('Wrong number soil layers.',1)
  ELSE IF( soilmodel_option == 2 .AND. nzsoil_ext /= nzsoilin_ext) THEN
    WRITE(6,'(2a,I3,a/2a,I2,a/,a,a)')                                   &
               'NCEP GFS grid #4 only provides 4 soil layers, However,',&
               ' You are trying to extract ',nzsoil_ext, ' layers.',    &
               ' for ARPS multi-layer OUSoil model ',                   &
               '(soilmodel_option = ',soilmodel_option,')',             &
               ' Please check the code ext2arps.f90 for NCEPGFS_grb2 grid.',&
               ' Terminating ...'
    CALL arpsstop('Wrong number of soil layers.',1)
  END IF

  IF(soilmodel_option == 1) THEN
    DO k=1,nzsoil_ext
      DO j=1,ny_ext
        DO i=1,nx_ext
          zpsoil_ext(i,j,k) = 0.5*(soildepth_ext(k)+soildepth_ext(k+1))
        END DO
      END DO
    END DO
  ELSE
    DO k=1,nzsoil_ext
      DO j=1,ny_ext
        DO i=1,nx_ext
          zpsoil_ext(i,j,k) = soildepth_ext(k)
        END DO
      END DO
    END DO
  END IF

!-----------------------------------------------------------------------
!
! Check GRIB file format
!
!-----------------------------------------------------------------------

  CALL getgrbfname(extdopt,extdfmt,dir_extd,extdname,                   &
                   extdinit,extdfcst,'    ',                            &
                   gribfile,grbflen,gribtime,grbtlen,istatus)
  IF (istatus /= 0) RETURN

  CALL chkgrb(gribfile(1:grbflen),grbflen,igrbfmt,istatus)
  IF (istatus /= 0) RETURN

  IF (igrbfmt == 2) THEN
    CALL getncepgfs_grb2(nx_ext,ny_ext,nz_ext,nzsoil_ext,nzsoilin_ext,  &
                      gribfile,grbflen, gribtime, soildepth_ext,        &
                      lon_0_360,iproj_ext,scale_ext,                    &
                      trlon_ext,latnot_ext,lat_ext,lon_ext,             &
                      p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,           &
                      qc_ext,                                           &
                      tsoil_ext,qsoil_ext,wetcanp_ext,                  &
                      snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext,        &
                      t_2m_ext, qv_2m_ext, u_10m_ext, v_10m_ext,        &
                      rain_ext,istatus)

    IF (istatus /= 0) RETURN
    istatus = 1

  ELSE

    ALLOCATE(var_grb2d(nx_ext,ny_ext,n2dvs,n2dlvt),        STAT = istatus)
    ALLOCATE(var_grb3d(nx_ext,ny_ext,nz_ext,n3dvs,n3dlvt), STAT = istatus)
    ALLOCATE(rcdata(nx_ext*ny_ext),          STAT = istatus)
    ALLOCATE(var_lev3d(nz_ext,n3dvs,n3dlvt), STAT = istatus)
!
!-----------------------------------------------------------------------
!
!  RDNMCGRB reads NMC GRIB data
!
!-----------------------------------------------------------------------
!
    gridtyp   = gfs04grid
    mproj_grb = gfs04proj

    n2dvars  = gfs04nvs2d
    n2dlvtps = gfs04nlvt2d

    DO k=1,n2dlvtps
      DO n=1,n2dvars
        var_id2d(n,k) = gfs04var_id2d(n,k)
      END DO
      levtyp2d(k) = gfs04levs2d(k)
    END DO

    n3dvars  = gfs04nvs3d
    n3dlvtps = gfs04nlvt3d

    DO m=1,n3dlvtps
      DO n=1,n3dvars
        var_id3d(n,m) = gfs04var_id3d(n,m)
      END DO
      levtyp3d(m) = gfs04levs3d(m)
    END DO

    CALL rdnmcgrb(nx_ext,ny_ext,nz_ext,gribfile,grbflen, gribtime,      &
                  gridesc, iproj_grb, gthin,                            &
                  ni_grb,nj_grb,np_grb, nk_grb,zk_grb, npeq,nit,        &
                  pi_grb,pj_grb,ipole, di_grb,dj_grb,                   &
                  latsw,lonsw, latne,lonne,                             &
                  latrot,lonrot,angrot,                                 &
                  latstr,lonstr,facstr,                                 &
                  lattru1,lattru2,lontrue,                              &
                  scanmode, iscan,jscan,kscan,                          &
                  ires,iearth,icomp,                                    &
                  jpenta,kpenta,mpenta,ispect,icoeff,                   &
                  xp_grb,yp_grb, xo_grb,yo_grb,zo_grb,                  &
                  rcdata,var_grb2d,var_grb3d,var_lev3d,lvldbg,istatus)

    IF (istatus /= 0)  THEN
      ISTATUS = -888
      GOTO 999
    END IF

    max_nr2d = 0
    DO n=1,n2dvars
      DO m=1,n2dlvtps
        max_nr2d = MAX( max_nr2d, var_nr2d(n,m) )
      END DO
    END DO

    max_nr3d = 0
    DO n=1,n3dvars
      max_nr3d = MAX( max_nr3d, var_nr3d(n,1) )
    END DO

    IF (max_nr3d /= 26) THEN
      WRITE(6,'(1x,I3,a,/,1x,a)')                                       &
        'The program is supposed to find 26 vertical levels, but it found ', max_nr3d,'.', &
        'Program terminating ...'
      istatus = -888
      GO TO 999
    ELSE
      max_nr3d = MIN(max_nr3d, 20)  ! Check the first 20 levels only
    END IF

    IF ( max_nr3d == 0 ) THEN
      WRITE (6,'(a)')                                                   &
          'No 3-D variable was found in the GRIB file',                 &
          'Dataset problem in GETNMCETA212.',                           &
          ' '
      !STOP
      istatus = -888
      GOTO 999
    END IF

    IF ( max_nr2d == 0 ) THEN
      WRITE (6,'(a)')                                                   &
          'No 2-D variables was found in the GRIB file'
    END IF

  !  WRITE (6,'(/a7,2x,6(i7))') 'Lev\VID',(var_id3d(n,1),n=1,n3dvars)
  !  DO k=1,max_nr3d
  !    write (6,'(/i5,4x,6(i7))') k,(var_lev3d(k,n,1),n=1,n3dvars)
  !   END DO

    DO k=1,max_nr3d
      DO n=2,n3dvars-1   ! Modified for tolerance of the EnKF dataset
        IF ( var_lev3d(k,1,1) /= var_lev3d(k,n,1) ) THEN
          WRITE (6,'(a,/,a,I0,a,I0,/,a,I0,a,I0,/,a,/)')                 &
          'Variables were not at the same level.',                      &
          'The K (',k,')th level for variable 1 is ',var_lev3d(k,1,1),  &
          'but it is ',var_lev3d(k,n,1),' for variable n = ',n,         &
          'Dataset problem in GETNCEPGFS GRIB1 format.'
          !STOP
          ISTATUS = -888
          GOTO 999
        END IF
      END DO
    END DO

    IF ( iproj_grb == 5 .AND. ipole == 0 ) THEN      ! Center in NP
      iproj_ext = 1
    ELSE IF ( iproj_grb == 5 .AND. ipole == 1 ) THEN  ! Center in SP
      iproj_ext = -1
    ELSE IF ( iproj_grb == 3 .AND. ipole == 0 ) THEN  ! Center in NP
      iproj_ext = 2
    ELSE IF ( iproj_grb == 3 .AND. ipole == 1 ) THEN  ! Center in SP
      iproj_ext = -2
    ELSE IF ( iproj_grb == 1 ) THEN
      iproj_ext = 3
    ELSE IF ( iproj_grb == 0 ) THEN
      iproj_ext = 4
    ELSE
      WRITE (6,'(a)')                                                   &
          'Unknown map projection. Set to non-projection.'
      iproj_ext = 0
    END IF

    scale_ext = 1.0
    latnot_ext(1) = lattru1
    latnot_ext(2) = lattru2
    trlon_ext = lontrue

    dx_ext = di_grb
    dy_ext = dj_grb

    DO j=1, ny_ext
      DO i=1, nx_ext
        lon_ext(i,j)= lonsw + (i-1) * dx_ext
        lat_ext(i,j)= latsw + (j-1) * dy_ext
      END DO
    END DO

    IF (myproc == 0) PRINT *,'LatSW = ',latsw,' LonSW = ',lonsw

    IF ( lon_0_360 ) THEN   ! Do not care whether median is in-between

      DO i = 1,nx_ext
        ii(i) = i
      END DO

    ELSE                    ! Median is in-between and 1:0 - 360:359
      IF (mod(nx_ext,2) /= 0) THEN
        WRITE(6,'(a/)') 'Wrong size of nx_ext in GETNCEPAVN3 for lon_0_360.'
        CALL arpsstop('Wrong nx_ext size.',1)
      END IF

      DO i = 1,nx_ext/2              ! map 1-360 to 361-720
        ii(i) = i + nx_ext/2
      END DO

      DO i = nx_ext/2+1, nx_ext      ! map 361-720 to 1-360
        ii(i) = i - nx_ext/2
      END DO

      WHERE (lon_ext >= 180) lon_ext = lon_ext - 360

      DO j = 1,ny_ext                ! swap lat_ext & lon_ext
        DO i = 1,nx_ext/2
          rtmp = lat_ext(ii(i),j)
          lat_ext(ii(i),j) = lat_ext(i,j)
          lat_ext(i,    j) = rtmp

          rtmp = lon_ext(ii(i),j)
          lon_ext(ii(i),j) = lon_ext(i,j)
          lon_ext(i,    j) = rtmp
        END DO
      END DO
    END IF

!
!-----------------------------------------------------------------------
!
!  Retrieve 2-D variables
!
!-----------------------------------------------------------------------
!

    WRITE(6,*) 'Filling 2D arrays from GRIB data.'
    DO j=1,ny_ext
      DO i=1,nx_ext
        IF (var_nr2d(1,1) == 0) THEN
          psfc_ext(i,j) = -999.0
        ELSE
          psfc_ext   (i,j) = var_grb2d(ii(i),j,1,1)
        END IF

        IF (var_nr2d(2,1) == 0) THEN
          trn_ext(i,j) = -999.0
        ELSE
          trn_ext   (i,j) = var_grb2d(ii(i),j,2,1)
        END IF

        wetcanp_ext(i,j) = -999.0
        !wetcanp_ext(i,j) = var_grb2d(ii(i),j,4,1)*1.e-3     ! in meter

        IF ( var_nr2d(4,1) == 0 ) THEN
          snowdpth_ext(i,j) = -999.
        ELSE
          ! Convert water equiv. of accum. snow depth (kg/m**2) to meters
          ! (where 1 meter liquid water is set equivqlent to 10 meters snow).
          !     0.01 = 10. (m snow/m liquid) / (1000 kg/m**3)
          snowdpth_ext(i,j) = 0.01 * var_grb2d(ii(i),j,4,1)  ! in meters
        END IF

        soiltyp_ext(i,j)= 0          ! We do not know the soil type
        IF (var_grb2d(ii(i),j,5,1) < 0.5) THEN   ! land cover = 0 for see
          soiltyp_ext(i,j) = 13
        END IF

        IF (var_grb2d(ii(i),j,6,1)  > 0.5) THEN  ! Ice cover = 1 for Ice
          soiltyp_ext(i,j) = 12
        END IF

        IF( var_nr2d(1,2) == 0 ) THEN
          t_2m_ext(i,j)= -999.0
        ELSE
          t_2m_ext(i,j)= var_grb2d(ii(i),j,1,2)
        END IF
        IF( var_nr2d(2,2) == 0 ) THEN
          qv_2m_ext(i,j)= -999.0
        ELSE
          qv_2m_ext(i,j)= var_grb2d(ii(i),j,2,2)
        END IF

        IF( var_nr2d(3,2) == 0 ) THEN
          u_10m_ext(i,j)= -999.0
        ELSE
          u_10m_ext(i,j)= var_grb2d(ii(i),j,3,2)
        END IF

        IF( var_nr2d(4,2) == 0 ) THEN
          v_10m_ext(i,j)= -999.0
        ELSE
          v_10m_ext(i,j)= var_grb2d(ii(i),j,4,2)
        END IF

        IF ( var_nr2d(1,3) == 0 ) THEN
          rain_ext(i,j)= -999.0
        ELSE
          rain_ext(i,j)= var_grb2d(ii(i),j,1,3)
        END IF
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Retrieve 3-D variables
!
!-----------------------------------------------------------------------
!
    WRITE(6,*) 'Filling 3D arrays from GRIB data.'
    DO k=1,nz_ext
      DO j=1,ny_ext
        DO i=1,nx_ext
          p_ext  (i,j,k) = 100.0 * REAL(var_lev3d(k,1,1)) ! Pressure
          hgt_ext(i,j,k) = var_grb3d(ii(i),j,k,1,1)
          t_ext  (i,j,k) = var_grb3d(ii(i),j,k,2,1)    ! Temperature (K)
          u_ext  (i,j,k) = var_grb3d(ii(i),j,k,3,1)    ! u wind (m/s)
          v_ext  (i,j,k) = var_grb3d(ii(i),j,k,4,1)    ! v wind (m/s)
          IF ( k > var_nr3d(6,1)) THEN
            qc_ext(i,j,k) = 0.0
          ELSE
            qc_ext (i,j,k) = var_grb3d(ii(i),j,k,6,1)    ! Cloud water
          END IF
          !qr_ext (i,j,k) = -999.
          !qi_ext (i,j,k) = -999.
          !qs_ext (i,j,k) = -999.
          !qh_ext (i,j,k) = -999.
        END DO
      END DO
    END DO

    WRITE(6,*) 'Computing saturated specific humidity.'

    CALL getqvs(nx_ext,ny_ext,nz_ext, 1,nx_ext,1,ny_ext,1,nz_ext,       &
                p_ext, t_ext, qv_ext )

    DO k=1,nz_ext
      DO j=1,ny_ext
        DO i=1,nx_ext
          IF ( k > var_nr3d(5,1)) THEN
            qv_ext(i,j,k) = 0.0
          ELSE
            qv_ext(i,j,k) = 0.01*var_grb3d(ii(i),j,k,5,1)*qv_ext(i,j,k)
          END IF
        END DO
      END DO
    END DO

    WRITE(6,*) 'Filling 3D Soil arrays from GRIB data.'

    tsoil_ext (:,:,:) = -999.0
    qsoil_ext (:,:,:) = -999.0

    IF(soilmodel_option == 1) THEN       ! ARPS: two-layer force-restore model)

      DO j = 1, ny_ext
        DO i = 1, nx_ext
          tsoil_ext(i,j,1) =     var_grb2d(ii(i),j,3,1)       ! sfc temp.
          tsoil_ext(i,j,2) = 0.1*var_grb3d(ii(i),j,1,1,2)               &
                           + 0.3*var_grb3d(ii(i),j,2,1,2)               &
                           + 0.6*var_grb3d(ii(i),j,3,1,2)

          qsoil_ext(i,j,1) =     var_grb3d(ii(i),j,1,2,2)
          qsoil_ext(i,j,2) = 0.1*var_grb3d(ii(i),j,1,2,2)               &
                           + 0.3*var_grb3d(ii(i),j,2,2,2)               &
                           + 0.6*var_grb3d(ii(i),j,3,2,2)

         IF (tsoil_ext(i,j,1) < 100. .OR. tsoil_ext(i,j,2) < 100.) THEN  ! Over water?
         !IF (soiltyp_ext(i,j) == 13) THEN   ! over see
           tsoil_ext(i,j,1) = var_grb2d(ii(i),j,3,1)  ! use ground (water surface) temperature
           tsoil_ext(i,j,2) = var_grb2d(ii(i),j,3,1)
           qsoil_ext(i,j,1) = 1.0
           qsoil_ext(i,j,2) = 1.0
         END IF
        END DO
      END DO

    ELSE                                ! ARPS: multi-layer OUSoil model
      DO k = 1,nzsoil_ext
        DO j = 1, ny_ext
          DO i = 1, nx_ext
            tsoil_ext(i,j,k) = var_grb3d(ii(i),j,k,1,2)
            qsoil_ext(i,j,k) = var_grb3d(ii(i),j,k,2,2)
            IF (tsoil_ext(i,j,k) < 100.) THEN  ! Over water?
            !IF (soiltyp_ext(i,j) == 13) THEN  ! over water
              tsoil_ext(i,j,k) = var_grb2d(ii(i),j,3,1)
              qsoil_ext(i,j,k) = 1.0
            END IF
          END DO
        END DO
      END DO

    END IF

    istatus = 1

    999 CONTINUE
    DEALLOCATE(var_grb2d,var_grb3d,rcdata,var_lev3d)

  END IF

  RETURN
END SUBROUTINE getncepgfs
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETNNARR221                ######
!######                                                      ######
!######                     Developed by                     ######
!######              Atmospheric Science Division            ######
!######        Lawrence Livermore National Laboratory        ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getnarr221(nx_ext,ny_ext,nz_ext,nzsoil_ext,                  &
                   nstyps_ext,dir_extd,extdname,extdopt,extdfmt,        &
                   extdinit,extdfcst,julfname,                          &
                   iproj_ext,scale_ext,                                 &
                   trlon_ext,latnot_ext,x0_ext,y0_ext,                  &
                   lat_ext,lon_ext,zpsoil_ext,                          &
                   p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,              &
                   qc_ext,qi_ext,                                       &
                   tsoil_ext,qsoil_ext,wetcanp_ext,                     &
                   snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext,           &
                   t_2m_ext, qv_2m_ext, u_10m_ext, v_10m_ext, rain_ext, &
                   istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads a NMC GRIB (Grid #221, 32km) ETA file from the North American
!  Regional Reanalysis (NARR) data for processing by ext2arps, a program
!  that converts external files to ARPS variables and format.
!  The data are in GRIB version 1 format.
!
!  NOTE:
!
!   The description about Eta 32km NARR output is available at
!   http://wwwt.emc.ncep.noaa.gov/mmb/rreanl/index.html
!
!   The data can be downloaded from
!   http://nomads.ncdc.noaa.gov/#narr_datasets
!
!   This subroutine assumes that data has been downloaded to local
!   disk and stored in "dir_extd"/"extdname"_YYYYMMDD_HHHH_000.grb
!
!   where:
!
!     o "dir_extd" is specified inside arps.input;
!     o "extdname" supposed to be "narr-a_221";
!     o CC is the model cycle, extracted from "extdinit" (00,06,12,18);
!     o HHHH is the forecast hour, extracted from "extdfcst";
!     o YYYY is the year
!     o MM is the month
!     o DD is the day
!   e.g. narr-a_221_20030701_0300_000.grb
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Tina Katopodes Chow (04/18/2005)
!  modified from getnmceta212 to create getnmceta221
!
!  MODIFICATION HISTORY:
!
!  2005/06/23 (Y. Wang)
!  Incorporated into ARPS package and renamed to getnarr221 as well did
!  some tests.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    dir_extd      Directory name for external file
!    extdname      Prefix string of external file name
!    extdopt       Option of external data sources
!    extdfmt       Option of external data format
!
!    extdinit      Initialized time in mm-dd-yyyy:hh:mm:ss format
!    extdfcst      Forecast hour in HHH:MM:SS format, does not used here
!    julfname      File name in yyjjjhhmm format
!
!  OUTPUT:
!
!    iproj_ext     Map projection number of external data
!    scale_ext     Scale factor of external data
!    trlon_ext     True longitude of external data (degrees E)
!    latnot_ext(2) True latitude(s) of external data (degrees N)
!    x0_ext        x coordinate of origin of external data
!    y0_ext        y coordinate of origin of external data
!    lat_ext       latitude of external data points (degrees N)
!    lon_ext       longitude of external data points (degrees E)
!    p_ext         pressure (Pascal)
!    hgt_ext       height (m)
!    t_ext         temperature (K)
!    qv_ext        specific humidity (kg/kg)
!    u_ext         u wind component (m/s)
!    v_ext         v wind component (m/s)
!    qc_ext        Cloud water mixing ratio (kg/kg)
!    qr_ext        Rain water mixing ratio (kg/kg)
!    qi_ext        Ice mixing ratio (kg/kg)
!    qs_ext        Snow mixing ratio (kg/kg)
!    qh_ext        Hail mixing ratio (kg/kg)
!
!    tsoil_ext     Soil temperature (K)
!    qsoil_ext     Soil moisture (m**3/m**3)
!    wetcanp_ext   Water content on canopy
!
!    trn_ext       External terrain (m)
!    psfc_ext      Surface pressure (Pa)
!
!    T_2m_ext      Temperature at 2m AGL
!    qv_2m_ext     Specific Humidity at 2m AGL
!    U_10m_ext     U at 10m AGL
!    V_10m_ext     V at 10m AGL
!
!    istatus       status indicator
!
!  WORK ARRAYS:
!
!    var_grb3d     Arrays to store the GRIB 3-D variables:
!                  var_grb3d(nxgrb,nygrb,nzgrb,1,1) - Temperature (K)
!                  var_grb3d(nxgrb,nygrb,nzgrb,2,1) - Specific humidity
!                                                     (kg/kg)
!                  var_grb3d(nxgrb,nygrb,nzgrb,3,1) - u wind (m/s)
!                  var_grb3d(nxgrb,nygrb,nzgrb,4,1) - v wind (m/s)
!                  var_grb3d(nxgrb,nygrb,nzgrb,5,1) - Geopotential
!                                                     height (gpm)
!                  var_grb3d(nxgrb,nygrb,nzgrb,6,1) - Pressure vertical
!                                                     velocity (Pa/s)
!                                                     (if applied)
!                  var_grb3d(nxgrb,nygrb,nzgrb,1,2) - soil temp. (K)
!                  var_grb3d(nxgrb,nygrb,nzgrb,2,2) - vol. soil moist.
!                                                     (m**3/m**3)
!
!    var_grb2d     Arrays to store the GRIB 2-D variables:
!                  var_grb2d(nxgrb,nygrb,1) - Surface pressure (Pa)
!                  var_grb2d(nxgrb,nygrb,2) - Geopotential height (gpm)
!                  var_grb2d(nxgrb,nygrb,3) - Surface temperature (K)
!                  var_grb2d(nxgrb,nygrb,4) - Plant canopy surface
!                                             water (kg/m**2)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'gribcst.inc'

  CHARACTER (LEN=*) :: dir_extd
  CHARACTER (LEN=*) :: extdname

  INTEGER :: extdopt
  INTEGER :: extdfmt

  CHARACTER (LEN=19) :: extdinit
  CHARACTER (LEN=9)  :: extdfcst  ! does not use
  CHARACTER (LEN=9)  :: julfname  ! does not use
!
!-----------------------------------------------------------------------
!
!  External grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj_ext
  REAL    :: scale_ext,trlon_ext
  REAL    :: latnot_ext(2)
  REAL    :: x0_ext,y0_ext
  REAL    :: dx_ext,dy_ext
!
!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx_ext,ny_ext,nz_ext
  INTEGER :: nzsoil_ext, nstyps_ext

  REAL :: lat_ext(nx_ext,ny_ext)
  REAL :: lon_ext(nx_ext,ny_ext)
  REAL :: p_ext  (nx_ext,ny_ext,nz_ext)   ! Pressure (Pascals)
  REAL :: hgt_ext(nx_ext,ny_ext,nz_ext)   ! Height (m)
  REAL :: t_ext  (nx_ext,ny_ext,nz_ext)   ! Temperature (K)
  REAL :: qv_ext (nx_ext,ny_ext,nz_ext)   ! Specific humidity (kg/kg)
  REAL :: u_ext  (nx_ext,ny_ext,nz_ext)   ! Eastward wind component
  REAL :: v_ext  (nx_ext,ny_ext,nz_ext)   ! Northward wind component
  REAL :: qc_ext (nx_ext,ny_ext,nz_ext)   ! Cloud H2O mixing ratio (kg/kg)
!  REAL :: qr_ext (nx_ext,ny_ext,nz_ext)   ! Rain  H2O mixing ratio (kg/kg)
  REAL :: qi_ext (nx_ext,ny_ext,nz_ext)   ! Ice   mixing ratio (kg/kg)
!  REAL :: qs_ext (nx_ext,ny_ext,nz_ext)   ! Snow  mixing ratio (kg/kg)
!  REAL :: qh_ext (nx_ext,ny_ext,nz_ext)   ! Hail  mixing ratio (kg/kg)
  REAL :: zpsoil_ext(nx_ext,ny_ext,nzsoil_ext) !Soil depths (m)

  REAL :: tsoil_ext (nx_ext,ny_ext,nzsoil_ext) ! Soil temperature (K)
  REAL :: qsoil_ext (nx_ext,ny_ext,nzsoil_ext) ! Soil moisture (m**3/m**3)
  REAL :: wetcanp_ext(nx_ext,ny_ext)      ! Canopy water amount
  REAL :: snowdpth_ext(nx_ext,ny_ext)     ! Snow depth (m)

  REAL :: trn_ext    (nx_ext,ny_ext)      ! External terrain (m)
  REAL :: psfc_ext   (nx_ext,ny_ext)      ! Surface pressure (Pa)

  INTEGER soiltyp_ext (nx_ext,ny_ext)     ! Soil type

  REAL :: t_2m_ext (nx_ext,ny_ext)
  REAL :: qv_2m_ext(nx_ext,ny_ext)
  REAL :: u_10m_ext(nx_ext,ny_ext)
  REAL :: v_10m_ext(nx_ext,ny_ext)
  REAL :: rain_ext(nx_ext,ny_ext)
!
!-----------------------------------------------------------------------
!
!  Other  external variable arrays
!
!-----------------------------------------------------------------------
!
  REAL :: x_ext(nx_ext)
  REAL :: y_ext(ny_ext)

  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Work arrays for storing grib data
!
!-----------------------------------------------------------------------
!
  REAL,    ALLOCATABLE :: var_grb2d(:,:,:,:)   ! GRIB variables
  REAL,    ALLOCATABLE :: var_grb3d(:,:,:,:,:) ! GRIB 3-D variables
  INTEGER, ALLOCATABLE :: var_lev3d(:,:,:)     ! Levels (hybrid) for
                                               ! each 3-D variable
  REAL,    ALLOCATABLE :: rcdata(:)            ! temporary data array
!
!-----------------------------------------------------------------------
!
!  Original grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj
  REAL    :: scale,trlon,x0,y0
  REAL    :: latnot(2)
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: gribfile
  CHARACTER (LEN=14)  :: gribtime
  INTEGER :: i,j,k,kk
  INTEGER :: grbflen, grbtlen, len1

  INTEGER :: m,n,nz1,max_nr2d,max_nr3d

  REAL :: govrd
  REAL :: soildepth_ext(5)  ! EMK 15 June 2002

  INTEGER :: chklev, lvscan

  INTEGER :: iret             ! Return flag
  INTEGER :: terrainflag
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'grid.inc'
!
!-----------------------------------------------------------------------
!
!  GRIB grid information
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=42) :: gridesc ! Grid description

  INTEGER :: iproj_grb    ! Map projection indicator
  INTEGER :: gthin        ! Indicator of whether the grid is "thinned"

  INTEGER :: ni_grb       ! Number of points along x-axis
  INTEGER :: nj_grb       ! Number of points along y-axis
  INTEGER :: np_grb       ! Total number of horizontal grid points

  INTEGER :: nk_grb       ! Number of vertical parameters
  REAL :: zk_grb(nz_ext)  ! Vertical coordinate parameters

  INTEGER :: npeq         ! Number of lat circles from pole to equator
  INTEGER :: nit(nz_ext)  ! Number of x-points for thinned grid

  REAL :: pi_grb          ! x-coordinate of pole point
  REAL :: pj_grb          ! y-coordinate of pole point
  INTEGER :: ipole        ! Projection center flag

  REAL :: di_grb          ! x-direction increment or grid length
  REAL :: dj_grb          ! y-direction increment or grid length

  REAL :: latsw           ! Latitude  of South West corner point
  REAL :: lonsw           ! Longitude of South West corner point
  REAL :: latne           ! Latitude  of North East corner point
  REAL :: lonne           ! Longitude of North East corner point

  REAL :: lattru1         ! Latitude (1st) at which projection is true
  REAL :: lattru2         ! Latitude (2nd) at which projection is true
  REAL :: lontrue         ! Longitude      at which projection is true

  REAL :: latrot          ! Latitude  of southern pole of rotation
  REAL :: lonrot          ! Longitude of southern pole of rotation
  REAL :: angrot          ! Angle of rotation

  REAL :: latstr          ! Latitude  of the pole of stretching
  REAL :: lonstr          ! Longitude of the pole of stretching
  REAL :: facstr          ! Stretching factor

  INTEGER :: scanmode     ! Scanning indicator
  INTEGER :: iscan        ! x-direction   scanning indicator
  INTEGER :: jscan        ! y-direction   scanning indicator
  INTEGER :: kscan        ! FORTRAN index scanning indicator

  INTEGER :: ires         ! Resolution direction increments indicator
  INTEGER :: iearth       ! Earth shape indicator: spherical or oblate?
  INTEGER :: icomp        ! (u,v) components decomposition indicator

  INTEGER :: jpenta       ! J-Pentagonal resolution parameter
  INTEGER :: kpenta       ! K-Pentagonal resolution parameter
  INTEGER :: mpenta       ! M-Pentagonal resolution parameter
  INTEGER :: ispect       ! Spectral representation type
  INTEGER :: icoeff       ! Spectral coefficient storage mode

  REAL :: xp_grb          ! X coordinate of sub-satellite point
  REAL :: yp_grb          ! Y coordinate of sub-satellite point
  REAL :: xo_grb          ! X coordinate of image sector origin
  REAL :: yo_grb          ! Y coordinate of image sector origin
  REAL :: zo_grb          ! Camera altitude from center of Earth

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ALLOCATE(var_grb2d(nx_ext,ny_ext,n2dvs,n2dlvt))
  ALLOCATE(var_grb3d(nx_ext,ny_ext,nz_ext,n3dvs,n3dlvt))
  ALLOCATE(rcdata(nx_ext*ny_ext))
  ALLOCATE(var_lev3d(nz_ext,n3dvs,n3dlvt))

  ! Define depths of soil layers
  DATA soildepth_ext/0.0,-0.05,-0.25,-0.70,-1.50/

  CALL getgrbfname(extdopt,extdfmt,dir_extd,extdname,                   &
                   extdinit,extdfcst,'    ',                            &
                   gribfile,grbflen,gribtime,grbtlen,istatus)
  IF (istatus /= 0 ) RETURN

!
!-----------------------------------------------------------------------
!
!  RDNMCGRB reads NMC GRIB data
!
!-----------------------------------------------------------------------
!
  govrd = g/rd

  !Get grid type and projection type
  gridtyp  = narr221grid
  mproj_grb = narr221proj

  !Get number and names of 2d variables
  n2dvars  = narr221nvs2d
  n2dlvtps = narr221nlvt2d

  DO k=1,n2dlvtps
    DO n=1,n2dvars
      var_id2d(n,k) = narr221var_id2d(n,k)
    END DO
    levtyp2d(k) = narr221levs2d(k)
  END DO

  !Get number and names of 3d variables
  n3dvars  = narr221nvs3d
  n3dlvtps = narr221nlvt3d

  DO m=1,n3dlvtps
    DO n=1,n3dvars
      var_id3d(n,m) = narr221var_id3d(n,m)
    END DO
    levtyp3d(m) = narr221levs3d(m)
  END DO

  !Get original map projection for restore purpose
  CALL getmapr(iproj,scale,latnot,trlon,x0,y0)

  len1 = LEN_TRIM(dir_extd)

  !Read fixed field file
  CALL rdnmcgrb(nx_ext,ny_ext,nz_ext,dir_extd(1:len1)//'AWIP32.fixed',  &
                len1+12,'1979110800f00',                                &
                gridesc, iproj_grb, gthin,                              &
                ni_grb,nj_grb,np_grb, nk_grb,zk_grb, npeq,nit,          &
                pi_grb,pj_grb,ipole, di_grb,dj_grb,                     &
                latsw,lonsw, latne,lonne,                               &
                latrot,lonrot,angrot,                                   &
                latstr,lonstr,facstr,                                   &
                lattru1,lattru2,lontrue,                                &
                scanmode, iscan,jscan,kscan,                            &
                ires,iearth,icomp,                                      &
                jpenta,kpenta,mpenta,ispect,icoeff,                     &
                xp_grb,yp_grb, xo_grb,yo_grb,zo_grb,                    &
                rcdata,var_grb2d,var_grb3d,var_lev3d,lvldbg,iret)

  IF (iret /= 0)  THEN
    ISTATUS = -888
    GOTO 999
  END IF

  max_nr2d = 0
  DO n=1,n2dvars
    DO m=1,n2dlvtps
      max_nr2d = MAX( max_nr2d, var_nr2d(n,m) )
    END DO
  END DO

  IF ( max_nr2d == 0 ) THEN
    WRITE (6,'(a)')                                                     &
        'No 2-D fixed field was found in the GRIB file, AWIP32.fixed.'
  END IF

  terrainflag = var_nr2d(2,1)

  !Read grib data
  print *,'Reading gribfile = ', TRIM(gribfile)
  CALL rdnmcgrb(nx_ext,ny_ext,nz_ext,gribfile,grbflen, gribtime,        &
                gridesc, iproj_grb, gthin,                              &
                ni_grb,nj_grb,np_grb, nk_grb,zk_grb, npeq,nit,          &
                pi_grb,pj_grb,ipole, di_grb,dj_grb,                     &
                latsw,lonsw, latne,lonne,                               &
                latrot,lonrot,angrot,                                   &
                latstr,lonstr,facstr,                                   &
                lattru1,lattru2,lontrue,                                &
                scanmode, iscan,jscan,kscan,                            &
                ires,iearth,icomp,                                      &
                jpenta,kpenta,mpenta,ispect,icoeff,                     &
                xp_grb,yp_grb, xo_grb,yo_grb,zo_grb,                    &
                rcdata,var_grb2d,var_grb3d,var_lev3d,lvldbg,iret)

  IF (iret /= 0)  THEN
    ISTATUS = -888
    GOTO 999
  END IF

  DO n=1,n2dvars
    DO m=1,n2dlvtps
      max_nr2d = MAX( max_nr2d, var_nr2d(n,m) )
    END DO
  END DO

  max_nr3d = 0
  DO n=1,n3dvars
    max_nr3d = MAX( max_nr3d, var_nr3d(n,1) )
  END DO

  IF ( max_nr3d == 0 ) THEN
    WRITE (6,'(a)')                                                     &
        'No 3-D variable was found in the GRIB file.',                  &
        'Dataset problem in GETNARR221.',                               &
        ' '
    ISTATUS = -888
    GOTO 999
  END IF

  IF ( max_nr2d == 0 ) THEN
    WRITE (6,'(a)')                                                     &
        'No 2-D variables was found in the GRIB file.'
  END IF

!  WRITE (6,'(/a7,2x,6(i7))')       &
!         'Lev\\VID',(var_id3d(n,1),n=1,n3dvars)

!  DO k=1,max_nr3d
!    WRITE (6,'(/i5,4x,6(i7))')     &
!         k,(var_lev3d(k,n,1),n=1,n3dvars)
!  END DO

  DO k=1,max_nr3d
    DO n=2,n3dvars
      IF ( var_lev3d(k,1,1) /= var_lev3d(k,n,1) ) THEN
        WRITE (6,'(a)')                                                 &
            'Variables were not at the same level.',                    &
            'Dataset problem in GETNARR221.',                           &
            ' '
        ISTATUS = -888
        GOTO 999
      END IF
    END DO
  END DO

  !Set projection and grid spacing
  IF ( iproj_grb == 5 .AND. ipole == 0 ) THEN      ! Center in NP
    iproj_ext = 1
  ELSE IF ( iproj_grb == 5 .AND. ipole == 1 ) THEN  ! Center in SP
    iproj_ext = -1
  ELSE IF ( iproj_grb == 3 .AND. ipole == 0 ) THEN  ! Center in NP
    iproj_ext = 2
  ELSE IF ( iproj_grb == 3 .AND. ipole == 1 ) THEN  ! Center in SP
    iproj_ext = -2
  ELSE IF ( iproj_grb == 1 ) THEN
    iproj_ext = 3
  ELSE IF ( iproj_grb == 0 ) THEN
    iproj_ext = 4
  ELSE
    WRITE (6,'(a)') 'Unknown map projection. Set to non-projection.'
    iproj_ext = 0
  END IF

  scale_ext = 1.0
  latnot_ext(1) = lattru1
  latnot_ext(2) = lattru2
  trlon_ext = lontrue

  dx_ext = di_grb
  dy_ext = dj_grb

  CALL setmapr(iproj_ext,scale_ext,latnot_ext,trlon_ext)
  CALL lltoxy(1,1,latsw,lonsw,x0_ext,y0_ext)

  DO i=1,nx_ext
    x_ext(i)=x0_ext+(i-1)*dx_ext
  END DO

  DO j=1,ny_ext
    y_ext(j)=y0_ext+(j-1)*dy_ext
  END DO

  CALL xytoll(nx_ext,ny_ext,x_ext,y_ext,lat_ext,lon_ext)

!---------------------------------------------------------------------
!  Define soil depths
!---------------------------------------------------------------------

  DO k=1,nzsoil_ext
    DO j=1,ny_ext
      DO i=1,nx_ext
        zpsoil_ext(i,j,k) = soildepth_ext(k)
      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Retrieve 2-D variables
!
!-----------------------------------------------------------------------
!

  DO j=1,ny_ext
    DO i=1,nx_ext

      IF ( var_nr2d(1,1) == 0 ) THEN   ! first level type, first variable
        psfc_ext   (i,j) = -999.0
      ELSE
        !psfc_ext   (i,j) = var_grb2d(i,j,1,1) * 100.0
        psfc_ext   (i,j) = var_grb2d(i,j,1,1) ! already with Pa unit (F.KONG)
!        IF (psfc_ext(i,j) < 30000.0) THEN
!          WRITE(6,'(1x,a,2(I3,a),F12.5,/,10x,a)')                       &
!            'WARNING: unrealistic surface pressure at (',i,',',j,       &
!            ') = ',psfc_ext(i,j),'It was reset to 100000.00.'
!          psfc_ext(i,j) = 100000.00
!        END IF
      END IF
      IF ( terrainflag == 0 ) THEN
        trn_ext    (i,j) = -999.0
        WRITE(6,'(/a/)') 'WARNING: no terrain found in grib file .... '
      ELSE
        trn_ext    (i,j) = var_grb2d(i,j,2,1)
      END IF

      IF ( var_nr3d(1,2) == 0 ) THEN  ! do not contain soil temperature
        DO k=1,nzsoil_ext
           tsoil_ext (i,j,k) = -999.0
           qsoil_ext (i,j,k) = -999.0
        END DO

      ELSE                            ! yes, it contain for NARR 221

        IF (soilmodel_option == 1) THEN  ! Old ARPS Force-Restore Soil Model
          tsoil_ext(i,j,1) = var_grb2d(i,j,3,1)      ! surface temperature

          IF ( nint(var_grb2d(i,j,5,1)) == 1 ) THEN  ! soil temp over land
            tsoil_ext(i,j,2) = 0.1 * var_grb3d(i,j,1,1,2) & !   0-10cm
                             + 0.3 * var_grb3d(i,j,2,1,2) & !  10-40cm
                             + 0.6 * var_grb3d(i,j,3,1,2)   ! 40-100cm
                                           ! should contain 100-200cm here?
            soiltyp_ext(i,j) = 0
          ELSE                                       ! sfc temp over sea
            tsoil_ext(i,j,2) = var_grb2d(i,j,3,1)
            soiltyp_ext(i,j) = 13            ! Set soil type to water
          END IF

          qsoil_ext(i,j,1) = var_grb3d(i,j,1,2,2)
          qsoil_ext(i,j,2) = 0.1 * var_grb3d(i,j,1,2,2)  & !   0-10cm
                           + 0.3 * var_grb3d(i,j,2,2,2)  & !  10-40cm
                           + 0.6 * var_grb3d(i,j,3,2,2)    ! 40-100cm

        ELSE                             ! OU Soil model

          IF ( nint(var_grb2d(i,j,5,1)) == 1 ) THEN    ! Land
            tsoil_ext (i,j,1) = var_grb2d(i,j,3,1)   ! Ground temperature
            qsoil_ext (i,j,1) = var_grb3d(i,j,1,2,2) ! Assumed to be same
                                                     ! as first below ground
                                                     ! level.
            DO k=2,nzsoil_ext
              ! "TSOIL" in GRIB is below ground, treated as separate
              ! variable from ground temperature.
              tsoil_ext (i,j,k) = var_grb3d(i,j,k-1,1,2) ! Not a bug;
              qsoil_ext (i,j,k) = var_grb3d(i,j,k-1,2,2)
            END DO

          ELSE                                         ! Water
            DO k=1,nzsoil_ext
              tsoil_ext (i,j,k) = var_grb2d(i,j,3,1) ! Water temperature
              qsoil_ext (i,j,k) = 1.                 ! 100% water
            END DO
          END IF ! Land or water?

        END IF                           ! soilmodel_option

      END IF

      IF ( var_nr2d(4,1) == 0 ) THEN
        wetcanp_ext(i,j) = -999.0
      ELSE                               !  Plant canopy surface water (kg/m**2)
        wetcanp_ext(i,j) = var_grb2d(i,j,4,1)*1.e-3     ! in meter
      END IF

      IF ( var_nr2d(6,1) == 0 ) THEN
        snowdpth_ext(i,j) = -999
      ELSE

        ! Convert water equiv. of accum. snow depth (kg/m**2) to meters
        ! (where 1 meter liquid water is set equivqlent to 10 meters snow).
        !      0.01 = 10. (m snow/m liquid) / (1000 kg/m**3)
        snowdpth_ext(i,j) = 0.01 * var_grb2d(i,j,6,1)  ! in meters
      END IF

      IF ( var_nr2d(7,1) == 0 ) THEN
        rain_ext(i,j)= -999.0
      ELSE
        rain_ext(i,j)= var_grb2d(i,j,7,1)
      END IF

      IF( var_nr2d(1,2) == 0 ) THEN
        t_2m_ext(i,j)= -999.0
      ELSE
        t_2m_ext(i,j)= var_grb2d(i,j,1,2)
!        IF (t_2m_ext(i,j) < 100.0) THEN
!          WRITE(6,'(1x,a,2(I3,a),F12.5,/,10x,a)')                       &
!            'WARNING: unrealistic 2m temperature at (',i,',',j,         &
!            ') = ',t_2m_ext(i,j),'It was reset to 273.15.'
!          t_2m_ext(i,j) = 273.15
!        END IF
      END IF
      IF( var_nr2d(2,2) == 0 ) THEN
        qv_2m_ext(i,j)= -999.0
      ELSE
        qv_2m_ext(i,j)= var_grb2d(i,j,2,2)
      END IF

      IF( var_nr2d(3,2) == 0 ) THEN
        u_10m_ext(i,j)= -999.0
      ELSE
        u_10m_ext(i,j)= var_grb2d(i,j,3,2)
      END IF

      IF( var_nr2d(4,2) == 0 ) THEN
        v_10m_ext(i,j)= -999.0
      ELSE
        v_10m_ext(i,j)= var_grb2d(i,j,4,2)
      END IF

    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Retrieve 3-D variables
!
!-----------------------------------------------------------------------
!
  nz1 = MIN(var_nr3d(1,1),nz_ext)

  IF ( var_lev3d(1,1,1) > var_lev3d(nz1,1,1) ) THEN  ! 1st level at sfc
    chklev = 1
    lvscan = 0
  ELSE
    chklev = -1
    lvscan = nz1+1
  END IF

  DO k=1,nz1
    kk = chklev * k + lvscan
    DO j=1,ny_ext
      DO i=1,nx_ext

        p_ext  (i,j,kk) = 100.0 * REAL(var_lev3d(k,1,1)) ! Pressure
        t_ext  (i,j,kk) = var_grb3d(i,j,k,1,1)    ! Temperature (K)
!        IF (t_ext(i,j,kk) < 100.0) THEN
!          WRITE(6,'(1x,a,3(I3,a),F12.5,/,10x,a)')                       &
!            'WARNING: unrealistic temperature at (',i,',',j,',',kk,     &
!            ') = ',t_ext(i,j,kk),'It was reset to 273.15.'
!          t_ext(i,j,kk) = 273.15 + 20.0
!        END IF
        qv_ext (i,j,kk) = var_grb3d(i,j,k,2,1)    ! Spec. humidity
                                                  ! (kg/kg)
        u_ext  (i,j,kk) = var_grb3d(i,j,k,3,1)    ! u wind (m/s)
        v_ext  (i,j,kk) = var_grb3d(i,j,k,4,1)    ! v wind (m/s)
        hgt_ext(i,j,kk) = var_grb3d(i,j,k,5,1)    ! height (m)
                         ! vertical velocity (pressure) - skipped.
        qc_ext (i,j,kk) = var_grb3d(i,j,k,7,1)    ! height (m)
        qi_ext (i,j,kk) = var_grb3d(i,j,k,8,1)    ! height (m)
!        qr_ext (i,j,kk) = -999.
!        qs_ext (i,j,kk) = -999.
!        qh_ext (i,j,kk) = -999.
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Rotate winds to be relative to true north.
!  The RUC data are sent as grid-relative.
!
!-----------------------------------------------------------------------
!
!NCEP NARR output is already earth relative - TINA

  istatus = 1

!
!-----------------------------------------------------------------------
!
!  Reset map projection to previous values
!
!-----------------------------------------------------------------------
!

  999  CONTINUE
  CALL setmapr(iproj,scale,latnot,trlon)
  CALL setorig(1,x0,y0)

  DEALLOCATE(var_grb2d,var_grb3d,rcdata,var_lev3d)

  RETURN
END SUBROUTINE getnarr221

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETECMF128                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getecmf128(nx_ext,ny_ext,nz_ext,nzsoil_ext,                  &
           nstyps_ext,dir_extd,extdname,extdopt,extdfmt,                &
           extdinit,extdfcst,julfname,                                  &
           iproj_ext,scale_ext,                                         &
           trlon_ext,latnot_ext,x0_ext,y0_ext,                          &
           lat_ext,lon_ext,zpsoil_ext,                                  &
           p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                      &
          !qc_ext,qr_ext,qi_ext,qs_ext,qh_ext,                          &
           tsoil_ext,qsoil_ext,wetcanp_ext,                             &
           snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext,                   &
           t_2m_ext, qv_2m_ext, u_10m_ext, v_10m_ext, rain_ext,         &
           istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads a ECMWF GRIB file for processing by ext2arps, a program that
!  converts external files to ARPS variables and format.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  04/06/2011
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    dir_extd      Directory name for external file
!    extdname      Prefix string of external file name
!    extdopt       Option of external data sources
!    extdfmt       Option of external data format

!    extdinit      Initialized time in mm-dd-yyyy:hh:mm:ss format
!    extdfcst      Forecast hour in HHH:MM:SS format
!    julfname      File name in yyjjjhhmm format
!
!  OUTPUT:
!
!    iproj_ext     Map projection number of external data
!    scale_ext     Scale factor of external data
!    trlon_ext     True longitude of external data (degrees E)
!    latnot_ext(2) True latitude(s) of external data (degrees N)
!    x0_ext        x coordinate of origin of external data
!    y0_ext        y coordinate of origin of external data
!    lat_ext       latitude of external data points (degrees N)
!    lon_ext       longitude of external data points (degrees E)
!    p_ext         pressure (Pascal)
!    hgt_ext       height (m)
!    t_ext         temperature (K)
!    qv_ext        specific humidity (kg/kg)
!    u_ext         u wind component (m/s)
!    v_ext         v wind component (m/s)
!    qc_ext        Cloud water mixing ratio (kg/kg)
!    qr_ext        Rain water mixing ratio (kg/kg)
!    qi_ext        Ice mixing ratio (kg/kg)
!    qs_ext        Snow mixing ratio (kg/kg)
!    qh_ext        Hail mixing ratio (kg/kg)
!
!    tsoil_ext     Soil temperature (K)
!    qsoil_ext     Soil moisture (m**3/m**3)
!    wetcanp_ext   Water content on canopy
!
!    trn_ext       External terrain (m)
!    psfc_ext      Surface pressure (Pa)
!
! F.KONG add four near surface variables
!    T_2m_ext      Temperature at 2m AGL
!    qv_2m_ext     Specific Humidity at 2m AGL
!    U_10m_ext     U at 10m AGL
!    V_10m_ext     V at 10m AGL
!    rain_ext
! end F.KONG mod
!
!    istatus       status indicator
!
!  WORK ARRAYS:
!
!    var_grb3d     Arrays to store the GRIB 3-D variables:
!                  var_grb3d(nxgrb,nygrb,nzgrb,1,1) - Temperature (K)
!                  var_grb3d(nxgrb,nygrb,nzgrb,2,1) - Specific humidity
!                                                     (kg/kg)
!                  var_grb3d(nxgrb,nygrb,nzgrb,3,1) - u wind (m/s)
!                  var_grb3d(nxgrb,nygrb,nzgrb,4,1) - v wind (m/s)
!                  var_grb3d(nxgrb,nygrb,nzgrb,5,1) - Geopotential
!                                                     height (gpm)
!                  var_grb3d(nxgrb,nygrb,nzgrb,6,1) - Pressure vertical
!                                                     velocity (Pa/s)
!                                                     (if applied)
!                  var_grb3d(nxgrb,nygrb,nzgrb,1,2) - soil temp. (K)
!                  var_grb3d(nxgrb,nygrb,nzgrb,2,2) - vol. soil moist.
!                                                     (m**3/m**3)
!
!    var_grb2d     Arrays to store the GRIB 2-D variables:
!                  var_grb2d(nxgrb,nygrb,1) - Surface pressure (Pa)
!                  var_grb2d(nxgrb,nygrb,2) - Geopotential height (gpm)
!                  var_grb2d(nxgrb,nygrb,3) - Surface temperature (K)
!                  var_grb2d(nxgrb,nygrb,4) - Plant canopy surface
!                                             water (kg/m**2)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'gribcst.inc'

  CHARACTER (LEN=*) :: dir_extd
  CHARACTER (LEN=*) :: extdname

  INTEGER :: extdopt
  INTEGER :: extdfmt

  CHARACTER (LEN=19) :: extdinit
  CHARACTER (LEN=9) :: extdfcst
  CHARACTER (LEN=9) :: julfname
!
!-----------------------------------------------------------------------
!
!  External grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj_ext
  REAL :: scale_ext,trlon_ext
  REAL :: latnot_ext(2)
  REAL :: x0_ext,y0_ext
  REAL :: dx_ext,dy_ext
!
!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx_ext,ny_ext,nz_ext
  INTEGER :: nzsoil_ext, nstyps_ext

  REAL :: lat_ext(nx_ext,ny_ext)
  REAL :: lon_ext(nx_ext,ny_ext)
  REAL :: p_ext  (nx_ext,ny_ext,nz_ext)   ! Pressure (Pascals)
  REAL :: hgt_ext(nx_ext,ny_ext,nz_ext)   ! Height (m)
  REAL :: t_ext  (nx_ext,ny_ext,nz_ext)   ! Temperature (K)
  REAL :: qv_ext (nx_ext,ny_ext,nz_ext)   ! Specific humidity (kg/kg)
  REAL :: u_ext  (nx_ext,ny_ext,nz_ext)   ! Eastward wind component
  REAL :: v_ext  (nx_ext,ny_ext,nz_ext)   ! Northward wind component
!  REAL :: qc_ext (nx_ext,ny_ext,nz_ext)   ! Cloud H2O mixing ratio (kg/kg)
!  REAL :: qr_ext (nx_ext,ny_ext,nz_ext)   ! Rain  H2O mixing ratio (kg/kg)
!  REAL :: qi_ext (nx_ext,ny_ext,nz_ext)   ! Ice   mixing ratio (kg/kg)
!  REAL :: qs_ext (nx_ext,ny_ext,nz_ext)   ! Snow  mixing ratio (kg/kg)
!  REAL :: qh_ext (nx_ext,ny_ext,nz_ext)   ! Hail  mixing ratio (kg/kg)
  REAL :: zpsoil_ext(nx_ext,ny_ext,nzsoil_ext) !Soil depths (m)

  REAL :: tsoil_ext (nx_ext,ny_ext,nzsoil_ext) ! Soil temperature (K)
  REAL :: qsoil_ext (nx_ext,ny_ext,nzsoil_ext) ! Soil moisture (m**3/m**3)
  REAL :: wetcanp_ext(nx_ext,ny_ext)      ! Canopy water amount
  REAL :: snowdpth_ext(nx_ext,ny_ext)     ! Snow depth (m)

  REAL :: trn_ext    (nx_ext,ny_ext)      ! External terrain (m)
  REAL :: psfc_ext   (nx_ext,ny_ext)      ! Surface pressure (Pa)

  INTEGER :: soiltyp_ext (nx_ext,ny_ext)     ! Soil type

  REAL :: t_2m_ext (nx_ext,ny_ext)
  REAL :: qv_2m_ext(nx_ext,ny_ext)
  REAL :: u_10m_ext(nx_ext,ny_ext)
  REAL :: v_10m_ext(nx_ext,ny_ext)
  REAL :: rain_ext (nx_ext,ny_ext)
!
!-----------------------------------------------------------------------
!
!  Other  external variable arrays
!
!-----------------------------------------------------------------------
!
  REAL :: x_ext(nx_ext)
  REAL :: y_ext(ny_ext)

  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Work arrays for storing grib data
!
!-----------------------------------------------------------------------
!
  REAL,    ALLOCATABLE :: var_grb2d(:,:,:,:)   ! GRIB variables
  REAL,    ALLOCATABLE :: var_grb3d(:,:,:,:,:) ! GRIB 3-D variables
  INTEGER, ALLOCATABLE :: var_lev3d(:,:,:)  ! Levels (hybrid) for
                                            ! each 3-D variable
  REAL, ALLOCATABLE :: rcdata(:)            ! temporary data array

  REAL, ALLOCATABLE :: qvs_ext(:,:,:)
!
!-----------------------------------------------------------------------
!
!  Original grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj
  REAL    :: scale,trlon,x0,y0
  REAL    :: latnot(2)
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256) :: gribfile
  CHARACTER (LEN=14)  :: gribtime
  INTEGER :: i,j,k,kk
  INTEGER :: grbflen, grbtlen

  INTEGER :: m,n,nz1,max_nr2d,max_nr3d

  REAL :: govrd
  INTEGER, PARAMETER :: nzsoilin_ext = 4
  REAL,    PARAMETER :: soildepth_ext(nzsoilin_ext+1) = (/0.0, 0.07, 0.28, 1.0, 2.55/)
           ! contains 4 soil layers 0-7cm, 7-28cm, 28-100cm, 100-255cm.

  INTEGER :: chklev, lvscan

  LOGICAL :: fexist
  INTEGER :: igrbfmt
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'grid.inc'
!
!-----------------------------------------------------------------------
!
!  GRIB grid information
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=42) :: gridesc ! Grid description

  INTEGER :: iproj_grb    ! Map projection indicator
  INTEGER :: gthin        ! Indicator of whether the grid is "thinned"

  INTEGER :: ni_grb       ! Number of points along x-axis
  INTEGER :: nj_grb       ! Number of points along y-axis
  INTEGER :: np_grb       ! Total number of horizontal grid points

  INTEGER :: nk_grb       ! Number of vertical parameters
  REAL    :: zk_grb(nz_ext)  ! Vertical coordinate parameters

  INTEGER :: npeq         ! Number of lat circles from pole to equator
  INTEGER :: nit(nz_ext)  ! Number of x-points for thinned grid

  REAL    :: pi_grb       ! x-coordinate of pole point
  REAL    :: pj_grb       ! y-coordinate of pole point
  INTEGER :: ipole        ! Projection center flag

  REAL :: di_grb          ! x-direction increment or grid length
  REAL :: dj_grb          ! y-direction increment or grid length

  REAL :: latsw           ! Latitude  of South West corner point
  REAL :: lonsw           ! Longitude of South West corner point
  REAL :: latne           ! Latitude  of North East corner point
  REAL :: lonne           ! Longitude of North East corner point

  REAL :: lattru1         ! Latitude (1st) at which projection is true
  REAL :: lattru2         ! Latitude (2nd) at which projection is true
  REAL :: lontrue         ! Longitude      at which projection is true

  REAL :: latrot          ! Latitude  of southern pole of rotation
  REAL :: lonrot          ! Longitude of southern pole of rotation
  REAL :: angrot          ! Angle of rotation

  REAL :: latstr          ! Latitude  of the pole of stretching
  REAL :: lonstr          ! Longitude of the pole of stretching
  REAL :: facstr          ! Stretching factor

  INTEGER :: scanmode     ! Scanning indicator
  INTEGER :: iscan        ! x-direction   scanning indicator
  INTEGER :: jscan        ! y-direction   scanning indicator
  INTEGER :: kscan        ! FORTRAN index scanning indicator

  INTEGER :: ires         ! Resolution direction increments indicator
  INTEGER :: iearth       ! Earth shape indicator: spherical or oblate?
  INTEGER :: icomp        ! (u,v) components decomposition indicator

  INTEGER :: jpenta       ! J-Pentagonal resolution parameter
  INTEGER :: kpenta       ! K-Pentagonal resolution parameter
  INTEGER :: mpenta       ! M-Pentagonal resolution parameter
  INTEGER :: ispect       ! Spectral representation type
  INTEGER :: icoeff       ! Spectral coefficient storage mode

  REAL :: xp_grb          ! X coordinate of sub-satellite point
  REAL :: yp_grb          ! Y coordinate of sub-satellite point
  REAL :: xo_grb          ! X coordinate of image sector origin
  REAL :: yo_grb          ! Y coordinate of image sector origin
  REAL :: zo_grb          ! Camera altitude from center of Earth
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  CALL getmapr(iproj,scale,latnot,trlon,x0,y0)

  CALL getgrbfname(extdopt,extdfmt,dir_extd,extdname,                   &
                   extdinit,extdfcst,'    ',                            &
                   gribfile,grbflen,gribtime,grbtlen,istatus)
  IF (istatus /=0 ) RETURN

  CALL chkgrb(gribfile(1:grbflen),grbflen,igrbfmt,istatus)
  IF (istatus /=0 ) RETURN

  IF (igrbfmt == 2) THEN

    WRITE(6,'(1x,a)') 'No GRIB 2 support for ECMWF at present'

    istatus = 0
    RETURN

  ELSE

    ALLOCATE(var_grb2d(nx_ext,ny_ext,n2dvs,n2dlvt))
    ALLOCATE(var_grb3d(nx_ext,ny_ext,nz_ext,n3dvs,n3dlvt))
    ALLOCATE(rcdata(nx_ext*ny_ext))
    ALLOCATE(var_lev3d(nz_ext,n3dvs,n3dlvt))
!
!-----------------------------------------------------------------------
!
!  RDNMCGRB reads ECMWF GRIB data using NMC grib facility
!
!-----------------------------------------------------------------------
!
    govrd = g/rd

    gridtyp   = ecmf128grid
    mproj_grb = ecmf128proj

    n2dvars  = ecmf128nvs2d
    n2dlvtps = ecmf128nlvt2d

    DO k=1,n2dlvtps
      DO n=1,n2dvars
        var_id2d(n,k) = ecmf128var_id2d(n,k)
      END DO
      levtyp2d(k) = ecmf128levs2d(k)
    END DO

    n3dvars  = ecmf128nvs3d
    n3dlvtps = ecmf128nlvt3d

    DO m=1,n3dlvtps
      DO n=1,n3dvars
        var_id3d(n,m) = ecmf128var_id3d(n,m)
      END DO
      levtyp3d(m) = ecmf128levs3d(m)
    END DO

    CALL rdnmcgrb(nx_ext,ny_ext,nz_ext,gribfile,grbflen, gribtime,      &
                  gridesc, iproj_grb, gthin,                            &
                  ni_grb,nj_grb,np_grb, nk_grb,zk_grb, npeq,nit,        &
                  pi_grb,pj_grb,ipole, di_grb,dj_grb,                   &
                  latsw,lonsw, latne,lonne,                             &
                  latrot,lonrot,angrot,                                 &
                  latstr,lonstr,facstr,                                 &
                  lattru1,lattru2,lontrue,                              &
                  scanmode, iscan,jscan,kscan,                          &
                  ires,iearth,icomp,                                    &
                  jpenta,kpenta,mpenta,ispect,icoeff,                   &
                  xp_grb,yp_grb, xo_grb,yo_grb,zo_grb,                  &
                  rcdata,var_grb2d,var_grb3d,var_lev3d,lvldbg,istatus)

    IF (istatus /= 0)  THEN
      ISTATUS = -888
      GOTO 999
    END IF

    max_nr2d = 0
    DO n=1,n2dvars
      DO m=1,n2dlvtps
        max_nr2d = MAX( max_nr2d, var_nr2d(n,m) )
      END DO
    END DO

    max_nr3d = 0
    DO n=1,n3dvars
      max_nr3d = MAX( max_nr3d, var_nr3d(n,1) )
    END DO

    IF ( max_nr3d == 0 ) THEN
      WRITE (6,'(a)')                                                   &
          'No 3-D variable was found in the GRIB file',                 &
          'Dataset problem in GETNMCECMF128.',                          &
          ' '
      !STOP
      ISTATUS = -888
      GOTO 999
    END IF

    IF ( max_nr2d == 0 ) THEN
      WRITE (6,'(a)')                                                   &
          'No 2-D variables was found in the GRIB file'
    END IF

  !  WRITE (6,'(/a7,2x,6(i7))') 'Lev\VID',(var_id3d(n,1),n=1,n3dvars)
  !  DO k=1,max_nr3d
  !    write (6,'(/i5,4x,6(i7))') k,(var_lev3d(k,n,1),n=1,n3dvars)
  !   END DO

    DO k=1,max_nr3d
      DO n=2,n3dvars
        IF ( var_lev3d(k,1,1) /= var_lev3d(k,n,1) ) THEN
          WRITE (6,'(a,/,a,I3,a,I7,/,a,I7,a,I3,/,a,/)')                   &
          'Variables were not at the same level.',                        &
          'The K (k= ',k,')th level for variable 1 is ',var_lev3d(k,1,1), &
          'but it is ',var_lev3d(k,n,1),' for variable n = ',n,           &
          'Dataset problem in GETNMCECMF128.'
          !STOP
          ISTATUS = -888
          GOTO 999
        END IF
      END DO
    END DO

    IF ( iproj_grb == 5 .AND. ipole == 0 ) THEN      ! Center in NP
      iproj_ext = 1
    ELSE IF ( iproj_grb == 5 .AND. ipole == 1 ) THEN  ! Center in SP
      iproj_ext = -1
    ELSE IF ( iproj_grb == 3 .AND. ipole == 0 ) THEN  ! Center in NP
      iproj_ext = 2
    ELSE IF ( iproj_grb == 3 .AND. ipole == 1 ) THEN  ! Center in SP
      iproj_ext = -2
    ELSE IF ( iproj_grb == 1 ) THEN
      iproj_ext = 3
    ELSE IF ( iproj_grb == 0 ) THEN
      iproj_ext = 4
    ELSE
      WRITE (6,'(a)')                                                     &
          'Unknown map projection. Set to non-projection.'
      iproj_ext = 0
    END IF

    scale_ext = 1.0
    latnot_ext(1) = lattru1
    latnot_ext(2) = lattru2
    trlon_ext = lontrue

    dx_ext = di_grb
    dy_ext = dj_grb
!
!-----------------------------------------------------------------------
!
!  Retrieve 2-D variables
!
!-----------------------------------------------------------------------
!

    DO j=1,ny_ext
      DO i=1,nx_ext

        IF ( var_nr2d(2,1) == 0 ) THEN
          trn_ext    (i,j) = -999.0
        ELSE
          trn_ext    (i,j) = var_grb2d(i,j,2,1) / g  ! m**2 s**-2
        END IF

        IF ( var_nr2d(3,1) == 0 ) THEN
          psfc_ext   (i,j) = -999.0
        ELSE
          psfc_ext   (i,j) = var_grb2d(i,j,3,1) ! already with Pa unit
        END IF

        !IF ( var_nr2d(4,1) == 0 ) THEN
          wetcanp_ext(i,j) = -999.0
        !ELSE
        !  wetcanp_ext(i,j) = var_grb2d(i,j,4,1)*1.e-3     ! in meter
        !END IF

        !IF ( var_nr2d(6,1) == 0 ) THEN
          snowdpth_ext(i,j) = -999.0
        !ELSE
        !
        !END IF

        IF( var_nr2d(4,1) == 0 ) THEN
          u_10m_ext(i,j)= -999.0
        ELSE
          u_10m_ext(i,j)= var_grb2d(i,j,4,1)
        END IF

        IF( var_nr2d(5,1) == 0 ) THEN
          v_10m_ext(i,j)= -999.0
        ELSE
          v_10m_ext(i,j)= var_grb2d(i,j,5,1)
        END IF

        IF( var_nr2d(6,1) == 0 ) THEN
          t_2m_ext(i,j)= -999.0
        ELSE
          t_2m_ext(i,j)= var_grb2d(i,j,6,1)
        END IF
        !IF( var_nr2d(2,2) == 0 ) THEN
          qv_2m_ext(i,j)= -999.0
        !ELSE
        !  qv_2m_ext(i,j)= var_grb2d(i,j,2,2)
        !END IF

        !IF ( var_nr2d(8,1) == 0 ) THEN
          rain_ext(i,j)= -999.0
        !ELSE
        !  rain_ext(i,j)= var_grb2d(i,j,8,1)
        !END IF

      END DO
    END DO

!
!-----------------------------------------------------------------------
!
!  Retrieve 3-D variables
!
!-----------------------------------------------------------------------
!
    nz1 = MIN(var_nr3d(1,1),nz_ext)

    IF ( var_lev3d(1,1,1) > var_lev3d(nz1,1,1) ) THEN  ! 1st level at sfc
      chklev = 1
      lvscan = 0
    ELSE
      chklev = -1
      lvscan = nz1+1
    END IF

    DO k=1,nz1
      kk = chklev * k + lvscan
      DO j=1,ny_ext
        DO i=1,nx_ext
          p_ext  (i,j,kk) = 100.0 * REAL(var_lev3d(k,1,1)) ! Pressure
          t_ext  (i,j,kk) = var_grb3d(i,j,k,1,1)    ! Temperature (K)
          u_ext  (i,j,kk) = var_grb3d(i,j,k,2,1)    ! u wind (m/s)
          v_ext  (i,j,kk) = var_grb3d(i,j,k,3,1)    ! v wind (m/s)
          hgt_ext(i,j,kk) = var_grb3d(i,j,k,4,1)    ! height (m)
          qv_ext (i,j,kk) = var_grb3d(i,j,k,5,1)    ! Rel. humidity (%)

!          qc_ext (i,j,kk) = -999.
!          qr_ext (i,j,kk) = -999.
!          qi_ext (i,j,kk) = -999.
!          qs_ext (i,j,kk) = -999.
!          qh_ext (i,j,kk) = -999.
        END DO
      END DO
    END DO

!-----------------------------------------------------------------------
!
! Convert relative humidity to specific humidity
!
!-----------------------------------------------------------------------

    ALLOCATE(qvs_ext(nx_ext,ny_ext,nz_ext), STAT = istatus)
    CALL check_alloc_status(istatus, "getecmf128:qvs_ext")

    CALL getqvs(nx_ext,ny_ext,nz_ext, 1,nx_ext,1,ny_ext,1,nz_ext,       &
                p_ext, t_ext, qvs_ext )

    DO k=1,nz_ext
      DO j=1,ny_ext
        DO i=1,nx_ext
          qv_ext(i,j,k) = 0.01*qv_ext(i,j,k)*qvs_ext(i,j,k)
        END DO
      END DO
    END DO
    DEALLOCATE(qvs_ext)
!
!-----------------------------------------------------------------------
!
!  Retrieve 3-D soil variables
!
!-----------------------------------------------------------------------
!
    IF ( var_nr2d(1,2) == 0 ) THEN
      DO k=1,nzsoil_ext
        DO j = 1,ny_ext
          DO i = 1,nx_ext
            tsoil_ext (i,j,k) = -999.0
            qsoil_ext (i,j,k) = -999.0
          END DO
        END DO
      END DO

    ELSE

      IF (soilmodel_option == 1) THEN ! Old ARPS Force-Restore Soil Model

        DO j = 1,ny_ext
          DO i = 1,nx_ext

            IF ( var_grb2d(i,j,1,1) > 170.0) THEN  ! SST abvailable, over ocean
            !IF ( nint(var_grb2d(i,j,7,1)) == 0 ) THEN  ! soil temp over ocean
              tsoil_ext(i,j,1) = var_grb2d(i,j,1,1)
              tsoil_ext(i,j,2) = var_grb2d(i,j,1,1)
              qsoil_ext(i,j,1) = 1.0
              qsoil_ext(i,j,2) = 1.0
              soiltyp_ext(i,j) = 13 ! Set soil type to water
            ELSE                                  ! sfc temp over land
              tsoil_ext(i,j,1) = var_grb2d(i,j,8,1)
              tsoil_ext(i,j,2) = 0.1 * var_grb2d(i,j,5,2) & !   0-7 cm
                               + 0.3 * var_grb2d(i,j,6,2) & !   7-28cm
                               + 0.6 * var_grb2d(i,j,7,2)   ! 28-100cm
              qsoil_ext(i,j,1) = var_grb2d(i,j,1,2)
              qsoil_ext(i,j,2) = 0.1 * var_grb2d(i,j,1,2)  & !   0-7 cm
                               + 0.3 * var_grb2d(i,j,2,2)  & !   7-28cm
                               + 0.6 * var_grb2d(i,j,3,2)    ! 28-100cm

              soiltyp_ext(i,j) = 0
            END IF

          END DO
        END DO

      ELSE ! OU Soil model

        DO j = 1,ny_ext
          DO i = 1,nx_ext

            IF ( var_grb2d(i,j,1,1) < 170.0) THEN  ! SST not abvailable, over land
            !IF ( nint(var_grb2d(i,j,7,1)) == 1 ) THEN  ! Land
              tsoil_ext (i,j,1) = var_grb2d(i,j,8,1)   ! Ground temperature
              qsoil_ext (i,j,1) = var_grb2d(i,j,1,2)   ! Assumed to be same
                                                       ! as first below ground
                                                       ! level.
              DO k=2,nzsoil_ext
                ! "TSOIL" in GRIB is below ground, treated as separate
                ! variable from ground temperature.
                tsoil_ext (i,j,k) = var_grb2d(i,j,k+3,2)
                qsoil_ext (i,j,k) = var_grb2d(i,j,k-1,2)
              END DO
              soiltyp_ext(i,j) = 0

            ELSE  ! Water
              DO k=1,nzsoil_ext
                tsoil_ext (i,j,k) = var_grb2d(i,j,1,1) ! Water temperature
                qsoil_ext (i,j,k) = 1.                 ! 100% water
              END DO
              soiltyp_ext(i,j) = 13 ! Set soil type to water
            END IF ! Land or water?

          END DO
        END DO

      END IF  ! soilmodel_option

    END IF

    DEALLOCATE(var_grb2d,var_grb3d,rcdata,var_lev3d)
  END IF

  IF (lonsw > 180) lonsw = lonsw-360.
  CALL setmapr(iproj_ext,scale_ext,latnot_ext,trlon_ext)
  CALL lltoxy(1,1,latsw,lonsw,x0_ext,y0_ext)

  DO i=1,nx_ext
    x_ext(i)=x0_ext+(i-1)*dx_ext
  END DO

  DO j=1,ny_ext
    y_ext(j)=y0_ext+(j-1)*dy_ext
  END DO

  CALL xytoll(nx_ext,ny_ext,x_ext,y_ext,lat_ext,lon_ext)

!---------------------------------------------------------------------
!
!  Define soil depths
!
!---------------------------------------------------------------------

  IF(soilmodel_option == 1 .AND. nzsoil_ext /= 2) THEN
    WRITE(6,'(2a,I3,a/2a,I2,a/,a,a)')                                   &
               'ECMWF only provides 4 soil layers, However,',           &
               ' You are trying to extract ',nzsoil_ext, ' layers.',    &
               ' for ARPS two-layer force-restore model ',              &
               '(soilmodel_option = ',soilmodel_option,')',             &
               ' Please check the code ext2arps.f90 for NAM212 grid.',  &
               ' Terminating ...'
    CALL arpsstop('Wrong number soil layers.',1)
  ELSE IF( soilmodel_option == 2 .AND. nzsoil_ext /= nzsoilin_ext+1) THEN
    WRITE(6,'(2a,I3,a/2a,I2,a/,a,a)')                                   &
               'ECMWF only provides 4 soil layers, However,',           &
               ' You are trying to extract ',nzsoil_ext, ' layers.',    &
               ' for ARPS multi-layer OUSoil model ',                   &
               '(soilmodel_option = ',soilmodel_option,')',             &
               ' Please check the code ext2arps.f90 for NAM212 grid.',  &
               ' Terminating ...'
    CALL arpsstop('Wrong number of soil layers.',1)
  END IF

  IF(soilmodel_option == 1) THEN
    DO j=1,ny_ext
      DO i=1,nx_ext
        zpsoil_ext(i,j,1) = 0.05  ! any values should work because
        zpsoil_ext(i,j,2) = 1.0   ! the ARPS system will ignore these values
      END DO
    END DO
  ELSE
    DO j=1,ny_ext
      DO i=1,nx_ext
        DO k=2,nzsoil_ext
          zpsoil_ext(i,j,k) = - (soildepth_ext(k)-soildepth_ext(k-1))/2 - soildepth_ext(k-1)
                                  ! The middle point for each layer
        END DO
        zpsoil_ext(i,j,1) = 0.0
      END DO
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Rotate winds to be relative to true north.
!  not needed since it is Latitude/Longitude grid
!-----------------------------------------------------------------------
!
  istatus = 1
!
!-----------------------------------------------------------------------
!
!  Reset map projection to previous values
!
!-----------------------------------------------------------------------
!

  999  CONTINUE
  CALL setmapr(iproj,scale,latnot,trlon)
  CALL setorig(1,x0,y0)

  RETURN
END SUBROUTINE getecmf128

!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE GETNCEPGFS_GRB2               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

SUBROUTINE getncepgfs_grb2(nx_ext,ny_ext,nz_ext,nzsoil_ext,nzsoilin_ext,&
                        gribfile, grbflen, gribtime, soildepth_ext,     &
                        lon_0_360,iproj_ext,scale_ext,                  &
                        trlon_ext,latnot_ext,lat_ext,lon_ext,           &
                        p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,         &
                        qc_ext,                                         &
                        tsoil_ext,qsoil_ext,wetcanp_ext,                &
                        snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext,      &
                        t_2m_ext, qv_2m_ext, u_10m_ext, v_10m_ext,      &
                        rain_ext,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads a NCEP GFS 0.5 degree global data in GRIB2 for processing by ext2arps.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  08/02/2006
!
!  MODIFICATION HISTORY:
!
!  08/30/2010 Y. Wang
!  Merged to add GRIB format support.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    gribfile      external file name
!    gribtime
!    grbflen
!    soildepth_ext
!
!
!  OUTPUT:
!
!    iproj_ext     Map projection number of external data
!    scale_ext     Scale factor of external data
!    trlon_ext     True longitude of external data (degrees E)
!    latnot_ext(2) True latitude(s) of external data (degrees N)
!    lat_ext       latitude of external data points (degrees N)
!    lon_ext       longitude of external data points (degrees E)
!    p_ext         pressure (Pascal)
!    hgt_ext       height (m)
!    t_ext         temperature (K)
!    qv_ext        specific humidity (kg/kg)
!    u_ext         u wind component (m/s)
!    v_ext         v wind component (m/s)
!    qc_ext        Cloud water mixing ratio (kg/kg)
!    qr_ext        Rain water mixing ratio (kg/kg)
!    qi_ext        Ice mixing ratio (kg/kg)
!    qs_ext        Snow mixing ratio (kg/kg)
!    qh_ext        Hail mixing ratio (kg/kg)
!
!    tsoil_ext     Soil temperature
!    qsoil_ext     Soil moisture
!    wetcanp_ext   Water content on canopy
!
!    trn_ext       External terrain (m)
!    psfc_ext      Surface pressure (Pa)
!
!    T_2m_ext      Temperature at 2m AGL
!    qv_2m_ext     Specific Humidity at 2m AGL
!    U_10m_ext     U at 10m AGL
!    V_10m_ext     V at 10m AGL
!
!    istatus       status indicator
!
!  WORK ARRAYS:
!
!    var_grb3d     Arrays to store the GRIB 3-D variables:
!
!    var_grb2d     Arrays to store the GRIB 2-D variables:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx_ext,ny_ext,nz_ext,nzsoil_ext
  INTEGER, INTENT(IN) :: nzsoilin_ext

  CHARACTER (LEN=256) :: gribfile
  CHARACTER (LEN=14)  :: gribtime
  INTEGER, INTENT(IN) :: grbflen

  REAL,    INTENT(IN) :: soildepth_ext(nzsoilin_ext)
  LOGICAL, INTENT(IN) :: lon_0_360
!
!-----------------------------------------------------------------------
!
!  External grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj_ext
  REAL    :: scale_ext,trlon_ext
  REAL    :: latnot_ext(2)
  REAL    :: dx_ext,dy_ext
!
!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------
!
  REAL :: lat_ext(nx_ext,ny_ext)
  REAL :: lon_ext(nx_ext,ny_ext)
  REAL :: p_ext  (nx_ext,ny_ext,nz_ext)   ! Pressure (Pascals)
  REAL :: hgt_ext(nx_ext,ny_ext,nz_ext)   ! Height (m)
  REAL :: t_ext  (nx_ext,ny_ext,nz_ext)   ! Temperature (K)
  REAL :: qv_ext (nx_ext,ny_ext,nz_ext)   ! Specific humidity (kg/kg)
  REAL :: u_ext  (nx_ext,ny_ext,nz_ext)   ! Eastward wind component
  REAL :: v_ext  (nx_ext,ny_ext,nz_ext)   ! Northward wind component
  REAL :: qc_ext (nx_ext,ny_ext,nz_ext)   ! Cloud H2O mixing ratio (kg/kg)
!  REAL :: qr_ext (nx_ext,ny_ext,nz_ext)   ! Rain  H2O mixing ratio (kg/kg)
!  REAL :: qi_ext (nx_ext,ny_ext,nz_ext)   ! Ice   mixing ratio (kg/kg)
!  REAL :: qs_ext (nx_ext,ny_ext,nz_ext)   ! Snow  mixing ratio (kg/kg)
!  REAL :: qh_ext (nx_ext,ny_ext,nz_ext)   ! Hail  mixing ratio (kg/kg)

  REAL :: tsoil_ext  (nx_ext,ny_ext,nzsoil_ext)      ! Temperature at surface (K)
  REAL :: qsoil_ext  (nx_ext,ny_ext,nzsoil_ext)      ! Deep soil temperature (K)
  REAL :: wetcanp_ext(nx_ext,ny_ext)                 ! Canopy water amount
  REAL :: snowdpth_ext(nx_ext,ny_ext)                ! Snow depth (m)

  REAL :: trn_ext    (nx_ext,ny_ext)      ! External terrain (m)
  REAL :: psfc_ext   (nx_ext,ny_ext)      ! Surface pressure (Pa)
  REAL :: t_2m_ext (nx_ext,ny_ext)
  REAL :: qv_2m_ext(nx_ext,ny_ext)
  REAL :: u_10m_ext(nx_ext,ny_ext)
  REAL :: v_10m_ext(nx_ext,ny_ext)
  REAL :: rain_ext (nx_ext,ny_ext)

  INTEGER :: soiltyp_ext(nx_ext,ny_ext)

  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Work arrays for storing grib data
!
!-----------------------------------------------------------------------
!
  REAL,    ALLOCATABLE :: var_grb2d(:,:,:)   ! GRIB 2D variables
  REAL,    ALLOCATABLE :: var_grb3d(:,:,:,:) ! GRIB 3-D variables
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,kk

  INTEGER :: ii(nx_ext)
  REAL    :: rtmp
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  GRIB grid information
!
!-----------------------------------------------------------------------
!
  INTEGER :: ni_grb       ! Number of points along x-axis
  INTEGER :: nj_grb       ! Number of points along y-axis

  REAL :: di_grb          ! x-direction increment or grid length
  REAL :: dj_grb          ! y-direction increment or grid length

  REAL :: latsw           ! Latitude  of South West corner point
  REAL :: lonsw           ! Longitude of South West corner point
  REAL :: latne           ! Latitude  of North East corner point
  REAL :: lonne           ! Longitude of North East corner point

  REAL :: lattru1         ! Latitude (1st) at which projection is true
  REAL :: lattru2         ! Latitude (2nd) at which projection is true
  REAL :: lontrue         ! Longitude      at which projection is true

  INTEGER :: uvearth

!-----------------------------------------------------------------------
!
! Definitions for GFS 0.5 degree global data
!
!-----------------------------------------------------------------------
  INTEGER, PARAMETER :: n3dvs = 9

  INTEGER, PARAMETER :: GPH = 1        ! Geopotential Height
  INTEGER, PARAMETER :: UWD = 2        ! U-component of Wind
  INTEGER, PARAMETER :: VWD = 3        ! V-component of Wind
  INTEGER, PARAMETER :: TMP = 4        ! Temperature
  INTEGER, PARAMETER :: RHD = 5        ! Relative Humidity
  INTEGER, PARAMETER :: WWD = 6        ! Pressure vertical velocity(Pa/s)
  INTEGER, PARAMETER :: QCM = 7        ! Cloud Mixing Ratio
  INTEGER, PARAMETER :: TSL = 8        ! Soil temperature
  INTEGER, PARAMETER :: QSL = 9        ! Soil Moisture

  INTEGER, PARAMETER :: n2dvs = 17

  INTEGER, PARAMETER :: TSFC = 1       ! Temperature at ground
  INTEGER, PARAMETER :: PSFC = 2       ! Pressure at ground
  INTEGER, PARAMETER :: GSFC = 3       ! Geopotential Height at ground
  INTEGER, PARAMETER :: SSFC = 4       ! Snow depth at ground (kg m-2)
  INTEGER, PARAMETER :: T2M  = 5       ! Temperature at 2m above ground
  INTEGER, PARAMETER :: QV2M = 6       ! Specific Humidity at 2m above ground
  INTEGER, PARAMETER :: U10M = 7       ! U-Component of Wind at 10m above ground
  INTEGER, PARAMETER :: V10M = 8       ! V-Component of Wind at 10m above ground

  INTEGER, PARAMETER :: PMSL = 9       ! Pressure Reduced to MSL
  INTEGER, PARAMETER :: TNER = 10      ! Temperature at 0.995 Sigma level
  INTEGER, PARAMETER :: PTNR = 11      ! Potential Temperature at 0.995 Sigma level
  INTEGER, PARAMETER :: RHNR = 12      ! Relative Humidity at 0.995 Sigma level
  INTEGER, PARAMETER :: UNER = 13      ! U-Component of Wind at 0.995 Sigma level
  INTEGER, PARAMETER :: VNER = 14      ! V-Component of Wind at 0.995 Sigma level

  INTEGER, PARAMETER :: LCVR = 15      ! V-Component of Wind at 0.995 Sigma level
  INTEGER, PARAMETER :: ICVR = 16      ! V-Component of Wind at 0.995 Sigma level
  INTEGER, PARAMETER :: PWTR = 17      ! V-Component of Wind at 0.995 Sigma level

  INTEGER, PARAMETER :: totvs = n2dvs + n3dvs
  INTEGER, PARAMETER :: varids(4,totvs) = RESHAPE(         & ! order is IMPORTANT
                   ! Discipline  Category  Parameter Layer/Level_Type
                        (/   0,   3,   5, 100,    & ! 1 GPH
                             0,   2,   2, 100,    & ! 2 UWD
                             0,   2,   3, 100,    & ! 3 VWD
                             0,   0,   0, 100,    & ! 4 TMP
                             0,   1,   1, 100,    & ! 5 RHD
                             0,   2,   8, 100,    & ! 6 WWD
                             0,   1,  22, 100,    & ! 7 QCM
                             0,   0,   0, 106,    & ! 8 TSL
                             2,   0, 192, 106,    & ! 9 QSL

                             0,   0,   0,   1,    & ! 1  - TSFC
                             0,   3,   0,   1,    & ! 2  - PSFC
                             0,   3,   5,   1,    & ! 3  - GSFC
                             0,   1,  13,   1,    & ! 4  - SSFC
                             0,   0,   0, 103,    & ! 5  - T2M
                             0,   1,   0, 103,    & ! 6  - QV2M
                             0,   2,   2, 103,    & ! 7  - U10M
                             0,   2,   3, 103,    & ! 8  - U10M
                             0,   3,   1, 101,    & ! 9  - PMSL
                             0,   0,   0, 104,    & ! 10 - TNER
                             0,   0,   2, 104,    & ! 11 - PTNR
                             0,   1,   1, 104,    & ! 12 - RHNR
                             0,   2,   2, 104,    & ! 13 - UNER
                             0,   2,   3, 104,    & ! 14 - VNER
                             2,   0,   0,   1,    & ! 15 - LCVR
                            10,   2,   0,   1,    & ! 16 - ICVR
                             0,   1,   3, 200 /), & ! 17 - PWTR
                        (/ 4, totvs /) )

  REAL,   PARAMETER :: var3dlvl(26) = (/ 100000, 97500, 95000, 92500, 90000, 85000, &
              80000, 75000, 70000, 65000, 60000, 55000, 50000, 45000, 40000, 35000, &
              30000, 25000, 20000, 15000, 10000,  7000,  5000,  3000,  2000,  1000 /)

  INTEGER, PARAMETER :: var3dindx(totvs) = (/ GPH, UWD, VWD, TMP, RHD, WWD, QCM, TSL, QSL, &
                                              (0, i = n3dvs+1,totvs) /)
  INTEGER, PARAMETER :: var2dindx(totvs) = (/ (0, i = 1, n3dvs), TSFC, PSFC, GSFC, SSFC,   &
           T2M, QV2M, U10M, V10M, PMSL, TNER, PTNR, RHNR, UNER, VNER, LCVR, ICVR, PWTR/)

  REAL,    PARAMETER :: var2dlvl(n2dvs) = (/ 0.,0.,0.,0.,2.,2.,10.,10., &
                                             (9950., i=1,6), & ! sigma = 0.995
                                             0.,0.,0./)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  RDNMCGRB reads NMC GRIB2 data
!
!-----------------------------------------------------------------------
!
  ALLOCATE(var_grb2d(nx_ext,ny_ext,n2dvs),        STAT = istatus)
  ALLOCATE(var_grb3d(nx_ext,ny_ext,nz_ext,n3dvs), STAT = istatus)
  var_grb2d = -999.0
  var_grb3d = 0.0

  CALL rdnmcgrb2(nx_ext,ny_ext,nz_ext,gribfile,grbflen, gribtime,       &
           1, nx_ext, 1, ny_ext,                                        &
           iproj_ext,ni_grb,nj_grb,di_grb,dj_grb,                       &
           latsw,lonsw,lattru1,lattru2,lontrue,uvearth,                 &
           n2dvs, n3dvs, totvs, nzsoilin_ext,                           &
           varids, var2dindx, var3dindx, var2dlvl, var3dlvl,            &
           100*soildepth_ext,                                           &
           var_grb2d, var_grb3d, lvldbg,                                &
           istatus)

  IF (istatus /= 0) RETURN

  scale_ext     = 1.0
  latnot_ext(1) = lattru1
  latnot_ext(2) = lattru2
  trlon_ext     = lontrue

  dx_ext = di_grb
  dy_ext = dj_grb

  DO j=1, ny_ext
    DO i=1, nx_ext
      lon_ext(i,j)= lonsw + (i-1) * dx_ext
      lat_ext(i,j)= latsw + (j-1) * dy_ext
    END DO
  END DO

  IF (myproc == 0) PRINT *,'LatSW = ',latsw,' LonSW = ',lonsw

  IF ( lon_0_360 ) THEN   ! Do not care whether median is in-between

    DO i = 1,nx_ext
      ii(i) = i
    END DO

  ELSE                    ! Median is in-between and 1:0 - 360:359
    IF (mod(nx_ext,2) /= 0) THEN
      WRITE(6,'(a/)') 'Wrong size of nx_ext in GETNCEPAVN3 for lon_0_360.'
      CALL arpsstop('Wrong nx_ext size.',1)
    END IF

    DO i = 1,nx_ext/2              ! map 1-180 to 181-360
      ii(i) = i + nx_ext/2
    END DO

    DO i = nx_ext/2+1, nx_ext      ! map 181-360 to 1-180
      ii(i) = i - nx_ext/2
    END DO

    WHERE (lon_ext >= 180) lon_ext = lon_ext - 360

    DO j = 1,ny_ext                ! swap lat_ext & lon_ext
      DO i = 1,nx_ext/2
        rtmp = lat_ext(ii(i),j)
        lat_ext(ii(i),j) = lat_ext(i,j)
        lat_ext(i,    j) = rtmp

        rtmp = lon_ext(ii(i),j)
        lon_ext(ii(i),j) = lon_ext(i,j)
        lon_ext(i,    j) = rtmp
      END DO
    END DO
  END IF

!
!-----------------------------------------------------------------------
!
!  Retrieve 2-D variables
!
!-----------------------------------------------------------------------
!

  WRITE(6,*) 'Filling 2D arrays from GRIB2 data.'
  DO j=1,ny_ext
    DO i=1,nx_ext
      psfc_ext   (i,j) = var_grb2d(ii(i),j,PSFC)
      trn_ext    (i,j) = var_grb2d(ii(i),j,GSFC)

      wetcanp_ext(i,j) = -999.0

!      wetcanp_ext(i,j) = var_grb2d(ii(i),j,4,1)*1.e-3     ! in meter

      IF ( var_grb2d(1,1,SSFC) > 0 ) THEN
!      Convert water equiv. of accum. snow depth (kg/m**2) to meters
!      (where 1 meter liquid water is set equivqlent to 10 meters snow).
!          0.01 = 10. (m snow/m liquid) / (1000 kg/m**3)
        snowdpth_ext(i,j) = 0.01 * var_grb2d(ii(i),j,SSFC)  ! in meters
      ELSE
        snowdpth_ext(i,j) = -999.
      END IF

      t_2m_ext(i,j)  = var_grb2d(ii(i),j,T2M)
      qv_2m_ext(i,j) = var_grb2d(ii(i),j,QV2M)
      u_10m_ext(i,j) = var_grb2d(ii(i),j,U10M)
      v_10m_ext(i,j) = var_grb2d(ii(i),j,V10M)

      soiltyp_ext(i,j)= 0          ! We do not know the soil type
      IF (var_grb2d(ii(i),j,LCVR) < 0.5) THEN  ! over see
        soiltyp_ext(i,j) = 13
      END IF
      IF (var_grb2d(ii(i),j,ICVR) > 0.5) THEN  ! Ice
        soiltyp_ext(i,j) = 12
      END IF

      rain_ext(i,j) = var_grb2d(ii(i),j,PWTR)

    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Retrieve 3-D variables
!
!-----------------------------------------------------------------------
!
  WRITE(6,*) 'Filling 3D arrays from GRIB2 data.'
  DO k=1,nz_ext
    DO j=1,ny_ext
      DO i=1,nx_ext
        p_ext  (i,j,k) = var3dlvl(k)             ! Pressure
        hgt_ext(i,j,k) = var_grb3d(ii(i),j,k,GPH)
        u_ext  (i,j,k) = var_grb3d(ii(i),j,k,UWD)    ! u wind (m/s)
        v_ext  (i,j,k) = var_grb3d(ii(i),j,k,VWD)    ! v wind (m/s)
        t_ext  (i,j,k) = var_grb3d(ii(i),j,k,TMP)    ! Temperature (K)
        qc_ext (i,j,k) = var_grb3d(ii(i),j,k,QCM)
!        qr_ext (i,j,k) = -999.
!        qi_ext (i,j,k) = -999.
!        qs_ext (i,j,k) = -999.
!        qh_ext (i,j,k) = -999.
      END DO
    END DO
  END DO

  WRITE(6,*) 'Computing saturated specific humidity.'

  CALL getqvs(nx_ext,ny_ext,nz_ext, 1,nx_ext,1,ny_ext,1,nz_ext,         &
              p_ext, t_ext, qv_ext )

  DO k=1,nz_ext
    DO j=1,ny_ext
      DO i=1,nx_ext
        qv_ext(i,j,k) = 0.01*var_grb3d(ii(i),j,k,RHD)*qv_ext(i,j,k)
      END DO
    END DO
  END DO

  WRITE(6,*) 'Filling 3D Soil arrays from GRIB2 data.'

  tsoil_ext (:,:,:) = -999.0
  qsoil_ext (:,:,:) = -999.0

  IF(soilmodel_option == 1) THEN       ! ARPS: two-layer force-restore model)

    DO j = 1, ny_ext
      DO i = 1, nx_ext
        tsoil_ext(i,j,1) =     var_grb2d(ii(i),j,TSFC)       ! sfc temp.
        tsoil_ext(i,j,2) = 0.1*var_grb3d(ii(i),j,1,TSL)                 &
                         + 0.3*var_grb3d(ii(i),j,2,TSL)                 &
                         + 0.6*var_grb3d(ii(i),j,3,TSL)

        qsoil_ext(i,j,1) =     var_grb3d(ii(i),j,1,QSL)
        qsoil_ext(i,j,2) = 0.1*var_grb3d(ii(i),j,1,QSL)                 &
                         + 0.3*var_grb3d(ii(i),j,2,QSL)                 &
                         + 0.6*var_grb3d(ii(i),j,3,QSL)

       IF (tsoil_ext(i,j,1) < 100. .OR. tsoil_ext(i,j,2) < 100) THEN  ! Over water?
       !IF (soiltyp_ext(i,j) == 13) THEN
         tsoil_ext(i,j,1) = var_grb2d(ii(i),j,TSFC)
         tsoil_ext(i,j,2) = var_grb2d(ii(i),j,TSFC)
         qsoil_ext(i,j,1) = 1.0
         qsoil_ext(i,j,2) = 1.0
       END IF
      END DO
    END DO

  ELSE                                ! ARPS: multi-layer OUSoil model
    DO k = 1,nzsoil_ext
      DO j = 1, ny_ext
        DO i = 1, nx_ext
          tsoil_ext(i,j,k) = var_grb3d(ii(i),j,k,TSL)
          qsoil_ext(i,j,k) = var_grb3d(ii(i),j,k,QSL)
          IF (tsoil_ext(i,j,k) < 100.) THEN  ! over water?
          !IF (soiltyp_ext(i,j) == 13) THEN
            tsoil_ext(i,j,k) = var_grb2d(ii(i),j,TSFC)
            qsoil_ext(i,j,k) = 1.0
          END IF
        END DO
      END DO
    END DO

  END IF

  DEALLOCATE(var_grb2d,var_grb3d, STAT = istatus)

  RETURN
END SUBROUTINE getncepgfs_grb2
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE GETMAM212_grb2                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE getnam212_grb2(nx_ext,ny_ext,nz_ext,nzsoil_ext,nzsoilin_ext,&
           gribfile,gribtime,grbflen,soildepth_ext,                    &
           dx_ext,dy_ext,                                              &
           iproj_ext,scale_ext,trlon_ext,latnot_ext,latsw_ext,lonsw_ext, &
           p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                     &
           qc_ext,                                                     &
           tsoil_ext,qsoil_ext,wetcanp_ext,                            &
           snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext,                  &
           t_2m_ext, qv_2m_ext, u_10m_ext, v_10m_ext, rain_ext,        &
           istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads NAM 40km GRIB2 data (grid #212)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  08/14/2008
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx_ext
!    ny_ext
!    nz_ext
!    nzsoil_ext
!    nstyps_ext
!    gribfile      File name to be read
!    gribtime
!
!  OUTPUT:
!
!    dx_ext
!    dy_ext
!    iproj_ext     Map projection number of external data
!    scale_ext     Scale factor of external data
!    trlon_ext     True longitude of external data (degrees E)
!    latnot_ext(2) True latitude(s) of external data (degrees N)
!    lawsw_ext     Latitude of the southwest corn of the domain
!    lonsw_ext     Longitude of the southwest corn of the domain
!    hgt_ext       height (m)
!    p_ext         pressure (Pascal)
!    t_ext         temperature (K)
!    qv_ext        specific humidity (kg/kg)
!    u_ext         u wind component (m/s)
!    v_ext         v wind component (m/s)
!    qc_ext        Cloud water mixing ratio (kg/kg)
!    qr_ext        Rain water mixing ratio (kg/kg)
!    qi_ext        Ice mixing ratio (kg/kg)
!    qs_ext        Snow mixing ratio (kg/kg)
!    qh_ext        Hail mixing ratio (kg/kg)
!
!    tsoil_ext     Soil temperature
!    qsoil_ext     Soil moisture
!    wetcanp_ext   Water content on canopy
!
!    trn_ext       External terrain (m)
!    psfc_ext      Surface pressure (Pa)
!
!    T_2m_ext      Temperature at 2m AGL
!    qv_2m_ext     Specific Humidity at 2m AGL
!    U_10m_ext     U at 10m AGL
!    V_10m_ext     V at 10m AGL
!    rain_ext      Accumulated rainfall contained in the data file
!
!    istatus       status indicator
!
!  WORK ARRAYS:
!
!    var_grb3d     Arrays to store the GRIB 3-D variables:
!    var_grb2d     Arrays to store the GRIB 2-D variables:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  CHARACTER (LEN=*), INTENT(IN) :: gribfile
  CHARACTER (LEN=*), INTENT(IN) :: gribtime
  INTEGER, INTENT(IN) :: grbflen

  INTEGER, INTENT(IN) :: nx_ext,ny_ext,nz_ext,nzsoil_ext

  INTEGER, INTENT(IN) :: nzsoilin_ext
  REAL,    INTENT(IN) :: soildepth_ext(nzsoilin_ext)

!  INTEGER, PARAMETER :: nzsoilin_ext = 4
!  REAL,    PARAMETER :: soildepth_ext(nzsoilin_ext) = (/0.0, 0.1, 0.4, 1.0/)
!
!-----------------------------------------------------------------------
!
!  External grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER, INTENT(OUT) :: iproj_ext
  REAL,    INTENT(OUT) :: scale_ext,trlon_ext
  REAL,    INTENT(OUT) :: latnot_ext(2)
  REAL,    INTENT(OUT) :: dx_ext,dy_ext
  REAL,    INTENT(OUT) :: latsw_ext, lonsw_ext
!
!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------
!
  REAL,    INTENT(OUT) :: p_ext  (nx_ext,ny_ext,nz_ext)   ! Pressure (Pascals)
  REAL,    INTENT(OUT) :: hgt_ext(nx_ext,ny_ext,nz_ext)   ! Height (m)
  REAL,    INTENT(OUT) :: t_ext  (nx_ext,ny_ext,nz_ext)   ! Temperature (K)
  REAL,    INTENT(OUT) :: qv_ext (nx_ext,ny_ext,nz_ext)   ! Specific humidity (kg/kg)
  REAL,    INTENT(OUT) :: u_ext  (nx_ext,ny_ext,nz_ext)   ! Eastward wind component
  REAL,    INTENT(OUT) :: v_ext  (nx_ext,ny_ext,nz_ext)   ! Northward wind component
  REAL,    INTENT(OUT) :: qc_ext (nx_ext,ny_ext,nz_ext)   ! Cloud H2O mixing ratio (kg/kg)
!  REAL,    INTENT(OUT) :: qr_ext (nx_ext,ny_ext,nz_ext)   ! Rain  H2O mixing ratio (kg/kg)
!  REAL,    INTENT(OUT) :: qi_ext (nx_ext,ny_ext,nz_ext)   ! Ice   mixing ratio (kg/kg)
!  REAL,    INTENT(OUT) :: qs_ext (nx_ext,ny_ext,nz_ext)   ! Snow  mixing ratio (kg/kg)
!  REAL,    INTENT(OUT) :: qh_ext (nx_ext,ny_ext,nz_ext)   ! Hail  mixing ratio (kg/kg)

  REAL,    INTENT(OUT) :: tsoil_ext  (nx_ext,ny_ext,nzsoil_ext)
  REAL,    INTENT(OUT) :: qsoil_ext  (nx_ext,ny_ext,nzsoil_ext)
  REAL,    INTENT(OUT) :: wetcanp_ext(nx_ext,ny_ext)      ! Canopy water amount
  REAL,    INTENT(OUT) :: snowdpth_ext(nx_ext,ny_ext)     ! Snow depth (m)

  REAL,    INTENT(OUT) :: trn_ext  (nx_ext,ny_ext)      ! External terrain (m)
  REAL,    INTENT(OUT) :: psfc_ext (nx_ext,ny_ext)      ! Surface pressure (Pa)
  REAL,    INTENT(OUT) :: t_2m_ext (nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: qv_2m_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: u_10m_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: v_10m_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: rain_ext (nx_ext,ny_ext)

  INTEGER, INTENT(OUT) :: soiltyp_ext(nx_ext,ny_ext)

  INTEGER, INTENT(OUT) :: istatus
!
!-----------------------------------------------------------------------
!
!  Work arrays for storing grib data
!
!-----------------------------------------------------------------------
!
  REAL,    ALLOCATABLE :: var_grb2d(:,:,:)   ! GRIB 2D variables
  REAL,    ALLOCATABLE :: var_grb3d(:,:,:,:) ! GRIB 3-D variables
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,kk

!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  GRIB grid information
!
!-----------------------------------------------------------------------
!
  INTEGER :: ni_grb       ! Number of points along x-axis
  INTEGER :: nj_grb       ! Number of points along y-axis

  REAL :: di_grb          ! x-direction increment or grid length
  REAL :: dj_grb          ! y-direction increment or grid length

  REAL :: latsw           ! Latitude  of South West corner point
  REAL :: lonsw           ! Longitude of South West corner point
!  REAL :: latne           ! Latitude  of North East corner point
!  REAL :: lonne           ! Longitude of North East corner point

  REAL :: lattru1         ! Latitude (1st) at which projection is true
  REAL :: lattru2         ! Latitude (2nd) at which projection is true
  REAL :: lontrue         ! Longitude      at which projection is true

  INTEGER :: uvearth

!-----------------------------------------------------------------------
!
! Definitions for grid #212 data
!
!-----------------------------------------------------------------------
  INTEGER, PARAMETER :: n3dvs = 9

  INTEGER, PARAMETER :: GPH = 1        ! Geopotential Height
  INTEGER, PARAMETER :: UWD = 2        ! U-component of Wind
  INTEGER, PARAMETER :: VWD = 3        ! V-component of Wind
  INTEGER, PARAMETER :: TMP = 4        ! Temperature
  INTEGER, PARAMETER :: SQV = 5        ! Specific Humidity
  INTEGER, PARAMETER :: WWD = 6        ! Pressure vertical velocity(Pa/s)
  INTEGER, PARAMETER :: QCM = 7        ! Cloud Mixing Ratio
  INTEGER, PARAMETER :: TSL = 8        ! Soil temperature
  INTEGER, PARAMETER :: QSL = 9        ! Soil Moisture

  INTEGER, PARAMETER :: n2dvs = 13

  INTEGER, PARAMETER :: TSFC = 1       ! Temperature at ground
  INTEGER, PARAMETER :: PSFC = 2       ! Pressure at ground
  INTEGER, PARAMETER :: GSFC = 3       ! Geopotential Height at ground
  INTEGER, PARAMETER :: SSFC = 4       ! Snow depth at ground (kg m-2)
  INTEGER, PARAMETER :: RAIN = 5       ! Total Precipitation (kg m-2)
  INTEGER, PARAMETER :: LAND = 6       ! Land Cover (1=land, 0=sea)
  INTEGER, PARAMETER :: ICEC = 7       ! Ice Cover (Proportion)
  INTEGER, PARAMETER :: SOIL = 8       ! Soil Type

  INTEGER, PARAMETER :: T2M  = 9       ! Temperature at 2m above ground
  INTEGER, PARAMETER :: QV2M = 10      ! Specific Humidity at 2m above ground
  INTEGER, PARAMETER :: U10M = 11      ! U-Component of Wind at 10m above ground
  INTEGER, PARAMETER :: V10M = 12      ! V-Component of Wind at 10m above ground

  INTEGER, PARAMETER :: PMSL = 13      ! Pressure Reduced to MSL

  INTEGER, PARAMETER :: totvs = n2dvs + n3dvs
  INTEGER, PARAMETER :: varids(4,totvs) = RESHAPE(         & ! order is IMPORTANT
                   ! Discipline  Category  Parameter Layer/Level_Type
                        (/   0,   3,   5, 100,    & ! 1 GPH
                             0,   2,   2, 100,    & ! 2 UWD
                             0,   2,   3, 100,    & ! 3 VWD
                             0,   0,   0, 100,    & ! 4 TMP
                             0,   1,   0, 100,    & ! 5 SQV
                             0,   2,   8, 100,    & ! 6 WWD
                             0,   1,  22, 100,    & ! 7 QCM
                             2,   0,   2, 106,    & ! 8 TSL
                             2,   0, 192, 106,    & ! 9 QSL

                             0,   0,   0,   1,    & ! 1  - TSFC
                             0,   3,   0,   1,    & ! 2  - PSFC
                             0,   3,   5,   1,    & ! 3  - GSFC
                             0,   1,  11,   1,    & ! 4  - SSFC
                             0,   1,   8,   1,    & ! 5  - RAIN
                             2,   0,   0,   1,    & ! 6  - LAND
                            10,   2,   0,   1,    & ! 7  - ICEC
                             2,   3,   0,   1,    & ! 8  - SOIL
                             0,   0,   0, 103,    & ! 9  - T2M
                             0,   1,   0, 103,    & ! 10 - QV2M
                             0,   2,   2, 103,    & ! 11 - U10M
                             0,   2,   3, 103,    & ! 12 - U10M
                             0,   3,   1, 101     & ! 13 - PMSL
                             /), (/ 4, totvs /) )

  REAL,   PARAMETER :: var3dlvl(39) = (/ 100000, 97500, 95000, 92500, 90000, 87500, &
              85000, 82500, 80000, 77500, 75000, 72500, 70000, 67500, 65000, 62500, &
              60000, 57500, 55000, 52500, 50000, 47500, 45000, 42500, 40000, 37500, &
              35000, 32500, 30000, 27500, 25000, 22500, 20000, 17500, 15000, 12500, &
              10000,  7500,  5000 /)

  INTEGER, PARAMETER :: var3dindx(totvs) = (/ GPH, UWD, VWD, TMP, SQV,  &
                        WWD, QCM, TSL, QSL, (0, i = n3dvs+1,totvs) /)

  INTEGER, PARAMETER :: var2dindx(totvs) = (/ (0, i = 1, n3dvs),        &
                        TSFC, PSFC, GSFC, SSFC, RAIN, LAND, ICEC, SOIL, &
                        T2M, QV2M, U10M, V10M, PMSL/)

  REAL,    PARAMETER :: var2dlvl(n2dvs) = (/ 0.,0.,0.,0.,0.,0.,0.,0.,   &
                                             2.,2.,10.,10.,0. /)
!
!
! This table converts WRF soil to those used in ARPS
!
! WRF uses 16 soil categories (FAO/STATSGO)
! ARPS uses 13-category soil  (Standard GRIB or GRIB2 definitions)
!
  INTEGER, PARAMETER :: soil_table(17) = (/ 1, 2, 3, 4, 4,              &
                                            5, 6, 7, 8, 9,              &
                                           10,11, 6,13, 1,              &
                                            2, 2/)

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
!  RDNMCGRB2 reads NMC GRIB2 data
!
!-----------------------------------------------------------------------
!
  ALLOCATE(var_grb2d(nx_ext,ny_ext,n2dvs),        STAT = istatus)
  ALLOCATE(var_grb3d(nx_ext,ny_ext,nz_ext,n3dvs), STAT = istatus)
  var_grb2d = -999.0
  var_grb3d = 0.0

  CALL rdnmcgrb2(nx_ext,ny_ext,nz_ext,gribfile,grbflen, gribtime,       &
           1, nx_ext, 1, ny_ext,                                        &
           iproj_ext,ni_grb,nj_grb,di_grb,dj_grb,                       &
           latsw,lonsw,lattru1,lattru2,lontrue,uvearth,                 &
           n2dvs, n3dvs, totvs, nzsoilin_ext,                           &
           varids, var2dindx, var3dindx, var2dlvl, var3dlvl,            &
           100*soildepth_ext,                                           &
           var_grb2d, var_grb3d, lvldbg,                                &
           istatus)

  IF (istatus /= 0) RETURN

  scale_ext     = 1.0
  latnot_ext(1) = lattru1
  latnot_ext(2) = lattru2
  trlon_ext     = lontrue

  dx_ext = di_grb
  dy_ext = dj_grb

  IF (myproc == 0) PRINT *,'LatSW = ',latsw,' LonSW = ',lonsw

  latsw_ext = latsw
  lonsw_ext = lonsw
!
!-----------------------------------------------------------------------
!
!  Retrieve 2-D variables
!
!-----------------------------------------------------------------------
!

  IF (myproc == 0 .AND. lvldbg > 0) WRITE(6,*) 'Filling 2D arrays from GRIB2 data.'
  DO j=1,ny_ext
    DO i=1,nx_ext
      psfc_ext   (i,j) = var_grb2d(i,j,PSFC)
      trn_ext    (i,j) = var_grb2d(i,j,GSFC)
      snowdpth_ext(i,j) = var_grb2d(i,j,SSFC)
      wetcanp_ext(i,j) = -999.0

      t_2m_ext (i,j) = var_grb2d(i,j,T2M)
      qv_2m_ext(i,j) = var_grb2d(i,j,QV2M)
      u_10m_ext(i,j) = var_grb2d(i,j,U10M)
      v_10m_ext(i,j) = var_grb2d(i,j,V10M)
      rain_ext (i,j) = var_grb2d(i,j,RAIN)

      !
      ! To make sure SOILTYPE and LANDCOVER consistent and
      ! convert SOILTYPE from to ARPS.
      !
      IF (var_grb2d(i,j,LAND) < 0.5) THEN      ! Water
        soiltyp_ext(i,j)= 13
      ELSE IF (var_grb2d(i,j,ICEC) > 0.5) THEN ! Ice
        soiltyp_ext(i,j) = 12
      !
      ! NAM grib2 files have soil type in them.  SREF doesn't.  Fall back
      ! to the grib1 method of just assigning 0 to the soil type.
      !
      ELSE IF (var_grb2d(i,j,SOIL) == -999.0) THEN
        soiltyp_ext(i,j) = 0
      ELSE IF (var_grb2d(i,j,SOIL) < 1. .OR. var_grb2d(i,j,SOIL) > 16.) THEN
        WRITE(6,'(1x,2(a,I4),a,I2,a,F5.0)')       &
                'WARNING: Soil type is unknown. var_grb2d(',i,',',j,    &
                ',SOIL=',SOIL,') = ',var_grb2d(i,j,SOIL)
      ELSE
        soiltyp_ext(i,j) = SOIL_TABLE(NINT(var_grb2d(i,j,SOIL)))
      END IF

    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Retrieve 3-D variables
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0 .AND. lvldbg > 0) WRITE(6,*) 'Filling 3D arrays from GRIB2 data.'
  DO k=1,nz_ext
    DO j=1,ny_ext
      DO i=1,nx_ext
        p_ext  (i,j,k) = var3dlvl(k)             ! Pressure
        hgt_ext(i,j,k) = var_grb3d(i,j,k,GPH)
        u_ext  (i,j,k) = var_grb3d(i,j,k,UWD)    ! u wind (m/s)
        v_ext  (i,j,k) = var_grb3d(i,j,k,VWD)    ! v wind (m/s)
        t_ext  (i,j,k) = var_grb3d(i,j,k,TMP)    ! Temperature (K)
        qv_ext (i,j,k) = var_grb3d(i,j,k,SQV)    ! Specific humidity (kg kg-1)
        qc_ext (i,j,k) = var_grb3d(i,j,k,QCM)
!        qr_ext (i,j,k) = -999.
!        qi_ext (i,j,k) = -999.
!        qs_ext (i,j,k) = -999.
!        qh_ext (i,j,k) = -999.
      END DO
    END DO
  END DO

  IF (myproc == 0 .AND. lvldbg > 0) WRITE(6,*) 'Filling 3D Soil arrays from GRIB2 data.'

  tsoil_ext (:,:,:) = -999.0
  qsoil_ext (:,:,:) = -999.0

  IF(soilmodel_option == 1) THEN       ! ARPS: two-layer force-restore model)

    DO j = 1, ny_ext
      DO i = 1, nx_ext
        tsoil_ext(i,j,1) = 0.5*( var_grb3d(i,j,1,TSL) + var_grb3d(i,j,2,TSL) )      ! sfc temp.
        tsoil_ext(i,j,2) = 0.5*( var_grb3d(i,j,3,TSL) + var_grb3d(i,j,4,TSL) )
        qsoil_ext(i,j,1) = 0.5*( var_grb3d(i,j,1,QSL) + var_grb3d(i,j,2,QSL) )      ! sfc moisture
        qsoil_ext(i,j,2) = 0.5*( var_grb3d(i,j,3,QSL) + var_grb3d(i,j,4,QSL) )
       IF (var_grb2d(i,j,LAND) < .5 ) THEN  ! Over water?
         tsoil_ext(i,j,1) = var_grb2d(i,j,TSFC)
         tsoil_ext(i,j,2) = var_grb2d(i,j,TSFC)
         qsoil_ext(i,j,1) = 1.0
         qsoil_ext(i,j,2) = 1.0
       END IF
      END DO
    END DO

  ELSE                                ! ARPS: multi-layer OUSoil model
    DO j = 1, ny_ext
      DO i = 1, nx_ext
        tsoil_ext(i,j,1) = var_grb2d(i,j,TSFC)
        qsoil_ext(i,j,1) = var_grb3d(i,j,1,QSL)
      END DO
    END DO

    DO k = 2,nzsoil_ext
      DO j = 1, ny_ext
        DO i = 1, nx_ext
          tsoil_ext(i,j,k) = var_grb3d(i,j,k-1,TSL)
          qsoil_ext(i,j,k) = var_grb3d(i,j,k-1,QSL)
          IF (var_grb2d(i,j,LAND) < .5) THEN  ! over water?
            tsoil_ext(i,j,k) = var_grb2d(i,j,TSFC)
            qsoil_ext(i,j,k) = 1.0
          END IF
        END DO
      END DO
    END DO

  END IF

  DEALLOCATE(var_grb2d,var_grb3d, STAT = istatus)

  RETURN
END SUBROUTINE getnam212_grb2
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE GETMAM218_grb2                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE getnam218_grb2(notiles,nx_ext,ny_ext,nz_ext,nzsoil_ext,      &
           nzsoilin_ext,nxgrb,nygrb,tile_x,tile_y,                      &
           gribfile,grbflen,gribtime,                                   &
           soildepth_ext,dx_ext,dy_ext,                                 &
           iproj_ext,scale_ext,trlon_ext,latnot_ext,latsw_ext,lonsw_ext,&
           p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                      &
           qc_ext,                                                      &
           tsoil_ext,qsoil_ext,wetcanp_ext,                             &
           snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext,                   &
           t_2m_ext, qv_2m_ext, u_10m_ext, v_10m_ext, rain_ext,         &
           istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads NAM 12km GRIB2 data (grid #218)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  08/21/2008
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx_ext
!    ny_ext
!    nz_ext
!    nzsoil_ext
!    nstyps_ext
!    gribfile      File name to be read
!    gribtime
!
!  OUTPUT:
!
!    dx_ext
!    dy_ext
!    iproj_ext     Map projection number of external data
!    scale_ext     Scale factor of external data
!    trlon_ext     True longitude of external data (degrees E)
!    latnot_ext(2) True latitude(s) of external data (degrees N)
!    lawsw_ext     Latitude of the southwest corn of the domain
!    lonsw_ext     Longitude of the southwest corn of the domain
!    hgt_ext       height (m)
!    p_ext         pressure (Pascal)
!    t_ext         temperature (K)
!    qv_ext        specific humidity (kg/kg)
!    u_ext         u wind component (m/s)
!    v_ext         v wind component (m/s)
!    qc_ext        Cloud water mixing ratio (kg/kg)
!    qr_ext        Rain water mixing ratio (kg/kg)
!    qi_ext        Ice mixing ratio (kg/kg)
!    qs_ext        Snow mixing ratio (kg/kg)
!    qh_ext        Hail mixing ratio (kg/kg)
!
!    tsoil_ext     Soil temperature
!    qsoil_ext     Soil moisture
!    wetcanp_ext   Water content on canopy
!
!    trn_ext       External terrain (m)
!    psfc_ext      Surface pressure (Pa)
!
!    T_2m_ext      Temperature at 2m AGL
!    qv_2m_ext     Specific Humidity at 2m AGL
!    U_10m_ext     U at 10m AGL
!    V_10m_ext     V at 10m AGL
!    rain_ext      Accumulated rainfall contained in the data file
!
!    istatus       status indicator
!
!  WORK ARRAYS:
!
!    var_grb3d     Arrays to store the GRIB 3-D variables:
!    var_grb2d     Arrays to store the GRIB 2-D variables:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  CHARACTER (LEN=*), INTENT(IN) :: gribfile
  CHARACTER (LEN=*), INTENT(IN) :: gribtime
  INTEGER, INTENT(IN) :: grbflen

  INTEGER, INTENT(IN) :: notiles
  INTEGER, INTENT(IN) :: tile_x, tile_y

  INTEGER, INTENT(IN) :: nxgrb, nygrb
  INTEGER, INTENT(IN) :: nx_ext,ny_ext,nz_ext,nzsoil_ext

  INTEGER, INTENT(IN) :: nzsoilin_ext
  REAL,    INTENT(IN) :: soildepth_ext(nzsoilin_ext)

!  INTEGER, PARAMETER :: nzsoilin_ext = 4
!  REAL,    PARAMETER :: soildepth_ext(nzsoilin_ext) = (/0.0, 0.1, 0.4, 1.0/)
!
!-----------------------------------------------------------------------
!
!  External grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER, INTENT(OUT) :: iproj_ext
  REAL,    INTENT(OUT) :: scale_ext,trlon_ext
  REAL,    INTENT(OUT) :: latnot_ext(2)
  REAL,    INTENT(OUT) :: dx_ext,dy_ext
  REAL,    INTENT(OUT) :: latsw_ext, lonsw_ext
!
!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------
!
  REAL,    INTENT(OUT) :: p_ext  (nx_ext,ny_ext,nz_ext)   ! Pressure (Pascals)
  REAL,    INTENT(OUT) :: hgt_ext(nx_ext,ny_ext,nz_ext)   ! Height (m)
  REAL,    INTENT(OUT) :: t_ext  (nx_ext,ny_ext,nz_ext)   ! Temperature (K)
  REAL,    INTENT(OUT) :: qv_ext (nx_ext,ny_ext,nz_ext)   ! Specific humidity (kg/kg)
  REAL,    INTENT(OUT) :: u_ext  (nx_ext,ny_ext,nz_ext)   ! Eastward wind component
  REAL,    INTENT(OUT) :: v_ext  (nx_ext,ny_ext,nz_ext)   ! Northward wind component
  REAL,    INTENT(OUT) :: qc_ext (nx_ext,ny_ext,nz_ext)   ! Cloud H2O mixing ratio (kg/kg)
!  REAL,    INTENT(OUT) :: qr_ext (nx_ext,ny_ext,nz_ext)   ! Rain  H2O mixing ratio (kg/kg)
!  REAL,    INTENT(OUT) :: qi_ext (nx_ext,ny_ext,nz_ext)   ! Ice   mixing ratio (kg/kg)
!  REAL,    INTENT(OUT) :: qs_ext (nx_ext,ny_ext,nz_ext)   ! Snow  mixing ratio (kg/kg)
!  REAL,    INTENT(OUT) :: qh_ext (nx_ext,ny_ext,nz_ext)   ! Hail  mixing ratio (kg/kg)

  REAL,    INTENT(OUT) :: tsoil_ext  (nx_ext,ny_ext,nzsoil_ext)
  REAL,    INTENT(OUT) :: qsoil_ext  (nx_ext,ny_ext,nzsoil_ext)
  REAL,    INTENT(OUT) :: wetcanp_ext(nx_ext,ny_ext)      ! Canopy water amount
  REAL,    INTENT(OUT) :: snowdpth_ext(nx_ext,ny_ext)     ! Snow depth (m)

  REAL,    INTENT(OUT) :: trn_ext  (nx_ext,ny_ext)      ! External terrain (m)
  REAL,    INTENT(OUT) :: psfc_ext (nx_ext,ny_ext)      ! Surface pressure (Pa)
  REAL,    INTENT(OUT) :: t_2m_ext (nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: qv_2m_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: u_10m_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: v_10m_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: rain_ext (nx_ext,ny_ext)

  INTEGER, INTENT(OUT) :: soiltyp_ext(nx_ext,ny_ext)

  INTEGER, INTENT(OUT) :: istatus
!
!-----------------------------------------------------------------------
!
!  Work arrays for storing grib data
!
!-----------------------------------------------------------------------
!
  REAL,    ALLOCATABLE :: var_grb2d(:,:,:)   ! GRIB 2D variables
  REAL,    ALLOCATABLE :: var_grb3d(:,:,:,:) ! GRIB 3-D variables

!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  GRIB grid information
!
!-----------------------------------------------------------------------
!
  INTEGER :: ni_grb       ! Number of points along x-axis
  INTEGER :: nj_grb       ! Number of points along y-axis

  REAL :: di_grb          ! x-direction increment or grid length
  REAL :: dj_grb          ! y-direction increment or grid length

  REAL :: latsw           ! Latitude  of South West corner point
  REAL :: lonsw           ! Longitude of South West corner point
!  REAL :: latne           ! Latitude  of North East corner point
!  REAL :: lonne           ! Longitude of North East corner point

  REAL :: lattru1         ! Latitude (1st) at which projection is true
  REAL :: lattru2         ! Latitude (2nd) at which projection is true
  REAL :: lontrue         ! Longitude      at which projection is true

  INTEGER :: uvearth

  INTEGER :: i,j,k,kk

!-----------------------------------------------------------------------
!
! Definitions for grid #218 data
!
!-----------------------------------------------------------------------
  INTEGER, PARAMETER :: n3dvs = 9

  INTEGER, PARAMETER :: GPH = 1        ! Geopotential Height
  INTEGER, PARAMETER :: UWD = 2        ! U-component of Wind
  INTEGER, PARAMETER :: VWD = 3        ! V-component of Wind
  INTEGER, PARAMETER :: TMP = 4        ! Temperature
  INTEGER, PARAMETER :: SRH = 5        ! Relative Humidity
  INTEGER, PARAMETER :: WWD = 6        ! Pressure vertical velocity(Pa/s)
  INTEGER, PARAMETER :: QCM = 7        ! Cloud Mixing Ratio
  INTEGER, PARAMETER :: TSL = 8        ! Soil temperature
  INTEGER, PARAMETER :: QSL = 9        ! Soil Moisture

  INTEGER, PARAMETER :: n2dvs = 14

  INTEGER, PARAMETER :: TSFC = 1       ! Temperature at ground
  INTEGER, PARAMETER :: PSFC = 2       ! Pressure at ground
  INTEGER, PARAMETER :: GSFC = 3       ! Geopotential Height at ground
  INTEGER, PARAMETER :: SSFC = 4       ! Snow depth at ground (kg m-2)
  INTEGER, PARAMETER :: RAIN = 5       ! Total Precipitation (kg m-2)
  INTEGER, PARAMETER :: LAND = 6       ! Land Cover (1=land, 0=sea)
  INTEGER, PARAMETER :: ICEC = 7       ! Ice Cover (Proportion)
  INTEGER, PARAMETER :: SOIL = 8       ! Soil Type
  INTEGER, PARAMETER :: CANP = 9       ! Canopy surface water (kg m-2)

  INTEGER, PARAMETER :: T2M  = 10      ! Temperature at 2m above ground
  INTEGER, PARAMETER :: QV2M = 11      ! Specific Humidity at 2m above ground
  INTEGER, PARAMETER :: U10M = 12      ! U-Component of Wind at 10m above ground
  INTEGER, PARAMETER :: V10M = 13      ! V-Component of Wind at 10m above ground

  INTEGER, PARAMETER :: PMSL = 14      ! Pressure Reduced to MSL

  INTEGER, PARAMETER :: totvs = n2dvs + n3dvs
  INTEGER, PARAMETER :: varids(4,totvs) = RESHAPE(         & ! order is IMPORTANT
                   ! Discipline  Category  Parameter Layer/Level_Type
                        (/   0,   3,   5, 100,    & ! 1 GPH
                             0,   2,   2, 100,    & ! 2 UWD
                             0,   2,   3, 100,    & ! 3 VWD
                             0,   0,   0, 100,    & ! 4 TMP
                             0,   1,   1, 100,    & ! 5 SRH
                             0,   2,   8, 100,    & ! 6 WWD
                             0,   1,  22, 100,    & ! 7 QCM
                             2,   0,   2, 106,    & ! 8 TSL
                             2,   0, 192, 106,    & ! 9 QSL

                             0,   0,   0,   1,    & ! 1  - TSFC
                             0,   3,   0,   1,    & ! 2  - PSFC
                             0,   3,   5,   1,    & ! 3  - GSFC
                             0,   1,  11,   1,    & ! 4  - SSFC
                             0,   1,   8,   1,    & ! 5  - RAIN
                             2,   0,   0,   1,    & ! 6  - LAND
                            10,   2,   0,   1,    & ! 7  - ICEC
                             2,   3,   0,   1,    & ! 8  - SOIL
                             2,   0, 196,   1,    & ! 9  - CANP
                             0,   0,   0, 103,    & ! 10 - T2M
                             0,   1,   1, 103,    & ! 11 - QV2M
                             0,   2,   2, 103,    & ! 12 - U10M
                             0,   2,   3, 103,    & ! 13 - U10M
                             0,   3,   1, 101     & ! 14 - PMSL
                             /), (/ 4, totvs /) )

  REAL,   PARAMETER :: var3dlvl(39) = (/ 100000, 97500, 95000, 92500, 90000, 87500, &
              85000, 82500, 80000, 77500, 75000, 72500, 70000, 67500, 65000, 62500, &
              60000, 57500, 55000, 52500, 50000, 47500, 45000, 42500, 40000, 37500, &
              35000, 32500, 30000, 27500, 25000, 22500, 20000, 17500, 15000, 12500, &
              10000,  7500,  5000 /)

  INTEGER, PARAMETER :: var3dindx(totvs) = (/ GPH, UWD, VWD, TMP, SRH,  &
                        WWD, QCM, TSL, QSL, (0, i = n3dvs+1,totvs) /)

  INTEGER, PARAMETER :: var2dindx(totvs) = (/ (0, i = 1, n3dvs),        &
                        TSFC, PSFC, GSFC, SSFC, RAIN, LAND, ICEC, SOIL, &
                        CANP,  T2M, QV2M, U10M, V10M, PMSL/)

  REAL,    PARAMETER :: var2dlvl(n2dvs) = (/ 0.,0.,0.,0.,0.,0.,0.,0.,   &
                                             0.,2.,2.,10.,10.,0. /)
!
!
! This table converts WRF soil to those used in ARPS
!
! WRF uses 16 soil categories (FAO/STATSGO)
! ARPS uses 13-category soil  (Standard GRIB or GRIB2 definitions)
!
  INTEGER, PARAMETER :: soil_table(17) = (/ 1, 2, 3, 4, 4,              &
                                            5, 6, 7, 8, 9,              &
                                           10,11, 6,13, 1,              &
                                            2, 2/)

!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: offset_x, offset_y

  INTEGER, PARAMETER :: nxmax = 69, nymax = 72
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
!  RDNMCGRB2 reads NMC GRIB2 data
!
!-----------------------------------------------------------------------
!
  ALLOCATE(var_grb2d(nxgrb,nygrb,n2dvs),        STAT = istatus)
  ALLOCATE(var_grb3d(nxgrb,nygrb,nz_ext,n3dvs), STAT = istatus)
  var_grb2d = -999.0
  var_grb3d = 0.0

  CALL rdnmcgrb2(nxgrb,nygrb,nz_ext,gribfile,grbflen, gribtime,         &
           1, nxgrb, 1, nygrb,                                          &
           iproj_ext,ni_grb,nj_grb,di_grb,dj_grb,                       &
           latsw,lonsw,lattru1,lattru2,lontrue,uvearth,                 &
           n2dvs, n3dvs, totvs, nzsoilin_ext,                           &
           varids, var2dindx, var3dindx, var2dlvl, var3dlvl,            &
           100*soildepth_ext,                                           &
           var_grb2d, var_grb3d, lvldbg,                                &
           istatus)

  IF (istatus /= 0) RETURN

  scale_ext     = 1.0
  latnot_ext(1) = lattru1
  latnot_ext(2) = lattru2
  trlon_ext     = lontrue

  dx_ext = di_grb
  dy_ext = dj_grb

  IF (myproc == 0) PRINT *,'LatSW = ',latsw,' LonSW = ',lonsw

  latsw_ext = latsw
  lonsw_ext = lonsw

  offset_x = (tile_x-1)*nxmax
  offset_y = (tile_y-1)*nymax
!
!-----------------------------------------------------------------------
!
!  Retrieve 2-D variables
!
!-----------------------------------------------------------------------
!

  IF (myproc == 0 .AND. lvldbg > 0) WRITE(6,*) 'Filling 2D arrays from GRIB2 data.'
  DO j=1,nygrb
    DO i=1,nxgrb
      psfc_ext    (offset_x+i,offset_y+j) = var_grb2d(i,j,PSFC)
      trn_ext     (offset_x+i,offset_y+j) = var_grb2d(i,j,GSFC)
      snowdpth_ext(offset_x+i,offset_y+j) = var_grb2d(i,j,SSFC)
      wetcanp_ext (offset_x+i,offset_y+j) = var_grb2d(i,j,CANP)*1.e-3     ! in meter

      t_2m_ext (offset_x+i,offset_y+j) = var_grb2d(i,j,T2M)
      qv_2m_ext(offset_x+i,offset_y+j) = var_grb2d(i,j,QV2M)
      u_10m_ext(offset_x+i,offset_y+j) = var_grb2d(i,j,U10M)
      v_10m_ext(offset_x+i,offset_y+j) = var_grb2d(i,j,V10M)
      rain_ext (offset_x+i,offset_y+j) = var_grb2d(i,j,RAIN)

      !
      ! To make sure SOILTYPE and LANDCOVER consistent and
      ! convert SOILTYPE from STATSGO to ARPS.
      !
      IF (var_grb2d(i,j,LAND) < 0.5) THEN      ! Water
        soiltyp_ext(offset_x+i,offset_y+j)= 13
      ELSE IF (var_grb2d(i,j,ICEC) > 0.5) THEN ! Ice
        soiltyp_ext(offset_x+i,offset_y+j) = 12
      ELSE IF (var_grb2d(i,j,SOIL) < 1. .OR. var_grb2d(i,j,SOIL) > 16.) THEN
        WRITE(6,'(1x,2(a,I4),a,I2,a,F3.0)')       &
                'WARNING: Soil type is unknown. var_grb2d(',offset_x+i,',',offset_y+j, &
                ',SOIL=',SOIL,') = ',var_grb2d(i,j,SOIL)
      ELSE
        soiltyp_ext(offset_x+i,offset_y+j) = SOIL_TABLE(NINT(var_grb2d(i,j,SOIL)))
      END IF

    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Retrieve 3-D variables
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0 .AND. lvldbg > 0) WRITE(6,*) 'Filling 3D arrays from GRIB2 data.'
  DO k=1,nz_ext
    DO j=1,nygrb
      DO i=1,nxgrb
        p_ext  (offset_x+i,offset_y+j,k) = var3dlvl(k)             ! Pressure
        hgt_ext(offset_x+i,offset_y+j,k) = var_grb3d(i,j,k,GPH)
        u_ext  (offset_x+i,offset_y+j,k) = var_grb3d(i,j,k,UWD)    ! u wind (m/s)
        v_ext  (offset_x+i,offset_y+j,k) = var_grb3d(i,j,k,VWD)    ! v wind (m/s)
        t_ext  (offset_x+i,offset_y+j,k) = var_grb3d(i,j,k,TMP)    ! Temperature (K)
        qv_ext (offset_x+i,offset_y+j,k) = var_grb3d(i,j,k,SRH)    ! Relative Humidity (%)
        qc_ext (offset_x+i,offset_y+j,k) = var_grb3d(i,j,k,QCM)
!        qr_ext (offset_x+i,offset_y+j,k) = -999.
!        qi_ext (offset_x+i,offset_y+j,k) = -999.
!        qs_ext (offset_x+i,offset_y+j,k) = -999.
!        qh_ext (offset_x+i,offset_y+j,k) = -999.
      END DO
    END DO
  END DO

  IF (myproc == 0 .AND. lvldbg > 0) WRITE(6,*) 'Filling 3D Soil arrays from GRIB2 data.'

  DO k = 1,nzsoil_ext
    DO j = 1, nygrb
      DO i = 1, nxgrb
        tsoil_ext (offset_x+i,offset_y+j,k) = -999.0
        qsoil_ext (offset_x+i,offset_y+j,k) = -999.0
      END DO
    END DO
  END DO

  IF(soilmodel_option == 1) THEN       ! ARPS: two-layer force-restore model)

    DO j = 1, nygrb
      DO i = 1, nxgrb
        ! sfc temp.
        tsoil_ext(offset_x+i,offset_y+j,1) = 0.5*( var_grb3d(i,j,1,TSL) + var_grb3d(i,j,2,TSL) )
        tsoil_ext(offset_x+i,offset_y+j,2) = 0.5*( var_grb3d(i,j,3,TSL) + var_grb3d(i,j,4,TSL) )
        ! sfc moisture
        qsoil_ext(offset_x+i,offset_y+j,1) = 0.5*( var_grb3d(i,j,1,QSL) + var_grb3d(i,j,2,QSL) )
        qsoil_ext(offset_x+i,offset_y+j,2) = 0.5*( var_grb3d(i,j,3,QSL) + var_grb3d(i,j,4,QSL) )
       IF (var_grb2d(i,j,LAND) < .5 ) THEN  ! Over water?
         tsoil_ext(offset_x+i,offset_y+j,1) = var_grb2d(i,j,TSFC)
         tsoil_ext(offset_x+i,offset_y+j,2) = var_grb2d(i,j,TSFC)
         qsoil_ext(offset_x+i,offset_y+j,1) = 1.0
         qsoil_ext(offset_x+i,offset_y+j,2) = 1.0
       END IF
      END DO
    END DO

  ELSE                                ! ARPS: multi-layer OUSoil model
    DO j = 1, nygrb
      DO i = 1, nxgrb
        tsoil_ext(offset_x+i,offset_y+j,1) = var_grb2d(i,j,TSFC)
        qsoil_ext(offset_x+i,offset_y+j,1) = var_grb3d(i,j,1,QSL)
      END DO
    END DO

    DO k = 2,nzsoil_ext
      DO j = 1, nygrb
        DO i = 1, nxgrb
          tsoil_ext(offset_x+i,offset_y+j,k) = var_grb3d(i,j,k-1,TSL)
          qsoil_ext(offset_x+i,offset_y+j,k) = var_grb3d(i,j,k-1,QSL)
          IF (var_grb2d(i,j,LAND) < .5) THEN  ! over water?
            tsoil_ext(offset_x+i,offset_y+j,k) = var_grb2d(i,j,TSFC)
            qsoil_ext(offset_x+i,offset_y+j,k) = 1.0
          END IF
        END DO
      END DO
    END DO

  END IF

  DEALLOCATE(var_grb2d,var_grb3d, STAT = istatus)

  RETURN
END SUBROUTINE getnam218_grb2
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE GETRUC236n_grb2                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE getruc236n_grb2(nx_ext,ny_ext,nz_ext,nzsoilin_ext,          &
           gribfile,gribtime,grbflen,soildepth_ext,                    &
           dx_ext,dy_ext,                                              &
           iproj_ext,scale_ext,trlon_ext,latnot_ext,latsw_ext,lonsw_ext, &
           p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                     &
           qc_ext,qr_ext,qi_ext,qs_ext,qh_ext,                         &
           tsfc_ext,tdeep_ext,wetsfc_ext,wetdp_ext,wetcanp_ext,        &
           snowdpth_ext,trn_ext,psfc_ext,t_2m_ext,qv_2m_ext,           &
           u_10m_ext,v_10m_ext,rain_ext,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads NAM 40km GRIB2 data (grid #212)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  12/31/2008
!
!  MODIFICATION HISTORY:
!
!  06/17/2009 (Y. Wang)
!  Added near surface fields.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx_ext
!    ny_ext
!    nz_ext
!    nzsoil_ext
!    nstyps_ext
!    gribfile      File name to be read
!    gribtime
!
!  OUTPUT:
!
!    dx_ext
!    dy_ext
!    iproj_ext     Map projection number of external data
!    scale_ext     Scale factor of external data
!    trlon_ext     True longitude of external data (degrees E)
!    latnot_ext(2) True latitude(s) of external data (degrees N)
!    lawsw_ext     Latitude of the southwest corn of the domain
!    lonsw_ext     Longitude of the southwest corn of the domain
!    hgt_ext       height (m)
!    p_ext         pressure (Pascal)
!    t_ext         temperature (K)
!    qv_ext        specific humidity (kg/kg)
!    u_ext         u wind component (m/s)
!    v_ext         v wind component (m/s)
!    qc_ext        Cloud water mixing ratio (kg/kg)
!    qr_ext        Rain water mixing ratio (kg/kg)
!    qi_ext        Ice mixing ratio (kg/kg)
!    qs_ext        Snow mixing ratio (kg/kg)
!    qh_ext        Hail mixing ratio (kg/kg)
!
!    tsoil_ext     Soil temperature
!    qsoil_ext     Soil moisture
!    wetcanp_ext   Water content on canopy
!
!    trn_ext       External terrain (m)
!    psfc_ext      Surface pressure (Pa)
!
!    T_2m_ext      Temperature at 2m AGL
!    qv_2m_ext     Specific Humidity at 2m AGL
!    U_10m_ext     U at 10m AGL
!    V_10m_ext     V at 10m AGL
!    rain_ext      Accumulated rainfall contained in the data file
!
!    istatus       status indicator
!
!  WORK ARRAYS:
!
!    var_grb3d     Arrays to store the GRIB 3-D variables:
!    var_grb2d     Arrays to store the GRIB 2-D variables:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  CHARACTER (LEN=*), INTENT(IN) :: gribfile
  CHARACTER (LEN=*), INTENT(IN) :: gribtime
  INTEGER, INTENT(IN) :: grbflen

  INTEGER, INTENT(IN) :: nx_ext,ny_ext,nz_ext

  INTEGER, INTENT(IN) :: nzsoilin_ext
  REAL,    INTENT(IN) :: soildepth_ext(nzsoilin_ext)
!
!-----------------------------------------------------------------------
!
!  External grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER, INTENT(OUT) :: iproj_ext
  REAL,    INTENT(OUT) :: scale_ext,trlon_ext
  REAL,    INTENT(OUT) :: latnot_ext(2)
  REAL,    INTENT(OUT) :: dx_ext,dy_ext
  REAL,    INTENT(OUT) :: latsw_ext, lonsw_ext
!
!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------
!
  REAL,    INTENT(OUT) :: p_ext  (nx_ext,ny_ext,nz_ext)   ! Pressure (Pascals)
  REAL,    INTENT(OUT) :: hgt_ext(nx_ext,ny_ext,nz_ext)   ! Height (m)
  REAL,    INTENT(OUT) :: t_ext  (nx_ext,ny_ext,nz_ext)   ! Temperature (K)
  REAL,    INTENT(OUT) :: qv_ext (nx_ext,ny_ext,nz_ext)   ! Specific humidity (kg/kg)
  REAL,    INTENT(OUT) :: u_ext  (nx_ext,ny_ext,nz_ext)   ! Eastward wind component
  REAL,    INTENT(OUT) :: v_ext  (nx_ext,ny_ext,nz_ext)   ! Northward wind component
  REAL,    INTENT(OUT) :: qc_ext (nx_ext,ny_ext,nz_ext)   ! Cloud H2O mixing ratio (kg/kg)
  REAL,    INTENT(OUT) :: qr_ext (nx_ext,ny_ext,nz_ext)   ! Rain  H2O mixing ratio (kg/kg)
  REAL,    INTENT(OUT) :: qi_ext (nx_ext,ny_ext,nz_ext)   ! Ice   mixing ratio (kg/kg)
  REAL,    INTENT(OUT) :: qs_ext (nx_ext,ny_ext,nz_ext)   ! Snow  mixing ratio (kg/kg)
  REAL,    INTENT(OUT) :: qh_ext (nx_ext,ny_ext,nz_ext)   ! Hail  mixing ratio (kg/kg)

  REAL,    INTENT(OUT) :: tsfc_ext  (nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: tdeep_ext  (nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: wetsfc_ext  (nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: wetdp_ext  (nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: wetcanp_ext(nx_ext,ny_ext)      ! Canopy water amount
  REAL,    INTENT(OUT) :: snowdpth_ext(nx_ext,ny_ext)     ! Snow depth (m)

  REAL,    INTENT(OUT) :: trn_ext  (nx_ext,ny_ext)      ! External terrain (m)
  REAL,    INTENT(OUT) :: psfc_ext (nx_ext,ny_ext)      ! Surface pressure (Pa)

  REAL,    INTENT(OUT) :: t_2m_ext  (nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: qv_2m_ext (nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: u_10m_ext (nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: v_10m_ext (nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: rain_ext  (nx_ext,ny_ext)

  INTEGER, INTENT(OUT) :: istatus
!
!-----------------------------------------------------------------------
!
!  Work arrays for storing grib data
!
!-----------------------------------------------------------------------
!
  REAL,    ALLOCATABLE :: var_grb2d(:,:,:)   ! GRIB 2D variables
  REAL,    ALLOCATABLE :: var_grb3d(:,:,:,:) ! GRIB 3-D variables
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,kk
  REAL    :: a, b, tvc_ext
  REAL    :: rovcp_p

!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  GRIB grid information
!
!-----------------------------------------------------------------------
!
  INTEGER :: ni_grb       ! Number of points along x-axis
  INTEGER :: nj_grb       ! Number of points along y-axis

  REAL :: di_grb          ! x-direction increment or grid length
  REAL :: dj_grb          ! y-direction increment or grid length

  REAL :: latsw           ! Latitude  of South West corner point
  REAL :: lonsw           ! Longitude of South West corner point
!  REAL :: latne           ! Latitude  of North East corner point
!  REAL :: lonne           ! Longitude of North East corner point

  REAL :: lattru1         ! Latitude (1st) at which projection is true
  REAL :: lattru2         ! Latitude (2nd) at which projection is true
  REAL :: lontrue         ! Longitude      at which projection is true

  INTEGER :: uvearth

!-----------------------------------------------------------------------
!
! Definitions for grid #236 data
!
!-----------------------------------------------------------------------
  INTEGER, PARAMETER :: n3dvs = 14

  INTEGER, PARAMETER :: PRE = 1        ! Pressure
  INTEGER, PARAMETER :: GPH = 2        ! Geopotential Height
  INTEGER, PARAMETER :: UWD = 3        ! U-component of Wind
  INTEGER, PARAMETER :: VWD = 4        ! V-component of Wind
  INTEGER, PARAMETER :: VPT = 5        ! Virtual Potential Temperature
  INTEGER, PARAMETER :: SQV = 6        ! Specific Humidity
  INTEGER, PARAMETER :: WWD = 7        ! Pressure vertical velocity(Pa/s)
  INTEGER, PARAMETER :: QCM = 8        ! Cloud Mixing Ratio
  INTEGER, PARAMETER :: QIM = 9        ! Ice Water Mixing Ratio
  INTEGER, PARAMETER :: QRM = 10       ! Rain Mixing Ratio
  INTEGER, PARAMETER :: QSM = 11       ! Snow Mixing Ratio
  INTEGER, PARAMETER :: QGM = 12       ! graupel Mixing Ratio
  INTEGER, PARAMETER :: TSL = 13       ! Soil temperature
  INTEGER, PARAMETER :: QSL = 14       ! Soil Moisture

  INTEGER, PARAMETER :: n2dvs = 14

  INTEGER, PARAMETER :: SSFC = 1       ! Snow depth at ground (m)
  INTEGER, PARAMETER :: RLRG = 2       ! Large-Scale Precipitation (kg m-2)
  INTEGER, PARAMETER :: RCVT = 3       ! Convective Precipitation (kg m-2)
  INTEGER, PARAMETER :: LAND = 4       ! Land Cover (1=land, 0=sea)
  INTEGER, PARAMETER :: ICEC = 5       ! Ice Cover (Proportion)
  INTEGER, PARAMETER :: SOIL = 6       ! Soil Type
  INTEGER, PARAMETER :: CANP = 7       ! Plant Canopy Surface Water (kg m-2)
  INTEGER, PARAMETER :: TMPS = 8       ! Soil Temperature (K)
  INTEGER, PARAMETER :: WETS = 9       ! Volumetric Soil Moisture Content(Fraction)

  INTEGER, PARAMETER :: T2M  = 10      ! Temperature at 2m above ground
  INTEGER, PARAMETER :: QV2M = 11      ! Specific Humidity at 2m above ground
  INTEGER, PARAMETER :: U10M = 12      ! U-Component of Wind at 10m above ground
  INTEGER, PARAMETER :: V10M = 13      ! V-Component of Wind at 10m above ground

  INTEGER, PARAMETER :: PMSL = 14      ! Pressure Reduced to MSL

  INTEGER, PARAMETER :: totvs = n2dvs + n3dvs
  INTEGER, PARAMETER :: varids(4,totvs) = RESHAPE(         & ! order is IMPORTANT
                   ! Discipline  Category  Parameter Layer/Level_Type
                        (/   0,   3,   0, 105,    & ! 1 PRE
                             0,   3,   5, 105,    & ! 2 GPH
                             0,   2,   2, 105,    & ! 3 UWD
                             0,   2,   3, 105,    & ! 4 VWD
                             0,   0,  15, 105,    & ! 5 VPT
                             0,   1,   2, 105,    & ! 6 SQV
                             0,   2,   8, 105,    & ! 7 WWD
                             0,   1,  22, 105,    & ! 8 QCM
                             0,   1,  23, 105,    & ! 9 QIM
                             0,   1,  24, 105,    & ! 10 QRM
                             0,   1,  25, 105,    & ! 11 QSM
                             0,   1,  32, 105,    & ! 12 QGM
                             2,   0,   2, 106,    & ! 13 TSL
                             2,   0, 192, 106,    & ! 14 QSL

                             0,   1,  11,   1,    & ! 1  - SSFC
                             0,   1,   9,   1,    & ! 2  - RLRG
                             0,   1,  10,   1,    & ! 3  - RCVT
                             2,   0,   0,   1,    & ! 4  - LAND
                            10,   2,   0,   1,    & ! 5  - ICEC
                             2,   3,   0,   1,    & ! 6  - SOIL
                             2,   0, 196,   1,    & ! 7  - CANP
                             2,   0,   2,   1,    & ! 8  - TMPS
                             2,   0, 192,   1,    & ! 9  - WETS
                             0,   0,   0, 103,    & ! 10 - T2M
                             0,   1,   0, 103,    & ! 11 - QV2M
                             0,   2,   2, 103,    & ! 12 - U10M
                             0,   2,   3, 103,    & ! 13 - U10M
                             0,   3,   1, 101     & ! 14 - PMSL
                             /), (/ 4, totvs /) )

  REAL,   PARAMETER :: var3dlvl(50) = (/ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, &
                    11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, &
                    26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, &
                    41, 42, 43, 44, 45, 46, 47, 48, 49, 50 /)

  INTEGER, PARAMETER :: var3dindx(totvs) = (/ PRE, GPH, UWD, VWD, VPT, SQV, WWD,&
                    QCM, QIM, QRM, QSM, QGM, TSL, QSL, (0, i = n3dvs+1,totvs) /)

  INTEGER, PARAMETER :: var2dindx(totvs) = (/ (0, i = 1, n3dvs),        &
                  SSFC, RLRG, RCVT, LAND, ICEC, SOIL, CANP, TMPS, WETS, &
                  T2M, QV2M, U10M, V10M, PMSL/)

  REAL,    PARAMETER :: var2dlvl(n2dvs) = (/ 0.,0.,0.,0.,0.,0.,0.,0.,0.,&
                                             2.,2.,10.,10.,0. /)
!
!
! This table converts WRF soil to those used in ARPS
!
! WRF uses 16 soil categories (FAO/STATSGO)
! ARPS uses 13-category soil  (Standard GRIB or GRIB2 definitions)
!
  INTEGER, PARAMETER :: soil_table(17) = (/ 1, 2, 3, 4, 4,              &
                                            5, 6, 7, 8, 9,              &
                                           10,11, 6,13, 1,              &
                                            2, 2/)

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
!  RDNMCGRB2 reads NMC GRIB2 data
!
!-----------------------------------------------------------------------
!
  ALLOCATE(var_grb2d(nx_ext,ny_ext,n2dvs),        STAT = istatus)
  ALLOCATE(var_grb3d(nx_ext,ny_ext,nz_ext,n3dvs), STAT = istatus)
  var_grb2d = -999.0
  var_grb3d = 0.0

  CALL rdnmcgrb2(nx_ext,ny_ext,nz_ext,gribfile,grbflen, gribtime,       &
           1, nx_ext, 1, ny_ext,                                        &
           iproj_ext,ni_grb,nj_grb,di_grb,dj_grb,                       &
           latsw,lonsw,lattru1,lattru2,lontrue,uvearth,                 &
           n2dvs, n3dvs, totvs, nzsoilin_ext,                           &
           varids, var2dindx, var3dindx, var2dlvl, var3dlvl,            &
           100*soildepth_ext,                                           &
           var_grb2d, var_grb3d, lvldbg,                                &
           istatus)

  IF (istatus /= 0) RETURN

  scale_ext     = 1.0
  latnot_ext(1) = lattru1
  latnot_ext(2) = lattru2
  trlon_ext     = lontrue

  dx_ext = di_grb
  dy_ext = dj_grb

  IF (myproc == 0) PRINT *,'LatSW = ',latsw,' LonSW = ',lonsw

  latsw_ext = latsw
  lonsw_ext = lonsw
!
!-----------------------------------------------------------------------
!
!  Retrieve 2-D variables
!
!-----------------------------------------------------------------------
!

  WRITE(6,*) 'Filling 2D arrays from GRIB2 data.'
  DO j=1,ny_ext
    DO i=1,nx_ext
      psfc_ext   (i,j) = -999.0
      trn_ext    (i,j) = -999.0
      snowdpth_ext(i,j) = var_grb2d(i,j,SSFC)
      wetcanp_ext(i,j) = var_grb2d(i,j,CANP)*1.0E-3

      t_2m_ext (i,j) = var_grb2d(i,j,T2M)
      qv_2m_ext(i,j) = var_grb2d(i,j,QV2M)
      u_10m_ext(i,j) = var_grb2d(i,j,U10M)
      v_10m_ext(i,j) = var_grb2d(i,j,V10M)
      rain_ext (i,j) = var_grb2d(i,j,RLRG) + var_grb2d(i,j,RCVT)

      !
      ! To make sure SOILTYPE and LANDCOVER consistent and
      ! convert SOILTYPE from to ARPS.
      !
!      IF (var_grb2d(i,j,LAND) < 0.5) THEN      ! Water
!        soiltyp_ext(i,j)= 13
!      ELSE IF (var_grb2d(i,j,ICEC) > 0.5) THEN ! Ice
!        soiltyp_ext(i,j) = 12
!
!!-----------------------------------------------------------------------
!!
!!     NAM grib2 files have soil type in them.  SREF doesn't.  Fall back
!!     to the grib1 method of just assigning 0 to the soil type.
!!
!!-----------------------------------------------------------------------
!!
!!      ELSE IF (var_grb2d(i,j,SOIL) == -999.0) THEN
!!        soiltyp_ext(i,j) = 0
!!      ELSE IF (var_grb2d(i,j,SOIL) < 1. .OR. var_grb2d(i,j,SOIL) > 16.) THEN
!!        WRITE(6,'(1x,2(a,I4),a,I2,a,F5.0)')       &
!!                'WARNING: Soil type is unknown. var_grb2d(',i,',',j,    &
!!                ',SOIL=',SOIL,') = ',var_grb2d(i,j,SOIL)
!!      ELSE
!!        soiltyp_ext(i,j) = SOIL_TABLE(NINT(var_grb2d(i,j,SOIL)))
!      ELSE    ! I don't know the definition of soil type for grid #236
!      	soiltyp_ext(i,j) = 0
!      END IF

    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Retrieve 3-D variables
!
!-----------------------------------------------------------------------
!
  WRITE(6,*) 'Filling 3D arrays from GRIB2 data.'
  rovcp_p = 0.285714    ! rd/cp used in RUC

  DO k=1,nz_ext
    DO j=1,ny_ext
      DO i=1,nx_ext
        p_ext  (i,j,k) = var_grb3d(i,j,k,PRE)    ! Pressure
        hgt_ext(i,j,k) = var_grb3d(i,j,k,GPH)
        u_ext  (i,j,k) = var_grb3d(i,j,k,UWD)    ! u wind (m/s)
        v_ext  (i,j,k) = var_grb3d(i,j,k,VWD)    ! v wind (m/s)

        a = REAL(100000)/var_grb3d(i,j,k,PRE)
        a = a**rovcp_p
        tvc_ext = var_grb3d(i,j,k,VPT)/a ! Virtual Temperature
        b = 0.61*var_grb3d(i,j,k,SQV)
        b = REAL(1) + b
        t_ext(i,j,k) = tvc_ext/b ! Temperature (K)

        a = var_grb3d(i,j,k,SQV)*var_grb3d(i,j,k,PRE)
        a = a/(0.622 - var_grb3d(i,j,k,SQV))
        qv_ext(i,j,k) = a*0.622/var_grb3d(i,j,k,PRE) ! SpecificHumidity

        qc_ext (i,j,k) = var_grb3d(i,j,k,QCM)
        qr_ext (i,j,k) = var_grb3d(i,j,k,QRM)
        qi_ext (i,j,k) = var_grb3d(i,j,k,QIM)
        qs_ext (i,j,k) = var_grb3d(i,j,k,QSM)
        qh_ext (i,j,k) = var_grb3d(i,j,k,QGM)
      END DO
    END DO
  END DO

  WRITE(6,*) 'Filling 3D Soil arrays from GRIB2 data.'

!  tsoil_ext (:,:,:) = -999.0
!  qsoil_ext (:,:,:) = -999.0

  IF(soilmodel_option == 1) THEN       ! ARPS: two-layer force-restore model)

    DO j = 1, ny_ext
      DO i = 1, nx_ext

        tsfc_ext   (i,j) = var_grb2d(i,j,TMPS)
        tdeep_ext  (i,j) = 0.1 * var_grb3d(i,j,2,TSL)  & !   5cm
                         + 0.2 * var_grb3d(i,j,3,TSL)  & !  20cm
                         + 0.4 * var_grb3d(i,j,4,TSL)  & !  40cm
                         + 0.3 * var_grb3d(i,j,5,TSL)    ! 160cm

        wetsfc_ext (i,j) = var_grb2d(i,j,WETS)
        wetdp_ext  (i,j) = 0.1 * var_grb3d(i,j,2,QSL)  & !   5cm
                         + 0.2 * var_grb3d(i,j,3,QSL)  & !  20cm
                         + 0.4 * var_grb3d(i,j,4,QSL)  & !  40cm
                         + 0.3 * var_grb3d(i,j,5,QSL)    ! 160cm

!        tsoil_ext(i,j,1) = 0.5*( var_grb3d(i,j,1,TSL) + var_grb3d(i,j,2,TSL) )      ! sfc temp.
!        tsoil_ext(i,j,2) = 0.5*( var_grb3d(i,j,3,TSL) + var_grb3d(i,j,4,TSL) )
!        qsoil_ext(i,j,1) = 0.5*( var_grb3d(i,j,1,QSL) + var_grb3d(i,j,2,QSL) )      ! sfc moisture
!        qsoil_ext(i,j,2) = 0.5*( var_grb3d(i,j,3,QSL) + var_grb3d(i,j,4,QSL) )
!       IF (var_grb2d(i,j,LAND) < .5 ) THEN  ! Over water?
!         tsoil_ext(i,j,1) = var_grb2d(i,j,TSFC)
!         tsoil_ext(i,j,2) = var_grb2d(i,j,TSFC)
!         qsoil_ext(i,j,1) = 1.0
!         qsoil_ext(i,j,2) = 1.0
!       END IF
      END DO
    END DO

  ELSE                      ! ARPS: multi-layer OUSoil model
    WRITE(6,'(1x,a)') 'To be implemented.'
    CALL arpsstop('Soilmodel_option not implemented.',1)
!    DO j = 1, ny_ext
!      DO i = 1, nx_ext
!        tsoil_ext(i,j,1) = var_grb2d(i,j,TSFC)
!        qsoil_ext(i,j,1) = var_grb3d(i,j,1,QSL)
!      END DO
!    END DO
!
!    DO k = 2,nzsoil_ext
!      DO j = 1, ny_ext
!        DO i = 1, nx_ext
!          tsoil_ext(i,j,k) = var_grb3d(i,j,k-1,TSL)
!          qsoil_ext(i,j,k) = var_grb3d(i,j,k-1,QSL)
!          IF (var_grb2d(i,j,LAND) < .5) THEN  ! over water?
!            tsoil_ext(i,j,k) = var_grb2d(i,j,TSFC)
!            qsoil_ext(i,j,k) = 1.0
!          END IF
!        END DO
!      END DO
!    END DO
!
  END IF

  DEALLOCATE(var_grb2d,var_grb3d, STAT = istatus)

  RETURN
END SUBROUTINE getruc236n_grb2
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE GETRUC130n_grb2                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE getruc130n_grb2(nx_ext,ny_ext,nz_ext,nzsoil_ext,nzsoilin_ext,&
           gribfile,gribtime,grbflen,soildepth_ext,                    &
           dx_ext,dy_ext,                                              &
           iproj_ext,scale_ext,trlon_ext,latnot_ext,latsw_ext,lonsw_ext,&
           p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                     &
           qc_ext,qr_ext,qi_ext,qs_ext,qh_ext,                         &
           tsoil_ext,qsoil_ext,wetcanp_ext,                            &
           snowdpth_ext,trn_ext,psfc_ext,t_2m_ext,qv_2m_ext,           &
           u_10m_ext,v_10m_ext,rain_ext,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads RUC 13km GRIB2 data (grid #130)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster, based on getruc236n_grb2
!  06/29/2011
!
!  MODIFICATION HISTORY:
!
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx_ext
!    ny_ext
!    nz_ext
!    nzsoil_ext
!    nstyps_ext
!    gribfile      File name to be read
!    gribtime
!
!  OUTPUT:
!
!    dx_ext
!    dy_ext
!    iproj_ext     Map projection number of external data
!    scale_ext     Scale factor of external data
!    trlon_ext     True longitude of external data (degrees E)
!    latnot_ext(2) True latitude(s) of external data (degrees N)
!    lawsw_ext     Latitude of the southwest corn of the domain
!    lonsw_ext     Longitude of the southwest corn of the domain
!    hgt_ext       height (m)
!    p_ext         pressure (Pascal)
!    t_ext         temperature (K)
!    qv_ext        specific humidity (kg/kg)
!    u_ext         u wind component (m/s)
!    v_ext         v wind component (m/s)
!    qc_ext        Cloud water mixing ratio (kg/kg)
!    qr_ext        Rain water mixing ratio (kg/kg)
!    qi_ext        Ice mixing ratio (kg/kg)
!    qs_ext        Snow mixing ratio (kg/kg)
!    qh_ext        Hail mixing ratio (kg/kg)
!
!    tsoil_ext     Soil temperature
!    qsoil_ext     Soil moisture
!    wetcanp_ext   Water content on canopy
!
!    trn_ext       External terrain (m)
!    psfc_ext      Surface pressure (Pa)
!
!    T_2m_ext      Temperature at 2m AGL
!    qv_2m_ext     Specific Humidity at 2m AGL
!    U_10m_ext     U at 10m AGL
!    V_10m_ext     V at 10m AGL
!    rain_ext      Accumulated rainfall contained in the data file
!
!    istatus       status indicator
!
!  WORK ARRAYS:
!
!    var_grb3d     Arrays to store the GRIB 3-D variables:
!    var_grb2d     Arrays to store the GRIB 2-D variables:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  CHARACTER (LEN=*), INTENT(IN) :: gribfile
  CHARACTER (LEN=*), INTENT(IN) :: gribtime
  INTEGER, INTENT(IN) :: grbflen

  INTEGER, INTENT(IN) :: nx_ext,ny_ext,nz_ext,nzsoil_ext

  INTEGER, INTENT(IN) :: nzsoilin_ext
  REAL,    INTENT(IN) :: soildepth_ext(nzsoilin_ext)
!
!-----------------------------------------------------------------------
!
!  External grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER, INTENT(OUT) :: iproj_ext
  REAL,    INTENT(OUT) :: scale_ext,trlon_ext
  REAL,    INTENT(OUT) :: latnot_ext(2)
  REAL,    INTENT(OUT) :: dx_ext,dy_ext
  REAL,    INTENT(OUT) :: latsw_ext, lonsw_ext
!
!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------
!
  REAL,    INTENT(OUT) :: p_ext  (nx_ext,ny_ext,nz_ext)   ! Pressure (Pascals)
  REAL,    INTENT(OUT) :: hgt_ext(nx_ext,ny_ext,nz_ext)   ! Height (m)
  REAL,    INTENT(OUT) :: t_ext  (nx_ext,ny_ext,nz_ext)   ! Temperature (K)
  REAL,    INTENT(OUT) :: qv_ext (nx_ext,ny_ext,nz_ext)   ! Specific humidity (kg/kg)
  REAL,    INTENT(OUT) :: u_ext  (nx_ext,ny_ext,nz_ext)   ! Eastward wind component
  REAL,    INTENT(OUT) :: v_ext  (nx_ext,ny_ext,nz_ext)   ! Northward wind component
  REAL,    INTENT(OUT) :: qc_ext (nx_ext,ny_ext,nz_ext)   ! Cloud H2O mixing ratio (kg/kg)
  REAL,    INTENT(OUT) :: qr_ext (nx_ext,ny_ext,nz_ext)   ! Rain  H2O mixing ratio (kg/kg)
  REAL,    INTENT(OUT) :: qi_ext (nx_ext,ny_ext,nz_ext)   ! Ice   mixing ratio (kg/kg)
  REAL,    INTENT(OUT) :: qs_ext (nx_ext,ny_ext,nz_ext)   ! Snow  mixing ratio (kg/kg)
  REAL,    INTENT(OUT) :: qh_ext (nx_ext,ny_ext,nz_ext)   ! Hail  mixing ratio (kg/kg)

  REAL,    INTENT(OUT) :: tsoil_ext  (nx_ext,ny_ext,nzsoil_ext)
  REAL,    INTENT(OUT) :: qsoil_ext  (nx_ext,ny_ext,nzsoil_ext)
  REAL,    INTENT(OUT) :: wetcanp_ext(nx_ext,ny_ext)      ! Canopy water amount
  REAL,    INTENT(OUT) :: snowdpth_ext(nx_ext,ny_ext)     ! Snow depth (m)

  REAL,    INTENT(OUT) :: trn_ext  (nx_ext,ny_ext)      ! External terrain (m)
  REAL,    INTENT(OUT) :: psfc_ext (nx_ext,ny_ext)      ! Surface pressure (Pa)

  REAL,    INTENT(OUT) :: t_2m_ext  (nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: qv_2m_ext (nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: u_10m_ext (nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: v_10m_ext (nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: rain_ext  (nx_ext,ny_ext)

  INTEGER, INTENT(OUT) :: istatus
!
!-----------------------------------------------------------------------
!
!  Work arrays for storing grib data
!
!-----------------------------------------------------------------------
!
  REAL,    ALLOCATABLE :: var_grb2d(:,:,:)   ! GRIB 2D variables
  REAL,    ALLOCATABLE :: var_grb3d(:,:,:,:) ! GRIB 3-D variables
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,kk
  REAL    :: a, b, tvc_ext
  REAL    :: rovcp_p

!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  GRIB grid information
!
!-----------------------------------------------------------------------
!
  INTEGER :: ni_grb       ! Number of points along x-axis
  INTEGER :: nj_grb       ! Number of points along y-axis

  REAL :: di_grb          ! x-direction increment or grid length
  REAL :: dj_grb          ! y-direction increment or grid length

  REAL :: latsw           ! Latitude  of South West corner point
  REAL :: lonsw           ! Longitude of South West corner point

  REAL :: lattru1         ! Latitude (1st) at which projection is true
  REAL :: lattru2         ! Latitude (2nd) at which projection is true
  REAL :: lontrue         ! Longitude      at which projection is true

  INTEGER :: uvearth

!-----------------------------------------------------------------------
!
! Definitions for grid #130 data
!
!-----------------------------------------------------------------------
  INTEGER, PARAMETER :: n3dvs = 14

  INTEGER, PARAMETER :: PRE = 1        ! Pressure
  INTEGER, PARAMETER :: GPH = 2        ! Geopotential Height
  INTEGER, PARAMETER :: UWD = 3        ! U-component of Wind
  INTEGER, PARAMETER :: VWD = 4        ! V-component of Wind
  INTEGER, PARAMETER :: TMP = 5        ! Virtual Potential Temperature
  INTEGER, PARAMETER :: SQV = 6        ! Specific Humidity
  INTEGER, PARAMETER :: WWD = 7        ! Pressure vertical velocity(Pa/s)
  INTEGER, PARAMETER :: QCM = 8        ! Cloud Mixing Ratio
  INTEGER, PARAMETER :: QIM = 9        ! Ice Water Mixing Ratio
  INTEGER, PARAMETER :: QRM = 10       ! Rain Mixing Ratio
  INTEGER, PARAMETER :: QSM = 11       ! Snow Mixing Ratio
  INTEGER, PARAMETER :: QGM = 12       ! graupel Mixing Ratio
  INTEGER, PARAMETER :: TSL = 13       ! Soil temperature
  INTEGER, PARAMETER :: QSL = 14       ! Soil Moisture

  INTEGER, PARAMETER :: n2dvs = 14

  INTEGER, PARAMETER :: HSFC = 1       ! Surface height (m)
  INTEGER, PARAMETER :: RLRG = 2       ! Large-Scale Precipitation (kg m-2)
  INTEGER, PARAMETER :: RCVT = 3       ! Convective Precipitation (kg m-2)
  INTEGER, PARAMETER :: SSFC = 4       ! Snow Cover (m)
  INTEGER, PARAMETER :: ICEC = 5       ! Ice Cover (Proportion)
  INTEGER, PARAMETER :: SOIL = 6       ! Soil Type
  INTEGER, PARAMETER :: CANP = 7       ! Plant Canopy Surface Water (kg m-2)
  INTEGER, PARAMETER :: TMPS = 8       ! Soil Temperature (K)
  INTEGER, PARAMETER :: WETS = 9       ! Volumetric Soil Moisture Content(Fraction)

  INTEGER, PARAMETER :: T2M  = 10      ! Temperature at 2m above ground
  INTEGER, PARAMETER :: QV2M = 11      ! Specific Humidity at 2m above ground
  INTEGER, PARAMETER :: U10M = 12      ! U-Component of Wind at 10m above ground
  INTEGER, PARAMETER :: V10M = 13      ! V-Component of Wind at 10m above ground

  INTEGER, PARAMETER :: PSFC = 14      ! Surface pressure (Pa)

  INTEGER, PARAMETER :: totvs = n2dvs + n3dvs
  INTEGER, PARAMETER :: varids(4,totvs) = RESHAPE(         & ! order is IMPORTANT
                   ! Discipline  Category  Parameter Layer/Level_Type
                        (/   0,   3,   0, 104,    & ! 1 PRE
                             0,   3,   5, 104,    & ! 2 GPH
                             0,   2,   2, 104,    & ! 3 UWD
                             0,   2,   3, 104,    & ! 4 VWD
                             0,   0,   0, 104,    & ! 5 TMP
                             0,   1,   0, 104,    & ! 6 SPFH
                             0,   2,   8, 104,    & ! 7 VVEL
                             0,   1,  22, 104,    & ! 8 QCM
                             0,   6,   0, 104,    & ! 9 QIM
                             0,   1,  24, 104,    & ! 10 QRM
                             0,   1,  25, 104,    & ! 11 QSM
                             0,   1,  32, 104,    & ! 12 QGM
                             2,   0,   2, 106,    & ! 13 TSL
                             2,   0, 192, 106,    & ! 14 QSL

                             0,   3,   5,   1,    & ! 1  - HSFC
                             0,   1,   9,   1,    & ! 2  - RLRG
                             0,   1,  10,   1,    & ! 3  - RCVT
                             0,   1, 201,   1,    & ! 4  - SSFC
                            10,   2,   0,   1,    & ! 5  - ICEC
                             2,   3,   0,   1,    & ! 6  - SOIL
                             2,   0, 196,   1,    & ! 7  - CANP
                             2,   0,   2,   1,    & ! 8  - TMPS
                             2,   0, 192,   1,    & ! 9  - WETS
                             0,   0,   0, 103,    & ! 10 - T2M
                             0,   1,   0, 103,    & ! 11 - QV2M
                             0,   2,   2, 103,    & ! 12 - U10M
                             0,   2,   3, 103,    & ! 13 - U10M
                             0,   3,   0,   1     & ! 14 - PSFC
                             /), (/ 4, totvs /) )

  REAL,   PARAMETER :: var3dlvl(50) = (/ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, &
                    11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, &
                    26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, &
                    41, 42, 43, 44, 45, 46, 47, 48, 49, 50 /)

  INTEGER, PARAMETER :: var3dindx(totvs) = (/ PRE, GPH, UWD, VWD, TMP, SQV, WWD,&
                    QCM, QIM, QRM, QSM, QGM, TSL, QSL, (0, i = n3dvs+1,totvs) /)

  INTEGER, PARAMETER :: var2dindx(totvs) = (/ (0, i = 1, n3dvs),        &
                  HSFC, RLRG, RCVT, SSFC, ICEC, SOIL, CANP, TMPS, WETS, &
                  T2M, QV2M, U10M, V10M, PSFC/)

  REAL,    PARAMETER :: var2dlvl(n2dvs) = (/ 0.,0.,0.,0.,0.,0.,0.,0.,0.,&
                                             2.,2.,10.,10.,0. /)
!
! This table converts WRF soil to those used in ARPS
!
! WRF uses 16 soil categories (FAO/STATSGO)
! ARPS uses 13-category soil  (Standard GRIB or GRIB2 definitions)
!
  INTEGER, PARAMETER :: soil_table(17) = (/ 1, 2, 3, 4, 4,              &
                                            5, 6, 7, 8, 9,              &
                                           10,11, 6,13, 1,              &
                                            2, 2/)

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
!  RDNMCGRB2 reads NMC GRIB2 data
!
!-----------------------------------------------------------------------
!
  ALLOCATE(var_grb2d(nx_ext,ny_ext,n2dvs),        STAT = istatus)
  ALLOCATE(var_grb3d(nx_ext,ny_ext,nz_ext,n3dvs), STAT = istatus)
  var_grb2d = -999.0
  var_grb3d = 0.0

  CALL rdnmcgrb2(nx_ext,ny_ext,nz_ext,gribfile,grbflen, gribtime,       &
           1, nx_ext, 1, ny_ext,                                        &
           iproj_ext,ni_grb,nj_grb,di_grb,dj_grb,                       &
           latsw,lonsw,lattru1,lattru2,lontrue,uvearth,                 &
           n2dvs, n3dvs, totvs, nzsoilin_ext,                           &
           varids, var2dindx, var3dindx, var2dlvl, var3dlvl,            &
           100*soildepth_ext,                                           &
           var_grb2d, var_grb3d, lvldbg,                                &
           istatus)

  IF (istatus /= 0) RETURN

  IF (myproc == 0) PRINT *,'LatTrue1 = ',lattru1,' LatTrue2 = ',lattru2
  IF (myproc == 0) PRINT *,'LonTrue  = ',lontrue
  scale_ext     = 1.0
  latnot_ext(1) = lattru1
  latnot_ext(2) = lattru2
  trlon_ext     = lontrue

  IF (myproc == 0) PRINT *,'DX = ',di_grb,' DY= ',dj_grb
  dx_ext = di_grb
  dy_ext = dj_grb

  IF (myproc == 0) PRINT *,'LatSW = ',latsw,' LonSW = ',lonsw

  latsw_ext = latsw
  lonsw_ext = lonsw
!
!-----------------------------------------------------------------------
!
!  Retrieve 2-D variables
!
!-----------------------------------------------------------------------
!

  WRITE(6,*) 'Filling 2D arrays from GRIB2 data.'
  DO j=1,ny_ext
    DO i=1,nx_ext
      psfc_ext   (i,j) = var_grb2d(i,j,PSFC)
      trn_ext    (i,j) = var_grb2d(i,j,HSFC)
      snowdpth_ext(i,j) = var_grb2d(i,j,SSFC)

      t_2m_ext (i,j) = var_grb2d(i,j,T2M)
      qv_2m_ext(i,j) = var_grb2d(i,j,QV2M)
      u_10m_ext(i,j) = var_grb2d(i,j,U10M)
      v_10m_ext(i,j) = var_grb2d(i,j,V10M)
      rain_ext (i,j) = var_grb2d(i,j,RLRG) + var_grb2d(i,j,RCVT)

      !
      ! To make sure SOILTYPE and LANDCOVER consistent and
      ! convert SOILTYPE from to ARPS.
      !
!      IF (var_grb2d(i,j,LAND) < 0.5) THEN      ! Water
!        soiltyp_ext(i,j)= 13
!      ELSE IF (var_grb2d(i,j,ICEC) > 0.5) THEN ! Ice
!        soiltyp_ext(i,j) = 12
!
!!-----------------------------------------------------------------------
!!
!!     NAM grib2 files have soil type in them.  SREF doesn't.  Fall back
!!     to the grib1 method of just assigning 0 to the soil type.
!!
!!-----------------------------------------------------------------------
!!
!!      ELSE IF (var_grb2d(i,j,SOIL) == -999.0) THEN
!!        soiltyp_ext(i,j) = 0
!!      ELSE IF (var_grb2d(i,j,SOIL) < 1. .OR. var_grb2d(i,j,SOIL) > 16.) THEN
!!        WRITE(6,'(1x,2(a,I4),a,I2,a,F5.0)')       &
!!                'WARNING: Soil type is unknown. var_grb2d(',i,',',j,    &
!!                ',SOIL=',SOIL,') = ',var_grb2d(i,j,SOIL)
!!      ELSE
!!        soiltyp_ext(i,j) = SOIL_TABLE(NINT(var_grb2d(i,j,SOIL)))
!      ELSE    ! I don't know the definition of soil type for grid #130
!      	soiltyp_ext(i,j) = 0
!      END IF

    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Retrieve 3-D variables
!
!-----------------------------------------------------------------------
!
  WRITE(6,*) 'Filling 3D arrays from GRIB2 data.'
  rovcp_p = 0.285714    ! rd/cp used in RUC

  DO k=1,nz_ext
    DO j=1,ny_ext
      DO i=1,nx_ext
        p_ext  (i,j,k) = var_grb3d(i,j,k,PRE)    ! Pressure
        hgt_ext(i,j,k) = var_grb3d(i,j,k,GPH)
        u_ext  (i,j,k) = var_grb3d(i,j,k,UWD)    ! u wind (m/s)
        v_ext  (i,j,k) = var_grb3d(i,j,k,VWD)    ! v wind (m/s)
        t_ext(i,j,k) = var_grb3d(i,j,k,TMP)      ! Temperature (K)
        qv_ext(i,j,k) = var_grb3d(i,j,k,SQV)     ! SpecificHumidity
        qc_ext (i,j,k) = var_grb3d(i,j,k,QCM)
        qr_ext (i,j,k) = var_grb3d(i,j,k,QRM)
        qi_ext (i,j,k) = var_grb3d(i,j,k,QIM)
        qs_ext (i,j,k) = var_grb3d(i,j,k,QSM)
        qh_ext (i,j,k) = var_grb3d(i,j,k,QGM)
      END DO
    END DO
  END DO

  WRITE(6,*) 'Filling 3D Soil arrays from GRIB2 data.'

!  tsoil_ext (:,:,:) = -999.0
!  qsoil_ext (:,:,:) = -999.0

  IF(soilmodel_option == 1) THEN       ! ARPS: two-layer force-restore model)

    DO j = 1, ny_ext
      DO i = 1, nx_ext

        tsoil_ext(i,j,1) = var_grb2d(i,j,TMPS)
        tsoil_ext(i,j,2) = 0.1 * var_grb3d(i,j,2,TSL)  & !   5cm
                         + 0.2 * var_grb3d(i,j,3,TSL)  & !  20cm
                         + 0.4 * var_grb3d(i,j,4,TSL)  & !  40cm
                         + 0.3 * var_grb3d(i,j,5,TSL)    ! 160cm

        qsoil_ext(i,j,1) = var_grb2d(i,j,WETS)
        qsoil_ext(i,j,2) = 0.1 * var_grb3d(i,j,2,QSL)  & !   5cm
                         + 0.2 * var_grb3d(i,j,3,QSL)  & !  20cm
                         + 0.4 * var_grb3d(i,j,4,QSL)  & !  40cm
                         + 0.3 * var_grb3d(i,j,5,QSL)    ! 160cm

!        tsoil_ext(i,j,1) = 0.5*( var_grb3d(i,j,1,TSL) + var_grb3d(i,j,2,TSL) )      ! sfc temp.
!        tsoil_ext(i,j,2) = 0.5*( var_grb3d(i,j,3,TSL) + var_grb3d(i,j,4,TSL) )
!        qsoil_ext(i,j,1) = 0.5*( var_grb3d(i,j,1,QSL) + var_grb3d(i,j,2,QSL) )      ! sfc moisture
!        qsoil_ext(i,j,2) = 0.5*( var_grb3d(i,j,3,QSL) + var_grb3d(i,j,4,QSL) )
!       IF (var_grb2d(i,j,LAND) < .5 ) THEN  ! Over water?
!         tsoil_ext(i,j,1) = var_grb2d(i,j,TSFC)
!         tsoil_ext(i,j,2) = var_grb2d(i,j,TSFC)
!         qsoil_ext(i,j,1) = 1.0
!         qsoil_ext(i,j,2) = 1.0
!       END IF
      END DO
    END DO

  ELSE                      ! ARPS: multi-layer OUSoil model
    WRITE(6,'(1x,a)') 'To be implemented.'
    CALL arpsstop('Soilmodel_option not implemented.',1)
!    DO j = 1, ny_ext
!      DO i = 1, nx_ext
!        tsoil_ext(i,j,1) = var_grb2d(i,j,TSFC)
!        qsoil_ext(i,j,1) = var_grb3d(i,j,1,QSL)
!      END DO
!    END DO
!
!    DO k = 2,nzsoil_ext
!      DO j = 1, ny_ext
!        DO i = 1, nx_ext
!          tsoil_ext(i,j,k) = var_grb3d(i,j,k-1,TSL)
!          qsoil_ext(i,j,k) = var_grb3d(i,j,k-1,QSL)
!          IF (var_grb2d(i,j,LAND) < .5) THEN  ! over water?
!            tsoil_ext(i,j,k) = var_grb2d(i,j,TSFC)
!            qsoil_ext(i,j,k) = 1.0
!          END IF
!        END DO
!      END DO
!    END DO
!
  END IF

  DEALLOCATE(var_grb2d,var_grb3d, STAT = istatus)

  RETURN
END SUBROUTINE getruc130n_grb2
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE GETMAM236p_grb2                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE getruc236p_grb2(nx_ext,ny_ext,nz_ext,nzsoilin_ext,           &
           gribfile,gribtime,grbflen,soildepth_ext,                     &
           dx_ext,dy_ext,                                               &
           iproj_ext,scale_ext,trlon_ext,latnot_ext,latsw_ext,lonsw_ext, &
           p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                      &
           tsfc_ext,tdeep_ext,wetsfc_ext,wetdp_ext,wetcanp_ext,         &
           snowdpth_ext,trn_ext,psfc_ext,t_2m_ext,qv_2m_ext,            &
           u_10m_ext,v_10m_ext,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads RUC2 isobaric GRIB2 data.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  04/07/2009
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx_ext
!    ny_ext
!    nz_ext
!    nzsoil_ext
!    nstyps_ext
!    gribfile      File name to be read
!    gribtime
!
!  OUTPUT:
!
!    dx_ext
!    dy_ext
!    iproj_ext     Map projection number of external data
!    scale_ext     Scale factor of external data
!    trlon_ext     True longitude of external data (degrees E)
!    latnot_ext(2) True latitude(s) of external data (degrees N)
!    lawsw_ext     Latitude of the southwest corn of the domain
!    lonsw_ext     Longitude of the southwest corn of the domain
!    hgt_ext       height (m)
!    p_ext         pressure (Pascal)
!    t_ext         temperature (K)
!    qv_ext        specific humidity (kg/kg)
!    u_ext         u wind component (m/s)
!    v_ext         v wind component (m/s)
!
!    tsoil_ext     Soil temperature
!    qsoil_ext     Soil moisture
!    wetcanp_ext   Water content on canopy
!
!    trn_ext       External terrain (m)
!    psfc_ext      Surface pressure (Pa)
!
!    T_2m_ext      Temperature at 2m AGL
!    qv_2m_ext     Specific Humidity at 2m AGL
!    U_10m_ext     U at 10m AGL
!    V_10m_ext     V at 10m AGL
!    rain_ext      Accumulated rainfall contained in the data file
!
!    istatus       status indicator
!
!  WORK ARRAYS:
!
!    var_grb3d     Arrays to store the GRIB 3-D variables:
!    var_grb2d     Arrays to store the GRIB 2-D variables:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  CHARACTER (LEN=*), INTENT(IN) :: gribfile
  CHARACTER (LEN=*), INTENT(IN) :: gribtime
  INTEGER, INTENT(IN) :: grbflen

  INTEGER, INTENT(IN) :: nx_ext,ny_ext,nz_ext

  INTEGER, INTENT(IN) :: nzsoilin_ext
  REAL,    INTENT(IN) :: soildepth_ext(nzsoilin_ext)
!
!-----------------------------------------------------------------------
!
!  External grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER, INTENT(OUT) :: iproj_ext
  REAL,    INTENT(OUT) :: scale_ext,trlon_ext
  REAL,    INTENT(OUT) :: latnot_ext(2)
  REAL,    INTENT(OUT) :: dx_ext,dy_ext
  REAL,    INTENT(OUT) :: latsw_ext, lonsw_ext
!
!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------
!
  REAL,    INTENT(OUT) :: p_ext  (nx_ext,ny_ext,nz_ext)   ! Pressure (Pascals)
  REAL,    INTENT(OUT) :: hgt_ext(nx_ext,ny_ext,nz_ext)   ! Height (m)
  REAL,    INTENT(OUT) :: t_ext  (nx_ext,ny_ext,nz_ext)   ! Temperature (K)
  REAL,    INTENT(OUT) :: qv_ext (nx_ext,ny_ext,nz_ext)   ! Specific humidity (kg/kg)
  REAL,    INTENT(OUT) :: u_ext  (nx_ext,ny_ext,nz_ext)   ! Eastward wind component
  REAL,    INTENT(OUT) :: v_ext  (nx_ext,ny_ext,nz_ext)   ! Northward wind component

  REAL,    INTENT(OUT) :: tsfc_ext  (nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: tdeep_ext  (nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: wetsfc_ext  (nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: wetdp_ext  (nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: wetcanp_ext(nx_ext,ny_ext)      ! Canopy water amount
  REAL,    INTENT(OUT) :: snowdpth_ext(nx_ext,ny_ext)     ! Snow depth (m)

  REAL,    INTENT(OUT) :: trn_ext  (nx_ext,ny_ext)      ! External terrain (m)
  REAL,    INTENT(OUT) :: psfc_ext (nx_ext,ny_ext)      ! Surface pressure (Pa)
  REAL,    INTENT(OUT) :: t_2m_ext (nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: qv_2m_ext (nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: u_10m_ext (nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: v_10m_ext (nx_ext,ny_ext)

  INTEGER, INTENT(OUT) :: istatus
!
!-----------------------------------------------------------------------
!
!  Work arrays for storing grib data
!
!-----------------------------------------------------------------------
!
  REAL,    ALLOCATABLE :: var_grb2d(:,:,:)   ! GRIB 2D variables
  REAL,    ALLOCATABLE :: var_grb3d(:,:,:,:) ! GRIB 3-D variables
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  REAL    :: qvsat

!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  GRIB grid information
!
!-----------------------------------------------------------------------
!
  INTEGER :: ni_grb       ! Number of points along x-axis
  INTEGER :: nj_grb       ! Number of points along y-axis

  REAL :: di_grb          ! x-direction increment or grid length
  REAL :: dj_grb          ! y-direction increment or grid length

  REAL :: latsw           ! Latitude  of South West corner point
  REAL :: lonsw           ! Longitude of South West corner point
!  REAL :: latne           ! Latitude  of North East corner point
!  REAL :: lonne           ! Longitude of North East corner point

  REAL :: lattru1         ! Latitude (1st) at which projection is true
  REAL :: lattru2         ! Latitude (2nd) at which projection is true
  REAL :: lontrue         ! Longitude      at which projection is true

  INTEGER :: uvearth

!-----------------------------------------------------------------------
!
! Definitions for grid #236 data
!
!-----------------------------------------------------------------------
  INTEGER, PARAMETER :: n3dvs = 6

  INTEGER, PARAMETER :: GPH = 1        ! Geopotential Height (gpm)
  INTEGER, PARAMETER :: TMP = 2        ! Temperature (K)
  INTEGER, PARAMETER :: UWD = 3        ! U-component of Wind (m s-1)
  INTEGER, PARAMETER :: VWD = 4        ! V-component of Wind (m s-1)
  INTEGER, PARAMETER :: RHD = 5        ! Relative Humidity (%)
  INTEGER, PARAMETER :: WWD = 6        ! Vertical Velocity(Pa/s)

  INTEGER, PARAMETER :: n2dvs = 10

  INTEGER, PARAMETER :: PSFC = 1       ! Pressrue (Pa)
  INTEGER, PARAMETER :: SGPH = 2       ! Geopotential Height (gpm)
  INTEGER, PARAMETER :: SSFC = 3       ! Snow depth at ground (m)
  INTEGER, PARAMETER :: RLRG = 4       ! Large-Scale Precipitation (kg m-2)
  INTEGER, PARAMETER :: RCVT = 5       ! Convective Precipitation (kg m-2)

  INTEGER, PARAMETER :: T2M  = 6       ! Temperature at 2m above ground
  INTEGER, PARAMETER :: RH2M = 7       ! Relative Humidity at 2m above ground
  INTEGER, PARAMETER :: U10M = 8       ! U-Component of Wind at 10m above ground
  INTEGER, PARAMETER :: V10M = 9       ! V-Component of Wind at 10m above ground
  INTEGER, PARAMETER :: QV2M = 10      ! Specific Humidity at 2m above ground

  INTEGER, PARAMETER :: totvs = n2dvs + n3dvs
  INTEGER, PARAMETER :: varids(4,totvs) = RESHAPE(         & ! order is IMPORTANT
                   ! Discipline  Category  Parameter Layer/Level_Type
                        (/   0,   3,   5, 100,    & ! 1 GPH
                             0,   0,   0, 100,    & ! 2 TMP
                             0,   2,   2, 100,    & ! 3 UWD
                             0,   2,   3, 100,    & ! 4 VWD
                             0,   1,   1, 100,    & ! 5 RHD
                             0,   2,   8, 100,    & ! 6 WWD

                             0,   3,   0,   1,    & ! 1 - PSFC
                             0,   3,   5,   1,    & ! 2 - SGPH
                             0,   1,  11,   1,    & ! 3 - SSFC
                             0,   1,   9,   1,    & ! 4 - RLRG
                             0,   1,  10,   1,    & ! 5 - RCVT
                             0,   0,   0, 103,    & ! 6 - T2M
                             0,   1,   1, 103,    & ! 7 - RH2M
                             0,   2,   2, 103,    & ! 8 - U10M
                             0,   2,   3, 103,    & ! 9 - U10M
                             0,   1,   0, 103     & ! 10 - QV2M
                             /), (/ 4, totvs /) )

  REAL,   PARAMETER :: var3dlvl(37) = (/ 100000., 97500., 95000., 92500., 90000., 87500., 85000.,  &
                  82500., 80000., 77500., 75000., 72500., 70000., 67500., 65000., 62500., 60000.,  &
                  57500., 55000., 52500., 50000., 47500., 45000., 42500., 40000., 37500., 35000.,  &
                  32500., 30000., 27500., 25000., 22500., 20000., 17500., 15000., 12500., 10000. /)

  INTEGER, PARAMETER :: var3dindx(totvs) = (/ GPH, TMP, UWD, VWD, RHD, WWD, &
                                              (0, i = n3dvs+1,totvs) /)

  INTEGER, PARAMETER :: var2dindx(totvs) = (/ (0, i = 1, n3dvs),            &
                      PSFC, SGPH, SSFC, RLRG, RCVT, T2M, RH2M, U10M, V10M, QV2M /)

  REAL,    PARAMETER :: var2dlvl(n2dvs) = (/ 0.,0.,0.,0.,0.,2.,2.,10.,10.,2. /)

  REAL :: f_qvsatl
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
!  RDNMCGRB2 reads NMC GRIB2 data
!
!-----------------------------------------------------------------------
!
  ALLOCATE(var_grb2d(nx_ext,ny_ext,n2dvs),        STAT = istatus)
  ALLOCATE(var_grb3d(nx_ext,ny_ext,nz_ext,n3dvs), STAT = istatus)
  var_grb2d = -999.0
  var_grb3d = 0.0

  CALL rdnmcgrb2(nx_ext,ny_ext,nz_ext,gribfile,grbflen, gribtime,       &
           1, nx_ext, 1, ny_ext,                                        &
           iproj_ext,ni_grb,nj_grb,di_grb,dj_grb,                       &
           latsw,lonsw,lattru1,lattru2,lontrue,uvearth,                 &
           n2dvs, n3dvs, totvs, nzsoilin_ext,                           &
           varids, var2dindx, var3dindx, var2dlvl, var3dlvl,            &
           100*soildepth_ext,                                           &
           var_grb2d, var_grb3d, lvldbg,                                &
           istatus)

  IF (istatus /= 0) RETURN

  scale_ext     = 1.0
  latnot_ext(1) = lattru1
  latnot_ext(2) = lattru2
  trlon_ext     = lontrue

  dx_ext = di_grb
  dy_ext = dj_grb

  IF (myproc == 0) PRINT *,'LatSW = ',latsw,' LonSW = ',lonsw

  latsw_ext = latsw
  lonsw_ext = lonsw
!
!-----------------------------------------------------------------------
!
!  Retrieve 2-D variables
!
!-----------------------------------------------------------------------
!

  WRITE(6,*) 'Filling 2D arrays from GRIB2 data.'
  DO j=1,ny_ext
    DO i=1,nx_ext
      psfc_ext   (i,j)  = var_grb2d(i,j,PSFC)
      trn_ext    (i,j)  = var_grb2d(i,j,SGPH)
      snowdpth_ext(i,j) = var_grb2d(i,j,SSFC)

      t_2m_ext (i,j) = var_grb2d(i,j,T2M)
      qv_2m_ext(i,j) = var_grb2d(i,j,QV2M)
      u_10m_ext(i,j) = var_grb2d(i,j,U10M)
      v_10m_ext(i,j) = var_grb2d(i,j,V10M)

      wetcanp_ext(i,j) = -999.0

    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Retrieve 3-D variables
!
!-----------------------------------------------------------------------
!
  WRITE(6,*) 'Filling 3D arrays from GRIB2 data.'

  DO k=1,nz_ext
    DO j=1,ny_ext
      DO i=1,nx_ext
        p_ext  (i,j,k) = var3dlvl(k)             ! Pressure
        t_ext  (i,j,k) = var_grb3d(i,j,k,TMP)
        hgt_ext(i,j,k) = var_grb3d(i,j,k,GPH)
        u_ext  (i,j,k) = var_grb3d(i,j,k,UWD)    ! u wind (m/s)
        v_ext  (i,j,k) = var_grb3d(i,j,k,VWD)    ! v wind (m/s)

        IF( (p_ext(i,j,k) > 0.0) .AND. (t_ext(i,j,k) > 0.0) ) THEN
          qvsat = f_qvsatl( p_ext(i,j,k), t_ext(i,j,k) )
          qv_ext(i,j,k) = var_grb3d(i,j,k,RHD) * qvsat * 0.01

          IF((k == 1) .OR. (p_ext(i,j,k) > psfc_ext(i,j)))THEN
            qv_ext(i,j,k)= var_grb2d(i,j,RH2M) * qvsat * 0.01
          END IF
        ELSE
          qv_ext(i,j,k)= 0.0
        END IF

      END DO
    END DO
  END DO

  WRITE(6,*) 'Filling 3D Soil arrays from GRIB2 data.'

!  tsoil_ext (:,:,:) = -999.0
!  qsoil_ext (:,:,:) = -999.0

  IF(soilmodel_option == 1) THEN       ! ARPS: two-layer force-restore model)

    DO j = 1, ny_ext
      DO i = 1, nx_ext

        tsfc_ext   (i,j) = -999.
        tdeep_ext  (i,j) = -999.

        wetsfc_ext (i,j) = -999.
        wetdp_ext  (i,j) = -999.

      END DO
    END DO

  ELSE                      ! ARPS: multi-layer OUSoil model
    WRITE(6,'(1x,a)') 'To be implemented.'
    CALL arpsstop('Soilmodel_option not implemented.',1)
  END IF

  DEALLOCATE(var_grb2d,var_grb3d, STAT = istatus)

  RETURN
END SUBROUTINE getruc236p_grb2
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE GETMAM212_grb2                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE getnam132_grb2(nx_ext,ny_ext,nz_ext,nzsoil_ext,nzsoilin_ext, &
           gribfile,gribtime,grbflen,soildepth_ext,                     &
           dx_ext,dy_ext,iboxs,iboxe,jboxs,jboxe,                       &
           iproj_ext,scale_ext,trlon_ext,latnot_ext,latsw_ext,lonsw_ext, &
           lat_ext,lon_ext,p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,      &
           qc_ext,                                                      &
           tsoil_ext,qsoil_ext,wetcanp_ext,                             &
           snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext,                   &
           t_2m_ext, qv_2m_ext, u_10m_ext, v_10m_ext, rain_ext,         &
           istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads NAM 16km GRIB2 data (grid #132)
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  01/24/2013
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx_ext
!    ny_ext
!    nz_ext
!    nzsoil_ext
!    nstyps_ext
!    gribfile      File name to be read
!    gribtime
!
!  OUTPUT:
!
!    dx_ext
!    dy_ext
!    iproj_ext     Map projection number of external data
!    scale_ext     Scale factor of external data
!    trlon_ext     True longitude of external data (degrees E)
!    latnot_ext(2) True latitude(s) of external data (degrees N)
!    lawsw_ext     Latitude of the southwest corn of the domain
!    lonsw_ext     Longitude of the southwest corn of the domain
!    hgt_ext       height (m)
!    p_ext         pressure (Pascal)
!    t_ext         temperature (K)
!    qv_ext        specific humidity (kg/kg)
!    u_ext         u wind component (m/s)
!    v_ext         v wind component (m/s)
!    qc_ext        Cloud water mixing ratio (kg/kg)
!    qr_ext        Rain water mixing ratio (kg/kg)
!    qi_ext        Ice mixing ratio (kg/kg)
!    qs_ext        Snow mixing ratio (kg/kg)
!    qh_ext        Hail mixing ratio (kg/kg)
!
!    tsoil_ext     Soil temperature
!    qsoil_ext     Soil moisture
!    wetcanp_ext   Water content on canopy
!
!    trn_ext       External terrain (m)
!    psfc_ext      Surface pressure (Pa)
!
!    T_2m_ext      Temperature at 2m AGL
!    qv_2m_ext     Specific Humidity at 2m AGL
!    U_10m_ext     U at 10m AGL
!    V_10m_ext     V at 10m AGL
!    rain_ext      Accumulated rainfall contained in the data file
!
!    istatus       status indicator
!
!  WORK ARRAYS:
!
!    var_grb3d     Arrays to store the GRIB 3-D variables:
!    var_grb2d     Arrays to store the GRIB 2-D variables:
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  CHARACTER (LEN=*), INTENT(IN) :: gribfile
  CHARACTER (LEN=*), INTENT(IN) :: gribtime
  INTEGER, INTENT(IN) :: grbflen

  INTEGER, INTENT(IN) :: nx_ext,ny_ext,nz_ext,nzsoil_ext

  INTEGER, INTENT(IN) :: nzsoilin_ext
  REAL,    INTENT(IN) :: soildepth_ext(nzsoilin_ext)

!  INTEGER, PARAMETER :: nzsoilin_ext = 4
!  REAL,    PARAMETER :: soildepth_ext(nzsoilin_ext) = (/0.0, 0.1, 0.4, 1.0/)

  INTEGER, INTENT(IN) :: iboxs, iboxe, jboxs, jboxe
!
!-----------------------------------------------------------------------
!
!  External grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER, INTENT(OUT) :: iproj_ext
  REAL,    INTENT(OUT) :: scale_ext,trlon_ext
  REAL,    INTENT(OUT) :: latnot_ext(2)
  REAL,    INTENT(OUT) :: dx_ext,dy_ext
  REAL,    INTENT(OUT) :: latsw_ext, lonsw_ext
!
!-----------------------------------------------------------------------
!
!  Output external variable arrays
!
!-----------------------------------------------------------------------
!
  REAL,    INTENT(OUT) :: lat_ext  (nx_ext,ny_ext)   ! Latitude (-90 to 90)
  REAL,    INTENT(OUT) :: lon_ext  (nx_ext,ny_ext)   ! East Longitude (0 to 360)

  REAL,    INTENT(OUT) :: p_ext  (nx_ext,ny_ext,nz_ext)   ! Pressure (Pascals)
  REAL,    INTENT(OUT) :: hgt_ext(nx_ext,ny_ext,nz_ext)   ! Height (m)
  REAL,    INTENT(OUT) :: t_ext  (nx_ext,ny_ext,nz_ext)   ! Temperature (K)
  REAL,    INTENT(OUT) :: qv_ext (nx_ext,ny_ext,nz_ext)   ! Specific humidity (kg/kg)
  REAL,    INTENT(OUT) :: u_ext  (nx_ext,ny_ext,nz_ext)   ! Eastward wind component
  REAL,    INTENT(OUT) :: v_ext  (nx_ext,ny_ext,nz_ext)   ! Northward wind component
  REAL,    INTENT(OUT) :: qc_ext (nx_ext,ny_ext,nz_ext)   ! Cloud H2O mixing ratio (kg/kg)
!  REAL,    INTENT(OUT) :: qr_ext (nx_ext,ny_ext,nz_ext)   ! Rain  H2O mixing ratio (kg/kg)
!  REAL,    INTENT(OUT) :: qi_ext (nx_ext,ny_ext,nz_ext)   ! Ice   mixing ratio (kg/kg)
!  REAL,    INTENT(OUT) :: qs_ext (nx_ext,ny_ext,nz_ext)   ! Snow  mixing ratio (kg/kg)
!  REAL,    INTENT(OUT) :: qh_ext (nx_ext,ny_ext,nz_ext)   ! Hail  mixing ratio (kg/kg)

  REAL,    INTENT(OUT) :: tsoil_ext  (nx_ext,ny_ext,nzsoil_ext)
  REAL,    INTENT(OUT) :: qsoil_ext  (nx_ext,ny_ext,nzsoil_ext)
  REAL,    INTENT(OUT) :: wetcanp_ext(nx_ext,ny_ext)      ! Canopy water amount
  REAL,    INTENT(OUT) :: snowdpth_ext(nx_ext,ny_ext)     ! Snow depth (m)

  REAL,    INTENT(OUT) :: trn_ext  (nx_ext,ny_ext)      ! External terrain (m)
  REAL,    INTENT(OUT) :: psfc_ext (nx_ext,ny_ext)      ! Surface pressure (Pa)
  REAL,    INTENT(OUT) :: t_2m_ext (nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: qv_2m_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: u_10m_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: v_10m_ext(nx_ext,ny_ext)
  REAL,    INTENT(OUT) :: rain_ext (nx_ext,ny_ext)

  INTEGER, INTENT(OUT) :: soiltyp_ext(nx_ext,ny_ext)

  INTEGER, INTENT(OUT) :: istatus
!
!-----------------------------------------------------------------------
!
!  Work arrays for storing grib data
!
!-----------------------------------------------------------------------
!
  REAL,    ALLOCATABLE :: var_grb2d(:,:,:)   ! GRIB 2D variables
  REAL,    ALLOCATABLE :: var_grb3d(:,:,:,:) ! GRIB 3-D variables

  REAL,    ALLOCATABLE :: qvs_ext(:,:,:)
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,kk

!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  GRIB grid information
!
!-----------------------------------------------------------------------
!
  INTEGER :: ni_grb       ! Number of points along x-axis
  INTEGER :: nj_grb       ! Number of points along y-axis

  REAL :: di_grb          ! x-direction increment or grid length
  REAL :: dj_grb          ! y-direction increment or grid length

  REAL :: latsw           ! Latitude  of South West corner point
  REAL :: lonsw           ! Longitude of South West corner point
!  REAL :: latne           ! Latitude  of North East corner point
!  REAL :: lonne           ! Longitude of North East corner point

  REAL :: lattru1         ! Latitude (1st) at which projection is true
  REAL :: lattru2         ! Latitude (2nd) at which projection is true
  REAL :: lontrue         ! Longitude      at which projection is true

  INTEGER :: uvearth

!-----------------------------------------------------------------------
!
! Definitions for grid #212 data
!
!-----------------------------------------------------------------------
  INTEGER, PARAMETER :: n3dvs = 10

  INTEGER, PARAMETER :: GPH = 1        ! Geopotential Height
  INTEGER, PARAMETER :: UWD = 2        ! U-component of Wind
  INTEGER, PARAMETER :: VWD = 3        ! V-component of Wind
  INTEGER, PARAMETER :: TMP = 4        ! Temperature
  INTEGER, PARAMETER :: SRH = 5        ! Specific Humidity
  INTEGER, PARAMETER :: WWD = 6        ! Pressure vertical velocity(Pa/s)
  INTEGER, PARAMETER :: QCM = 7        ! Cloud Mixing Ratio
  INTEGER, PARAMETER :: QIM = 8        ! Cloud Ice Ratio
  INTEGER, PARAMETER :: TSL = 9        ! Soil temperature
  INTEGER, PARAMETER :: QSL = 10       ! Soil Moisture

  INTEGER, PARAMETER :: n2dvs = 14

  INTEGER, PARAMETER :: TSFC = 1       ! Temperature at ground
  INTEGER, PARAMETER :: PSFC = 2       ! Pressure at ground
  INTEGER, PARAMETER :: GSFC = 3       ! Geopotential Height at ground
  INTEGER, PARAMETER :: SSFC = 4       ! Water equivalent of accumulated snow depth (kg m-2)
  INTEGER, PARAMETER :: RAIN = 5       ! Total Precipitation (kg m-2)
  INTEGER, PARAMETER :: LAND = 6       ! Land Cover (1=land, 0=sea)
  INTEGER, PARAMETER :: PCAN = 7       ! Plant Canopy Surface Water (kg m-2)
  INTEGER, PARAMETER :: WTMP = 8       ! Water Temperature (K)

  INTEGER, PARAMETER :: T2M  = 9       ! Temperature at 2m above ground
  INTEGER, PARAMETER :: QV2M = 10      ! Specific Humidity at 2m above ground
  INTEGER, PARAMETER :: U10M = 11      ! U-Component of Wind at 10m above ground
  INTEGER, PARAMETER :: V10M = 12      ! V-Component of Wind at 10m above ground

  INTEGER, PARAMETER :: LATI = 13       ! Latitude
  INTEGER, PARAMETER :: LONG = 14       ! East longitude

  INTEGER, PARAMETER :: totvs = n2dvs + n3dvs
  INTEGER, PARAMETER :: varids(4,totvs) = RESHAPE(         & ! order is IMPORTANT
                   ! Discipline  Category  Parameter Layer/Level_Type
                        (/   0,   3,   5, 100,    & !  1 GPH
                             0,   2,   2, 100,    & !  2 UWD
                             0,   2,   3, 100,    & !  3 VWD
                             0,   0,   0, 100,    & !  4 TMP
                             0,   1,   1, 100,    & !  5 SRH
                             0,   2,   8, 100,    & !  6 WWD
                             0,   1,  22, 100,    & !  7 QCM
                             0,   6,   0, 100,    & !  8 QIM
                             2,   0,   2, 106,    & !  9 TSL
                             2,   0, 192, 106,    & ! 10 QSL

                             0,   0,   0,   1,    & ! 1  - TSFC
                             0,   3,   0,   1,    & ! 2  - PSFC
                             0,   3,   5,   1,    & ! 3  - GSFC
                             0,   1,  13,   1,    & ! 4  - SSFC
                             0,   1,   8,   1,    & ! 5  - RAIN
                             2,   0,   0,   1,    & ! 6  - LAND
                             2,   0, 196,   1,    & ! 7  - PCAN
                            10,   3,   0,   1,    & ! 8  - WTMP
                             0,   0,   0, 103,    & ! 9  - T2M
                             0,   1,   0, 103,    & ! 10 - QV2M
                             0,   2,   2, 103,    & ! 11 - U10M
                             0,   2,   3, 103,    & ! 12 - U10M
                             0, 191, 192,   1,    & ! 13 - LATI
                             0, 191, 193,   1     & ! 14 - LONG
                             /), (/ 4, totvs /) )

  REAL,   PARAMETER :: var3dlvl(39) = (/ 100000, 97500, 95000, 92500, 90000, 87500, &
              85000, 82500, 80000, 77500, 75000, 72500, 70000, 67500, 65000, 62500, &
              60000, 57500, 55000, 52500, 50000, 47500, 45000, 42500, 40000, 37500, &
              35000, 32500, 30000, 27500, 25000, 22500, 20000, 17500, 15000, 12500, &
              10000,  7500,  5000 /)

  INTEGER, PARAMETER :: var3dindx(totvs) = (/ GPH, UWD, VWD, TMP, SRH,  &
                        WWD, QCM, QIM, TSL, QSL, (0, i = n3dvs+1,totvs) /)

  INTEGER, PARAMETER :: var2dindx(totvs) = (/ (0, i = 1, n3dvs),        &
                        TSFC, PSFC, GSFC, SSFC, RAIN, LAND, PCAN, WTMP, &
                        T2M, QV2M, U10M, V10M, LATI, LONG /)

  REAL,    PARAMETER :: var2dlvl(n2dvs) = (/ 0.,0.,0.,0.,0.,0.,0.,0.,   &
                                             2.,2.,10.,10., 0.,0. /)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ALLOCATE(var_grb2d(nx_ext,ny_ext,n2dvs),        STAT = istatus)
  ALLOCATE(var_grb3d(nx_ext,ny_ext,nz_ext,n3dvs), STAT = istatus)
  var_grb2d = -999.0
  var_grb3d = 0.0
!
!-----------------------------------------------------------------------
!
!  RDNMCGRB2 reads NMC GRIB2 data
!
!-----------------------------------------------------------------------
!
  CALL rdnmcgrb2(nx_ext,ny_ext,nz_ext,gribfile,grbflen, gribtime,       &
           iboxs,iboxe,jboxs,jboxe,                                     &
           iproj_ext,ni_grb,nj_grb,di_grb,dj_grb,                       &
           latsw,lonsw,lattru1,lattru2,lontrue,uvearth,                 &
           n2dvs, n3dvs, totvs, nzsoilin_ext,                           &
           varids, var2dindx, var3dindx, var2dlvl, var3dlvl,            &
           100*soildepth_ext,                                           &
           var_grb2d, var_grb3d, lvldbg,                                &
           istatus)

  IF (istatus /= 0) RETURN

  scale_ext     = 1.0
  latnot_ext(1) = lattru1
  latnot_ext(2) = lattru2
  trlon_ext     = lontrue

  dx_ext = di_grb
  dy_ext = dj_grb

  IF (myproc == 0) PRINT *,'LatSW = ',latsw,' LonSW = ',lonsw

  latsw_ext = latsw
  lonsw_ext = lonsw
!
!-----------------------------------------------------------------------
!
!  Retrieve 2-D variables
!
!-----------------------------------------------------------------------
!

  IF (myproc == 0 .AND. lvldbg > 0) WRITE(6,*) 'Filling 2D arrays from GRIB2 data.'
  DO j=1,ny_ext
    DO i=1,nx_ext
      psfc_ext   (i,j)  = var_grb2d(i,j,PSFC)
      trn_ext    (i,j)  = var_grb2d(i,j,GSFC)
      wetcanp_ext(i,j)  = var_grb2d(i,j,PCAN)*1.e-3     ! in meter

      ! Convert water equiv. of accum. snow depth (kg/m**2) to meters
      ! (where 1 meter liquid water is set equivqlent to 10 meters snow).
      !     0.01 = 10. (m snow/m liquid) / (1000 kg/m**3)
      snowdpth_ext(i,j) = 0.01 * var_grb2d(i,j,SSFC)

      t_2m_ext (i,j) = var_grb2d(i,j,T2M)
      qv_2m_ext(i,j) = var_grb2d(i,j,QV2M)
      u_10m_ext(i,j) = var_grb2d(i,j,U10M)
      v_10m_ext(i,j) = var_grb2d(i,j,V10M)

      rain_ext (i,j) = var_grb2d(i,j,RAIN)

      lat_ext(i,j)   = var_grb2d(i,j,LATI)
      lon_ext(i,j)   = var_grb2d(i,j,LONG)

      !
      ! To make sure SOILTYPE and LANDCOVER consistent and
      ! convert SOILTYPE from to ARPS.
      !
      IF (var_grb2d(i,j,LAND) < 0.5) THEN      ! Water
        soiltyp_ext(i,j)= 13
        !
        ! NAM grib2 files have soil type in them.  SREF doesn't.  Fall back
        ! to the grib1 method of just assigning 0 to the soil type.
        !
      ELSE
        soiltyp_ext(i,j) = 0
      END IF

    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  Retrieve 3-D variables
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0 .AND. lvldbg > 0) WRITE(6,*) 'Filling 3D arrays from GRIB2 data.'
  DO k=1,nz_ext
    DO j=1,ny_ext
      DO i=1,nx_ext
        p_ext  (i,j,k) = var3dlvl(k)             ! Pressure
        hgt_ext(i,j,k) = var_grb3d(i,j,k,GPH)
        u_ext  (i,j,k) = var_grb3d(i,j,k,UWD)    ! u wind (m/s)
        v_ext  (i,j,k) = var_grb3d(i,j,k,VWD)    ! v wind (m/s)
        t_ext  (i,j,k) = var_grb3d(i,j,k,TMP)    ! Temperature (K)
        qv_ext (i,j,k) = var_grb3d(i,j,k,SRH)    ! Relative humidity (%)
        qc_ext (i,j,k) = var_grb3d(i,j,k,QCM)
!        qr_ext (i,j,k) = -999.
!        qi_ext (i,j,k) = var_grb3d(i,j,k,QIM)
!        qs_ext (i,j,k) = -999.
!        qh_ext (i,j,k) = -999.
      END DO
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
! Convert relative humidity to specific humidity
!
!-----------------------------------------------------------------------
!
  ALLOCATE(qvs_ext(nx_ext,ny_ext,nz_ext), STAT = istatus)
  CALL check_alloc_status(istatus, "getnam132_grb2:qvs_ext")

  CALL getqvs(nx_ext,ny_ext,nz_ext, 1,nx_ext,1,ny_ext,1,nz_ext,         &
              p_ext, t_ext, qvs_ext )

  DO k=1,nz_ext
    DO j=1,ny_ext
      DO i=1,nx_ext
        qv_ext(i,j,k) = 0.01*qv_ext(i,j,k)*qvs_ext(i,j,k)
      END DO
    END DO
  END DO

  DEALLOCATE(qvs_ext)
!
!-----------------------------------------------------------------------
!
! Retrieve 3D soil arrays
!
!-----------------------------------------------------------------------

  IF (myproc == 0 .AND. lvldbg > 0) WRITE(6,*) 'Filling 3D Soil arrays from GRIB2 data.'

  tsoil_ext (:,:,:) = -999.0
  qsoil_ext (:,:,:) = -999.0

  IF(soilmodel_option == 1) THEN       ! ARPS: two-layer force-restore model)

    DO j = 1, ny_ext
      DO i = 1, nx_ext
        tsoil_ext(i,j,1) = 0.5*( var_grb3d(i,j,1,TSL) + var_grb3d(i,j,2,TSL) )      ! sfc temp.
        tsoil_ext(i,j,2) = 0.5*( var_grb3d(i,j,3,TSL) + var_grb3d(i,j,4,TSL) )
        qsoil_ext(i,j,1) = 0.5*( var_grb3d(i,j,1,QSL) + var_grb3d(i,j,2,QSL) )      ! sfc moisture
        qsoil_ext(i,j,2) = 0.5*( var_grb3d(i,j,3,QSL) + var_grb3d(i,j,4,QSL) )
        IF (var_grb2d(i,j,LAND) < .5 ) THEN  ! Over water?
          tsoil_ext(i,j,1) = var_grb2d(i,j,WTMP)
          tsoil_ext(i,j,2) = var_grb2d(i,j,WTMP)
          qsoil_ext(i,j,1) = 1.0
          qsoil_ext(i,j,2) = 1.0
        END IF
      END DO
    END DO

  ELSE                                ! ARPS: multi-layer OUSoil model

    DO j = 1, ny_ext
      DO i = 1, nx_ext
        tsoil_ext(i,j,1) = var_grb2d(i,j,TSFC)
        qsoil_ext(i,j,1) = var_grb3d(i,j,1,QSL)
        IF (var_grb2d(i,j,LAND) < .5) THEN  ! over water?
          tsoil_ext(i,j,1) = var_grb2d(i,j,WTMP)
        END IF
      END DO
    END DO

    DO k = 2,nzsoil_ext
      DO j = 1, ny_ext
        DO i = 1, nx_ext
          tsoil_ext(i,j,k) = var_grb3d(i,j,k-1,TSL)
          qsoil_ext(i,j,k) = var_grb3d(i,j,k-1,QSL)
          !IF (var_grb2d(i,j,LAND) < .5) THEN  ! over water?
          !  tsoil_ext(i,j,k) = var_grb2d(i,j,TSFC)
          !  qsoil_ext(i,j,k) = 1.0
          !END IF
        END DO
      END DO
    END DO

  END IF

  DEALLOCATE(var_grb2d,var_grb3d, STAT = istatus)

  RETURN
END SUBROUTINE getnam132_grb2

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE getgrbfname                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE getgrbfname(extdopt,extdfmt,dir_extd,extdname,               &
                       extdinit,extdfcst,appstr,                        &
                       gribfile,grbflen,gribtime,grbtlen,istatus)
!
!-----------------------------------------------------------------------
!
! PURPOSE:
!
!   Construct external data file name based on extdopt/extdfmt and
!   check its existence.
!
!-----------------------------------------------------------------------
!
! AUTHOR:
!   11/15/2007 (Yunheng Wang)
!   Modified based on common block from the origial subroutines
!   "getxxxx" etc above.
!
!-----------------------------------------------------------------------
!
! MODIFICATIONS:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  CHARACTER (LEN=*)   :: dir_extd
  CHARACTER (LEN=*)   :: extdname

  INTEGER, INTENT(IN) :: extdopt
  INTEGER, INTENT(IN) :: extdfmt
  CHARACTER(LEN=4), INTENT(IN) :: appstr

  CHARACTER(LEN=19),  INTENT(IN)  :: extdinit
  CHARACTER(LEN=9) ,  INTENT(IN)  :: extdfcst

  CHARACTER(LEN=256), INTENT(OUT) :: gribfile
  CHARACTER(LEN=14),  INTENT(OUT) :: gribtime

  INTEGER, INTENT(OUT) :: grbtlen
  INTEGER, INTENT(OUT) :: grbflen
  INTEGER, INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

  INTEGER :: lendir, lenext, fcstlen
  LOGICAL :: fexist

  INTEGER :: iyr,imo,iday,myr, jldy
  INTEGER :: ihr,imin,isec
  INTEGER :: ifhr,ifmin,ifsec

  CHARACTER(LEN=256) :: gribfile_new
  INTEGER            :: grbflen_new

  CHARACTER(LEN=32)  :: fmtstr

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  READ (extdinit,'(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)')                  &
                  iyr,  imo, iday,  ihr, imin, isec

  CALL julday(iyr,imo,iday,jldy)

  myr=MOD(iyr,100)
  ifhr=0
  ifmin=0
  ifsec=0

  IF (extdfcst(3:3) == ':') THEN
    READ(extdfcst,'(i2)') ifhr
  ELSE
    READ(extdfcst,'(i3)') ifhr
  END IF

  fcstlen = 2
  IF (ifhr > 99) fcstlen = 3

  WRITE(fmtstr,'(a,2(I1,a))') '(i4.4,i2.2,i2.2,i2.2,a1,i',fcstlen,'.',fcstlen,')'
  WRITE (gribtime,fmtstr) iyr, imo,iday, ihr,'f',ifhr
  grbtlen = 11 + fcstlen

!-----------------------------------------------------------------------
!
! File name and length
!
!-----------------------------------------------------------------------

  lendir = LEN_TRIM(dir_extd)

  IF( lendir == 0 .OR. dir_extd(1:lendir) == ' ' ) THEN
    dir_extd = './'
    lendir   = 2
  END IF

  IF( dir_extd(lendir:lendir) /= '/' ) THEN
    lendir = lendir + 1
    dir_extd(lendir:lendir) = '/'
  END IF

  IF (extdfmt == 0) THEN               ! NCEP name conventions

    WRITE(fmtstr,'(a,2(I1,a))') '(2a,a,I2.2,a,I',fcstlen,'.',fcstlen,',a)'

    SELECT CASE (extdopt)
    CASE (2)
      WRITE(gribfile,fmtstr) TRIM(dir_extd),                            &
             TRIM(extdname),'.t',ihr,'z.awip3d',ifhr,'.tm00'
    CASE (11)   ! Added by Eddie Natenberg for 13km RUC B grid
      IF (ifhr .gt. 0) THEN
           WRITE(gribfile,fmtstr) dir_extd(1:lendir),                   &
                  TRIM(extdname),'.t',ihr,'z.bgrb13f',ifhr
      ELSE
           WRITE(gribfile,fmtstr) dir_extd(1:lendir),                   &
                  TRIM(extdname),'.t',ihr,'z.bgrb13anl'
      END IF

    CASE (13)
      WRITE(gribfile,fmtstr) dir_extd(1:lendir),                        &
             TRIM(extdname),'.t',ihr,'Z.pgrbf',ifhr
    CASE (16)
      WRITE(gribfile,fmtstr) dir_extd(1:lendir),                        &
             TRIM(extdname),'.t',ihr,'z.awip218',ifhr
    CASE (17)
      WRITE(gribfile,'(3a,i4.4,2i2.2,a1,i2.2,a)') dir_extd(1:lendir),   &
             TRIM(extdname),'_',iyr,imo,iday,'_',ihr,'00_000.grb'
    CASE (20)
      WRITE(gribfile,fmtstr) dir_extd(1:lendir),                        &
             TRIM(extdname),'.t',ihr,'z.pgbf',ifhr
    CASE DEFAULT
      WRITE(6,'(1x,a,I4,a,/)')                                          &
        'ERROR: NCEP file name for extdopt = ',extdopt,                 &
              ' is still not implemented in GETGRBFNAME.'
      istatus = -2
      RETURN
    END SELECT

    grbflen = LEN_TRIM(gribfile)

  ELSE IF (extdfmt == 100) THEN   ! NOMADS data

    SELECT CASE (extdopt)
    CASE (13,16)
      WRITE(gribfile,'(3a,i4.4,2i2.2,a1,i2.2,a,i3.3,a)')                &
             dir_extd(1:lendir),TRIM(extdname),'_',                     &
             iyr,imo,iday,'_',ihr,'00_',ifhr,'.grb'
    CASE DEFAULT
      WRITE(6,'(1x,a,I4,a,/)')                                          &
        'ERROR: NOMADS file name for extdopt = ',extdopt,               &
              ' is still not implemented in GETGRBFNAME.'
      istatus = -2
      RETURN
    END SELECT

    grbflen = LEN_TRIM(gribfile)

  ELSE

    lenext = LEN_TRIM( extdname )

    gribfile = dir_extd(1:lendir)//extdname(1:lenext)                   &
                                 //'.'//gribtime(1:grbtlen)
    grbflen = lendir + lenext + grbtlen + 1

  END IF

  IF (LEN_TRIM(appstr) > 0) THEN
    gribfile_new = gribfile
    WRITE(gribfile,"(a,a)") TRIM(gribfile_new),TRIM(appstr)
    grbflen = grbflen + LEN_TRIM(appstr)
  END IF

!-----------------------------------------------------------------------
!
! Check file for existence
!
!-----------------------------------------------------------------------

  INQUIRE(FILE=gribfile(1:grbflen),EXIST=fexist)
  IF ( fexist ) RETURN

  IF ( ifhr < 100 .AND. MOD(extdfmt,100) > 0 ) THEN ! Check ARPS convention with 3 digit forecast time

    fcstlen = 3
    !WRITE(fmtstr,'(a,2(I0,a))') '(a,i4.4,i2.2,i2.2,i2.2,a1,i',fcstlen,'.',fcstlen,',a)'

    lenext = LEN_TRIM( extdname )

    WRITE( gribfile_new,'(a,i4.4,3i2.2,a1,i3.3,a)' )                    &
                         dir_extd(1:lendir)//extdname(1:lenext)//'.',   &
                         iyr, imo,iday, ihr,'f',ifhr,                   &
                         appstr
    grbflen_new = lendir + lenext + 11 + fcstlen + 1 + LEN_TRIM(appstr)

    INQUIRE(FILE=gribfile_new(1:grbflen_new), EXIST=fexist)
    IF (fexist) THEN
      grbtlen = 11 + fcstlen
      WRITE (gribtime,'(i4.4,i2.2,i2.2,i2.2,a1,i3.3)') iyr, imo,iday, ihr,'f',ifhr
      gribfile = gribfile_new
      grbflen  = grbflen_new
      RETURN
    END IF

  END IF

  gribfile_new = gribfile(1:grbflen)//'.grib2'
  grbflen_new  = grbflen + 6
  INQUIRE(FILE=gribfile_new(1:grbflen_new),EXIST=fexist)
  IF ( fexist) THEN
    gribfile = gribfile_new
    grbflen  = grbflen_new
    RETURN
  END IF

  INQUIRE(FILE=gribfile(1:grbflen)//'.Z', EXIST = fexist )
  IF( fexist ) THEN
    CALL unixcmd( 'cp '//trim(gribfile(1:grbflen))//'.Z tem_ext_file.Z' )
    INQUIRE(FILE='tem_ext_file', EXIST = fexist )
    IF( fexist ) call unixcmd( '/bin/rm tem_ext_file' )
    CALL uncmprs( 'tem_ext_file.Z' )
    gribfile = 'tem_ext_file'
    grbflen  = 12
    RETURN
  END IF

  INQUIRE(FILE=trim(gribfile(1:grbflen))//'.gz', EXIST = fexist )
  IF( fexist ) THEN
    CALL unixcmd( 'cp '//trim(gribfile(1:grbflen))//'.gz tem_ext_file.gz' )
    INQUIRE(FILE='tem_ext_file', EXIST = fexist )
    IF( fexist ) call unixcmd( '/bin/rm tem_ext_file' )
    CALL uncmprs( 'tem_ext_file.gz' )
    gribfile = 'tem_ext_file'
    grbflen  = 12
    RETURN
  END IF

  WRITE(6,'(/1x,a/,10x,a/)')                                            &
        'WARNING: GRIB file "'//gribfile(1:grbflen)//'" or its',        &
        'compressed version not found with GETGRBFNAME.'

  istatus = -888
  RETURN
END SUBROUTINE getgrbfname
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE getgrbdims                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getgrbdims(dir_extd,extdname,extdopt,extdfmt,                &
                      extdinit,extdfcst,ni_grb,nj_grb,iret)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine extracts dimension size from either GRIB file or GRIB2 file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  08/17/2007 Initialized.
!
!  MODIFICATIONS:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    dir_extd      GRIB file directory
!    extdname      Use to construct the GRIB file name
!    extdinit
!    extdfcst
!    extdopt
!    extdfmt
!
!  OUTPUT:
!
!    ni_grb        Grid size in west-east direction
!    nj_grb        Grid size in south-north direction
!
!    iret          Return flag
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  CHARACTER(LEN=*)              :: dir_extd
  CHARACTER(LEN=*), INTENT(IN)  :: extdname
  CHARACTER(LEN=19),INTENT(IN)  :: extdinit
  CHARACTER(LEN=9)              :: extdfcst

  INTEGER,          INTENT(IN)  :: extdopt
  INTEGER,          INTENT(IN)  :: extdfmt

  INTEGER,          INTENT(OUT) :: ni_grb       ! Number of points along x-axis
  INTEGER,          INTENT(OUT) :: nj_grb       ! Number of points along y-axis

  INTEGER,          INTENT(OUT) :: iret         ! Return flag

  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. Local variables
!
!-----------------------------------------------------------------------
!
  CHARACTER(LEN=256) :: gribfile
  CHARACTER(LEN=14)  :: gribtime

  INTEGER            :: grbflen, grbtlen, grbfmt
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  CALL getgrbfname(extdopt,extdfmt,dir_extd,extdname,                   &
                   extdinit,extdfcst,'    ',                            &
                   gribfile,grbflen,gribtime,grbtlen,iret)
  IF (iret /=0 ) RETURN

  IF (myproc == 0) THEN

    CALL chkgrb(gribfile(1:grbflen),grbflen,grbfmt,iret)

    IF (iret /=0 ) RETURN

    IF (grbfmt == 2) THEN
      CALL rdgrb2dims(gribfile,grbflen,ni_grb,nj_grb,iret)
    ELSE
      CALL rdgrbdims (gribfile,grbflen,ni_grb,nj_grb,iret)
    END IF

  END IF         ! IF (myproc == 0)

  CALL mpupdatei(ni_grb,1)
  CALL mpupdatei(nj_grb,1)

  CALL mpupdatei(iret, 1)

  RETURN
END SUBROUTINE getgrbdims
