!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETWRFDATA                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getwrfdata(nx_ext,ny_ext,nz_ext,nzsoil_ext,nstyps_ext,nscalar,  &
           dir_extd,extdname,extdopt,extdfmt,                           &
           extdinit,extdfcst,julfname,                                  &
           soilmodel_option,nzsoil_extin,                               &
           iproj_ext,scale_ext,                                         &
           trlon_ext,latnot_ext,x0_ext,y0_ext,                          &
           lat_ext,lon_ext,latu_ext,lonu_ext,latv_ext,lonv_ext,         &
           zpsoil_ext,                                                  &
           p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                      &
           qscalar_ext,qnames_wrf,                                      &
           tsoil_ext,qsoil_ext,wetcanp_ext,                             &
           snowdpth_ext,trn_ext,psfc_ext,soiltyp_ext,                   &
           t_2m_ext, qv_2m_ext, u_10m_ext, v_10m_ext, rain_ext,         &
           istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads a WRF ARW file in netCDF format for processing by
!  ext2arps, a program that converts external files to ARPS variables
!  and format.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  11/22/1911
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
!    T_2m_ext      Temperature at 2m AGL
!    qv_2m_ext     Specific Humidity at 2m AGL
!    U_10m_ext     U at 10m AGL
!    V_10m_ext     V at 10m AGL
!    rain_ext
!
!    istatus       status indicator
!
!  WORK ARRAYS:
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  CHARACTER (LEN=*) :: dir_extd
  CHARACTER (LEN=*) :: extdname

  INTEGER :: extdopt
  INTEGER :: extdfmt

  CHARACTER (LEN=19) :: extdinit
  CHARACTER (LEN=9)  :: extdfcst
  CHARACTER (LEN=9)  :: julfname

  INTEGER, INTENT(IN) :: soilmodel_option, nzsoil_extin

  INTEGER, INTENT(IN) :: nscalar
  CHARACTER(LEN=10) :: qnames_wrf(nscalar)
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
  REAL :: latu_ext(nx_ext,ny_ext)
  REAL :: lonu_ext(nx_ext,ny_ext)
  REAL :: latv_ext(nx_ext,ny_ext)
  REAL :: lonv_ext(nx_ext,ny_ext)

  REAL :: p_ext  (nx_ext,ny_ext,nz_ext)   ! Pressure (Pascals)
  REAL :: hgt_ext(nx_ext,ny_ext,nz_ext)   ! Height (m)
  REAL :: t_ext  (nx_ext,ny_ext,nz_ext)   ! Temperature (K)
  REAL :: qv_ext (nx_ext,ny_ext,nz_ext)   ! Specific humidity (kg/kg)
  REAL :: u_ext  (nx_ext,ny_ext,nz_ext)   ! Eastward wind component
  REAL :: v_ext  (nx_ext,ny_ext,nz_ext)   ! Northward wind component
  REAL :: qscalar_ext (nx_ext,ny_ext,nz_ext,nscalar)

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

  INTEGER :: istatus

!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  !INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Work arrays for storing grib data
!
!-----------------------------------------------------------------------
!
  REAL :: t0

  REAL,    ALLOCATABLE :: var1d(:)
  REAL,    ALLOCATABLE :: var2d(:)
  REAL,    ALLOCATABLE :: var3d(:)
  INTEGER, ALLOCATABLE :: var2di(:,:)

  REAL, ALLOCATABLE :: utmp(:,:), vtmp(:,:)

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
  INTEGER :: i,j,k,nq
  INTEGER :: grbflen, grbtlen
  INTEGER :: ncid

!
!-----------------------------------------------------------------------
!
!  GRIB grid information
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj_grb    ! Map projection indicator

  INTEGER :: ni_grb       ! Number of points along x-axis
  INTEGER :: nj_grb       ! Number of points along y-axis

  INTEGER :: nk_grb       ! Number of vertical parameters

  INTEGER :: nxs_ext, nys_ext, nzs_ext
  INTEGER :: itime

  REAL :: latsw           ! Latitude  of South West corner point
  REAL :: lonsw           ! Longitude of South West corner point
  REAL :: latctr, lonctr

  REAL :: lattru1         ! Latitude (1st) at which projection is true
  REAL :: lattru2         ! Latitude (2nd) at which projection is true
  REAL :: lontrue         ! Longitude      at which projection is true

  INTEGER :: sfcphys, mpphys

! WRF uses 16 soil categories and 24-category (USGS) vegetation
! ARPS uses 13-category soil  and 14-category (ND) vegetation
!
! The following two tables were provided by Jerry Brotzge on Oct. 20, 2003
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

  CALL getmapr(iproj,scale,latnot,trlon,x0,y0)

!-----------------------------------------------------------------------
!
! Reading data file
!
!-----------------------------------------------------------------------
  CALL get_wrf_fname(extdopt,extdfmt,dir_extd,extdname,                 &
                   extdinit,extdfcst,'    ',                            &
                   gribfile,grbflen,gribtime,grbtlen,istatus)
  IF (istatus /=0 ) RETURN


!-----------------------------------------------------------------------
!
! Open file and get metadata
!
!-----------------------------------------------------------------------

  IF (myproc == 0) THEN
    CALL open_ncd_file(gribfile,ncid)

    CALL get_ncd_dom_ti_integer(ncid,'WEST-EAST_GRID_DIMENSION',  ni_grb,istatus)
    CALL get_ncd_dom_ti_integer(ncid,'SOUTH-NORTH_GRID_DIMENSION',nj_grb,istatus)
    CALL get_ncd_dom_ti_integer(ncid,'BOTTOM-TOP_GRID_DIMENSION', nk_grb,istatus)

    CALL get_ncd_dom_ti_integer(ncid,'MP_PHYSICS',        mpphys, istatus)
    CALL get_ncd_dom_ti_integer(ncid,'SF_SURFACE_PHYSICS',sfcphys,istatus)

    CALL get_ncd_dom_ti_real   (ncid,'DX',  dx_ext,    istatus)
    CALL get_ncd_dom_ti_real   (ncid,'DY',  dy_ext,    istatus)

    CALL get_ncd_dom_ti_real   (ncid,'CEN_LAT',      latctr, istatus)
    CALL get_ncd_dom_ti_real   (ncid,'CEN_LON',      lonctr, istatus)

    CALL get_ncd_dom_ti_real   (ncid,'TRUELAT1',     lattru1, istatus)
    CALL get_ncd_dom_ti_real   (ncid,'TRUELAT2',     lattru2, istatus)
    CALL get_ncd_dom_ti_real   (ncid,'STAND_LON',    lontrue,  istatus)

    CALL get_ncd_dom_ti_integer(ncid,'MAP_PROJ',     iproj_grb,  istatus)

  END IF
  CALL mpupdatei(ni_grb,1)
  CALL mpupdatei(nj_grb,1)
  CALL mpupdatei(nk_grb,1)

  CALL mpupdatei(iproj_grb,1)
  CALL mpupdater(lattru1,1)
  CALL mpupdater(lattru2,1)
  CALL mpupdater(lontrue,1)

  CALL mpupdater(dx_ext,1)
  CALL mpupdater(dy_ext,1)

  IF (ni_grb /= nx_ext .OR. nj_grb /= ny_ext .OR. nk_grb /= nz_ext) THEN
    WRITE(*,'(1x,a,3I4,/,1x,a,3I4)')                                    &
           'Dimensions found in data file: ',ni_grb,nj_grb,nk_grb,      &
           '                     expected: ',nx_ext,ny_ext,nz_ext
    CALL arpsstop('ERROR: unmatched dimensions.',1)
  END IF
  nxs_ext = nx_ext - 1
  nys_ext = ny_ext - 1
  nzs_ext = nz_ext - 1

  IF(iproj_grb == 0) THEN        ! No projection
    iproj_ext = 0
  ELSE IF(iproj_grb == 1) THEN   ! LAMBERT CONFORMAL
    iproj_ext = 2
  ELSE IF(iproj_grb == 2) THEN   ! POLAR STEREOGRAPHIC
    iproj_ext = 1
  ELSE IF(iproj_grb == 3) THEN   ! MERCATOR
    iproj_ext = 3
  ELSE
    WRITE(6,*) 'Unknown map projection, ', iproj_grb
    istatus = -555
    CALL arpsstop('WRONG WRF map projection parameter.',1)
  END IF

  IF(lattru1 < 0.0) iproj_ext = -1*iproj_ext

  scale_ext = 1.0
  latnot_ext(1) = lattru1
  latnot_ext(2) = lattru2
  trlon_ext = lontrue

  itime = 1

  ALLOCATE(var2d (nx_ext*ny_ext),      STAT = istatus)

!-----------------------------------------------------------------------
!
! Get model latitude and Longitude
!
!-----------------------------------------------------------------------

  CALL get_wrf_2d(ncid,itime,1,(/'XLAT   '/),nxs_ext,nys_ext,var2d,     &
                  nx_ext,ny_ext,lat_ext,istatus)

  CALL get_wrf_2d(ncid,itime,1,(/'XLONG  '/),nxs_ext,nys_ext,var2d,     &
                  nx_ext,ny_ext,lon_ext,istatus)

  CALL get_wrf_2d(ncid,itime,1,(/'XLAT_U '/),nx_ext,nys_ext,var2d,      &
                  nx_ext,ny_ext,latu_ext,istatus)

  CALL get_wrf_2d(ncid,itime,1,(/'XLONG_U'/),nx_ext,nys_ext,var2d,      &
                  nx_ext,ny_ext,lonu_ext,istatus)

  CALL get_wrf_2d(ncid,itime,1,(/'XLAT_V '/),nxs_ext,ny_ext,var2d,      &
                  nx_ext,ny_ext,latv_ext,istatus)

  CALL get_wrf_2d(ncid,itime,1,(/'XLONG_V'/),nxs_ext,ny_ext,var2d,      &
                  nx_ext,ny_ext,lonv_ext,istatus)

  DO j = 1, ny_ext
    lon_ext (nx_ext,j) = lonu_ext(nx_ext,j) + lonu_ext(nx_ext,j) - lon_ext (nx_ext-1,j)
    lonv_ext(nx_ext,j) = lonu_ext(nx_ext,j) + lonu_ext(nx_ext,j) - lonv_ext(nx_ext-1,j)
  END DO

  DO i = 1, nx_ext
    lat_ext (i,ny_ext) = latv_ext(i,ny_ext) + latv_ext(i,ny_ext) - lat_ext (i,ny_ext-1)
    latu_ext(i,ny_ext) = latv_ext(i,ny_ext) + latv_ext(i,ny_ext) - latu_ext(i,ny_ext-1)
  END DO
!
!-----------------------------------------------------------------------
!
!  Retrieve 2-D variables one by one
!
!-----------------------------------------------------------------------
!
  CALL get_wrf_2d(ncid,itime,1,(/'PSFC   '/),nxs_ext,nys_ext,var2d,     &
                  nx_ext,ny_ext,psfc_ext,istatus)

  CALL get_wrf_2d(ncid,itime,1,(/'HGT    '/),nxs_ext,nys_ext,var2d,     &
                  nx_ext,ny_ext,trn_ext,istatus)

  CALL get_wrf_2d(ncid,itime,1,(/'CANWAT '/),nxs_ext,nys_ext,var2d,     &
                  nx_ext,ny_ext,wetcanp_ext,istatus)

  CALL get_wrf_2d(ncid,itime,1,(/'SNOWH  '/),nxs_ext,nys_ext,var2d,     &
                  nx_ext,ny_ext,snowdpth_ext,istatus)

  CALL get_wrf_2d(ncid,itime,1,(/'T2     '/),nxs_ext,nys_ext,var2d,     &
                  nx_ext,ny_ext,t_2m_ext,istatus)

  CALL get_wrf_2d(ncid,itime,1,(/'Q2     '/),nxs_ext,nys_ext,var2d,     &
                  nx_ext,ny_ext,qv_2m_ext,istatus)

  CALL get_wrf_2d(ncid,itime,1,(/'U10    '/),nxs_ext,nys_ext,var2d,     &
                  nx_ext,ny_ext,u_10m_ext,istatus)

  CALL get_wrf_2d(ncid,itime,1,(/'V10    '/),nxs_ext,nys_ext,var2d,     &
                  nx_ext,ny_ext,v_10m_ext,istatus)

  CALL get_wrf_2d(ncid,itime,2,(/'RAINC  ','RAINNC '/),nxs_ext,nys_ext,var2d,  &
                  nx_ext,ny_ext,rain_ext,istatus)

!
!-----------------------------------------------------------------------
!
!  Retrieve 3-D variables
!
!-----------------------------------------------------------------------
!
  ALLOCATE(var3d (nx_ext*ny_ext*nz_ext), STAT = istatus)

  t0 = 300.0

  CALL get_wrf_3d(ncid,itime,2,(/'P         ','PB        '/),           &
                  nxs_ext,nys_ext,nzs_ext,var3d,                        &
                  nx_ext,ny_ext,nz_ext,p_ext,istatus)

  CALL get_wrf_3d(ncid,itime,2,(/'PH        ','PHB       '/),           &
                  nxs_ext,nys_ext,nz_ext,var3d,                         &
                  nx_ext,ny_ext,nz_ext,hgt_ext,istatus)

  hgt_ext(:,:,:) = hgt_ext(:,:,:)/g  ! convert to height (m) from geopotential height (m2 s-2)

  CALL get_wrf_3d(ncid,itime,1,(/'T         '/),                        &
                  nxs_ext,nys_ext,nzs_ext,var3d,                        &
                  nx_ext,ny_ext,nz_ext,t_ext,istatus)

  !
  ! Convert potential temperature to temperature
  !
  t_ext(:,:,:) = (t_ext(:,:,:) + t0) * ((p_ext(:,:,:)/p0)**rddcp)


  CALL get_wrf_3d(ncid,itime,1,(/'U         '/),                        &
                  nx_ext,nys_ext,nzs_ext,var3d,                         &
                  nx_ext,ny_ext,nz_ext,u_ext,istatus)

  CALL get_wrf_3d(ncid,itime,1,(/'V         '/),                        &
                  nxs_ext,ny_ext,nzs_ext,var3d,                         &
                  nx_ext,ny_ext,nz_ext,v_ext,istatus)

  CALL get_wrf_3d(ncid,itime,1,(/'QVAPOR    '/),                        &
                  nxs_ext,nys_ext,nzs_ext,var3d,                        &
                  nx_ext,ny_ext,nz_ext,qv_ext,istatus)

  qscalar_ext(:,:,:,:) = -999.90

  DO nq = 1, nscalar
    CALL get_wrf_3d(ncid,itime,1,qnames_wrf(nq),                        &
                    nxs_ext,nys_ext,nzs_ext,var3d,                      &
                    nx_ext,ny_ext,nz_ext,qscalar_ext(:,:,:,nq),istatus)
  END DO

!---------------------------------------------------------------------
!
!  get soil types
!
!---------------------------------------------------------------------

  ALLOCATE(var2di(nxs_ext,nys_ext),      STAT = istatus)

  IF (myproc == 0) CALL get_ncd_2di(ncid,itime,'ISLTYP',nxs_ext,nys_ext,var2di,istatus)
  CALL mpupdatei(var2di,nxs_ext*nys_ext)

  DO j = 1, nys_ext
    DO i = 1,nxs_ext
      soiltyp_ext(i,j) = soil_table(var2di(i,j))  ! convert to ARPS categories
    END DO
  END DO
  CALL iedgfill(soiltyp_ext,1,nx_ext,1,nxs_ext,1,ny_ext,1,nys_ext,1,1,1,1)

  DEALLOCATE( var2di)

!---------------------------------------------------------------------
!
!  get soil temperature and moisture
!
!---------------------------------------------------------------------

  IF (soilmodel_option == 1) THEN    ! old ARPS Force-Restore Soil model

     CALL get_wrf_soil(ncid,itime,nxs_ext,nys_ext,nzsoil_extin,         &
            nx_ext,ny_ext,nzsoil_ext,zpsoil_ext,tsoil_ext,qsoil_ext,    &
            var3d,istatus)
  ELSE

    !Get soil depths

    ALLOCATE(var1d (nzsoil_extin), STAT = istatus)

!    IF (lvldbg > 25) WRITE(*,'(1x,a)') 'Reading variable ZS ...'
    IF (myproc == 0) CALL get_ncd_1d(ncid,itime,'ZS',nzsoil_extin,var1d,istatus)
    CALL mpupdater(var1d,nzsoil_extin)

    DO k = 1, nzsoil_extin
      zpsoil_ext(:,:,k+1) = var1d(k)
    END DO
    zpsoil_ext(:,:,1) = 0.0

!    IF (lvldbg > 25) WRITE(*,'(1x,a)') 'Reading variable TSLB ...'
    CALL get_wrf_3d(ncid,itime,1,(/'TSLB      '/),                      &
                    nxs_ext,nys_ext,nzsoil_extin,var3d,                 &
                    nx_ext,ny_ext,nzsoil_extin,tsoil_ext,istatus)

!    IF (lvldbg > 25) WRITE(*,'(1x,a)') 'Reading variable SMOIS ...'
    CALL get_wrf_3d(ncid,itime,1,(/'SMOIS     '/),                      &
                    nxs_ext,nys_ext,nzsoil_extin,var3d,                 &
                    nx_ext,ny_ext,nzsoil_extin,qsoil_ext,istatus)

    DEALLOCATE(var1d)

    DO k = nzsoil_extin,1,-1
      DO j = 1,ny_ext
        DO i = 1,nx_ext
          tsoil_ext(i,j,k+1) = tsoil_ext(i,j,k)
          qsoil_ext(i,j,k+1) = qsoil_ext(i,j,k)
        END DO
      END DO
    END DO

!    IF (lvldbg > 25) WRITE(*,'(1x,a)') 'Reading variable TSK ...'
    CALL get_wrf_2d(ncid,itime,1,(/'TSK    '/),nxs_ext,nys_ext,var2d,   &
                  nx_ext,ny_ext,tsoil_ext(:,:,1),istatus)

!    IF (lvldbg > 25) WRITE(*,'(1x,a)') 'Reading variable SST ...'
    CALL get_wrf_2d(ncid,itime,1,(/'SST    '/),nxs_ext,nys_ext,var2d,   &
                  nx_ext,ny_ext,qsoil_ext(:,:,1),istatus)
                  ! qsoil_ext is used as temporary array for SST

    DO j = 1, ny_ext
      DO i = 1, nx_ext
        IF (soiltyp_ext(i,j) == 13 .OR. soiltyp_ext(i,j) == 12) THEN
          tsoil_ext(i,j,1) = qsoil_ext(i,j,1)
        END IF
        qsoil_ext(i,j,1) = qsoil_ext(i,j,2)
        ! set surface value to be the same as the highest soil layer
      END DO
    END DO
  END IF

!
!-----------------------------------------------------------------------
!
!  Close data file
!
!-----------------------------------------------------------------------
!

  IF (myproc == 0) CALL close_ncd_file(ncid)

  DEALLOCATE( var2d )
  DEALLOCATE( var3d )

!  IF (myproc == 0) THEN
!    CALL get_ncd_scalar(ncid,itime,'LAT_LL_D',latsw,istatus)
!    CALL get_ncd_scalar(ncid,itime,'LON_LL_D',lonsw,istatus)
!  END IF
!  CALL mpupdater(latsw,1)
!  CALL mpupdater(lonsw,1)

   latsw = latv_ext(1,1)
   lonsw = lonu_ext(1,1)
!
!-----------------------------------------------------------------------
!
!  Rotate winds to be relative to true north.
!  The WRF data are sent as grid-relative.
!
!-----------------------------------------------------------------------
!
  ALLOCATE(utmp(nx_ext,ny_ext), STAT = istatus)
  ALLOCATE(vtmp(nx_ext,ny_ext), STAT = istatus)

  CALL setmapr(iproj_ext,scale_ext,latnot_ext,trlon_ext)
  CALL lltoxy(1,1,latsw,lonsw,x0_ext,y0_ext)

  DO k=1,nz_ext
    CALL uvmptoe(nx_ext,ny_ext,u_ext(:,:,k),v_ext(:,:,k),               &
                 lon_ext,utmp,vtmp)
    u_ext(:,:,k) = utmp(:,:)
    v_ext(:,:,k) = vtmp(:,:)
  END DO

  DEALLOCATE(utmp)
  DEALLOCATE(vtmp)

  istatus = 1
!
!-----------------------------------------------------------------------
!
!  Reset map projection to previous values
!
!-----------------------------------------------------------------------
!

  CALL setmapr(iproj,scale,latnot,trlon)
  CALL setorig(1,x0,y0)

  RETURN
END SUBROUTINE getwrfdata

SUBROUTINE get_wrf_2d(ncid,itime,numvar,varnames,nxs,nys,tem2d,         &
                      nx,ny,var2d,istatus)
!#######################################################################
!
! Read a 2D field from a WRF netCDF file.
!
!-----------------------------------------------------------------------
!
! AUTHOR: Y. Wang (11/23/2011)
!
!-----------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ncid
  INTEGER, INTENT(IN) :: itime, numvar
  INTEGER, INTENT(IN) :: nxs, nys
  INTEGER, INTENT(IN) :: nx,ny

  CHARACTER(LEN=7), INTENT(IN) :: varnames(numvar)

  REAL,    INTENT(OUT) :: tem2d(nxs,nys)
  REAL,    INTENT(OUT) :: var2d(nx,ny)

  INTEGER, INTENT(OUT) :: istatus

  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------

  INTEGER :: i,j
  INTEGER :: ivar

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (myproc == 0) CALL get_ncd_2d(ncid,itime,varnames(1),nxs,nys,tem2d,istatus)

  CALL mpupdatei(istatus,1)
  IF (istatus /= 0) THEN
    WRITE(*,'(1x,3a)') 'WARNING: Variable (',varnames(1),') is not found.'
    !CALL arpsstop('WARNING: variable reading problem in get_wrf_2d.',1)
    RETURN
  ELSE
    CALL mpupdater(tem2d,nxs*nys)
    DO j = 1,nys
      DO i = 1,nxs
        var2d(i,j) = tem2d(i,j)
      END DO
    END DO
  END IF

  DO ivar = 2, numvar

    IF (myproc == 0) CALL get_ncd_2d(ncid,itime,varnames(ivar),nxs,nys,tem2d,istatus)

    CALL mpupdatei(istatus,1)
    IF (istatus /= 0) THEN
      WRITE(*,'(1x,3a)') 'WARNING: Variable (',varnames(ivar),') is not found.'
      !CALL arpsstop('WARNING: variable reading problem in get_wrf_2d.',1)
      RETURN
    ELSE
      CALL mpupdater(tem2d,nxs*nys)
      DO j = 1,nys
        DO i = 1,nxs
          var2d(i,j) = var2d(i,j) + tem2d(i,j)
        END DO
      END DO
    END IF

  END DO

  CALL edgfill(var2d,1,nx,1,nxs,1,ny,1,nys,1,1,1,1)

  RETURN
END SUBROUTINE get_wrf_2d

SUBROUTINE get_wrf_3d(ncid,itime,numvar,varnames,                       &
                      nxs,nys,nzs,tem3d,                                &
                      nx,ny,nz,var3d,istatus)
!#######################################################################
!
! Read a 3D field from a WRF netCDF file.
!
!-----------------------------------------------------------------------
!
! AUTHOR: Y. Wang (11/23/2011)
!
!-----------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ncid
  INTEGER, INTENT(IN) :: itime, numvar
  INTEGER, INTENT(IN) :: nxs,nys,nzs
  INTEGER, INTENT(IN) :: nx,ny,nz

  CHARACTER(LEN=10), INTENT(IN) :: varnames(numvar)

  REAL,    INTENT(OUT) :: tem3d(nxs,nys,nzs)
  REAL,    INTENT(OUT) :: var3d(nx,ny,nz)

  INTEGER, INTENT(OUT) :: istatus

  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------

  INTEGER :: i,j,k
  INTEGER :: ivar

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  IF (myproc == 0) CALL get_ncd_3d(ncid,itime,varnames(1),nxs,nys,nzs,tem3d,istatus)

  CALL mpupdatei(istatus,1)
  IF (istatus /= 0) THEN
    !WRITE(*,'(1x,3a)') 'WARNING: Variable (',varnames(1),') is not found.'
    !CALL arpsstop('WARNING: variable reading problem in get_wrf_3d.',1)
    RETURN
  ELSE
    WRITE(6,*) 'Read in variable '//TRIM(varnames(1))
    CALL mpupdater(tem3d,nxs*nys*nzs)
    DO k = 1,nzs
      DO j = 1,nys
        DO i = 1,nxs
          var3d(i,j,k) = tem3d(i,j,k)
        END DO
      END DO
    END DO
  END IF

  DO ivar = 2, numvar

    IF (myproc == 0) CALL get_ncd_3d(ncid,itime,varnames(ivar),nxs,nys,nzs,tem3d,istatus)

    CALL mpupdatei(istatus,1)
    IF (istatus /= 0) THEN
      !WRITE(*,'(1x,3a)') 'WARNING: Variable (',varnames(ivar),') is not found.'
      !CALL arpsstop('WARNING: variable reading problem in get_wrf_3d.',1)
      RETURN
    ELSE
      WRITE(6,*) 'Read in variable '//TRIM(varnames(ivar))
      CALL mpupdater(tem3d,nxs*nys*nzs)
      DO k = 1, nzs
        DO j = 1,nys
          DO i = 1,nxs
            var3d(i,j,k) = var3d(i,j,k) + tem3d(i,j,k)
            !WRITE(*,*) i,j,k,var3d(i,j,k)
          END DO
        END DO
      END DO
    END IF

  END DO

  CALL edgfill(var3d,1,nx,1,nxs,1,ny,1,nys,1,nz,1,nzs)

  RETURN
END SUBROUTINE get_wrf_3d

SUBROUTINE get_wrf_soil(ncid,itime,nxs,nys,nzsoils,                     &
                       nx,ny,nzsoil,zpsoil,tsoil,qsoil,                 &
                       tem1,istatus)

!#######################################################################
!
! Read WRF soil variables and interpolate to ARPS 2 soil layers.
!
!-----------------------------------------------------------------------
!
! AUTHOR: Y. Wang (11/23/2011)
!
!-----------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ncid
  INTEGER, INTENT(IN) :: itime
  INTEGER, INTENT(IN) :: nxs,nys,nzsoils
  INTEGER, INTENT(IN) :: nx,ny,nzsoil
  REAL,    INTENT(OUT) :: zpsoil(nx,ny,nzsoil)
  REAL,    INTENT(OUT) :: tsoil(nx,ny,nzsoil), qsoil(nx,ny,nzsoil)
  REAL,    INTENT(OUT) :: tem1(nxs,nys,nzsoils)
  INTEGER, INTENT(OUT) :: istatus

  INCLUDE 'mp.inc'

!-----------------------------------------------------------------------
  CHARACTER(LEN=7) :: varname

  REAL, ALLOCATABLE :: tsoil_in(:,:,:), qsoil_in(:,:,:)
  REAL, ALLOCATABLE :: zpsoil_in(:)

  REAL    :: w1, w2, dist
  INTEGER :: m, i, j, k

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


  ALLOCATE(zpsoil_in(nzsoils),       STAT = istatus)
  ALLOCATE(tsoil_in (nx,ny,nzsoils), STAT = istatus)
  ALLOCATE(qsoil_in (nx,ny,nzsoils), STAT = istatus)

  IF (myproc == 0) CALL get_ncd_1d(ncid,itime,'ZS',nzsoils,zpsoil_in,istatus)
  CALL mpupdater(zpsoil_in,nzsoils)

  CALL get_wrf_3d(ncid,itime,1,(/'TSLB   '/),                           &
                  nxs,nys,nzsoils,tem1,                                 &
                  nx,ny,nzsoils,tsoil_in,istatus)

  CALL get_wrf_3d(ncid,itime,1,(/'SMOIS  '/),                           &
                  nxs,nys,nzsoils,tem1,                                 &
                  nx,ny,nzsoils,qsoil_in,istatus)

!
!-----------------------------------------------------------------------
!
!  Loop through all ARPS grid points
!
!-----------------------------------------------------------------------
!
  zpsoil(:,:,1) = 0.05
  zpsoil(:,:,2) = 0.55

  DO k = 1, nzsoil
    IF(zpsoil(2,2,k) < zpsoil_in(1)) THEN            ! extrapolation
      dist = zpsoil_in(1) - zpsoil(2,2,k)

      DO j=1,ny
        DO i=1,nx
          w1 = (tsoil_in(i,j,1) - tsoil_in(i,j,2)) /                    &
               (zpsoil_in(2)- zpsoil_in(1))
          w2 = (qsoil_in(i,j,1) - qsoil_in(i,j,2)) /                    &
               (zpsoil_in(2)- zpsoil_in(1))
          tsoil(i,j,k) = tsoil_in(i,j,1) + w1*dist
          qsoil(i,j,k) = qsoil_in(i,j,1) + w2*dist
        END DO
      END DO

    ELSE IF(zpsoil(2,2,k) >= zpsoil_in(nzsoils)) THEN ! constant layers

      DO j=1,ny
        DO i=1,nx
          tsoil(i,j,k) = tsoil_in(i,j,nzsoils)
          qsoil(i,j,k) = qsoil_in(i,j,nzsoils)
        END DO
      END DO

    ELSE                         ! linear interpolation

      DO m = 1, nzsoils
        IF(zpsoil_in(m) > zpsoil(2,2,k)) EXIT
      END DO

      DO j=1,ny
        DO i=1,nx
          w1 = (zpsoil_in(m) - zpsoil(i,j,k)) /                         &
               (zpsoil_in(m) - zpsoil_in(m-1))
          w2 = (zpsoil(i,j,k) - zpsoil_in(m-1) ) /                      &
               (zpsoil_in(m) - zpsoil_in(m-1))
          tsoil(i,j,k) = tsoil_in(i,j,m-1)*w1 + tsoil_in(i,j,m)*w2
          qsoil(j,j,k) = qsoil_in(i,j,m-1)*w1 + qsoil_in(i,j,m)*w2

        END DO
      END DO

    END IF
  END DO

!-----------------------------------------------------------------------

  DEALLOCATE(zpsoil_in,tsoil_in, qsoil_in)

  RETURN
END SUBROUTINE get_wrf_soil
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE getwrffname                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_wrf_fname(extdopt,extdfmt,dir_extd,extdname,             &
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
!   11/22/2011 (Yunheng Wang)
!   Modified based on getgrbfname.
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
  CHARACTER(LEN=4),   INTENT(IN)  :: appstr

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
  INTEGER :: grid_id

  INTEGER :: iyr,imo,iday,myr, jldy
  INTEGER :: ihr,imin,isec
  INTEGER :: ifhr,ifmin,ifsec

  INTEGER :: abstsec
  INTEGER :: iyear,imonth,idayy,ihour,iminute,isecond

  CHARACTER(LEN=256) :: gribfile_new
  INTEGER            :: grbflen_new

  CHARACTER(LEN=32)  :: fmtstr

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = 0

  grid_id = 1

  READ (extdinit,'(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)')                  &
                  iyr,  imo, iday,  ihr, imin, isec

  !CALL julday(iyr,imo,iday,jldy)

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

    WRITE(fmtstr,'(a,2(I1,a))') '(2a,I2.2,a,I4.4,5(a,I2.2))'

    CALL ctim2abss( iyr,imo,iday,ihr,imin,isec, abstsec )
    abstsec = abstsec + ifhr*3600 + ifmin*60 + ifsec
    CALL abss2ctim( abstsec, iyear, imonth, idayy, ihour, iminute, isecond )

    SELECT CASE (extdopt)
    CASE (23)
      WRITE(gribfile,fmtstr) dir_extd(1:lendir),'wrfout_d',grid_id,'_', &
               iyear,'-',imonth,'-',idayy,'_',ihour,':',iminute,':',isecond
    CASE DEFAULT
      WRITE(6,'(1x,a,I4,a,/)')                                          &
        'ERROR: Expecting WRF file name for extdopt = ',extdopt,        &
              ' is not valid in GETWRFFNAME.'
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
END SUBROUTINE get_wrf_fname

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE getwrfdims                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getwrfdims(dir_extd,extdname,extdopt,extdfmt,                &
                      extdinit,extdfcst,nx_wrf,ny_wrf,nz_wrf,nzsoil_wrf,&
                      qnames_wrf,iret)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This subroutine extracts dimension size from WRF history file in netCDF format.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yunheng Wang
!  11/22/2011 Initialized.
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
!    nx_wrf        Grid size in west-east direction
!    ny_wrf        Grid size in south-north direction
!    nz_wrf        Grid size in bottom-top direction
!    nzsoil_wrf
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

  INTEGER,          INTENT(OUT) :: nx_wrf       ! Number of points along x-axis
  INTEGER,          INTENT(OUT) :: ny_wrf       ! Number of points along y-axis
  INTEGER,          INTENT(OUT) :: nz_wrf       ! Number of points along y-axis
  INTEGER,          INTENT(OUT) :: nzsoil_wrf       ! Number of points along y-axis

  CHARACTER(LEN=10), INTENT(OUT) :: qnames_wrf(20)

  INTEGER,          INTENT(OUT) :: iret         ! Return flag

  INCLUDE 'mp.inc'
  INCLUDE 'globcst.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. Local variables
!
!-----------------------------------------------------------------------
!
  CHARACTER(LEN=256) :: wrffile
  CHARACTER(LEN=14)  :: gribtime

  INTEGER            :: grbflen, grbtlen
  INTEGER            :: nid

  INTEGER            :: mphys
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  CALL get_wrf_fname(extdopt,extdfmt,dir_extd,extdname,                 &
                   extdinit,extdfcst,'    ',                            &
                   wrffile,grbflen,gribtime,grbtlen,iret)
  IF (iret /=0 ) RETURN

  IF (myproc == 0) THEN

    CALL open_ncd_file(wrffile,nid)

    CALL get_ncd_dimension(nid,nx_wrf,ny_wrf,nz_wrf,nzsoil_wrf,iret)

    CALL get_ncd_dom_ti_integer(nid,'MP_PHYSICS',mphys,iret)

    CALL close_ncd_file(nid)

  END IF         ! IF (myproc == 0)
  CALL mpupdatei(iret, 1)

  CALL mpupdatei(nx_wrf,1)
  CALL mpupdatei(ny_wrf,1)
  CALL mpupdatei(nz_wrf,1)
  CALL mpupdatei(nzsoil_wrf,1)

  CALL mpupdatei(mphys,1)

!-----------------------------------------------------------------------
!
! Set up moist and scalar variables to be read based on mp_physics
!
!-----------------------------------------------------------------------

  P_QC = 0;   P_NC = 0
  P_QR = 0;   P_NR = 0
  P_QI = 0;   P_NI = 0
  P_QS = 0;   P_NS = 0
  P_QG = 0;   P_NG = 0
  P_QH = 0;
  nscalar  = 0
  nscalarq = 0
  qnames(:)= ' '; qdescp(:)= ' '

  SELECT CASE (mphys)
  CASE (0)                           ! passiveqv
    nscalar = 0
  CASE (1,3)                         ! kesslerscheme, wsm3scheme
    P_QC      = 1
    P_QR      = 2
    nscalar   = 2
    nscalarq  = 2
  CASE (2,6)                         ! linscheme, wsm6scheme
    P_QC      = 1
    P_QR      = 2
    P_QI      = 3
    P_QS      = 4
    P_QG      = 5
    nscalar    = 5
    nscalarq   = 5

  CASE (4,99)                        ! wsm5scheme,ncepcloud5
    P_QC      = 1
    P_QR      = 2
    P_QI      = 3
    P_QS      = 4
    nscalar    = 4
    nscalarq   = 4
  CASE (5)                           ! etampnew
    P_QC      = 1
    P_QR      = 2
    P_QS      = 3
    nscalar    = 3
    nscalarq   = 3

  CASE (8)                           ! thompson
    P_QC      = 1
    P_QR      = 2
    P_QI      = 3
    P_QS      = 4
    P_QG      = 5

    P_NR      = 6
    P_NI      = 7

    nscalarq  = 5;    nscalar   = 7

  CASE (10)                            ! Morrison DB
    P_QC = 1
    P_QR = 2
    P_QI = 3
    P_QS = 4
    P_QG = 5

    P_NR = 6
    P_NI = 7
    P_NS = 8
    P_NG = 9

    nscalarq = 5;     nscalar  = 9

  CASE (14)                            ! DM-5
    P_QC = 1
    P_QR = 2
    P_QI = 3
    P_QS = 4

    P_NC = 5
    P_NR = 6
    !P_QNN = 3               ! Do not know how to handle it in the ARPS grid?

    nscalarq = 4;    nscalar  = 6

  CASE (16)                            ! DM-6
    P_QC = 1
    P_QR = 2
    P_QI = 3
    P_QS = 4
    P_QG = 5

    P_NC = 6
    P_NR = 7
    !P_QNN = 3

    nscalarq = 5;    nscalar  = 7

  CASE (98)                            ! Old Thompson (07)
    P_QC = 1
    P_QR = 2
    P_QI = 3
    P_QS = 4
    P_QG = 5

    P_NI = 6

    nscalarq = 5;    nscalar  = 6

  CASE DEFAULT
    iret = -1
    WRITE(6,'(/,1x,a,I2,a,/)') 'ERROR: Wrong parameter - mp_physics = ',&
                                mphys,'.'
  END SELECT

  IF (P_QC > 0) THEN
    qnames(P_QC) = 'qc'; qdescp(P_QC) = 'Cloud water mixing ratio (kg/kg)'
  END IF
  IF (P_QR > 0) THEN
    qnames(P_QR) = 'qr'; qdescp(P_QR) = 'Rain  water mixing ratio (kg/kg)'
  END IF
  IF (P_QI > 0) THEN
    qnames(P_QI) = 'qi'; qdescp(P_QI) = 'Cloud ice   mixing ratio (kg/kg)'
  END IF
  IF (P_QS > 0) THEN
    qnames(P_QS) = 'qs'; qdescp(P_QS) = 'Snow mixing ratio (kg/kg)'
  END IF
  IF (P_QG > 0) THEN
    qnames(P_QG) = 'qg'; qdescp(P_QG) = 'Graupel mixing ratio (kg/kg)'
  END IF

  ! Number of concentrations
  IF (P_NC > 0) THEN
    qnames(P_NC) = 'nc'; qdescp(P_NC) = 'Cloud water concentrations (#/m3)'
  END IF
  IF (P_NR > 0) THEN
    qnames(P_NR) = 'nr'; qdescp(P_NR) = 'Rain water concentrations (#/m3)'
  END IF
  IF (P_NI > 0) THEN
    qnames(P_NI) = 'ni'; qdescp(P_NI) = 'Cloud ice concentrations (#/m3)'
  END IF
  IF (P_NS > 0) THEN
    qnames(P_NS) = 'ns'; qdescp(P_NS) = 'Snow concentrations (#/m3)'
  END IF
  IF (P_NG > 0) THEN
    qnames(P_QG) = 'qg'; qdescp(P_NG) = 'Graupel concentrations (#/m3)'
  END IF

!-----------------------------------------------------------------------
!
! WRF variable names to be read.
!
!-----------------------------------------------------------------------

  IF (P_QC > 0)   qnames_WRF(P_QC) = 'QCLOUD    '
  IF (P_QR > 0)   qnames_WRF(P_QR) = 'QRAIN     '
  IF (P_QI > 0)   qnames_WRF(P_QI) = 'QICE      '
  IF (P_QS > 0)   qnames_WRF(P_QS) = 'QSNOW     '
  IF (P_QG > 0)   qnames_WRF(P_QG) = 'QGRAUP    '

  IF (P_NC > 0)   qnames_WRF(P_NC) = 'QNCLOUD   '
  IF (P_NR > 0)   qnames_WRF(P_NR) = 'QNRAIN    '
  IF (P_NI > 0)   qnames_WRF(P_NI) = 'QNICE     '
  IF (P_NS > 0)   qnames_WRF(P_NS) = 'QNSNOW    '
  IF (P_NG > 0)   qnames_WRF(P_NG) = 'QNGRAUPEL '

  RETURN
END SUBROUTINE getwrfdims
