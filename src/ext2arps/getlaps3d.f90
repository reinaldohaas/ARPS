!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETLAPS                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getlaps(nx_ext,ny_ext,nz_ext,dir_extd,                       &
           extdinit,extdfcst,julfname,i4time,                           &
           iproj_ext,scale_ext,                                         &
           trlon_ext,latnot_ext,x0_ext,y0_ext,                          &
           lat_ext,lon_ext,                                             &
           p_ext,hgt_ext,t_ext,qv_ext,u_ext,v_ext,                      &
           istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  OLAPS version.
!
!  Reads OLAPS file for processing by ext2arps, a program
!  that converts external files to ARPS variables and format.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  May, 1994
!
!  MODIFICATION HISTORY:
!  September, 1994 (KB)
!  Upgrade for LLNL effort.
!
!  November, 1994 (KB)
!  Beefed up documentation.
!  OLAPS Version.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    dir_extd      Directory name for external file
!    extdinit      Initialized time in mm-dd-yyyy:hh:mm:ss format
!    extdfcst      Forecast hour in HHH:MM:SS format
!    julfname      File name in yyjjjhhmm format
!    i4time        Absolute time in seconds (from 1 Jan 1970) of desired data
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
!    istatus       status indicator
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  CHARACTER (LEN=*) :: dir_extd
  CHARACTER (LEN=19) :: extdinit
  CHARACTER (LEN=9) :: extdfcst
  CHARACTER (LEN=9) :: julfname
  INTEGER :: i4time
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
!  External grid variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iproj_ext
  REAL :: scale_ext,trlon_ext
  REAL :: latnot_ext(2)
  REAL :: x0_ext,y0_ext
  REAL :: olatsw,olonsw
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
  REAL :: p_ext(nx_ext,ny_ext,nz_ext)     ! Pressure (Pascals)
  REAL :: hgt_ext(nx_ext,ny_ext,nz_ext)   ! Height (m)
  REAL :: t_ext(nx_ext,ny_ext,nz_ext)     ! Temperature (K)
  REAL :: qv_ext(nx_ext,ny_ext,nz_ext)    ! Specific humidity (kg/kg)
  REAL :: u_ext(nx_ext,ny_ext,nz_ext)     ! Eastward wind component
  REAL :: v_ext(nx_ext,ny_ext,nz_ext)     ! Northward wind component
!  REAL :: qc_ext(nx_ext,ny_ext,nz_ext)    ! Cloud H2O mixing ratio (kg/kg)
!  REAL :: qr_ext(nx_ext,ny_ext,nz_ext)    ! Rain  H2O mixing ratio (kg/kg)
!  REAL :: qi_ext(nx_ext,ny_ext,nz_ext)    ! Ice   mixing ratio (kg/kg)
!  REAL :: qs_ext(nx_ext,ny_ext,nz_ext)    ! Snow  mixing ratio (kg/kg)
!  REAL :: qh_ext(nx_ext,ny_ext,nz_ext)    ! Hail  mixing ratio (kg/kg)

  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  OLAPS file grid variables
!
!-----------------------------------------------------------------------
!
  REAL :: dxlaps,dylaps
  PARAMETER (dxlaps=10000., dylaps=10000.)
!
  REAL :: xlaps(nx_ext),ylaps(ny_ext)
  REAL :: lat(nx_ext,ny_ext),lon(nx_ext,ny_ext)
  REAL :: terr(nx_ext,ny_ext)
!
!-----------------------------------------------------------------------
!
!  OLAPS Forecast variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=50) :: dir_lsx,dir_lw3,dir_lq3,dir_lt1,dir_static
  CHARACTER (LEN=31) :: ext,ext_s

  INTEGER :: nz_ext_mx
  PARAMETER( nz_ext_mx = 41)
  INTEGER :: lvl(nz_ext_mx)
  CHARACTER (LEN=3) :: var_t(nz_ext_mx),var_q(nz_ext_mx)
  CHARACTER (LEN=3) :: var_u(nz_ext_mx),var_v(nz_ext_mx)
  CHARACTER (LEN=4) :: lvl_coord(nz_ext_mx)

  CHARACTER (LEN=4) :: sfc_coord
  CHARACTER (LEN=10) :: units(nz_ext)
  CHARACTER (LEN=125) :: comment(nz_ext)
!
  REAL :: sfct(nx_ext,ny_ext)
  REAL :: sfcmr(nx_ext,ny_ext)
  REAL :: sfcp_pa(nx_ext,ny_ext)
  REAL :: psnd(nz_ext)
  REAL :: tsnd(nz_ext)
  REAL :: htsnd(nz_ext)
  REAL :: mrsnd(nz_ext)
!
!-----------------------------------------------------------------------
!
!  OLAPS initializations for reading subroutines
!
!-----------------------------------------------------------------------
!
  DATA lvl /1100,1075,1050,1025,1000,975,950,925,900,                   &
                  875, 850, 825, 800,775,750,725,700,                   &
                  675, 650, 625, 600,575,550,525,500,                   &
                  475, 450, 425, 400,375,350,325,300,                   &
                  275, 250, 225, 200,175,150,125,100/
!
  DATA lvl_coord                                                        &
          /'HPA ','HPA ','HPA ','HPA ','HPA ','HPA ','HPA ','HPA ',     &
           'HPA ','HPA ','HPA ','HPA ','HPA ','HPA ','HPA ','HPA ',     &
           'HPA ','HPA ','HPA ','HPA ','HPA ','HPA ','HPA ','HPA ',     &
           'HPA ','HPA ','HPA ','HPA ','HPA ','HPA ','HPA ','HPA ',     &
           'HPA ','HPA ','HPA ','HPA ','HPA ','HPA ','HPA ','HPA ',     &
           'HPA '/
  DATA var_t /'t  ','t  ','t  ','t  ','t  ','t  ','t  ','t  ',          &
              't  ','t  ','t  ','t  ','t  ','t  ','t  ','t  ',          &
              't  ','t  ','t  ','t  ','t  ','t  ','t  ','t  ',          &
              't  ','t  ','t  ','t  ','t  ','t  ','t  ','t  ',          &
              't  ','t  ','t  ','t  ','t  ','t  ','t  ','t  ',          &
              't  '/
  DATA var_q /'sh ','sh ','sh ','sh ','sh ','sh ','sh ','sh ',          &
              'sh ','sh ','sh ','sh ','sh ','sh ','sh ','sh ',          &
              'sh ','sh ','sh ','sh ','sh ','sh ','sh ','sh ',          &
              'sh ','sh ','sh ','sh ','sh ','sh ','sh ','sh ',          &
              'sh ','sh ','sh ','sh ','sh ','sh ','sh ','sh ',          &
              'sh '/
  DATA var_u /'u  ','u  ','u  ','u  ','u  ','u  ','u  ','u  ',          &
              'u  ','u  ','u  ','u  ','u  ','u  ','u  ','u  ',          &
              'u  ','u  ','u  ','u  ','u  ','u  ','u  ','u  ',          &
              'u  ','u  ','u  ','u  ','u  ','u  ','u  ','u  ',          &
              'u  ','u  ','u  ','u  ','u  ','u  ','u  ','u  ',          &
              'u  '/
  DATA var_v /'v  ','v  ','v  ','v  ','v  ','v  ','v  ','v  ',          &
              'v  ','v  ','v  ','v  ','v  ','v  ','v  ','v  ',          &
              'v  ','v  ','v  ','v  ','v  ','v  ','v  ','v  ',          &
              'v  ','v  ','v  ','v  ','v  ','v  ','v  ','v  ',          &
              'v  ','v  ','v  ','v  ','v  ','v  ','v  ','v  ',          &
              'v  '/
!
!-----------------------------------------------------------------------
!
!  Misc internal variables
!
!-----------------------------------------------------------------------
!
  REAL :: smfact
  PARAMETER (smfact=0.5)
  INTEGER :: i,j,k,ldir,kdim,lsfc
  REAL :: swlat,swlon,xsw,ysw,grid_spacing
  REAL :: qvlast,t_sfc,w_sfc
!
!-----------------------------------------------------------------------
!
!  Functions
!
!-----------------------------------------------------------------------
!
  REAL :: qvtomxr
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'lapsparms.cmn'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  iproj_ext=1          ! polar stereographic
  scale_ext=1.0        ! report lengths in m
  trlon_ext=-105.      ! orientation of OLAPS data grids
  latnot_ext(1)=40.
  latnot_ext(2)=40.
  x0_ext=0.
  y0_ext=0.
  swlat=32.5353
  swlon=-103.7612
!
  IF(dir_extd(1:5) == '     ') THEN
    dir_extd='/vortex/vortex95/lapsprd/'
    WRITE(6,'(a,a)')                                                    &
        ' Use the default directory of OLAPS analyses: ',               &
        ' /vortex/vortex95/lapsprd'
  ELSE
    WRITE(6,'(a,a)')                                                    &
        ' Use the input directory of OLAPS analyses: ',dir_extd
  END IF
!
  ldir=LEN(dir_extd)
  CALL strlnth(dir_extd,ldir)
  IF(dir_extd(ldir:ldir) /= '/') THEN
    ldir=ldir+1
    dir_extd(ldir:ldir)='/'
  END IF
!
  WRITE(dir_lsx,'(a,a)') dir_extd(1:ldir),'lsx/'
  WRITE(dir_lq3,'(a,a)') dir_extd(1:ldir),'lq3/'
  WRITE(dir_lw3,'(a,a)') dir_extd(1:ldir),'lw3/'
  WRITE(dir_lt1,'(a,a)') dir_extd(1:ldir),'lt1/'
!
  WRITE(6,'(a)')                                                        &
      ' The time in format yydddhhff of OLAPS initial file is ',julfname
!
  CALL cv_asc_i4time(julfname,i4time)
!
  DO k=1,nz_ext
    DO j=1,ny_ext
      DO i=1,nx_ext
        p_ext(i,j,k)=100.*FLOAT(lvl(k))
      END DO
    END DO
  END DO
!
  CALL getmapr(iproj,scale,latnot,trlon,x0,y0)
  CALL setmapr(iproj_ext,scale_ext,latnot_ext,trlon_ext)
  CALL setorig(1,x0_ext,y0_ext)
!
!-----------------------------------------------------------------------
!
!  Fake that get_laps_parms was called.
!  This is to insure proper reading of OLAPS terrain
!  file.
!
!-----------------------------------------------------------------------
!
  standard_latitude=latnot_ext(1)
  standard_longitude=trlon_ext
  iflag_lapsparms_cmn = 1
!
  dir_static='/vortex/vortex95/static/'
!
  ext_s='vortex95'
!
!-----------------------------------------------------------------------
!
!  Get OLAPS static variables
!
!-----------------------------------------------------------------------
!
  CALL rd_laps_static(dir_static,                                       &
               ext_s,nx_ext,ny_ext,1,                                   &
               'LAT',units,comment,                                     &
               lat_ext,grid_spacing,istatus)
  CALL rd_laps_static(dir_static,                                       &
               ext_s,nx_ext,ny_ext,1,                                   &
               'LON',units,comment,                                     &
               lon_ext,grid_spacing,istatus)
  CALL rd_laps_static(dir_static,                                       &
               ext_s,nx_ext,ny_ext,1,                                   &
               'AVG',units,comment,                                     &
               terr,grid_spacing,istatus)
!
!-----------------------------------------------------------------------
!
!  Reset map projection to previous values
!
!-----------------------------------------------------------------------
!
  CALL setmapr(iproj,scale,latnot,trlon)
  CALL setorig(1,x0,y0)
!
!-----------------------------------------------------------------------
!
!  Get file name information
!
!-----------------------------------------------------------------------
!
  ldir=LEN(dir_extd)
  CALL strlnth( dir_extd, ldir)
  IF(dir_extd(ldir:ldir) /= '/') THEN
    ldir=ldir+1
    dir_extd(ldir:ldir)='/'
  END IF
!
!-----------------------------------------------------------------------
!
!  Read LAPS surface pressure data
!
!-----------------------------------------------------------------------
!
  ext='lsx'
  kdim=1
  lsfc=0
  sfc_coord='AGL'
  CALL read_laps_data(i4time,dir_lsx,ext,                               &
                      nx_ext,ny_ext,1,kdim,                             &
                      'ps ',lsfc,sfc_coord,units,comment,               &
                      sfcp_pa,istatus)
  WRITE(6,'(a)')    ' Sfc Pressure Read'
  WRITE(6,'(a,i6)') ' istatus= ',istatus
  WRITE(6,'(a,a)')  ' units= ',units(1)
  IF(istatus /= 1) GO TO 598
!
!-----------------------------------------------------------------------
!
!  Read LAPS surface temperature data
!
!-----------------------------------------------------------------------
!
  CALL read_laps_data(i4time,dir_lsx,ext,nx_ext,ny_ext,1,kdim,          &
                      't  ',lsfc,sfc_coord,units,comment,               &
                      sfct,istatus)
  WRITE(6,'(a)')    ' Sfc Temp Read'
  WRITE(6,'(a,i6)') ' istatus= ',istatus
  IF(istatus /= 1) GO TO 598
!
!-----------------------------------------------------------------------
!
!  Read LAPS surface specific humidity
!
!-----------------------------------------------------------------------
!
  CALL read_laps_data(i4time,dir_lsx,ext,nx_ext,ny_ext,1,kdim,          &
                      'mr ',lsfc,sfc_coord,units,comment,               &
                      sfcmr,istatus)
  WRITE(6,'(a)')    ' Sfc specific humidity read'
  WRITE(6,'(a,i6)') ' istatus= ',istatus
  WRITE(6,'(a,a)')  ' units= ',units(1)
  IF(istatus /= 1) GO TO 598
!
!-----------------------------------------------------------------------
!
!  Read LAPS Temperature data
!
!-----------------------------------------------------------------------
!
  ext='lt1'
  kdim=nz_ext
  CALL read_laps_data(i4time,dir_lt1,ext,                               &
                      nx_ext,ny_ext,nz_ext,kdim,                        &
                      var_t,lvl,lvl_coord,units,comment,                &
                      t_ext, istatus)
  WRITE(6,'(a)')    ' Temperature Read'
  WRITE(6,'(a,i6)') ' istatus= ',istatus
  WRITE(6,'(a,a)')  ' units= ',units(1)
  IF(istatus /= 1) GO TO 598
!
!-----------------------------------------------------------------------
!
!  Read LAPS specific humidity field
!
!-----------------------------------------------------------------------
!
  ext='lq3'
  kdim=nz_ext
  CALL read_laps_data(i4time,dir_lq3,ext,                               &
                      nx_ext,ny_ext,nz_ext,kdim,                        &
                      var_q,lvl,lvl_coord,units,comment,                &
                      qv_ext,istatus)
  WRITE(6,'(a)')    ' Specific humidity read'
  WRITE(6,'(a,i6)') ' istatus= ',istatus
  WRITE(6,'(a,a)')  ' units= ',units(1)
  IF(istatus /= 1) GO TO 598
!
!-----------------------------------------------------------------------
!
!  Fix-up the q field so that q is the same as the first
!  q above ground.
!
!-----------------------------------------------------------------------
!
  DO j=1,ny_ext
    DO i=1,nx_ext
      qvlast=0.
      DO k=nz_ext,1,-1
        IF(qv_ext(i,j,k) > 1. .OR. qv_ext(i,j,k) < 0.) THEN
          qv_ext(i,j,k)=qvlast
        ELSE
          qvlast=qv_ext(i,j,k)
        END IF
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Integrate temperature from surface to get
!  height of pressure levels.
!  Convert temperature to C for hydrostatic integration
!  Convert specific humidity to mixing ratio
!
!-----------------------------------------------------------------------
!
  DO j=1,ny_ext
    DO i=1,nx_ext
      t_sfc=sfct(i,j)-273.15
      w_sfc=sfcmr(i,j)
      DO k=1,nz_ext
        psnd(k)=p_ext(i,j,k)
        tsnd(k)=t_ext(i,j,k)-273.15
        IF(qv_ext(i,j,k) > 0.) THEN
          mrsnd(k)=1000.*qvtomxr(qv_ext(i,j,k))
        ELSE
          mrsnd(k)=0.
        END IF
      END DO
      CALL getht(nz_ext,nz_ext,                                         &
             tsnd,mrsnd,psnd,                                           &
             terr(i,j),sfcp_pa(i,j),t_sfc,w_sfc,                        &
             htsnd)
      DO k=1,nz_ext
        hgt_ext(i,j,k)=htsnd(k)
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Read u component
!
!-----------------------------------------------------------------------
!
  ext='lw3'
  kdim=nz_ext
  CALL read_laps_data(i4time,dir_lw3,ext,                               &
                      nx_ext,ny_ext,nz_ext,nz_ext,                      &
                      var_u,lvl,lvl_coord,units,comment,                &
                      u_ext,istatus)
  PRINT *, ' U Component read'
  PRINT *, ' istatus= ',istatus
  PRINT *, ' units= ',units(1)
  IF(istatus /= 1) GO TO 598
!
!-----------------------------------------------------------------------
!
!  Read v component
!
!-----------------------------------------------------------------------
!
  ext='lw3'
  kdim=nz_ext
  CALL read_laps_data(i4time,dir_lw3,ext,                               &
                      nx_ext,ny_ext,nz_ext,nz_ext,                      &
                      var_v,lvl,lvl_coord,units,comment,                &
                      v_ext,istatus)
  PRINT *, ' V Component read'
  PRINT *, ' istatus= ',istatus
  PRINT *, ' units= ',units(1)
  IF(istatus /= 1) GO TO 598
!
!-----------------------------------------------------------------------
!
!  Fill qc,qr,qi,qs,qh arrays with missing value.
!
!-----------------------------------------------------------------------
!
!  DO k=1,nz_ext
!    DO j=1,ny_ext
!      DO i=1,nx_ext
!        qc_ext(i,j,k)=-999.
!        qr_ext(i,j,k)=-999.
!        qi_ext(i,j,k)=-999.
!        qs_ext(i,j,k)=-999.
!        qh_ext(i,j,k)=-999.
!      END DO
!    END DO
!  END DO
!
  istatus=1
  RETURN
!
!-----------------------------------------------------------------------
!
!  Error destination
!
!-----------------------------------------------------------------------
!
  598 CONTINUE
  WRITE(6,'(a)')  ' Error reading last field, returning'
  RETURN
END SUBROUTINE getlaps
!

SUBROUTINE getht(maxlev,nlevel,                                         &
           t,w,p_pa,elev,p_sfc_pa,t_sfc,w_sfc,ht)
  IMPLICIT NONE
!
!  Integrate hypsometric equation to get hydrostatic
!  heights from temperatures on pressure levels.
!
!  Input variables:   t  temperature (C)
!                     w  mixing ratio (g/kg)
!                     p_pa pressure in Pascals
!                     elev   elevation of surface in meters
!                     p_sfc_pa  surface pressure in Pascals
!                     t_sfc     surface temperature in C
!                     w_sfc     surface mixing ratio in g/kg
!
!  Output variable:
!                     ht        heights on pressure levels in meters
!
!
!  Keith Brewster Feb, 1994
!
!
!  Input variables
!
  INTEGER :: maxlev,nlevel
  REAL :: t(maxlev),w(maxlev),p_pa(maxlev)
  REAL :: elev,p_sfc_pa,t_sfc,w_sfc
!
!  Output
!
  REAL :: ht(maxlev)
!
!  Parameters
!
  REAL :: g,rd,rdovg
  PARAMETER (g=9.80,rd=287.0)
  PARAMETER (rdovg=(rd/g))
!
!  Functions
!
  REAL :: tctotv
!
!  Misc internal variables
!
  INTEGER :: k,ksfc
  REAL :: tvsfc,tvhi,tvlo,tvbar,thick
!
  DO k=1,nlevel-1
    IF(p_pa(k) < p_sfc_pa) EXIT
  END DO
  ksfc=k
!
  tvsfc=tctotv(t_sfc,w_sfc)
  tvhi =tctotv(t(ksfc),w(ksfc))
  tvbar=0.5*(tvsfc+tvhi)
  ht(ksfc)=elev+rdovg*tvbar*ALOG(p_sfc_pa/p_pa(ksfc))
  IF(ksfc > 1) THEN
    tvlo =tctotv(t(ksfc-1),w(ksfc-1))
    tvbar=0.5*(tvsfc+tvlo)
    ht(ksfc-1)=elev+rdovg*tvbar*ALOG(p_sfc_pa/p_pa(ksfc-1))
  END IF
  DO k=ksfc-2,1,-1
    tvlo=tctotv(t(k),w(k))
    tvhi=tctotv(t(k+1),w(k+1))
    tvbar=0.5*(tvlo+tvhi)
    ht(k)=ht(k+1)-rdovg*tvbar*ALOG(p_pa(k)/p_pa(k+1))
  END DO
  DO k=ksfc+1,nlevel
    tvlo=tctotv(t(k-1),w(k-1))
    tvhi=tctotv(t(k),w(k))
    tvbar=0.5*(tvlo+tvhi)
    ht(k)=ht(k-1)+rdovg*tvbar*ALOG(p_pa(k-1)/p_pa(k))
  END DO
  RETURN
END SUBROUTINE getht
!

  FUNCTION tctotv(tt,ww)
  IMPLICIT NONE
!
!  Virtual Temperature
!
!  Given T in Celcius and mixing ratio in g/kg
!  find the virtual temperature.
!
  REAL :: tctotv,tt,ww
  tctotv=(tt+273.15)*(1.+0.0006*ww)
  RETURN
  END FUNCTION tctotv
!

  FUNCTION qvtomxr(qv)
!
!  This is an approximation to the formula
!
!      qv = w/(1+w)
!
!  Or w = qv + qv*w
!       = qv + qv*(qv + qv*w)
!  So, when q is small, w is approx equal to q
!  so
!     w = qv + qv*(qv + qv*qv)
!
!  Keith Brewster, CAPS
!  Feb 1994
!
!
  IMPLICIT NONE
  REAL :: qvtomxr,qv
!
  qvtomxr=qv + qv*(qv + qv*qv)
!  print *, ' qv, mixr = ',qv,qvtomxr
!
  RETURN
  END FUNCTION qvtomxr
