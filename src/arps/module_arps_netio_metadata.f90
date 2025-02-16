MODULE arps_netio_metadata

!-----------------------------------------------------------------------
!
!  Define meta types
!
!-----------------------------------------------------------------------

  TYPE :: VARMETA               ! Structure to hold CF variable attributes
     CHARACTER(80) :: standard_name
     CHARACTER(80) :: long_name
     CHARACTER(40) :: units
     CHARACTER(1)  :: stagger
  END TYPE VARMETA

  TYPE :: ARPSVARS            ! All of ARPS variables are defined here
    TYPE(VARMETA) :: Time
    TYPE(VARMETA) :: x_stag, y_stag, z_stag, zp, zpsoil      ! grdout = 1
    TYPE(VARMETA) :: ubar, vbar, wbar, ptbar, pbar, qvbar    ! grdbas = 1
    TYPE(VARMETA) :: soiltyp, stypfrct,vegtyp,lai,roufns,veg ! landout= 1
    TYPE(VARMETA) ::    u,    v,    w,    pt, p    ! totout = 1, varout = 1
    TYPE(VARMETA) :: uprt, vprt, wprt, ptprt, pprt ! totout = 0, varout = 1
    TYPE(VARMETA) :: qvprt,qv, qscalar(20)                   ! mstout = 1
    TYPE(VARMETA) :: prcrate1, prcrate2, prcrate3, prcrate4  ! prcout = 1
    TYPE(VARMETA) :: raing, rainc
    TYPE(VARMETA) :: tke                                     ! tkeout = 1
    TYPE(VARMETA) :: kmh, kmv                                ! trbout = 1
    TYPE(VARMETA) :: tsoil, qsoil, wetcanp, snowdpth         ! sfcout = 1
    TYPE(VARMETA) :: radfrc, radsw, rnflx, radswnet, radlwin ! radflg = 1
    TYPE(VARMETA) :: usflx, vsflx, ptsflx, qvsflx            ! flxflg = 1
  END TYPE ARPSVARS

  TYPE TRNVARS
    TYPE(VARMETA) :: hterrain
  END TYPE TRNVARS

  TYPE SFCVARS
    TYPE(VARMETA) :: soiltyp
    TYPE(VARMETA) :: stypfrct
    TYPE(VARMETA) :: vegtyp
    TYPE(VARMETA) :: lai
    TYPE(VARMETA) :: roufns
    TYPE(VARMETA) :: veg
    TYPE(VARMETA) :: ndvi
  END TYPE SFCVARS

  TYPE SOILVARS
    TYPE(VARMETA) :: zpsoil
    TYPE(VARMETA) :: tsoil
    TYPE(VARMETA) :: qsoil
    TYPE(VARMETA) :: wetcanp
    TYPE(VARMETA) :: snowdpth
    TYPE(VARMETA) :: soiltyp
  END TYPE SOILVARS

  TYPE BDYVARS
    TYPE(VARMETA) :: ctime
    TYPE(VARMETA) :: u, v, w, pt, p
    TYPE(VARMETA) :: qv, qscalar(20)
  END TYPE BDYVARS

!-----------------------------------------------------------------------
!
! Define meta variables
!
!-----------------------------------------------------------------------

  TYPE (ARPSVARS), SAVE  :: arpsmeta
  TYPE (TRNVARS),  SAVE  :: trnmeta
  TYPE (SFCVARS),  SAVE  :: sfcmeta
  TYPE (SOILVARS), SAVE  :: soilmeta
  TYPE (BDYVARS),  SAVE  :: bdymeta

  LOGICAL,         SAVE  :: initialized = .false.

  PRIVATE :: init_varmeta

CONTAINS

  SUBROUTINE init_varmeta(P_QC, P_QR, P_QI, P_QS, P_QG, P_QH,           &
                          P_NC, P_NR, P_NI, P_NS, P_NG, P_NH,           &
                                P_ZR, P_ZI, P_ZS, P_ZG, P_ZH)

!#######################################################################
!
! Assign values to all defined variables
!
! Note: It should NOT be called outside of this module.
!
!#######################################################################
!
    IMPLICIT NONE
  
   INTEGER, INTENT(IN) :: P_QC,P_QR,P_QI,P_QS,P_QG,P_QH
   INTEGER, INTENT(IN) :: P_NC,P_NR,P_NI,P_NS,P_NG,P_NH
   INTEGER, INTENT(IN) ::      P_ZR,P_ZI,P_ZS,P_ZG,P_ZH

!----------------------------------------------------------------------
! 
! ARPS history variables
!
!----------------------------------------------------------------------

    arpsmeta%x_stag = VARMETA('projection_x_coordinate',                &
                                                 'X coordinate','m','X')
    arpsmeta%y_stag = VARMETA('projection_y_coordinate',                &
                                                 'Y coordinate','m','Y')
    arpsmeta%z_stag = VARMETA('height',                                 &
                                                 'Z coordinate','m','Z')
    arpsmeta%zp     = VARMETA('geopotential_height',                    &
                             'Physical height coordinate (MSL)','m','Z')
    arpsmeta%zpsoil = VARMETA('geopotential_height',                    &
          'Physical height coordinate (soil), meters above msl','m',' ')

    arpsmeta%ubar  = VARMETA('base_x_wind',                             &
                                    'Base state u-velocity','m s-1','X') !!
    arpsmeta%vbar  = VARMETA('base_y_wind',                             &
                                    'Base state v-velocity','m s-1','Y') !!
    arpsmeta%wbar  = VARMETA('base_upward_air_velocity',                &
                                    'Base state w-velocity','m s-1','Z') !!
    arpsmeta%ptbar = VARMETA('base_air_potential_temperature',          &
                           'Base state potential temperature','K',  ' ') !!
    arpsmeta%pbar  = VARMETA('base_air_pressure',                       &
                                        'Base state pressure','Pa', ' ') !!
    arpsmeta%qvbar = VARMETA('base_specific_humidity',                  &
                     'Base state water vapor specific humidity','1',' ') !!

    arpsmeta%soiltyp  = VARMETA('soil_category',                        &
                                          'Soil type',      'index',' ') !!
    arpsmeta%stypfrct = VARMETA('soil_fraction',                        &
                                'Soil type fractional coverage','1',' ') !!
    arpsmeta%vegtyp   = VARMETA('vegetation_category',                  &
                                          'Vegetation type','index',' ') !!
    arpsmeta%lai      = VARMETA('leaf_area_index',                      &
                                            'Leaf Area Index',  '1',' ')
    arpsmeta%roufns   = VARMETA('surface_roughness_length',             &
                                            'Surface roughness','m',' ')
    arpsmeta%veg      = VARMETA('vegetation_area_fraction',             &
                                          'Vegetation fraction','1',' ')

    arpsmeta%Time  = VARMETA('time','Data valid time','s',' ')

    arpsmeta%uprt  = VARMETA('perturbation_x_wind',                     &
                                  'Perturbation u-velocity','m s-1','X') !!
    arpsmeta%vprt  = VARMETA('perturbation_y_wind',                     &
                                  'Perturbation u-velocity','m s-1','Y') !!
    arpsmeta%wprt  = VARMETA('perturbation_upward_air_velocity',        &
                                  'Perturbation w-velocity','m s-1','Z') !!
    arpsmeta%ptprt = VARMETA('perturbation_air_potential_temperature',  &
                           'Perturbation potential temperature','K',' ') !!
    arpsmeta%pprt  = VARMETA('perturbation_air_pressure',               &
                                       'Perturbation pressure','Pa',' ') !!
    arpsmeta%qvprt = VARMETA('perturbation_specific_humidity',          &
                   'Perturbation water vapor specific humidity','1',' ') !!

    arpsmeta%u  = VARMETA('x_wind',                                     &
                                               'U-velocity','m s-1','X')
    arpsmeta%v  = VARMETA('y_wind',                                     &
                                               'V-velocity','m s-1','Y')
    arpsmeta%w  = VARMETA('upward_air_velocity',                        &
                                               'W-velocity','m s-1','Z')
    arpsmeta%pt = VARMETA('air_potential_temperature',                  &
                                        'Potential temperature','K',' ')
    arpsmeta%p  = VARMETA('air_pressure',                               &
                                                    'Pressure','Pa',' ')

    arpsmeta%qv = VARMETA('specific_humidity',                          &
                                'Water vapor specific humidity','1',' ') 
    IF (P_QC > 0) arpsmeta%qscalar(P_QC) =                              &
                  VARMETA('mass_fraction_of_cloud_liquid_water_in_air', &
                                     'Cloud water mixing ratio','1',' ')
    IF (P_QR > 0) arpsmeta%qscalar(P_QR) =                              &
                  VARMETA('mass_fraction_of_rain_water_in_air',         &
                                      'Rain water mixing ratio','1',' ') !!
    IF (P_QI > 0) arpsmeta%qscalar(P_QI) =                              &
                  VARMETA('mass_fraction_of_cloud_ice_in_air',          &
                                       'Cloud ice mixing ratio','1',' ')
    IF (P_QS > 0) arpsmeta%qscalar(P_QS) =                              &
                  VARMETA('mass_fraction_of_cloud_snow_water_in_air',   &
                                            'Snow mixing ratio','1',' ') !!
    IF (P_QG > 0) arpsmeta%qscalar(P_QG) =                              &
                  VARMETA('mass_fraction_of_cloud_graupel_water_in_air',&
                                         'Graupel mixing ratio','1',' ') !!
    IF (P_QH > 0) arpsmeta%qscalar(P_QH) =                              &
                  VARMETA('mass_fraction_of_cloud_graupel_water_in_air',&
                                            'Hail mixing ratio','1',' ') !!

    IF (P_NC > 0) arpsmeta%qscalar(P_NC) =                              &
                  VARMETA('mass_fraction_of_cloud_liquid_water_in_air', &
                          'Cloud water concentration number','m-3',' ')
    IF (P_NR > 0) arpsmeta%qscalar(P_NR) =                              &
                  VARMETA('mass_fraction_of_rain_water_in_air',         &
                           'Rain water concentration number','m-3',' ')
    IF (P_NI > 0) arpsmeta%qscalar(P_NI) =                              &
                  VARMETA('mass_fraction_of_cloud_ice_in_air',          &
                            'Cloud ice concentration number','m-3',' ')
    IF (P_NS > 0) arpsmeta%qscalar(P_NS) =                              &
                  VARMETA('mass_fraction_of_cloud_snow_water_in_air',   &
                                 'Snow concentration number','m-3',' ')
    IF (P_NG > 0) arpsmeta%qscalar(P_NG) =                              &
                  VARMETA('mass_fraction_of_cloud_graupel_water_in_air',&
                              'Graupel concentration number','m-3',' ')
    IF (P_NH > 0) arpsmeta%qscalar(P_NH) =                              &
                  VARMETA('mass_fraction_of_cloud_graupel_water_in_air',&
                                 'Hail concentration number','m-3',' ')

    IF (P_ZR > 0) arpsmeta%qscalar(P_ZR) =                              &
                  VARMETA('mass_fraction_of_rain_water_in_air',         &
                                'Rain water reflectivity','m6 m-3',' ') !!
    IF (P_ZI > 0) arpsmeta%qscalar(P_ZI) =                              &
                  VARMETA('mass_fraction_of_cloud_ice_in_air',          &
                                'Cloud ice reflectivity','m6 m-3',' ') !!
    IF (P_ZS > 0) arpsmeta%qscalar(P_ZS) =                              &
                  VARMETA('mass_fraction_of_cloud_snow_water_in_air',   &
                                     'Snow reflectivity','m6 m-3',' ') !!
    IF (P_ZG > 0) arpsmeta%qscalar(P_ZG) =                              &
                  VARMETA('mass_fraction_of_cloud_graupel_water_in_air',&
                                  'Graupel reflectivity','m6 m-3',' ') !!
    IF (P_ZH > 0) arpsmeta%qscalar(P_ZH) =                              &
                  VARMETA('mass_fraction_of_cloud_graupel_water_in_air',&
                                     'Hail reflectivity','m6 m-3',' ') !!

    arpsmeta%prcrate1 = VARMETA('lwe_large_scale_precipitation_rate',   &
                            'Grid scale precipitation rate','m s-1',' ')
    arpsmeta%prcrate2 = VARMETA('lwe_convective_precipitation_rate',    &
                               'Cumulus precipitation rate','m s-1',' ')
    arpsmeta%prcrate3 = VARMETA('lwe_microphysics_precipitation_rate',  &
                          'Microphysics precipitation rate','m s-1',' ') !!
    arpsmeta%prcrate4 = VARMETA('lwe_precipitation_rate',               &
                                 'Total precipitation rate','m s-1',' ')

    arpsmeta%raing = VARMETA('thickness_of_large_scale_rainfall_amount',&
                                   'Grid supersaturation rain','mm',' ')  !
    arpsmeta%rainc = VARMETA('thickness_of_convective_rainfall_amount', &
                                     'Cumulus convective rain','mm',' ')  !

    arpsmeta%tke = VARMETA('turbulent_kinetic_energy',                  &
                                'Turbulent Kinetic Energy','m2 s-2',' ') !!
    arpsmeta%kmh = VARMETA('turbulent_mixing_coefficient_for_horizontal_momentum', &
                   'Hori. turb. mixing coef. for momentum','m2 s-1',' ') !!
    arpsmeta%kmv = VARMETA('turbulent_mixing_coefficient_for_vertical_momentum', &
                   'Vert. turb. mixing coef. for momentum','m2 s-1',' ') !!

    arpsmeta%tsoil   = VARMETA('soil_temperature',                      &
                                             'Soil temperature','K',' ')
    arpsmeta%qsoil   = VARMETA('soil_moisture_content',                 &
                                           'Soil moisture','kg m-2',' ')
    arpsmeta%wetcanp = VARMETA('canopy_water_amount',                   &
                                     'Canopy water amount','kg m-2',' ')
    arpsmeta%snowdpth= VARMETA('surface_snow_thickness',                &
                                                   'Snow depth','m',' ')

    arpsmeta%radfrc   = VARMETA('radiation_forcing',                    &
                                          'Radiation forcing','K/s',' ') !!
    arpsmeta%radsw    = VARMETA('surface_downwelling_shortwave_flux_in_air', &
                     'Solar radiation reaching the surface','W m-2',' ')
    arpsmeta%rnflx    = VARMETA('surface_net_downward_longwave_flux',   &
                   'Net radiation flux absorbed by surface','W m-2',' ') !
    arpsmeta%radswnet = VARMETA('surface_net_downward_shortwave_flux',  &
                                      'Net solar radiation','W m-2',' ') !
    arpsmeta%radlwin  = VARMETA('surface_downwelling_longwave_flux_in_air',&
                              'Incoming longwave radiation','W m-2',' ') !

    arpsmeta%usflx  = VARMETA('surface_u_momentum_flux',                &
                          'Surface flux of u-momentum','kg m-1 s-2',' ') !!
    arpsmeta%vsflx  = VARMETA('surface_v_momentum_flux',                &
                          'Surface flux of v-momentum','kg m-1 s-2',' ') !!
    arpsmeta%ptsflx = VARMETA('surface_downward_sensible_heat_flux',    &
                                 'Surface heat flux','K kg m-1 s-2',' ') !
    arpsmeta%qvsflx = VARMETA('surface_downward_water_flux',            &
                               'Surface moisture flux','kg m-2 s-1',' ') !

!----------------------------------------------------------------------
! 
! ARPS terrain variables
!
!----------------------------------------------------------------------
    trnmeta%hterrain = VARMETA('geopotential_height',                   &
                                               'Terrain height','m',' ') !

!----------------------------------------------------------------------
! 
! ARPS surface variables
!
!----------------------------------------------------------------------

    sfcmeta%soiltyp  = arpsmeta%soiltyp
    sfcmeta%stypfrct = arpsmeta%stypfrct
    sfcmeta%vegtyp   = arpsmeta%vegtyp
    sfcmeta%lai      = arpsmeta%lai
    sfcmeta%roufns   = arpsmeta%roufns
    sfcmeta%veg      = arpsmeta%veg
    sfcmeta%ndvi     = VARMETA('ndvi',                                  &
                     'Normalized differential vegetation index','1',' ') !!

!----------------------------------------------------------------------
! 
! ARPS soil variables
!
!----------------------------------------------------------------------

    soilmeta%zpsoil   = arpsmeta%zpsoil
    soilmeta%tsoil    = arpsmeta%tsoil
    soilmeta%qsoil    = arpsmeta%qsoil
    soilmeta%wetcanp  = arpsmeta%wetcanp
    soilmeta%snowdpth = arpsmeta%snowdpth
    soilmeta%soiltyp  = arpsmeta%soiltyp
    
!----------------------------------------------------------------------
! 
! ARPS lateral boundary variables
!
!----------------------------------------------------------------------

    bdymeta%ctime = VARMETA('ctime',                                    &
                                'Data valid time','YYYYMMDD.HHMMSS',' ') !!
    bdymeta%u     = arpsmeta%u
    bdymeta%v     = arpsmeta%v
    bdymeta%w     = arpsmeta%w
    bdymeta%pt    = arpsmeta%pt
    bdymeta%p     = arpsmeta%p
    bdymeta%qv    = arpsmeta%qv

    bdymeta%qscalar(:) = arpsmeta%qscalar(:)

    initialized = .TRUE.

    RETURN
  END SUBROUTINE init_varmeta
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE netwrt_general_att             ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
  SUBROUTINE netwrt_general_att(ncid,packed,title,nx,ny,dx,dy,mapproj,  &
                                sclfct,trulat1,trulat2,trulon,          &
                                ctrlat,ctrlon,istatus)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
!    Write general global attributes to ARPS data file
!
! Note: It will be called outside of this module. So keep the interface
!       unchanged as long as possilbe.
!
!-----------------------------------------------------------------------
    IMPLICIT NONE
    INTEGER,          INTENT(IN)  :: ncid
    CHARACTER(LEN=*), INTENT(IN)  :: title
    INTEGER,          INTENT(IN)  :: packed
    INTEGER,          INTENT(IN)  :: nx, ny    ! global domain size
    REAL,             INTENT(IN)  :: dx, dy
    INTEGER,          INTENT(IN)  :: mapproj
    REAL,             INTENT(IN)  :: sclfct
    REAL,             INTENT(IN)  :: trulat1, trulat2, trulon
    REAL,             INTENT(IN)  :: ctrlat, ctrlon
    INTEGER,          INTENT(OUT) :: istatus
  
!-----------------------------------------------------------------------
!
! Misc. local variables
!
!-----------------------------------------------------------------------

    CHARACTER(LEN=80) :: tmpstr,  history
    CHARACTER(LEN=10) :: datestr, timestr, zonestr

    CHARACTER(LEN=17), PARAMETER :: fmtstr = '005.30 NetCDF 3.0'
   
    INTEGER           :: lenstr

    INTEGER           :: iproj_orig
    REAL              :: scale_orig
    REAL              :: latnot_orig(2),trlon_orig
    REAL              :: x0_orig,y0_orig

    REAL              :: trulats(2)
    REAL              :: xctr, yctr, x0,y0
    REAL              :: xno, yno, offeast, offnorth

    INCLUDE 'netcdf.inc'
    INCLUDE 'globcst.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Begining of executable code ....
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    IF (.NOT. initialized) THEN

      CALL init_varmeta(P_QC,P_QR,P_QI,P_QS,P_QG,P_QH,                  &
                        P_NC,P_NR,P_NI,P_NS,P_NG,P_NH,                  &
                        P_ZR,P_ZI,P_ZS,P_ZG,P_ZH)

    END IF

!-----------------------------------------------------------------------
!
! Define global attributes
!
!-----------------------------------------------------------------------

    istatus = NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'Conventions',9,'ARPS/ADAS')
                                           ! For Lead IDV visualization

    lenstr  = LEN_TRIM(title)
    istatus = NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'Title',lenstr,title)

    CALL DATE_AND_TIME(datestr,timestr,zonestr)
    WRITE(tmpstr,'(13a)') datestr(1:4),'-',datestr(5:6),'-',            &
                          datestr(7:8),'_',timestr(1:2),':',            &
                          timestr(3:4),':',timestr(5:10),' ',           &
                          zonestr(1:3)
    WRITE(history,'(2a)')  'Created from ARPS NetCDF API at ', TRIM(tmpstr)
    lenstr  = LEN_TRIM(history)
    istatus = NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'History',lenstr,history)

    lenstr  = LEN_TRIM(fmtstr)
    istatus = NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'FMTVER',lenstr,fmtstr)

    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'PACKED',NF_INT,1,packed)

    istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'DX',NF_FLOAT,1,dx)
    istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'DY',NF_FLOAT,1,dy)

    !
    ! Map projection
    !
  
    SELECT CASE (mapproj)
    CASE (0)
      WRITE(tmpstr,'(a)') 'no_mapping'
    CASE (1,-1)
      WRITE(tmpstr,'(a)') 'polar_stereographic'
    CASE (2,-2)
      WRITE(tmpstr,'(a)') 'lambert_conformal_conic'
    CASE (3,-3)
      WRITE(tmpstr,'(a)') 'mercator'
    CASE (4)
      WRITE(tmpstr,'(a)') 'lat_lon_projection'
    CASE (5)
      WRITE(tmpstr,'(a)') 'flat_earth_projection'
    CASE DEFAULT
      WRITE(tmpstr,'(a)') 'unknown_map_projection'
    END SELECT

    lenstr  = LEN_TRIM(tmpstr)
    istatus = NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,'grid_mapping_name',lenstr,tmpstr)

    istatus = NF_PUT_ATT_INT (ncid,NF_GLOBAL,'MAPPROJ',   NF_INT,  1,mapproj)
    istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'SCLFCT',    NF_FLOAT,1,sclfct)
    istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'TRUELAT1',  NF_FLOAT,1,trulat1)
    istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'TRUELAT2',  NF_FLOAT,1,trulat2)
    istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'TRUELON',   NF_FLOAT,1,trulon)
    istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'CTRLAT',    NF_FLOAT,1,ctrlat)
    istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'CTRLON',    NF_FLOAT,1,ctrlon)

    !-------------------------------------------------------------------
    !
    ! Used only by IDV
    ! 
    !-------------------------------------------------------------------
    !
    ! save original map projection
    !
    CALL getmapr(iproj_orig,scale_orig,latnot_orig,                     &
                 trlon_orig,x0_orig,y0_orig)

    !
    ! set ARPS map projection
    !
    trulats(1) = trulat1
    trulats(2) = trulat2
    CALL setmapr(mapproj,sclfct,trulats,trulon)
    CALL lltoxy(1,1,ctrlat,ctrlon,xctr,yctr)
    x0 = xctr - 0.5*(nx-3)*dx
    y0 = yctr - 0.5*(ny-3)*dy

    CALL lltoxy(1,1,trulat1,trulon,xno,yno)
    offeast  = xno - x0      ! false easting
    offnorth = yno - y0      ! false northing

    !
    ! Restore orginal map projection
    !
    CALL setmapr(iproj_orig,scale_orig,latnot_orig,trlon_orig)
    CALL setorig(1,x0_orig,y0_orig)

    istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'latitude_of_projection_origin', &
                              NF_FLOAT,1,trulat1)
    istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'longitude_of_central_meridian', &
                              NF_FLOAT,1,trulon)
    istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'standard_parallel',             &
                              NF_FLOAT,2,trulats)
    istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'false_easting',                 &
                              NF_FLOAT,1,offeast)
    istatus = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,'false_northing',                &
                              NF_FLOAT,1,offnorth)

    RETURN
  END SUBROUTINE netwrt_general_att
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE net_define_var_meta            ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
  SUBROUTINE net_define_var_meta(ncid,varid,vartype,vmeta)
!
!-----------------------------------------------------------------------
!
! PURPOSE:
!
!   Write ARPS variable meta data (NetCDF variable attributes)
!  
! Note: It will be called outside of this module. So keep the interface
!       unchanged as long as possilbe.
!
!-----------------------------------------------------------------------
  
    IMPLICIT NONE
  
    INTEGER,          INTENT(IN) :: ncid
    INTEGER,          INTENT(IN) :: varid
    CHARACTER(LEN=*), INTENT(IN) :: vartype
    TYPE(VARMETA),    INTENT(IN) :: vmeta
  
!-----------------------------------------------------------------------
!
! Include file
!
!-----------------------------------------------------------------------
  
    INCLUDE 'netcdf.inc'
  
!-----------------------------------------------------------------------
!
! Local variables
!
!-----------------------------------------------------------------------
  
    INTEGER :: lenstr, istatus
  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code ... ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
    lenstr = LEN_TRIM(vmeta%standard_name)
    istatus = NF_PUT_ATT_TEXT(ncid,varid,'standard_name',lenstr,vmeta%standard_name)
  
    lenstr = LEN_TRIM(vmeta%long_name)
    istatus = NF_PUT_ATT_TEXT(ncid,varid,'long_name',lenstr,vmeta%long_name)
  
    lenstr = LEN_TRIM(vmeta%units)
    istatus = NF_PUT_ATT_TEXT(ncid,varid,'units',lenstr,vmeta%units)
  
    istatus = NF_PUT_ATT_TEXT(ncid,varid,'stagger',1,vmeta%stagger)
  
    IF ( TRIM(vartype) == 'INT' ) THEN
      istatus = NF_PUT_ATT_INT (ncid,varid,'_FillValue',NF_INT,1,0)
    END IF

  !  IF ( TRIM(vartype) == 'INT' ) THEN
  !    istatus = NF_PUT_ATT_INT (ncid,varid,'scale_factor',NF_INT,1,1)
  !    istatus = NF_PUT_ATT_INT (ncid,varid,'add_offset',NF_INT,1,0)
  !  ELSE
  !    istatus = NF_PUT_ATT_REAL (ncid,varid,'scale_factor',NF_INT,1,1.0)
  !    istatus = NF_PUT_ATT_REAL (ncid,varid,'add_offset',NF_INT,1,0.0)
  !  END IF
  
    RETURN
  END SUBROUTINE net_define_var_meta

END MODULE arps_netio_metadata
