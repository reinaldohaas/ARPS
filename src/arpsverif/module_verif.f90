MODULE verif

  IMPLICIT NONE

  INTEGER,PARAMETER :: flag_p = 1, flag_sfc = 1 ! EMK
  INTEGER,PARAMETER :: nvar_p=5, nvar_sfc=5
  INTEGER,PARAMETER :: nvar_p_stats=nvar_p+2, nvar_sfc_stats=nvar_sfc+2
  INTEGER,PARAMETER :: nlevel=14
  
  REAL, PARAMETER :: pressure(nlevel) =		&
      (/ 1000.0E2, 950.0E2, 900.0E2, 850.0E2, 800.0E2, 700.0E2, 600.0E2, &
          500.0E2, 400.0E2, 300.0E2, 250.0E2, 200.0E2, 150.0E2, 100.0E2 /)

  ! Note that p/h, t, rh, u, v are in the same position for the _sfc and _p
  ! arrays. This may not be true for future variables.
  INTEGER, PARAMETER :: id_p=1, id_h=1, id_t=2, id_rh=3, id_u=4, id_v=5, &
      id_wspd=6, id_wvec=7

  CHARACTER*4,PARAMETER :: &
      varid_p(nvar_p_stats) = &
          (/ 'gpht', 't   ', 'rh  ', 'u   ', 'v   ', 'wspd', 'wvec' /), &
      varid_sfc(nvar_sfc_stats) = &
          (/ 'mslp', 't   ', 'rh  ', 'u   ', 'v   ', 'wspd', 'wvec' /)

  ! variable names

  CHARACTER*16,PARAMETER :: varname_p(nvar_p_stats) =  &
      (/ 'Gptl Height     ', &
	 'Temperature     ', &
	 'Rel Humidity    ', &
	 'U-wind          ', &
	 'V-wind          ', &
	 'Wind Speed      ', &
	 'Wind Vector     ' /)
  
  CHARACTER*16,PARAMETER :: varname_sfc(nvar_sfc_stats) =  &
      (/ 'MSL Pressure    ', &
	 'Sfc Temp        ', &
	 'Sfc Rel Hum     ', &
	 'Sfc U-wind      ', &
	 'Sfc V-wind      ', &
	 'Sfc Wind Speed  ', &
	 'Sfc Wind Vector ' /)

  ! units

  CHARACTER*3,PARAMETER :: &
      varunit_p(nvar_p_stats) = &
          (/ 'm  ', 'K  ', '%  ', 'm/s', 'm/s', 'm/s', 'm/s' /), &
      varunit_sfc(nvar_sfc_stats) = &
          (/ 'Pa ', 'K  ', '%  ', 'm/s', 'm/s', 'm/s', 'm/s' /)

END MODULE verif
