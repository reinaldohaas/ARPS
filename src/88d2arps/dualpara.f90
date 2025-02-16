!########################################################################
!########################################################################
!#########                                                      #########
!#########                  module dualpara                     #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

MODULE DUALPARA

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Declare some constants used for calculaion of dual polarization
! parameters such as Zhh, Zdr, and Kdp. (It can be expanded...)
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 12/3/2004
!
!-----------------------------------------------------------------------
! Declare parameters.
!-----------------------------------------------------------------------

  IMPLICIT NONE
  SAVE

  REAL, PARAMETER :: pi = 3.141592   ! pi

  REAL, PARAMETER :: lambda = 107.   ! wave length of radar (mm)

  REAL,PARAMETER :: Ki2 = 0.176 ! Dielectric factor for ice.
  REAL,PARAMETER :: Kw2 = 0.93 ! Dielectric factor for water.

  REAL,PARAMETER :: alphaa = 4.28e-4   ! backscattering amplitude constant
                                       ! along major axis for rain
  REAL,PARAMETER :: beta_ra = 3.04
  REAL,PARAMETER :: alphab = 4.28e-4   ! backscattering amplitude constant
                                       ! along minor axis for rain
  REAL,PARAMETER :: beta_rb = 2.77
  REAL,PARAMETER :: alphak = 1.3e-5   ! differential forward scattering
                                       ! amplitude for rain
  REAL,PARAMETER :: alphasa = 3.01e-5   ! backscattering amplitude constant
                                        ! along major axis for snow
  REAL,PARAMETER :: alphasb = 2.89e-5   ! backscattering amplitude constant
                                        ! along minor axis for snow
  REAL,PARAMETER :: alphask = 8.53e-7   ! differential forward scattering
                                        ! amplitude for snow
  REAL,PARAMETER :: betasa = 2.246      ! backscattering amplitude power
                                        ! along major axis for snow
  REAL,PARAMETER :: betasb = 2.260      ! backscattering amplitude power
                                        ! along minor axis for snow
  REAL,PARAMETER :: betask = 1.59       ! differential forward scattering
                                        ! power for snow
  REAL,PARAMETER :: alphaa_ds = 1.94e-5 ! for dry snow at horz plane
  REAL,PARAMETER :: alphab_ds = 1.91e-5 ! for dry snow at vert plane
  REAL,PARAMETER :: alphaa_dh = 1.91e-4 ! for dry hail at horz plane
  REAL,PARAMETER :: alphab_dh = 1.65e-4 ! for dry hail at vert plane

  REAL,PARAMETER :: alphak_ds = 0.03e-5 ! alphaa_ds - alphab_ds
  REAL,PARAMETER :: alphak_dh = 0.26e-4 ! alphaa_dh - alphab_dh

  REAL,PARAMETER :: fos = 0.5    ! Maximum fraction of rain-snow mixture
  REAL,PARAMETER :: foh = 0.5    ! Maximum fraction of rain-hail mixture

  REAL,PARAMETER :: rho_0r = 1.0      ! rho_0 for rain
  REAL,PARAMETER :: rho_0s = 1.0      ! rho_0 for snow
  REAL,PARAMETER :: rho_0h = 0.97     ! rho_0 for hail
  REAL,PARAMETER :: rho_0rsi = 0.82   ! lower limit of rho_0rs (rain-snow mixture)
  REAL,PARAMETER :: rho_0rsf = 0.95   ! upper limit of rho_0rs (rain-snow mixture)
  REAL,PARAMETER :: rho_0rhi = 0.85   ! lower limit of rho_0rh (rain-hail mixture)
  REAL,PARAMETER :: rho_0rhf = 0.95   ! upper limit of rho_0rh (rain-hail mixture)

  REAL,PARAMETER :: Coef = 5.6212976E+26
                    ! 4.*(lambda)**4.*gamma(7)/(pi**4.*Kw2)*(6*10**x/(pi*gamma(4)))^1.75
  REAL,PARAMETER :: kdpCoef = 1.1709e10    ! 180*lambda*1.e6*6/(pi**2)

  REAL,PARAMETER :: degKtoC=273.15 ! Conversion factor from degrees K to
                                   !   degrees C
  REAL, PARAMETER :: htmax = 275.65   ! maximum temp. for dry hail(+2.5C)
  REAL, PARAMETER :: htmin = 270.65   ! minimum temp. for wet hail(-2.5C)

  REAL, PARAMETER :: stmax = 273.15   ! maximum temp. for dry snow(+0.0C)
  REAL, PARAMETER :: stmin = 268.15   ! minimum temp. for wet snow(-5.0C)

  REAL,PARAMETER :: rhoi=917.  ! Density of ice (kg m**-3)

  REAL,PARAMETER :: sqpow=1.655     ! Power to which product rho * q
                                    !   is raised for snow.

  REAL,PARAMETER :: lg10mul=10.0 ! Log10 multiplier

  REAL,PARAMETER :: Zefact=720.0    ! Multiplier for Ze components.

  REAL,PARAMETER :: m3todBZ=1.0E+18 ! Conversion factor from m**3 to
                                    !   mm**6 m**-3.
  REAL,PARAMETER :: mm3todBZ=1.0E+9 ! Conversion factor from mm**3 to
                                    !   mm**6 m**-3.

! Parameters for temporary use for wet hail.
  REAL,PARAMETER :: approxpow=0.95 ! Approximation power for hail
                                   !   integral.
  REAL,PARAMETER :: rqhpowf=(7.0/4.0)*approxpow ! Power to which product
                                                !   rho * qh is raised.
  REAL,PARAMETER :: N0xpowf=3.0/4.0 ! Power to which N0r,N0s & N0h are
                                    !   raised.
  REAL,PARAMETER :: rhoxpowf=7.0/4.0 ! Power to which rhoh is raised.

! Missing value
  REAL,PARAMETER :: missing = -9999.0

  LOGICAL :: firstcall = .true.

!-----------------------------------------------------------------------
! Precalculated complete gamma function values
!-----------------------------------------------------------------------
  REAL,PARAMETER :: gamma7_04 = 776.1029053
  REAL,PARAMETER :: gamma7_08 = 836.7818
  REAL,PARAMETER :: gamma6_81 = 505.8403
  REAL,PARAMETER :: gamma6_38 = 232.4367523
  REAL,PARAMETER :: gamma6_54 = 309.3308
  REAL,PARAMETER :: gamma4_8 = 17.83786774
  REAL,PARAMETER :: gamma2_9 = 1.827355266
  REAL,PARAMETER :: gamma5_492 = 51.67281342
  REAL,PARAMETER :: gamma3_29 = 2.655858755
  REAL,PARAMETER :: gamma5_52 = 54.05897522
  REAL,PARAMETER :: gamma5_61 = 62.56661224
  REAL,PARAMETER :: gamma5_63 = 64.6460
  REAL,PARAMETER :: gamma2_055 = 1.024513245
  REAL,PARAMETER :: gamma4_18 = 7.556306362
  REAL,PARAMETER :: gamma2_59 = 1.418960810

!-----------------------------------------------------------------------
! Variables to can be changed by parameter retrieval
!-----------------------------------------------------------------------
  REAL :: N0r        ! Intercept parameter in 1/(m^4) for rain.
  REAL :: N0h        ! Intercept parameter in 1/(m^4) for hail.
  REAL :: N0s        ! Intercept parameter in 1/(m^4) for snow.

  REAL :: rhor=1000. ! Density of rain (kg m**-3)
  REAL :: rhoh       ! Density of hail (kg m**-3)
  REAL :: rhos       ! Density of snow (kg m**-3)

!-----------------------------------------------------------------------
! Constants
!-----------------------------------------------------------------------

  REAL :: Zerhf, Zervf, Zesf, Zesnegf, Zesposf, Zehnegf, Zehposf
  REAL :: Mfactor
  REAL :: ZerhfGamma, ZervfGamma, ZerhvfGamma
  REAL :: Zeshf, Zesvf
  REAL :: constKdpr,constKdpsZ,constKdpsR

!-----------------------------------------------------------------------
! Scattering matrix coefficient for snow
!
! phi=0.       (Mean orientation)
! sigmas=pi/9
! As=1/8*(3+4*cos(2*phi)*exp(-2*sigmas**2)+cos(4*phi)*exp(-8*sigmas**2))
! Bs=1/8*(3-4*cos(2*phi)*exp(-2*sigmas**2)+cos(4*phi)*exp(-8*sigmas**2))
! Cs=1/8*(1-cos(4*phi)*exp(-8*sigmas**2))
! Ds=1/8*(3+cos(4*phi)*exp(-8*sigmas**2))
! Cks=cos(2*phi)*exp(-2*sigmas**2)
!-----------------------------------------------------------------------

  REAL,PARAMETER :: sigmas = 0.3491
  REAL,PARAMETER :: As = 0.8140
  REAL,PARAMETER :: Bs = 0.0303
  REAL,PARAMETER :: Cs = 0.0778
  REAL,PARAMETER :: Ds = 0.4221
  REAL,PARAMETER :: Cks = 0.7837

!-----------------------------------------------------------------------
! Scattering matrix coefficient for hail
!
! phi=0.     (Mean orientation)
! sigmah=pi/3*(1-sf*fw)
! Ah=1/8*(3+4*cos(2*phi)*exp(-2*sigmah**2)+cos(4*phi)*exp(-8*sigmah**2))
! Bh=1/8*(3-4*cos(2*phi)*exp(-2*sigmah**2)+cos(4*phi)*exp(-8*sigmah**2))
! Ch=1/8*(1-cos(4*phi)*exp(-8*sigmah**2))
! Dh=1/8*(3+cos(4*phi)*exp(-8*sigmah**2))
! Ckh=cos(2*phi)*exp(-2*sigmah**2)
! 
! corresponding coefficient for dry hail: Ahd, Bhd, Chd, Dhd, Ckhd
!-----------------------------------------------------------------------

  REAL,PARAMETER :: sigmahd = 1.0472
  REAL,PARAMETER :: Ahd = 0.4308
  REAL,PARAMETER :: Bhd = 0.3192
  REAL,PARAMETER :: Chd = 0.1250
  REAL,PARAMETER :: Dhd = 0.3750
  REAL,PARAMETER :: Ckhd = 0.1116

  REAL,PARAMETER :: q_threshold = 2.e-4
  REAL :: sf
  REAL :: sigmah, Ah, Bh, Ch, Dh, Ckh

  TYPE T_obs_dual
       REAL :: T_log_ref, T_sum_ref_h, T_sum_ref_v
       REAL :: T_log_zdr, T_sum_ref_hv
  END TYPE T_obs_dual


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! SUBROUTINES AND FUNCTIONS
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  CONTAINS

  SUBROUTINE calcConstants()

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Precalculate commonly unsed constants to save computations.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/28/2005
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! For raindrop DSD from 0 to infinity
!-----------------------------------------------------------------------

  ZerhfGamma = mm3todBZ * (4. * lambda**4. * (alphaa)**2. * gamma7_08) &
      / (pi**5.77 * Kw2 * rhor ** 1.77 * (N0r*1.e-12) ** 0.77)
  ZervfGamma = mm3todBZ * (4. * lambda**4. * (alphab)**2. * gamma6_54) &
      / (pi**5.635 * Kw2 * rhor ** 1.635 * (N0r*1.e-12)**0.635)
  ZerhvfGamma = mm3todBZ * (4. * lambda**4. * alphaa*alphab * gamma6_81)&
      / (pi**5.703 * Kw2 * rhor ** 1.703 * (N0r*1.e-12)**0.703)

!-----------------------------------------------------------------------
! For raindrop DSD from 0 to Dmax at 8 mm
!-----------------------------------------------------------------------

  Zerhf = mm3todBZ * (4. * lambda**4.0 * (alphaa)**2.0) &
         / (pi**5.77 * Kw2 * rhor**1.77 * (N0r*1.e-12)**0.77)
  Zervf = mm3todBZ * (4. * lambda**4.0 * (alphab)**2.0) &
         / (pi**5.635 * Kw2 * rhor**1.635 * (N0r*1.e-12)**0.635)

!-----------------------------------------------------------------------
! For ice status hydrometeors
!
!    Zesf    The density of dry snow depends on the drop size.
!            rho = 0.07 D**(-1.1)    (g/cm**3)
!            Ref: Ryzhkov et al.(1998)
!    Zesnegf Preexisting ARPS formular for dry snot
!    Zesposf Preexisting ARPS formular for wet snow
!
!    Zehposf Preexisting ARPS formular for wet hail/graupel
!    Zehnegf Preexisting ARPS formular for dry hail/gsaupel
!-----------------------------------------------------------------------

  Zesf = 7.3e3 * gamma4_8 / (9 * Kw2 * gamma2_9)
  Zesnegf = ((m3todBZ * Zefact   * Ki2 * (rhos ** 0.25)) /  &
            ((pi ** 1.75) * Kw2 * (N0s ** 0.75) *  &
            (rhoi ** 2.)))
  Zesposf = ((m3todBZ * Zefact) /  &
            ((pi ** 1.75) * (N0s ** 0.75) *  (rhos ** 1.75)))
  Zehnegf = Ki2/Kw2 * (m3todBZ * Zefact) /  &
            (pi ** 1.75 * N0h ** 0.75 * rhoh ** 1.75)
  Zehposf = (((m3todBZ * Zefact) /  &
            ((pi ** 1.75) * (N0h ** 0.75) *  &
            (rhoh ** 1.75))) ** 0.95)
  Mfactor = (pi*0.07*1.e3/6.) * (6.e9/(0.07*pi*N0s*gamma2_9))**0.655

  Zeshf = (4. * (lambda)**4. * (alphasa)**2.              & 
          * (12.99)**1.669 * gamma5_492)                      &
          / (pi**4. * Kw2 * (N0s*1.e-3)** 0.669 * gamma3_29**1.669 )
  Zesvf = (4. * (lambda)**4. * (alphasb)**2. *            & 
          (12.99)**1.678 * gamma5_52)                         &
          / (pi**4. * Kw2 * (N0s*1.e-3)** 0.678 * gamma3_29**1.678 )

!-----------------------------------------------------------------------
! For Kdp constants
!-----------------------------------------------------------------------

  constKdpr = 180. * lambda  * alphak * gamma5_63 * 1.0e6 &
      / (pi**2.41 * (rhor) ** 1.41 * (N0r*1.e-12) ** 0.41)
  constKdpsR = 0.03 * pi * 0.0624 * N0s * gamma2_055 &
      / (lambda * (3.67e-2) ** 0.526)
  constKdpsZ = 180. * lambda  * alphask * gamma2_59 * 1.0e6 &
      / (pi**1.65 * (rhor) ** 0.65 * (N0r*1.e-24) ** (-0.35) )

  END SUBROUTINE calcConstants

  SUBROUTINE init_dsd()

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Setup default dsd values or reinialize default dsd values
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/20/2006
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Variables can be changed by parameter retrieval
!-----------------------------------------------------------------------
  N0r=8.0E+06 ! Intercept parameter in 1/(m^4) for rain.
  N0h=4.0E+04 ! Intercept parameter in 1/(m^4) for hail.
  N0s=3.0E+06 ! Intercept parameter in 1/(m^4) for snow.

  rhoh=913.  ! Density of hail (kg m**-3)
  rhos=100.  ! Density of snow (kg m**-3)

  END SUBROUTINE init_dsd

  SUBROUTINE model_dsd(n0rain,n0snow,n0hail,rhosnow,rhohail)

!-----------------------------------------------------------------------
!   
! PURPOSE: 
!
! Set dsd values to those used in the arps forecasts
! 
!-----------------------------------------------------------------------
!       
! AUTHOR:  Youngsun Jung, 1/20/2006
!       
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Define variables:
!-----------------------------------------------------------------------

  REAL :: n0rain,n0snow,n0hail,rhosnow,rhohail

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!         
! Beginning of executable code...
!         
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            
  N0r=n0rain
  N0s=n0snow
  N0h=n0hail

  rhos=rhosnow
  rhoh=rhohail

  END SUBROUTINE model_dsd

  SUBROUTINE coeff_hail(fw,qml)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Scattering matrix coefficient for hail
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/27/2007
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Define variables:
!-----------------------------------------------------------------------

  REAL, INTENT(IN) :: fw, qml

  IF(qml < q_threshold) THEN
     sf = 4*qml*1.e3
  ELSE
     sf = 0.8
  ENDIF

  sigmah=pi/3*(1-sf*fw)
  Ah=.125*(3+4*exp(-2*sigmah**2)+exp(-8*sigmah**2))
  Bh=.125*(3-4*exp(-2*sigmah**2)+exp(-8*sigmah**2))
  Ch=.125*(1-exp(-8*sigmah**2))
  Dh=.125*(3+exp(-8*sigmah**2))
  Ckh=exp(-2*sigmah**2)

  END SUBROUTINE coeff_hail

  TYPE(T_obs_dual) FUNCTION assign_Refl(var1,var2,var3,var4)
       REAL :: var1,var2,var3,var4,var5
       assign_Refl%T_sum_ref_h = var1
       assign_Refl%T_sum_ref_v = var2
       assign_Refl%T_log_zdr = var3
       assign_Refl%T_log_ref = var4
  END FUNCTION assign_Refl

  TYPE(T_obs_dual) FUNCTION init_Refl()
       init_Refl%T_sum_ref_h = 0.
       init_Refl%T_sum_ref_v = 0.
       init_Refl%T_log_zdr = missing
       init_Refl%T_log_ref = 0.
       init_Refl%T_sum_ref_hv = 0.
  END FUNCTION init_Refl

TYPE(T_obs_dual) FUNCTION rainIceRefl(qr,qs,qh,rho,flg)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine calculates the partial reflectivity factor
! of melting(wet) snow/hail at horizontal polarization
! and compute total reflectivity as a sum of those.
! The same formula used in shfactor is used with different
! alpha and beta coefficients that contain the effect of the fraction
! of water in the melting snow to take the melting layer into account.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/29/2007
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! External Function declaration
!-----------------------------------------------------------------------

  REAL, EXTERNAL :: snow_alpha_a, hail_alpha_a
  REAL, EXTERNAL :: snow_alpha_b, hail_alpha_b

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL :: qr,qs,qh,rho
  REAL :: rainIceRefl_hh,rainIceRefl_vv,rainIceRefl_hv,zdr
  REAL :: fracqrs,fracqrh,fracqs,fracqh,fms,fmh,fws,fwh,rhoms,rhomh
  REAL :: qrf,qsf,qhf
  REAL :: alphaa_ws,alphab_ws,alphaa_wh,alphab_wh
  REAL :: alphak_ws,alphak_wh
  REAL :: rainReflH,ZdrysnowH,ZwetsnowH,ZdryhailH,ZwethailH
  REAL :: rainReflV,ZdrysnowV,ZwetsnowV,ZdryhailV,ZwethailV
  REAL :: rainReflHV,ZdrysnowHV,ZwetsnowHV,ZdryhailHV,ZwethailHV
  REAL :: log_ref
  REAL :: rho_0rs,rho_0rh,temp

  INTEGER :: flg

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  qrf = 0.
  qsf = 0.
  qhf = 0.

  rainReflH = 0.
  rainReflV = 0.
  rainReflHV = 0.
  ZdrysnowH = 0.
  ZdrysnowV = 0.
  ZdrysnowHV = 0.
  ZwetsnowH = 0.
  ZwetsnowV = 0.
  ZwetsnowHV = 0.
  ZdryhailH = 0.
  ZdryhailV = 0.
  ZdryhailHV = 0.
  ZwethailH = 0.
  ZwethailV = 0.
  ZwethailHV = 0.

  rainIceRefl_hh = 0.
  rainIceRefl_vv = 0.
  rainIceRefl_hv = 0.
  zdr = missing
  log_ref = 0.

  rho_0rs = rho_0rsf
  rho_0rh = rho_0rhf

  if(qr < 0.0) qr =0.0
  if(qs < 0.0) qs =0.0
  if(qh < 0.0) qh =0.0

!-----------------------------------------------------------------------
! Calculate the fraction of water and ice.
!   qrf  pure rain water mixing ratio
!   qsf  dry snow mixing ratio
!   qhf  dry hail mixing ratio
!   fms  wet snow mixing ratio
!   fmh  wet hail mixing ratio
!   rhoms  density of wet snow (kg/m-3)
!   rhomh  density of wet hail (kg/m-3)
!-----------------------------------------------------------------------
  CALL fractionWater(qr,qs,fos,rhos,fracqrs,fracqs,fms,fws,rhoms)
  CALL fractionWater(qr,qh,foh,rhoh,fracqrh,fracqh,fmh,fwh,rhomh)

  qrf = qr - fracqrs - fracqrh
  qsf = qs - fracqs
  qhf = qh - fracqh

!-----------------------------------------------------------------------
! Calculate the matrix coefficient for hail (Ah,Bh,Ch,Ckh)
!-----------------------------------------------------------------------
  CALL coeff_hail(fwh,fmh)

!-----------------------------------------------------------------------
! Calculate alpha values
!-----------------------------------------------------------------------
  IF(fms > 0.) THEN
    alphaa_ws = snow_alpha_a(fws)
    alphab_ws = snow_alpha_b(fws)
    alphak_ws = alphaa_ws - alphab_ws
  ENDIF

  IF(fmh > 0.) THEN
    alphaa_wh = hail_alpha_a(fwh)
    alphab_wh = hail_alpha_b(fwh)
    alphak_wh = alphaa_wh - alphab_wh
  ENDIF

!-----------------------------------------------------------------------
! Calculate rho_0rs and rho_0rh
!-----------------------------------------------------------------------
  IF(flg > 2 .and. qsf > 0.) THEN
    temp=rho*qsf*1.e3
    if(temp > 1.) then
      rho_0rs = rho_0rsi
    elseif (1.e-2 < temp .AND. temp <= 1.) then
      rho_0rs = rho_0rsi - .5*log10(temp)*(rho_0rsf-rho_0rsi)
    endif
  ENDIF

  IF(flg > 2 .and. qhf > 0.) THEN
    temp=rho*qhf*1.e3
    if(temp > 1.) then
      rho_0rh = rho_0rhi
    elseif (1.e-2 < temp .AND. temp <= 1.) then
      rho_0rh = rho_0rhi - .5*log10(temp)*(rho_0rhf-rho_0rhi)
    endif
  ENDIF

!-----------------------------------------------------------------------
! Calculate reflectivity (Zhh and Zvv (and Zhv, if necessary))
!-----------------------------------------------------------------------
  IF(qsf > 0.) THEN
    CALL partialRefIce(Coef,N0s,As,Bs,Cs,alphaa_ds,       &
                       alphab_ds,rho,qsf,rhos,ZdrysnowH,ZdrysnowV)
    IF(flg > 2) THEN
      CALL partialRhoIce(Coef,N0s,Cs,Ds,alphaa_ds,        &
                       alphab_ds,rho,rho_0s,qsf,rhos,ZdrysnowHV)
    ENDIF
  ENDIF
  IF(fms > 0.) THEN
    CALL partialRefIce(Coef,N0s,As,Bs,Cs,alphaa_ws,       &
                       alphab_ws,rho,fms,rhoms,ZwetsnowH,ZwetsnowV)
    IF(flg > 2) THEN
      CALL partialRhoIce(Coef,N0s,Cs,Ds,alphaa_ws,        &
                       alphab_ws,rho,rho_0rs,fms,rhoms,ZwetsnowHV)
    ENDIF
  ENDIF

  IF(qhf > 0.) THEN
    CALL partialRefIce(Coef,N0h,Ahd,Bhd,Chd,alphaa_dh,    &
                       alphab_dh,rho,qhf,rhoh,ZdryhailH,ZdryhailV)
    IF(flg > 2) THEN
      CALL partialRhoIce(Coef,N0h,Chd,Dhd,alphaa_dh,      &
                       alphab_dh,rho,rho_0h,qhf,rhoh,ZdryhailHV)
    ENDIF
  ENDIF
  IF(fmh > 0.) THEN
    CALL partialRefIce(Coef,N0h,Ah,Bh,Ch,alphaa_wh,       &
                       alphab_wh,rho,fmh,rhomh,ZwethailH,ZwethailV)
    IF(flg > 2) THEN
      CALL partialRhoIce(Coef,N0h,Ch,Dh,alphaa_wh,        &
                       alphab_wh,rho,rho_0rh,fmh,rhomh,ZwethailHV)
    ENDIF
  ENDIF

  IF(qrf > 0. ) THEN
    rainReflH = ZerhfGamma*(qrf*rho)**1.77
    IF(flg > 2) THEN
      rainReflHV = ZerhvfGamma*(qrf*rho)**1.703
    ENDIF
  ENDIF

  rainIceRefl_hh=rainReflH+ZdrysnowH+ZwetsnowH+ZdryhailH+ZwethailH
  log_ref = 10.*LOG10(MAX(1.0,rainIceRefl_hh))

  IF(flg == 1) THEN
    rainIceRefl = assign_Refl(rainIceRefl_hh,rainIceRefl_vv,zdr,log_ref)

  ELSE IF(flg > 1) THEN
!-----------------------------------------------------------------------
! Calculate differential reflectivity (Zdr)
!-----------------------------------------------------------------------
    IF(qrf > 0. ) THEN
      rainReflV = ZervfGamma*(qrf*rho)**1.635
      rainReflV = MIN(rainReflV,rainReflH)
    ENDIF
    rainIceRefl_vv=rainReflV+ZdrysnowV+ZwetsnowV+ZdryhailV+ZwethailV

    if(rainIceRefl_vv > 0.) then
      zdr = 10.*LOG10(MAX(1.0,rainIceRefl_hh/rainIceRefl_vv))
    endif

    rainIceRefl = assign_Refl(rainIceRefl_hh,rainIceRefl_vv,zdr,log_ref)
  
    IF(flg > 2) THEN 
      rainIceRefl_hv=rainReflHV+ZdrysnowHV+ZwetsnowHV                  &
                   +ZdryhailHV+ZwethailHV 
      rainIceRefl%T_sum_ref_hv = rainIceRefl_hv
    ENDIF
  ENDIF

END FUNCTION rainIceRefl

REAL FUNCTION rainIceKdp(qr,qs,qh,rho)
  
!-----------------------------------------------------------------------
! 
! PURPOSE:
!
! This subroutine calculates specific differential phase. 
! 
!-----------------------------------------------------------------------
! 
! AUTHOR:  Youngsun Jung, 1/29/2007
! 
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------
  
  IMPLICIT NONE 
  
!-----------------------------------------------------------------------
! External Function declaration 
!-----------------------------------------------------------------------
  
  REAL, EXTERNAL :: snow_alpha_a, hail_alpha_a
  REAL, EXTERNAL :: snow_alpha_b, hail_alpha_b

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL :: qr,qs,qh,rho,kdp
  REAL :: fracqrs,fracqrh,fracqs,fracqh,fms,fmh,fws,fwh,rhoms,rhomh
  REAL :: qrf,qsf,qhf
  REAL :: alphaa_ws,alphab_ws,alphaa_wh,alphab_wh
  REAL :: alphak_ws,alphak_wh
  REAL :: rainKdp,drysnowKdp,wetsnowKdp,dryhailKdp,wethailKdp

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  qrf = 0.
  qsf = 0.
  qhf = 0.

  drysnowKdp = 0.
  wetsnowKdp = 0.
  dryhailKdp = 0.
  wethailKdp = 0.
  rainKdp = 0.
  rainIceKdp = 0.

  if(qr < 0.0) qr =0.0
  if(qs < 0.0) qs =0.0
  if(qh < 0.0) qh =0.0

!-----------------------------------------------------------------------
! Calculate the faction of water and ice.
! For variable names, see FUNCTION rainIceRefl
!-----------------------------------------------------------------------
  CALL fractionWater(qr,qs,fos,rhos,fracqrs,fracqs,fms,fws,rhoms)
  CALL fractionWater(qr,qh,foh,rhoh,fracqrh,fracqh,fmh,fwh,rhomh)

  qrf = qr - fracqrs - fracqrh
  qsf = qs - fracqs
  qhf = qh - fracqh

!-----------------------------------------------------------------------
! Calculate the matrix coefficient for hail (Ah,Bh,Ch,Ckh)
!-----------------------------------------------------------------------
  CALL coeff_hail(fwh,fmh)

!-----------------------------------------------------------------------
! Calculate alpha values
!-----------------------------------------------------------------------
  IF(fms > 0.) THEN
    alphaa_ws = snow_alpha_a(fws)
    alphab_ws = snow_alpha_b(fws)
    alphak_ws = alphaa_ws - alphab_ws
  ENDIF

  IF(fmh > 0.) THEN
    alphaa_wh = hail_alpha_a(fwh)
    alphab_wh = hail_alpha_b(fwh)
    alphak_wh = alphaa_wh - alphab_wh
  ENDIF

!-----------------------------------------------------------------------
! Calculate specific differential phase (Kdp)
!-----------------------------------------------------------------------
  IF(qsf > 0.) THEN
    CALL partialKdpIce(kdpCoef,Cks,alphak_ds,rho,qsf,rhos,drysnowKdp)
  ENDIF

  IF(fms > 0.) THEN
    CALL partialKdpIce(kdpCoef,Cks,alphak_ws,rho,fms,rhoms,wetsnowKdp)
  ENDIF

  IF(qhf > 0.) THEN
    CALL partialKdpIce(kdpCoef,Ckhd,alphak_dh,rho,qhf,rhoh,dryhailKdp)
  ENDIF

  IF(fmh > 0.) THEN
    CALL partialKdpIce(kdpCoef,Ckh,alphak_wh,rho,fmh,rhomh,wethailKdp)
  ENDIF

  IF(qrf > 0.) THEN
    rainKdp = constKdpr*(qrf*rho)**1.41
  ENDIF

  rainIceKdp=rainKdp+drysnowKdp+wetsnowKdp+dryhailKdp+wethailKdp

END FUNCTION rainIceKdp

!########################################################################
!########################################################################
!#########                                                      #########
!#########              FUNCTION calculate_obs                  #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

TYPE(T_obs_dual) FUNCTION calculate_obs(rho,qr,qs,qh,t,flg,drysnowopt)

!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 2/27/2007
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! External Function declaration
!-----------------------------------------------------------------------

  REAL, EXTERNAL :: rhfactorCG, rvfactorCG
  REAL, EXTERNAL :: hfactorNeg, hfactorPos
  REAL, EXTERNAL :: sfactorPos, sfactorNeg, sfactor,shfactor,svfactor
  REAL, EXTERNAL :: rainSnowRefl_hh,rainSnowRefl_vv
  REAL, EXTERNAL :: combine

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------

  REAL, INTENT(IN) :: rho ! Air density (kg m**-3)
  REAL, INTENT(IN) :: qr  ! Rain mixing ratio (kg kg**-1)
  REAL, INTENT(IN) :: qs  ! Snow mixing ratio (kg kg**-1)
  REAL, INTENT(IN) :: qh  ! Hail mixing ratio (kg kg**-1)
  REAL, INTENT(IN) :: t   ! Temperature (K)
  INTEGER, INTENT(IN) :: flg   ! flag for ref(1) and zdr(2)
  INTEGER, INTENT(IN) :: drysnowopt

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL :: rhcomp,rvcomp,hcomp,shcomp,svcomp,sumhcomp,sumvcomp
  REAL :: tempo1, tempo2, tempo3

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  rhcomp = 0.0
  rvcomp = 0.0
  hcomp = 0.0
  shcomp = 0.0
  svcomp = 0.0
  sumhcomp = 0.0
  sumvcomp = 0.0
  tempo1 = 0.0
  tempo2 = 0.0
  tempo3 = 0.0
  calculate_obs = init_Refl()

!-----------------------------------------------------------------------
! Check for bad air density value.
!-----------------------------------------------------------------------

  IF (rho > 0.0) THEN        ! outer if starts

!-----------------------------------------------------------------------
! Calculate reflectivity contribution from rain.
!-----------------------------------------------------------------------

    IF (qr > 0.0 .and. drysnowopt < 4) THEN
      rhcomp = rhfactorCG(qr,rho,ZerhfGamma)
      IF(flg == 2) THEN
        rvcomp = rvfactorCG(qr,rho,ZervfGamma)
        rvcomp = MIN(rvcomp, rhcomp)
      ENDIF
    END IF

!-----------------------------------------------------------------------
! Calculate reflectivity contribution from hail.
!-----------------------------------------------------------------------

    IF (qh > 0.0 .and. drysnowopt < 4) THEN
      IF (t <= htmax .and. t >= htmin) THEN
         tempo1 = hfactorNeg(qh,rho,Zehnegf)
         tempo2 = hfactorPos(qh,rho,rqhpowf,Zehposf)
         hcomp = combine(htmax, htmin, t, tempo1, tempo2)
      ELSE IF (t > htmax) THEN
         tempo2 = hfactorPos(qh,rho,rqhpowf,Zehposf)
         hcomp = tempo2
      ELSE IF (t < htmin) THEN
         tempo1 = hfactorNeg(qh,rho,Zehnegf)
         hcomp = tempo1
      ENDIF
    END IF

    tempo1 = 0.0
    tempo2 = 0.0
!-----------------------------------------------------------------------
! Calculate reflectivity contribution from snow.
!-----------------------------------------------------------------------

    IF (qs > 0.0 .and. drysnowopt < 4) THEN
      IF (t <= stmax .and. t >= stmin) THEN
         tempo3 = sfactorPos(qs,rho,Zesposf)
         SELECT CASE (drysnowopt)
         CASE (1)
           tempo1 = sfactorNeg(qs,rho,Zesnegf)
           tempo2 = tempo1
         CASE (2)
           tempo1 = sfactor(Mfactor,qs,rho,sqpow,Zesf)
           tempo2 = tempo1
         CASE (3)
           tempo1 = shfactor(qs,rho,Zeshf)
           IF(flg == 2) THEN
             tempo2 = svfactor(qs,rho,Zesvf)
           ENDIF
         END SELECT
         shcomp = combine(stmax, stmin, t, tempo1, tempo3)
         IF(flg == 2) THEN
           svcomp = combine(stmax, stmin, t, tempo2, tempo3)
         ENDIF
      ELSE IF (t > stmax) THEN
         shcomp = sfactorPos(qs,rho,Zesposf)
         svcomp = shcomp
      ELSE IF (t < stmin) THEN
         SELECT CASE (drysnowopt)
         CASE (1)
           shcomp = sfactorNeg(qs,rho,Zesnegf)
           svcomp = shcomp
         CASE (2)
           shcomp = sfactor(Mfactor,qs,rho,sqpow,Zesf)
           svcomp = shcomp
         CASE (3)
           shcomp = shfactor(qs,rho,Zeshf)
           IF(flg == 2) THEN
             svcomp = svfactor(qs,rho,Zesvf)
           ENDIF
         END SELECT
      ENDIF
      IF(flg == 2) THEN
        svcomp = MIN(svcomp, shcomp)
      ENDIF
    END IF

    IF (drysnowopt < 4) THEN
      calculate_obs%T_sum_ref_h = rhcomp + shcomp + hcomp
      calculate_obs%T_sum_ref_v = rvcomp + svcomp + hcomp
      calculate_obs%T_log_ref =                                   &
             lg10mul * LOG10(MAX(1.0,calculate_obs%T_sum_ref_h))

      IF (flg == 2 .AND. sumvcomp > 0.0) THEN
        calculate_obs%T_log_zdr =                                 &
             lg10mul * LOG10(MAX(1.0,sumhcomp)/MAX(1.0,sumvcomp))
      ENDIF
    ENDIF

    IF ((qr > 0. .or. qs > 0. .or. qh > 0.) &
         .and. drysnowopt == 4) THEN
      calculate_obs = rainIceRefl(qr,qs,qh,rho,flg)
    END IF

  END IF                             ! outer if ends

END FUNCTION  calculate_obs

END MODULE DUALPARA
