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

  REAL :: lambda           ! wavelength of radar (mm)

  REAL,PARAMETER :: Kw2 = 0.93 ! Dielectric factor for water.

  REAL,PARAMETER :: alphaa = 4.28e-4   ! backscattering amplitude constant
                                       ! along major axis for rain
  REAL,PARAMETER :: beta_ra = 3.04
  REAL,PARAMETER :: alphab = 4.28e-4   ! backscattering amplitude constant
                                       ! along minor axis for rain
  REAL,PARAMETER :: beta_rb = 2.77
  REAL,PARAMETER :: alphak = 3.88e-4   ! differential forward scattering
                                       ! amplitude for rain
  REAL,PARAMETER :: alphask = 8.53e-7   ! differential forward scattering
                                        ! amplitude for snow
  REAL,PARAMETER :: alphaa_ds = 1.94e-5 ! for dry snow at horz plane
  REAL,PARAMETER :: alphab_ds = 1.91e-5 ! for dry snow at vert plane
  REAL,PARAMETER :: alphaa_dh = 1.91e-4 ! for dry hail at horz plane
  REAL,PARAMETER :: alphab_dh = 1.65e-4 ! for dry hail at vert plane
  REAL,PARAMETER :: alphaa_dg = 0.81e-4 ! for dry graupel at horz plane
  REAL,PARAMETER :: alphab_dg = 0.76e-4 ! for dry graupel at vert plane

  REAL,PARAMETER :: alphak_ds = 0.03e-5 ! alphaa_ds - alphab_ds
  REAL,PARAMETER :: alphak_dh = 0.26e-4 ! alphaa_dh - alphab_dh
  REAL,PARAMETER :: alphak_dg = 0.05e-4 ! alphaa_dh - alphab_dh

  REAL,PARAMETER :: rho_0r = 1.0      ! rho_0 for rain
  REAL,PARAMETER :: rho_0s = 1.0      ! rho_0 for snow
  REAL,PARAMETER :: rho_0h = 0.97     ! rho_0 for hail
  REAL,PARAMETER :: rho_0g = 0.40     ! rho_0 for hail
  REAL,PARAMETER :: rho_0rsi = 0.82   ! lower limit of rho_0rs (rain-snow mixture)
  REAL,PARAMETER :: rho_0rsf = 0.95   ! upper limit of rho_0rs (rain-snow mixture)
  REAL,PARAMETER :: rho_0rhi = 0.85   ! lower limit of rho_0rh (rain-hail mixture)
  REAL,PARAMETER :: rho_0rhf = 0.95   ! upper limit of rho_0rh (rain-hail mixture)
  REAL,PARAMETER :: rho_0rgi = 0.82   ! lower limit of rho_0rg (rain-graupel mixture)
  REAL,PARAMETER :: rho_0rgf = 0.95   ! upper limit of rho_0rg (rain-graupel mixture)

  REAL,PARAMETER :: Coef = 5.6212976E+26
                    ! 4.*(lambda)**4.*gamma(7)/(pi**4.*Kw2)*(6*10**x/(pi*gamma(4)))^1.75
  REAL,PARAMETER :: kdpCoef = 1.1709e10    ! 180*lambda*1.e6*6/(pi**2)

  REAL,PARAMETER :: degKtoC=273.15 ! Conversion factor from degrees K to
                                   !   degrees C

  REAL,PARAMETER :: rhoi=917.  ! Density of ice (kg m**-3)

  REAL,PARAMETER :: mm3todBZ=1.0E+9 ! Conversion factor from mm**3 to
                                    !   mm**6 m**-3.

! Missing value
  REAL,PARAMETER :: missing = -9999.0
  REAL :: grpl_miss
  REAL :: hl_miss 
  
  LOGICAL :: firstcall = .true.
  INTEGER :: grpl_ON
  INTEGER :: hl_ON 
  INTEGER :: qgh_opt

  INTEGER :: attn_ON

!-----------------------------------------------------------------------
! Precalculated complete gamma function values
!-----------------------------------------------------------------------
  REAL,PARAMETER :: gamma7_08 = 836.7818
  REAL,PARAMETER :: gamma6_81 = 505.8403
  REAL,PARAMETER :: gamma6_54 = 309.3308
  REAL,PARAMETER :: gamma5_63 = 64.6460
  REAL,PARAMETER :: gamma4_16 = 7.3619
  REAL,PARAMETER :: gamma3_97 = 5.7788

!-----------------------------------------------------------------------
! Variables to can be changed by parameter retrieval
!-----------------------------------------------------------------------
  REAL :: N0r        ! Intercept parameter in 1/(m^4) for rain
  REAL :: N0s        ! Intercept parameter in 1/(m^4) for snow
  REAL :: N0h        ! Intercept parameter in 1/(m^4) for hail
  REAL :: N0g        ! Intercept parameter in 1/(m^4) for hail

  REAL :: rhor=1000. ! Density of rain (kg m**-3)
  REAL :: rhoh       ! Density of hail (kg m**-3)
  REAL :: rhos       ! Density of snow (kg m**-3)
  REAL :: rhog       ! Density of graupel (kg m**-3)

!-----------------------------------------------------------------------
! Variables to can be changed for meling ice
!-----------------------------------------------------------------------
  REAL :: fos        ! Maximum fraction of rain-snow mixture
  REAL :: foh        ! Maximum fraction of rain-hail mixture
  REAL :: fog        ! Maximum fraction of rain-hail mixture

!-----------------------------------------------------------------------
! Constants
!-----------------------------------------------------------------------

  REAL :: ZerhfGamma, ZervfGamma, ZerhvfGamma
  REAL :: constKdpr

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

!-----------------------------------------------------------------------
! Scattering matrix coefficient for graupel
! 
! phi=0.     (Mean orientation)
! sigmag=pi/3*(1-sf*fw)
! Ag=1/8*(3+4*cos(2*phi)*exp(-2*sigmag**2)+cos(4*phi)*exp(-8*sigmag**2))
! Bg=1/8*(3-4*cos(2*phi)*exp(-2*sigmag**2)+cos(4*phi)*exp(-8*sigmag**2))
! Cg=1/8*(1-cos(4*phi)*exp(-8*sigmag**2))
! Dg=1/8*(3+cos(4*phi)*exp(-8*sigmag**2))
! Ckg=cos(2*phi)*exp(-2*sigmag**2)
! 
! corresponding coefficient for dry graupel: Agd, Bgd, Cgd, Dgd, Ckgd
!-----------------------------------------------------------------------
  
  REAL,PARAMETER :: sigmagd = 1.0472
  REAL,PARAMETER :: Agd = 0.4308
  REAL,PARAMETER :: Bgd = 0.3192
  REAL,PARAMETER :: Cgd = 0.1250
  REAL,PARAMETER :: Dgd = 0.3750
  REAL,PARAMETER :: Ckgd = 0.1116
  
  REAL :: sigmag, Ag, Bg, Cg, Dg, Ckg

!-----------------------------------------------------------------------
!  Declare new observation type
!-----------------------------------------------------------------------

! TYPE T_obs_dual
  REAL :: T_log_ref, T_sum_ref_h, T_sum_ref_v
  REAL :: T_log_zdr, T_sum_ref_hv, T_kdp
  REAL :: T_Ahh,     T_Avv
! END TYPE T_obs_dual

!-----------------------------------------------------------------------
!  Declare new DSD parameter data type
!-----------------------------------------------------------------------

! TYPE T_para_dsd
  REAL :: T_qr, T_qs, T_qh, T_qg
  REAL :: T_Ntr, T_Nts, T_Nth, T_Ntg
  REAL*8 :: T_alfr,T_alfs,T_alfh,T_alfg
! END TYPE T_para_dsd
! DTD: add additional shape parameter mu
  REAL*8 :: T_mur, T_mus, T_mug, T_muh
! DTD: add mixing ratios of liquid water on snow, graupel and hail (for the ZVD scheme)
  REAL :: T_qsw, T_qgw, T_qhw

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! SUBROUTINES AND FUNCTIONS
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  CONTAINS

  SUBROUTINE setgrplhl(graupel_ON, hail_ON)
!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Set graupel and hail parameters
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Dan Dawson, 02/11/2011
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER :: graupel_ON, hail_ON

  grpl_ON = graupel_ON
  hl_ON = hail_ON

  RETURN

  END SUBROUTINE setgrplhl

  SUBROUTINE init_fox()

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Setup default maximum fraction of water in the melting ice.
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
! Variables can vary depend on whether graupel/hail exists. 
!-----------------------------------------------------------------------
  fos = 0.25             ! Maximum fraction of rain-snow mixture
  foh = 0.2              ! Maximum fraction of rain-hail mixture
  fog = 0.25             ! Maximum fraction of rain-hail mixture

  END SUBROUTINE init_fox

  SUBROUTINE init_fox_no_grpl()

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Setup default maximum fraction of water in the melting ice 
! when graupel is suppressed.
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
  fos = 0.5              ! Maximum fraction of rain-snow mixture
  foh = 0.3              ! Maximum fraction of rain-hail mixture
  fog = 0.0              ! Maximum fraction of rain-hail mixture
      
  END SUBROUTINE init_fox_no_grpl

  SUBROUTINE init_fox_no_hail() 

!-----------------------------------------------------------------------
!
! PURPOSE:  
!
!  Setup default maximum fraction of water in the melting ice 
!  when hail is suprressed. 
!
!-----------------------------------------------------------------------
!
! AUTHOR: Bryan Putnam, 12/14/10
!
!-----------------------------------------------------------------------
! Force explicit declarations. 
!-----------------------------------------------------------------------

  IMPLICIT NONE 

!-----------------------------------------------------------------------
! Variables can be changed by parameter retrieval 
!-----------------------------------------------------------------------

  fos = 0.5             ! Maximum fraction of rain-snow mixture
  foh = 0.0             ! Maximum fraction of rain-hail mixture
  fog = 0.3             ! Maximum fraction of rain-hail mixture

  END SUBROUTINE init_fox_no_hail

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
  N0g=4.0E+05 ! Intercept parameter in 1/(m^4) for graupel.

  rhoh=913.  ! Density of hail (kg m**-3)
  rhos=100.  ! Density of snow (kg m**-3)
  rhog=400.  ! Density of graupel (kg m**-3)

  END SUBROUTINE init_dsd

  SUBROUTINE model_dsd(n0rain,n0snow,n0hail,n0grpl,rhosnow,rhohail,rhogrpl)

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

  REAL :: n0rain,n0snow,n0hail,n0grpl,rhosnow,rhohail,rhogrpl

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  N0r=n0rain
  N0s=n0snow
  N0h=n0hail
  N0g=n0grpl

  rhos=rhosnow
  rhoh=rhohail
  rhog=rhogrpl

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

! TYPE(T_obs_dual) FUNCTION assign_Refl(var1,var2,var3,var4)
!      REAL :: var1,var2,var3,var4,var5
!      assign_Refl%T_sum_ref_h = var1
!      assign_Refl%T_sum_ref_v = var2
!      assign_Refl%T_log_zdr = var3
!      assign_Refl%T_log_ref = var4
! END FUNCTION assign_Refl

! TYPE(T_obs_dual) FUNCTION init_Refl()
!      init_Refl%T_sum_ref_h = 0.
!      init_Refl%T_sum_ref_v = 0.
!      init_Refl%T_log_zdr = missing
!      init_Refl%T_log_ref = 0.
!      init_Refl%T_sum_ref_hv = 0.
!      init_Refl%T_kdp = 0.
!      init_Refl%T_Ahh = 0.
!      init_Refl%T_Avv = 0.
! END FUNCTION init_Refl

! TYPE(T_para_dsd) FUNCTION init_para_dsd()
!   init_para_dsd%T_qr = 0.0
!   init_para_dsd%T_qs = 0.0
!   init_para_dsd%T_qh = 0.0
!   init_para_dsd%T_qg = 0.0
!   init_para_dsd%T_Ntr = 0.0
!   init_para_dsd%T_Nts = 0.0
!   init_para_dsd%T_Nth = 0.0
!   init_para_dsd%T_Ntg = 0.0
!   init_para_dsd%T_alfr = 0.d0
!   init_para_dsd%T_alfs = 0.d0
!   init_para_dsd%T_alfh = 0.d0
!   init_para_dsd%T_alfg = 0.d0
! END FUNCTION init_para_dsd

! TYPE(T_para_dsd) FUNCTION assign_para_dsd_SM(var1,var2,var3,var4)
!   REAL :: var1,var2,var3,var4

!   assign_para_dsd_SM%T_qr = var1
!   assign_para_dsd_SM%T_qs = var2
!   assign_para_dsd_SM%T_qh = var3
!   assign_para_dsd_SM%T_qg = var4
!   assign_para_dsd_SM%T_Ntr = 0.0 
!   assign_para_dsd_SM%T_Nts = 0.0 
!   assign_para_dsd_SM%T_Nth = 0.0 
!   assign_para_dsd_SM%T_Ntg = 0.0 
!   assign_para_dsd_SM%T_alfr = 0.d0
!   assign_para_dsd_SM%T_alfs = 0.d0
!   assign_para_dsd_SM%T_alfh = 0.d0
!   assign_para_dsd_SM%T_alfg = 0.d0
! END FUNCTION assign_para_dsd_SM

! TYPE(T_para_dsd) FUNCTION assign_para_dsd_TM(var1,var2,var3,var4, &
!                           var5,var6,var7,var8,var9,var10,var11,var12)
!   REAL :: var1,var2,var3,var4,var5,var6,var7,var8
!   REAL*8 :: var9,var10,var11,var12

!   assign_para_dsd_TM%T_qr = var1
!   assign_para_dsd_TM%T_qs = var2
!   assign_para_dsd_TM%T_qh = var3
!   assign_para_dsd_TM%T_qg = var4
!   assign_para_dsd_TM%T_Ntr = var5
!   assign_para_dsd_TM%T_Nts = var6
!   assign_para_dsd_TM%T_Nth = var7
!   assign_para_dsd_TM%T_Ntg = var8
!   assign_para_dsd_TM%T_alfr = var9
!   assign_para_dsd_TM%T_alfs = var10
!   assign_para_dsd_TM%T_alfh = var11
!   assign_para_dsd_TM%T_alfg = var12
! END FUNCTION assign_para_dsd_TM

SUBROUTINE cal_N0(rhoa,q,Ntx,rhox,alpha,N0,mu)

!
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates intercept parameter and "effective" intercept parameter
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson
!  (02/06/2008)
!
!  MODIFICATION HISTORY:
!
!  (03/26/2008)
!  Recast N0 as a double precision variable, and used double precision for
!  all intermediate calculations.  The calling subroutine should
!  also define it as double precision.  For situations with large alpha,
!  N0 can become very large, and loss of precision can result.
!  Also tweaked the calculation of N0 a bit to avoid overflow, in keeping
!  With Jason Milbrandt's calculation of N0 just before evaporation in
!  the multi-moment code.
!
!  (02/14/2011)
!  Added additional shape parameter mu, i.e. for a gamma distribution of
!  the form N(D) = N0*D^alpha*exp(-lambda*D^3mu).  Setting mu to 1/3
!  retrieves the original formulation for the standard gamma distribution.
!  This value assumes that the mass diameter relationship for the hydrometeor
!  distribution in question is of the form m(D) = c*D^3 and won't be 
!  correct for other values of the exponent.
!
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!

  REAL, PARAMETER :: pi = 3.141592   ! pi
  REAL :: rhoa,q,Ntx
  REAL*8 :: N0,alpha,mu
!JYS  REAL*8 :: N0_eff
  REAL :: rhox
  REAL*8 :: gamma1, gamma4

  REAL*8 :: gamma

  DOUBLE PRECISION :: lamda

  gamma1 = gamma((1.d0+alpha)/(3.d0*mu))
  gamma4 = gamma((4.d0+alpha)/(3.d0*mu))

  IF(rhoa > 0.0 .and. q > 0.0) THEN
    lamda = ((gamma4/gamma1)*dble(pi/6.*rhox)*dble(Ntx)/(dble(rhoa)*  &
        dble(q)))**mu
  ELSE
    lamda = 0.d0
  END IF

  N0 = 3*mu*dble(Ntx)*lamda**(0.5d0*((1.d0+alpha)/(3*mu)))*                         &
              (1.d0/gamma1)*lamda**(0.5d0*((1.d0+alpha))/(3*mu))

!JYS  IF(lamda /= 0.d0) THEN
!JYS    N0_eff = N0*(((4.d0+alpha)/lamda)**alpha)*gamma4*(128.d0/3.d0)/   &
!JYS                ((4.d0+alpha)**(4.d0+alpha))
!JYS  ELSE
!JYS    N0_eff = 0.d0
!JYS  END IF

END SUBROUTINE cal_N0

SUBROUTINE cal_lamda(rhoa,q,Ntx,rhox,alpha,lamda,mu)

!
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates slope parameter lamda
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson!  (02/06/2008)
!
!  MODIFICATION HISTORY:
!  (03/31/2008)
!  Converted intermediate calculations and arrays alpha and lamda to
!  double precision.
!
!  (02/14/2011)
!  Added additional shape parameter mu, i.e. for a gamma distribution of
!  the form N(D) = N0*D^alpha*exp(-lambda*D^3mu).  Setting mu to 1/3
!  retrieves the original formulation for the standard gamma distribution.
!  This value assumes that the mass diameter relationship for the hydrometeor
!  distribution in question is of the form m(D) = c*D^3 and won't be 
!  correct for other values of the exponent.
!
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!

  REAL, PARAMETER :: pi = 3.141592   ! pi
  REAL :: rhoa,q
  REAL*8 :: lamda,alpha,mu
  REAL :: rhox
  REAL :: Ntx
  REAL*8 :: gamma1, gamma4

  REAL*8 :: gamma

  gamma1 = gamma((1.d0+alpha)/(3.d0*mu))
  gamma4 = gamma((4.d0+alpha)/(3.d0*mu))

  IF(rhoa > 0.0 .and. q > 0.0) THEN
    lamda = ((gamma4/gamma1)*dble(pi/6.*rhox)*dble(Ntx)/(dble(rhoa)*  &
          dble(q)))**mu
  ELSE
    lamda = 0.d0
  END IF

END SUBROUTINE cal_lamda

SUBROUTINE cal_Nt(rhoa,q,N0,cx,alpha,Ntx,mu)

!   
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates number concentration at scalar points
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson 
!  (02/06/2008) 
!   
!  MODIFICATION HISTORY:
!  
!  03/31/08 - converted intermediate calculations to double precision
!             as well as a few of the input arguments.
!
!  (02/15/2011)
!  Added additional shape parameter mu, i.e. for a gamma distribution of
!  the form N(D) = N0*D^alpha*exp(-lambda*D^3mu).  Setting mu to 1/3
!  retrieves the original formulation for the standard gamma distribution.
!  This value assumes that the mass diameter relationship for the hydrometeor
!  distribution in question is of the form m(D) = c*D^3 and won't be 
!  correct for other values of the exponent.
! 
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!     

  REAL :: rhoa,q
  REAL*8 :: N0,alpha,mu
  REAL :: cx
  REAL :: Ntx
  REAL*8 :: gamma1,gamma4

  REAL*8 :: gamma

  gamma1 = gamma((1.d0+alpha)/(3.d0*mu))
  gamma4 = gamma((4.d0+alpha)/(3.d0*mu))

  Ntx = sngl((N0*gamma1/(3.d0*mu))**(3.d0/(4.d0+alpha))*   &
                   ((gamma1/gamma4)*dble(rhoa)* &
                   dble(q)/dble(cx))**((1.d0+alpha)/(4.d0+alpha)))

END SUBROUTINE cal_Nt

SUBROUTINE solve_alpha(nx,ny,nz,rhoa,cx,q,Ntx,Z,alpha)
!
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates shape parameter alpha
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson
!  (02/06/2008)
!
!  MODIFICATION HISTORY:
!  (03/31/2008)
!  Changed alpha array to double precision
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!

  REAL, PARAMETER :: pi = 3.141592   ! pi
  INTEGER :: nx,ny,nz
  REAL    :: rhoa(nx,ny,nz),q(nx,ny,nz),Ntx(nx,ny,nz),Z(nx,ny,nz)
  REAL*8  :: alpha(nx,ny,nz)

  REAL*8  :: solveAlpha

  REAL :: cx

  INTEGER i,j,k

  DO k=2,nz-2
    DO j=2,ny-2
      DO i=2,nx-2
        IF(q(i,j,k) > 0.0 .and. Ntx(i,j,k) > 0.0 .and. Z(i,j,k) > 0.0) THEN

          alpha(i,j,k) = solveAlpha(q(i,j,k),Ntx(i,j,k),Z(i,j,k),cx,rhoa(i,j,k))

        ELSE
          alpha(i,j,k) = 0.d0
        END IF
      END DO
    END DO
  END DO

END SUBROUTINE solve_alpha


SUBROUTINE fractionWater(qr,qi,fo,density_ice,fracqr,fracqi,fm,fw,rhom)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine calculates the fractions of water, dry ice (snow or
! hail), the mixture. It also calculate the density of mixture.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 5/30/2006
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare variables.
!-----------------------------------------------------------------------

  REAL :: qr, qi, fo, density_ice
  REAL,INTENT(OUT) :: fracqr, fracqi, fm, fw, rhom

  REAL :: fr

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  fr = 0.
  fw = 0.
  fracqr = 0.
  fracqi = 0.
  fm = 0.
  rhom = 0.

!-----------------------------------------------------------------------
! Calculate the fraction of mleting ice (fr) based on the ratio between
! qr and qi. fo is the maximum allowable fraction of melting snow.
!-----------------------------------------------------------------------
  IF (qr > 0. .AND. qi > 0.) THEN
    fr = fo*(MIN(qi/qr,qr/qi))**.3
  ENDIF

!-----------------------------------------------------------------------
! Calculate the faction of water and ice.
! fracqr : the mass of water in the melting ice
! fracqi : the mass of ice in the melting ice
! fm     : total mass of melting ice
! fw     : the fraction of water within melting ice
! rhom   : density of mixture
!-----------------------------------------------------------------------
  fracqr = fr * qr
  fracqi = fr * qi
  fm = fracqr + fracqi

  IF (fm .EQ. 0. .AND. qr > 0.) THEN
    fw = 1.
  ELSE IF (fm > 0.) THEN
    fw = fracqr/fm
  ENDIF

  rhom = 1000.*fw**2. + (1.-fw**2.)*density_ice

END SUBROUTINE fractionWater

!########################################################################
!########################################################################
!#########                                                      #########
!#########              FUNCTION refl_rsa                       #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

!SUBROUTINE refl_rsa (rhoa,T_qr,T_qs,T_qh,T_qg,                 &
!           T_Ntr,T_Nts,T_Nth,T_Ntg,T_alfr,T_alfs,T_alfh,       &
!           T_alfg,T_log_ref,T_sum_ref_h,T_sum_ref_v, &
!           T_log_zdr,T_sum_ref_hv,T_kdp,T_Ahh,T_Avv)
SUBROUTINE refl_rsa (MFflg,rhoa)
!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine calculates specific differential phase.
! Compute radar observations by integrating radar scattering amplitudes
! over the drop size range. Radar scattering amplitues were calculated
! using T-matrix method and stored in the table.
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 2/27/2007
!
! MODIFIED: Dan Dawson, 2/08/2011
!           Changed from a function to a subroutine.
!           Put derived type variables explicitly in argument list.
!           Removed references to derived type and instead refer 
!           explicitly to the variables within
!
!           Dan Dawson, 2/14/2011
!           Modified DSD integration to allow for gamma distribution of
!           the form n(D) = N0*D^alpha*exp(-lambda*D^(3*mu)),
!           that is, an additional shape parameter, mu, as in the ZVD
!           scheme. Also included the case (through the MFflg variable)
!           where the melted water fraction on ice is explicitly predicted
!           (the ZVDM scheme variants).  
!    
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------
  USE rsa_table
  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------

  INTEGER :: MFflg
  REAL, PARAMETER :: pi4 = 97.409           ! pi^4
  REAL*8, DIMENSION (nd) :: Ndr, Nds, Ndh, Ndg, Ndrs, Ndrh, Ndrg

  REAL*8, EXTERNAL :: gamma
  INTEGER, EXTERNAL :: get_qgh_opt

  REAL :: rhoa,cx
  REAL :: qr,qs,qh,qg
  REAL :: ntr,nts,nth,ntg
  REAL*8 :: alfr,alfs,alfh,alfg
  !DTD add mu shape parameter
  REAL*8 :: mur,mus,mug,muh
  !DTD add liquid water fraction variables
  REAL*8 :: qsw,qgw,qhw
  REAL*8 :: db_N0,Ntw,Ntd
  REAL :: fracqrs,fracqrh,fracqrg
  REAL :: fracqs,fracqh,fracqg
  REAL :: fms,fmh,fmg
  REAL :: fws,fwh,fwg
  REAL :: rhoms,rhomh,rhomg
  REAL :: qrf,qsf,qhf,qgf
  REAL*8 :: lamr,lams,lamrs,lamh,lamrh,lamg,lamrg
  REAL :: tsar_h,tsas_h,tsah_h,tsag_h,tsars_h,tsarh_h,tsarg_h
  REAL :: tsar_v,tsas_v,tsah_v,tsag_v,tsars_v,tsarh_v,tsarg_v
  COMPLEX :: tsar_hv,tsas_hv,tsah_hv,tsag_hv,tsars_hv,tsarh_hv,tsarg_hv
  REAL :: tfsar,tfsas,tfsah,tfsag,tfsars,tfsarh,tfsarg
  REAL :: D,intv,far,lambda4
  REAL :: fa2,fb2,temph,tempv,temphv,tempk,temp
  COMPLEX :: fab,fba,fconj
  INTEGER :: i, idx
  REAL :: tempAhh,tempAvv
  REAL :: Ar_h,As_h,Ag_h,Ah_h,Ars_h,Arg_h,Arh_h
  REAL :: Ar_v,As_v,Ag_v,Ah_v,Ars_v,Arg_v,Arh_v

!  TYPE(T_para_dsd) :: var_dsd

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!-----------------------------------------------------------------------
! Initialization
!-----------------------------------------------------------------------

  ! DTD: Set attenuation to off

  attn_ON = 0

  IF(firstcall) THEN
    qgh_opt = get_qgh_opt(grpl_ON,hl_ON)

    SELECT CASE (qgh_opt)
     CASE (1)
       fos = 0.5; foh = 0.0; fog = 0.0
     CASE (2)
       CALL init_fox_no_grpl()
     CASE (3)
       CALL init_fox_no_hail()
     CASE (4)
       CALL init_fox()
    END SELECT

    firstcall = .false. 
  END IF

  lambda4 = lambda**4.

! qr = var_dsd%T_qr
! qs = var_dsd%T_qs
! qh = var_dsd%T_qh
! qg = var_dsd%T_qg
! ntr = var_dsd%T_Ntr
! nts = var_dsd%T_Nts
! nth = var_dsd%T_Nth
! ntg = var_dsd%T_Ntg
! alfr = var_dsd%T_alfr
! alfs = var_dsd%T_alfs
! alfh = var_dsd%T_alfh
! alfg = var_dsd%T_alfg

  qr = T_qr
  qs = T_qs
  qh = T_qh
  qg = T_qg
  ntr = T_Ntr
  nts = T_Nts
  nth = T_Nth
  ntg = T_Ntg
  alfr = T_alfr
  alfs = T_alfs
  alfh = T_alfh
  alfg = T_alfg

  mur = T_mur  !
  mus = T_mus  ! New shape parameters
  mug = T_mug  !
  muh = T_muh  !

  qsw = T_qsw  !
  qgw = T_qgw  ! New variables holding liquid water fraction
  qhw = T_qhw  ! mixing ratios

  
  
  qrf = 0.; qsf = 0.; qhf = 0.; qgf = 0.
  Ndr = 0.; Nds = 0.; Ndh = 0.; Ndg = 0.
  Ndrs = 0.; Ndrh = 0.; Ndrg = 0.
  temph = 0.; tempv = 0.; temphv = 0.; temp = 0.; tempk = 0.
  tempAhh = 0.; tempAvv = 0.
  fracqrs = 0.; fracqs = 0.; fms = 0.; fws = 0.; rhoms = 100.
  fracqrh = 0.; fracqh = 0.; fmh = 0.; fwh = 0.; rhomh = 913.
  fracqrg = 0.; fracqg = 0.; fmg = 0.; fwg = 0.; rhomg = 400.
  !refl_rsa = init_Refl()
  T_sum_ref_h = 0.
  T_sum_ref_v = 0.
  T_log_zdr = missing
  T_log_ref = 0.
  T_sum_ref_hv = 0.
  T_kdp = 0.
  T_Ahh = 0.
  T_Avv = 0.

  if(qr < 0.0) qr =0.0
  if(qs < 0.0) qs =0.0
  if(qh < 0.0) qh =0.0
  if(qg < 0.0) qg =0.0

!-----------------------------------------------------------------------
! Calculate the fraction of water and ice.
! For the variable definition, see "FUNCTION rainIceRefl".
! DTD: In the case of the ZVDM scheme variants, the liquid water mixing
! ratio on ice is explicitly predicted
!-----------------------------------------------------------------------

  IF (MFflg == 0) THEN       ! Melted fraction not predicted: diagnose it!
                             ! Or maybe we just want to compute it this way anyway,
                             ! for comparison purposes.

    CALL fractionWater(qr,qs,fos,rhos,fracqrs,fracqs,fms,fws,rhoms)
    IF(hl_ON == 1)  &
    CALL fractionWater(qr,qh,foh,rhoh,fracqrh,fracqh,fmh,fwh,rhomh)
    IF(grpl_ON == 1) &
    CALL fractionWater(qr,qg,fog,rhog,fracqrg,fracqg,fmg,fwg,rhomg)

    qrf = qr - fracqrs - fracqrh - fracqrg
    if(qrf < 0.0) qrf = 0.0

    qsf = qs - fracqs
    if(qsf < 0.0) qsf = 0.0
    qhf = qh - fracqh
    if(qhf < 0.0) qhf = 0.0
    qgf = qg - fracqg
    if(qgf < 0.0) qgf = 0.0
  ELSE                          ! Melted fraction explicitly predicted
   
    ! Note, the handling of this case is somewhat different than the above
    ! Above, there is an explicit mixing ratio of a rain/snow mixture = fms
    ! while qsf represents the mixing ratio of the nonmelting dry snow.  In the ZVD
    ! case, it is assumed that if melting is occuring, all of the snow
    ! at a given point is melting: thus, fms = qs, and qsf = 0.
    ! Similarly for the other ice hydrometeors.
    ! For rain, since the water fraction on ice is carried separately,
    ! qrf = qr

    qrf = qr

    IF(qsw > 0.0) THEN
      qsf = 0.0
      fms = qs
      fracqrs = qsw
      fracqs = qs-qsw
      IF(qs > 0.0) THEN
        fws = qsw/qs
      ELSE 
        fws = 0.0
      END IF
      ! For now just assume rhoms is the same as original dry snow (i.e. constant)
      rhoms = rhos        
    ELSE
      qsf = qs
      fms = 0.0
      fracqrs = 0.0
      fracqs = qs
      fws = 0.0
      rhoms = rhos
    END IF
    IF(qgw > 0.0) THEN
      qgf = 0.0
      fmg = qg
      fracqrg = qgw
      fracqg = qg-qgw
      IF(qg > 0.0) THEN
        fwg = qgw/qg
      ELSE
        fwg = 0.0
      END IF
      ! The (variable) density of graupel should already take into account melting
      rhomg = rhog
    ELSE
      qgf = qg
      fmg = 0.0
      fracqrg = 0.0
      fracqg = qg
      fwg = 0.0
      rhomg = rhog
    END IF
    IF(qhf > 0.0) THEN
      qhf = 0.0
      fmh = qh
      fracqrh = qhw
      fracqh = qh-qhw
      IF(qh > 0.0) THEN
        fwh = qhw/qh
      ELSE
        fwh = 0.0
      END IF
      ! The (variable) density of hail should already take into account melting
      rhomh = rhoh
    ELSE
      qhf = qh
      fmh = 0.0
      fracqrh = 0.0
      fracqh = qh
      fwh = 0.0
      rhomh = rhoh
    END IF
  END IF

!-----------------------------------------------------------------------
! Calculate the matrix coefficient for hail (Ah,Bh,Ch,Ckh)
!-----------------------------------------------------------------------
  IF(hl_ON == 1) CALL coeff_hail(fwh,fmh)
  IF(grpl_ON == 1) CALL coeff_hail(fwg,fmg)

!-----------------------------------------------------------------------
! Calculate slopes of the invers exponential curves of N(d)
!-----------------------------------------------------------------------
  lamr = 0.d0; lams = 0.d0; lamh = 0.d0; lamrs = 0.d0; lamrh = 0.d0;
  db_N0 = 0.d0

  if(qrf > 0.) then
    intv = (dsr(2) - dsr(1))*1.e-3
    Ntw = 0.
    if(ntr > 0.0) then
      CALL cal_N0(rhoa,qrf,ntr,rhor,alfr,db_N0,mur)
      CALL cal_lamda(rhoa,qrf,ntr,rhor,alfr,lamr,mur)
      Ntw = ntr
    else
      !print*,'Here rain!'
      db_N0 = dble(N0r)
      !print*,'N0r = ',N0r
      !lamr = ((pi*rhor*N0r) / (qrf*rhoa))**0.25
      cx = rhor*(pi/6.)
      !print*,'rhor = ',rhor 
      CALL cal_Nt(rhoa,qrf,db_N0,cx,alfr,ntr,mur)  
      CALL cal_lamda(rhoa,qrf,ntr,rhor,alfr,lamr,mur)
      Ntw = ntr
    endif
    !Ntw = db_N0/lamr
    Do i = 1,nd
      Ndr(i) = db_N0*(dsr(i)*1.e-3)**alfr*exp(-lamr*(dsr(i)*1.e-3)**(3.0*mur))*intv
      Ntw = Ntw - Ndr(i)
      if(Ntw <= 0) Ndr(i) = 0.
    ENDDO
  endif

  db_N0 = 0.d0

  if(qsf > 0.) then
    intv = (dss(2) - dss(1))*1.e-3
    Ntd = 0.
    if(nts > 0.0) then
      CALL cal_N0(rhoa,qsf,nts,rhos,alfs,db_N0,mus)
      CALL cal_lamda(rhoa,qsf,nts,rhos,alfs,lams,mus)
      Ntd = nts
    else
      !print*,'Here snow!'
      db_N0 = dble(N0s)
      cx = rhos*(pi/6.)
      !lams = ((pi*rhos*N0s) / (qsf*rhoa))**0.25
      CALL cal_Nt(rhoa,qsf,db_N0,cx,alfs,nts,mus)
      CALL cal_lamda(rhoa,qsf,nts,rhos,alfs,lams,mus)
      Ntd = nts
    endif
    !Ntd = db_N0/lams
    Do i = 1,nd
      Nds(i) = db_N0*(dss(i)*1.e-3)**alfs*exp(-lams*(dss(i)*1.e-3)**(3.0*mus))*intv
      Ntd = Ntd - Nds(i)
      if(Ntd <= 0) Nds(i) = 0.
    ENDDO
  endif

  db_N0 = 0.d0

  if(fms > 0.) then
    intv = (dss(2) - dss(1))*1.e-3
    Ntw = 0.
    if(nts > 0.0) then
      CALL cal_N0(rhoa,fms,nts,rhoms,alfs,db_N0,mus)
      CALL cal_lamda(rhoa,fms,nts,rhoms,alfs,lamrs,mus)
      Ntw = nts
    else
      !print*,'Here snow!'
      db_N0 = dble(N0s)
      cx = rhoms*(pi/6.)
      !lamrs = ((pi*rhoms*N0s) / (fms*rhoa))**0.25
      CALL cal_Nt(rhoa,fms,db_N0,cx,alfs,nts,mus)
      CALL cal_lamda(rhoa,fms,nts,rhoms,alfs,lamrs,mus)
      Ntw = nts
    endif
    !Ntw = db_N0/lamrs
    Do i = 1,nd
      Ndrs(i) = db_N0*(dss(i)*1.e-3)**alfs*exp(-lamrs*(dss(i)*1.e-3)**(3.0*mus))*intv
      Ntw = Ntw - Ndrs(i)
      if(Ntw <= 0) Ndrs(i) = 0.
    ENDDO
  endif

  IF(hl_ON == 1) THEN  
    db_N0 = 0.d0

    if(qhf > 0.) then
      intv = (dsh(2) - dsh(1))*1.e-3
      Ntd = 0.
      if(nth > 0.0) then
        CALL cal_N0(rhoa,qhf,nth,rhoh,alfh,db_N0,muh)
        CALL cal_lamda(rhoa,qhf,nth,rhoh,alfh,lamh,muh)
        Ntd = nth
      else
        !print*,'Here hail!'
        db_N0 = dble(N0h)
        cx = rhoh*(pi/6.)
        !lamh = ((pi*rhoh*N0h) / (qhf*rhoa))**0.25
        CALL cal_Nt(rhoa,qhf,db_N0,cx,alfh,nth,muh)
        CALL cal_lamda(rhoa,qhf,nth,rhoh,alfh,lamh,muh)
        Ntd = nth
      endif
      !Ntd = db_N0/lamh
      Do i = 1,nd
        Ndh(i) = db_N0*(dsh(i)*1.e-3)**alfh*exp(-lamh*(dsh(i)*1.e-3)**(3.0*muh))*intv
        Ntd = Ntd - Ndh(i)
        if(Ntd <= 0) Ndh(i) = 0.
      ENDDO
    endif

    db_N0 = 0.d0

    if(fmh > 0.) then
      intv = (dsh(2) - dsh(1))*1.e-3
      Ntw = 0.
      if(nth > 0.0) then
        CALL cal_N0(rhoa,fmh,nth,rhomh,alfh,db_N0,muh)
        CALL cal_lamda(rhoa,fmh,nth,rhomh,alfh,lamrh,muh)
        Ntw = nth
      else
        !print*,'Here hail!'
        db_N0 = dble(N0h)
        cx = rhomh*(pi/6.)
        !lamrh = ((pi*rhomh*N0h) / (fmh*rhoa))**0.25
        CALL cal_Nt(rhoa,fmh,db_N0,cx,alfh,nth,muh)
        CALL cal_lamda(rhoa,fmh,nth,rhomh,alfh,lamrh,muh)
        Ntw = nth
      endif
      !Ntw = db_N0/lamrh
      Do i = 1,nd
        Ndrh(i) = db_N0*(dsh(i)*1.e-3)**alfh*exp(-lamrh*(dsh(i)*1.e-3)**(3.0*muh))*intv
        Ntw = Ntw - Ndrh(i)
        if(Ntw <= 0) Ndrh(i) = 0.
      ENDDO
    endif
  ENDIF 

  IF(grpl_ON == 1) THEN
    db_N0 = 0.d0

    if(qgf > 0.) then
      intv = (dsg(2) - dsg(1))*1.e-3
      Ntd = 0.
      if(ntg > 0.0) then
        CALL cal_N0(rhoa,qgf,ntg,rhog,alfg,db_N0,mug)
        CALL cal_lamda(rhoa,qgf,ntg,rhog,alfg,lamg,mug)
        Ntd = ntg
      else
        !print*,'Here graupel!'
        db_N0 = dble(N0g)
        !print*,'n0g = ',db_N0
        !print*,'rhog = ',rhog
        !print*,'qg = ',qgf
        !print*,'alphag = ',alfg
        !print*,'mug = ',mug
        cx = rhog*(pi/6.)
        !lamg = ((pi*rhog*N0g) / (qgf*rhoa))**0.25
        CALL cal_Nt(rhoa,qgf,db_N0,cx,alfg,ntg,mug)
        CALL cal_lamda(rhoa,fmg,ntg,rhog,alfg,lamg,mug)
        Ntd = ntg
      endif
      !Ntd = db_N0/lamg
      Do i = 1,nd
        Ndg(i) = db_N0*(dsg(i)*1.e-3)**alfg*exp(-lamg*(dsg(i)*1.e-3)**(3.0*mug))*intv
        Ntd = Ntd - Ndg(i)
        if(Ntd <= 0) Ndg(i) = 0.
      ENDDO
    endif

    db_N0 = 0.d0

    if(fmg > 0.) then
      intv = (dsg(2) - dsg(1))*1.e-3
      Ntw = 0.
      if(ntg > 0.0) then
        CALL cal_N0(rhoa,fmg,ntg,rhomg,alfg,db_N0,mug)
        CALL cal_lamda(rhoa,fmg,ntg,rhomg,alfg,lamrg,mug)
        Ntw = ntg
      else
        !print*,'Here graupel!'
        db_N0 = dble(N0g)
        !print*,'n0g = ',db_N0
        !print*,'rhomg = ',rhomg
        !print*,'fmg = ',fmg
        !print*,'alphag = ',alfg
        !print*,'mug = ',mug
        cx = rhomg*(pi/6.)
        !lamrg = ((pi*rhomg*N0g) / (fmg*rhoa))**0.25
        CALL cal_Nt(rhoa,fmg,db_N0,cx,alfg,ntg,mug)
        !print*,'Here after cal_Nt!'
        !print*,'Ntg = ',ntg
        CALL cal_lamda(rhoa,fmg,ntg,rhomg,alfg,lamrg,mug)
        !print*,'Here after cal_lamda'
        !print*,'lamrg = ',lamrg
        Ntw = ntg
      endif
      !Ntw = db_N0/lamrg
      Do i = 1,nd
        Ndrg(i) = db_N0*(dsg(i)*1.e-3)**alfg*exp(-lamrg*(dsg(i)*1.e-3)**(3.0*mug))*intv
        Ntw = Ntw - Ndrg(i)
        if(Ntw <= 0) Ndrg(i) = 0.
      ENDDO
    endif
  ENDIF

!-----------------------------------------------------------------------
! Calculate radar observations.
!-----------------------------------------------------------------------
  tsar_h=0.; tsas_h=0.; tsah_h=0.; tsag_h=0.
  tsars_h=0.; tsarh_h=0.; tsarg_h=0.
  tsar_v=0.; tsas_v=0.; tsah_v=0.; tsag_v=0.
  tsars_v=0.; tsarh_v=0.; tsarg_v=0.
  tsar_hv=0.; tsas_hv=0.; tsah_hv=0.; tsag_hv = 0.
  tsars_hv=0.; tsarh_hv=0.; tsarg_hv=0.
  tfsar=0.; tfsas=0.; tfsah=0.; tfsag=0.
  tfsars=0.; tfsarh=0.; tfsarg=0.
  fa2=0.; fb2=0.; fab=0.; far=0.
  Ar_h=0.; As_h=0.; Ag_h=0.; Ah_h=0.
  Ars_h=0.; Arg_h=0.; Arh_h=0.
  Ar_v=0.; As_v=0.; Ag_v=0.; Ah_v=0.
  Ars_v=0.; Arg_v=0.; Arh_v=0.

  !print*,'Here rain 2!'
  if(qrf > 0.) then
    do i=1,nd
      fa2 = ABS(far_b(i))**2
      fb2 = ABS(fbr_b(i))**2
      fab = far_b(i)*CONJG(fbr_b(i))
      fba = fbr_b(i)*CONJG(far_b(i))
      tsar_h = tsar_h + fa2*Ndr(i)
      tsar_v = tsar_v + fb2*Ndr(i)
      tsar_hv = tsar_hv + fab*Ndr(i)
      far=REAL(far_f(i) - fbr_f(i))
      tfsar = tfsar + far*Ndr(i)
      IF(attn_ON == 2) THEN
        Ar_h = Ar_h + AIMAG(far_f(i)*Ndr(i))
        Ar_v = Ar_v + AIMAG(fbr_f(i)*Ndr(i))
      ENDIF
    enddo
  endif

  fa2=0.; fb2=0.; fab=0.; far=0.

  !print*,'Here snow 2!'
  if(qsf > 0.) then
    do i=1,nd
      fa2 = ABS(fas_b(i,1))**2
      fb2 = ABS(fbs_b(i,1))**2
      fab = fas_b(i,1)*CONJG(fbs_b(i,1))
      fba = fbs_b(i,1)*CONJG(fas_b(i,1))
      far=(As*fa2 + Bs*fb2 + 2*Cs*REAL(fab))
      tsas_h = tsas_h + far*Nds(i)
      far=(Bs*fa2 + As*fb2 + 2*Cs*REAL(fab))
      tsas_v = tsas_v + far*Nds(i)
      fconj=(Cs*(fa2+fb2)+As*fab+Bs*fba)
      tsas_hv = tsas_hv + fconj*Nds(i)
      far=Cks*REAL(fas_f(i,1) - fbs_f(i,1))
      tfsas = tfsas + far*Nds(i)
      IF(attn_ON == 2) THEN
        As_h = As_h + AIMAG(fas_f(i,1)*Nds(i))
        As_v = As_v + AIMAG(fbs_f(i,1)*Nds(i))
      ENDIF
    enddo
  endif

  fa2=0.; fb2=0.; fab=0.; far=0.

  !print*,'Here snow 3!'
  if(fms > 0.) then
    idx = INT(fws * 20 + 0.5) + 1
    do i=1,nd
      fa2 = ABS(fas_b(i,idx))**2
      fb2 = ABS(fbs_b(i,idx))**2
      fab = fas_b(i,idx)*CONJG(fbs_b(i,idx))
      fba = fbs_b(i,idx)*CONJG(fas_b(i,idx))
      far=(As*fa2 + Bs*fb2 + 2*Cs*REAL(fab))
      tsars_h = tsars_h + far*Ndrs(i)
      far=(Bs*fa2 + As*fb2 + 2*Cs*REAL(fab))
      tsars_v = tsars_v + far*Ndrs(i)
      fconj=(Cs*(fa2+fb2)+As*fab+Bs*fba)
      tsars_hv = tsars_hv + fconj*Ndrs(i)
      far=Cks*REAL(fas_f(i,idx)-fbs_f(i,idx))
      tfsars = tfsars + far*Ndrs(i)
      IF(attn_ON == 2) THEN
        Ars_h = Ars_h + AIMAG(fas_f(i,idx)*Ndrs(i))
        Ars_v = Ars_v + AIMAG(fbs_f(i,idx)*Ndrs(i))
      ENDIF
    enddo
  endif

  IF(hl_ON == 1) THEN
    fa2=0.; fb2=0.; fab=0.; far=0.

    !print*,'Here hail 2!'
    if(qhf > 0.) then
      do i=1,nd
        fa2 = ABS(fah_b(i,1))**2
        fb2 = ABS(fbh_b(i,1))**2
        fab = fah_b(i,1)*CONJG(fbh_b(i,1))
        fba = fbh_b(i,1)*CONJG(fah_b(i,1))
        far=(Ahd*fa2 + Bhd*fb2 + 2*Chd*REAL(fab))
        tsah_h = tsah_h + far*Ndh(i)
        far=(Bhd*fa2 + Ahd*fb2 + 2*Chd*REAL(fab))
        tsah_v = tsah_v + far*Ndh(i)
        fconj=(Chd*(fa2+fb2)+Ahd*fab+Bhd*fba)
        tsah_hv = tsah_hv + fconj*Ndh(i)
        far=Ckhd*REAL(fah_f(i,1) - fbh_f(i,1))
        tfsah = tfsah + far*Ndh(i)
        IF(attn_ON == 2) THEN
          Ah_h = Ah_h + AIMAG(fah_f(i,1)*Ndh(i))
          Ah_v = Ah_v + AIMAG(fbh_f(i,1)*Ndh(i))
        ENDIF
      enddo
    endif

    fa2=0.; fb2=0.; fab=0.; far=0.
    !print*,'Here hail 3!'
    if(fmh > 0.) then
      idx = INT(fwh * 20 + 0.5) + 1
      do i=1,nd
        fa2 = ABS(fah_b(i,idx))**2
        fb2 = ABS(fbh_b(i,idx))**2
        fab = fah_b(i,idx)*CONJG(fbh_b(i,idx))
        fba = fbh_b(i,idx)*CONJG(fah_b(i,idx))
        far=(Ah*fa2 + Bh*fb2 + 2*Ch*REAL(fab))
        tsarh_h = tsarh_h + far*Ndrh(i)
        far=(Bh*fa2 + Ah*fb2 + 2*Ch*REAL(fab))
        tsarh_v = tsarh_v + far*Ndrh(i)
        fconj=(Ch*(fa2+fb2)+Ah*fab+Bh*fba)
        tsarh_hv = tsarh_hv + fconj*Ndrh(i)
        far=Ckh*REAL(fah_f(i,idx)-fbh_f(i,idx))
        tfsarh = tfsarh + far*Ndrh(i)
        IF(attn_ON == 2) THEN
          Arh_h = Arh_h + AIMAG(fah_f(i,idx)*Ndrh(i))
          Arh_v = Arh_v + AIMAG(fbh_f(i,idx)*Ndrh(i))
        ENDIF
      enddo
    endif
  ENDIF 

  IF(grpl_ON == 1) THEN
    fa2=0.; fb2=0.; fab=0.; far=0.

    !print*,'Here graupel 2!'
    if(qgf > 0.) then
      do i=1,nd
        fa2 = ABS(fag_b(i,1))**2
        fb2 = ABS(fbg_b(i,1))**2
        fab = fag_b(i,1)*CONJG(fbg_b(i,1))
        fba = fbg_b(i,1)*CONJG(fag_b(i,1))
        far=(Agd*fa2 + Bgd*fb2 + 2*Cgd*REAL(fab))
        tsag_h = tsag_h + far*Ndg(i)
        far=(Bgd*fa2 + Agd*fb2 + 2*Cgd*REAL(fab))
        tsag_v = tsag_v + far*Ndg(i)
        fconj=(Cgd*(fa2+fb2)+Agd*fab+Bgd*fba)
        tsag_hv = tsag_hv + fconj*Ndg(i)
        far=Ckgd*REAL(fag_f(i,1) - fbg_f(i,1))
        tfsag = tfsag + far*Ndg(i)
        IF(attn_ON == 2) THEN
          Ag_h = Ag_h + AIMAG(fag_f(i,1)*Ndg(i))
          Ag_v = Ag_v + AIMAG(fbg_f(i,1)*Ndg(i))
        ENDIF
      enddo
    endif

    fa2=0.; fb2=0.; fab=0.; far=0.

    !print*,'Here graupel 3!'
    if(fmg > 0.) then
      idx = INT(fwg * 20 + 0.5) + 1
      do i=1,nd
        fa2 = ABS(fag_b(i,idx))**2
        fb2 = ABS(fbg_b(i,idx))**2
        fab = fag_b(i,idx)*CONJG(fbg_b(i,idx))
        fba = fbg_b(i,idx)*CONJG(fag_b(i,idx))
        far=(Ag*fa2 + Bg*fb2 + 2*Cg*REAL(fab))
        tsarg_h = tsarg_h + far*Ndrg(i)
        far=(Bg*fa2 + Ag*fb2 + 2*Cg*REAL(fab))
        tsarg_v = tsarg_v + far*Ndrg(i)
        fconj=(Cg*(fa2+fb2)+Ag*fab+Bg*fba)
        tsarg_hv = tsarg_hv + fconj*Ndrg(i)
        far=Ckg*REAL(fag_f(i,idx)-fbg_f(i,idx))
        tfsarg = tfsarg + far*Ndrg(i)
        IF(attn_ON == 2) THEN
          Arg_h = Arg_h + AIMAG(fag_f(i,idx)*Ndrg(i))
          Arg_v = Arg_v + AIMAG(fbg_f(i,idx)*Ndrg(i))
        ENDIF
      enddo
    endif
  ENDIF

  temph = 4*lambda4/(pi4*Kw2)*(tsar_h+tsas_h+tsah_h+tsag_h+tsars_h+tsarh_h+tsarg_h)
  tempv = 4*lambda4/(pi4*Kw2)*(tsar_v+tsas_v+tsah_v+tsag_v+tsars_v+tsarh_v+tsarg_v)
  !refl_rsa%T_sum_ref_h = temph
  !refl_rsa%T_sum_ref_v = tempv
  T_sum_ref_h = temph
  T_sum_ref_v = tempv

  !print*,'temph',temph
  !print*,'tempv',tempv

  IF(attn_ON == 2) THEN
    tempAhh = 2.*lambda*(Ar_h+As_h+Ag_h+Ah_h+Ars_h+Arg_h+Arh_h)*1.e-3
    tempAvv = 2.*lambda*(Ar_v+As_v+Ag_v+Ah_v+Ars_v+Arg_v+Arh_v)*1.e-3
    !refl_rsa%T_Ahh = tempAhh
    !refl_rsa%T_Avv = tempAvv
    T_Ahh = tempAhh
    T_Avv = tempAvv
  ENDIF

  temphv = 4*lambda4/(pi4*Kw2)*ABS(tsar_hv+tsas_hv+tsah_hv+tsag_hv+tsars_hv+tsarh_hv+tsarg_hv) 
  !refl_rsa%T_sum_ref_hv = temphv
  T_sum_ref_hv = temphv

  if(temph > 0.) T_log_ref = 10*log10(temph)

  if(tempv > 0.) then
!    refl_rsa%T_log_zdr = 10.*LOG10(MAX(1.0,temph/tempv))
    T_log_zdr = 10.*LOG10(MAX(1.0,temph/tempv))
    !print*,'T_log_zdr',T_log_zdr
  endif

!JYS  if(tempk < 0.) tempk = 0.0
  tempk = 180.*lambda/pi*(tfsar+tfsas+tfsah+tfsag+tfsars+tfsarh+tfsarg)*1.e-3
  !refl_rsa%T_kdp = tempk
  T_kdp = tempk
END SUBROUTINE refl_rsa

END MODULE DUALPARA

SUBROUTINE refl_rsa_array(ZVDflg,MFflg,nx,ny,nz,rsafndir,wavelen,rhoa_arr,qr,qs,qg,qh,          &
           Ntr,Nts,Ntg,Nth,alphar,alphas,alphag,alphah,qsw,qgw,qhw,rhogrpl,rhohail,                &
           logZ,sumZh,sumZv,logZdr,sumZhv,Kdp,Ahh,Avv)
!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This is a wrapper subroutine for the subroutine refl_rsa, to allow iteration
! over a 3D array, and to accomodate direct passing of microphysics variables
! through the argument list
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Dan Dawson, 2/10/2011
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------
  USE rsa_table
  USE dualpara
  IMPLICIT NONE

  INTEGER :: ZVDflg,MFflg
  INTEGER :: nx,ny,nz
  CHARACTER (LEN=256) :: rsafndir
  REAL, INTENT(IN) :: wavelen
  REAL, INTENT(IN) :: rhoa_arr(nx,ny,nz)
  REAL, INTENT(IN) :: qr(nx,ny,nz),qs(nx,ny,nz),qg(nx,ny,nz),qh(nx,ny,nz)
  REAL, INTENT(IN) :: Ntr(nx,ny,nz),Nts(nx,ny,nz),Ntg(nx,ny,nz),Nth(nx,ny,nz)
  REAL, INTENT(IN) :: alphar(nx,ny,nz),alphas(nx,ny,nz),alphag(nx,ny,nz),alphah(nx,ny,nz)
  REAL, INTENT(IN) :: qsw(nx,ny,nz),qgw(nx,ny,nz),qhw(nx,ny,nz)
  REAL, INTENT(IN) :: rhogrpl(nx,ny,nz),rhohail(nx,ny,nz)
  REAL, INTENT(OUT) :: logZ(nx,ny,nz),sumZh(nx,ny,nz),sumZv(nx,ny,nz),logZdr(nx,ny,nz)
  REAL, INTENT(OUT) :: sumZhv(nx,ny,nz),Kdp(nx,ny,nz),Ahh(nx,ny,nz),Avv(nx,ny,nz)

  REAL :: rhoa
  INTEGER :: i,j,k

! Initialize radar wavelength

  lambda = wavelen

! Initialize DSD parameters

  !CALL init_dsd()

! DTD: Read RSA table here

  CALL read_table(rsafndir)

  ! DEBUG: far_b before loop

 !DO k=1,100
 !  print*,'far_b = ',far_b(k)
 !END DO

  ! Assume that mu is constant for each category, and assign it here based on the value of ZVDflg

  IF(ZVDflg == 1) THEN ! We are using the ZVD scheme, which has an additional shape parameter mu for
                       ! the gamma distribution (setting it to 1/3 simplifies it to the standard gamma)
    T_mur = 1.0d0
    T_mus = 1.0d0
    T_mug = 1.0d0/3.0d0
    T_muh = 1.0d0/3.0d0

  ELSE                 ! Milbrandt and Yau scheme, LFO, etc. (based on standard gamma or exponential, no mu)

    T_mur = 1.0d0/3.0d0
    T_mus = 1.0d0/3.0d0
    T_mug = 1.0d0/3.0d0
    T_muh = 1.0d0/3.0d0

  END IF

  DO k=1,nz
    DO j=1,ny
      DO i=1,nx

        ! Assign input variables
        T_qr = qr(i,j,k)
        T_qs = qs(i,j,k) 
        T_qh = qh(i,j,k)
        T_qg = qg(i,j,k)
        T_Ntr = Ntr(i,j,k)
        T_Nts = Nts(i,j,k)
        T_Ntg = Ntg(i,j,k)
        T_Nth = Nth(i,j,k)
        T_alfr = alphar(i,j,k)
        T_alfs = alphas(i,j,k)
        T_alfg = alphag(i,j,k)
        T_alfh = alphah(i,j,k)

        IF(MFflg == 1) THEN ! Explicit melted water fraction
          T_qsw = qsw(i,j,k)
          T_qgw = qgw(i,j,k)
          T_qhw = qhw(i,j,k)
        ELSE ! No explicit melted water fraction and/or use fractionWater subroutine
          T_qsw = 0.0
          T_qgw = 0.0
          T_qhw = 0.0
        END IF

        rhog = rhogrpl(i,j,k)
        rhoh = rhohail(i,j,k)

        !print*,'rhog,rhoh',rhog,rhoh

        rhoa = rhoa_arr(i,j,k)

        ! Compute the polarimetric variables
        CALL refl_rsa(MFflg,rhoa)

        !print*,'i,j,k = ',i,j,k

        ! Store result in arrays for output
        logZ(i,j,k) = T_log_ref
        sumZh(i,j,k) = T_sum_ref_h
        sumZv(i,j,k) = T_sum_ref_v
        logZdr(i,j,k) = T_log_zdr
        sumZhv(i,j,k) = T_sum_ref_hv
        Kdp(i,j,k) = T_kdp
        Ahh(i,j,k) = T_Ahh
        Avv(i,j,k) = T_Avv
        
      END DO
    END DO
  END DO

  RETURN

END SUBROUTINE refl_rsa_array

SUBROUTINE solve_alpha_iter(nx,ny,nz,rhoa,mu,q,Ntx,Z,rhox,alpha)
!
!-----------------------------------------------------------------------
!  PURPOSE:  Calculates shape parameter alpha
!-----------------------------------------------------------------------
!
!  AUTHOR: Dan Dawson
!  (06/09/2011)
!
!  Calculates the shape parameter (alpha) for a gamma distribution
!  of the form N(D) = N0*D^alpha*exp(-lambda*D^3mu), using an iterative
!  approach.  When mu=1/3 the standard gamma distribution is recovered.
!  When mu = 1.0, you get a gamma distribution in terms of volume.  Right
!  now the subroutine only allows a fixed mu=1/3 or 1.
!
!  MODIFICATION HISTORY:
!-----------------------------------------------------------------------
!  Variable Declarations:
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: pi = 3.141592   ! pi
  INTEGER, INTENT(IN) :: nx,ny,nz
  REAL, INTENT(IN) :: mu
  REAL, INTENT(IN)    :: rhoa(nx,ny,nz),q(nx,ny,nz),Ntx(nx,ny,nz),Z(nx,ny,nz)
  REAL, INTENT(IN)    :: rhox(nx,ny,nz)
  REAL, INTENT(OUT)  :: alpha(nx,ny,nz)

  REAL :: temp
  REAL :: nu

  INTEGER :: i,j,k,iter

  REAL, PARAMETER :: rnumin = -0.8
  REAL, PARAMETER :: rnumax = 8.0
  REAL, PARAMETER :: alpmin = 0.0
  REAL, PARAMETER :: alpmax = 15.0

  alpha = 0.0

  IF(mu == 1.0) THEN
    alpha = rnumin
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          IF(q(i,j,k) > 0.0 .and. Ntx(i,j,k) > 0.0 .and. Z(i,j,k) > 0.0) THEN
            ! Compute mean volume
            temp = rhoa(i,j,k)*q(i,j,k)/(1000.*max(1.0e-9,Ntx(i,j,k)))
            ! Compute initial value of nu before iteration
            nu = 36.*(alpha(i,j,k)+2.0)*Ntx(i,j,k)*temp**2/(Z(i,j,k)*pi**2)-1.
            ! Perform the iteration on nu
            DO iter=1,10 ! Ten iterations should be enough, according to Ted
              IF(abs(nu - alpha(i,j,k)) .lt. 0.01) EXIT
              alpha(i,j,k) = max(rnumin, min(rnumax,nu))
              nu = 36.*(alpha(i,j,k)+2.0)*Ntx(i,j,k)*temp**2/(Z(i,j,k)*pi**2)-1.
              nu = max(rnumin,min(rnumax,nu))
            END DO
            ! Convert alpha to "standard gamma" alpha (i.e. alpha = 3*nu+2)
            alpha(i,j,k) = max(alpmin,min(alpmax,3.*alpha(i,j,k)+2.))
          ELSE
            alpha(i,j,k) = 0.0
          END IF   
        END DO
      END DO
    END DO
  ELSE IF (mu == 1./3.) THEN
    alpha = alpmin
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          IF(q(i,j,k) > 0.0 .and. Ntx(i,j,k) > 0.0 .and. Z(i,j,k) > 0.0) THEN
            temp = Z(i,j,k)*(pi/6.*rhox(i,j,k))**2*Ntx(i,j,k)/((rhoa(i,j,k)*q(i,j,k))**2)
            ! Compute initial value of nu before iteration (note, nu is actually alpha)
            nu = (6.+alpha(i,j,k))*(5.+alpha(i,j,k))*(4.+alpha(i,j,k))/  &
                 ((3.+alpha(i,j,k))*(2.+alpha(i,j,k))*temp) - 1.0
            DO iter=1,10
              IF(abs(nu - alpha(i,j,k)) .lt. 0.01) EXIT
              nu = (6.+alpha(i,j,k))*(5.+alpha(i,j,k))*(4.+alpha(i,j,k))/  &
                 ((3.+alpha(i,j,k))*(2.+alpha(i,j,k))*temp) - 1.0
              alpha(i,j,k) = max(alpmin,min(alpmax,nu))
            END DO
          ELSE
            alpha(i,j,k) = 0.0
          END IF
        END DO
      END DO
    END DO
  ELSE
    print*,'Sorry, only mu = 1/3 or mu = 1 accepted for now. Try again later!'
    RETURN
  END IF

  RETURN

END SUBROUTINE solve_alpha_iter

INTEGER FUNCTION get_qgh_opt (graupel_ON, hail_ON)

  INTEGER :: graupel_ON,hail_ON

  get_qgh_opt = 0

  IF(graupel_ON == 0 .and. hail_ON == 0) THEN
    get_qgh_opt = 1
  ELSE IF(graupel_ON == 0 .and. hail_ON == 1) THEN
    get_qgh_opt = 2
  ELSE IF(graupel_ON == 1 .and. hail_ON == 0) THEN
    get_qgh_opt = 3
  ELSE IF(graupel_ON == 1 .and. hail_ON == 1) THEN 
    get_qgh_opt = 4
  ENDIF 

END FUNCTION get_qgh_opt

FUNCTION gamma(xx)

!  Modified from "Numerical Recipes"

  IMPLICIT NONE

! PASSING PARAMETERS:
  DOUBLE PRECISION, INTENT(IN) :: xx

! LOCAL PARAMETERS:
  DOUBLE PRECISION  :: gamma
  INTEGER  :: j
  DOUBLE PRECISION  :: ser,stp,tmp,x,y,cof(6)


  SAVE cof,stp
  DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,               &
       24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,  &
       -.5395239384953d-5,2.5066282746310005d0/
  x=xx
  y=x
  tmp=x+5.5d0
  tmp=(x+0.5d0)*log(tmp)-tmp
  ser=1.000000000190015d0
! do j=1,6   !original
  do j=1,4
!!do j=1,3   !gives result to within ~ 3 %
     y=y+1.d0
     ser=ser+cof(j)/y
  enddo
  gamma=tmp+log(stp*ser/x)
  gamma= exp(gamma)

END FUNCTION gamma

FUNCTION solveAlpha(Q,N,Z,Cx,rho)

 IMPLICIT NONE

! PASSING PARAMETERS:
  real, INTENT(IN) :: Q, N, Z, Cx, rho

! LOCAL PARAMETERS:
  real*8 :: solveAlpha
  real   :: a,g,a1,g1,g2,tmp1
  integer :: i
  real, parameter :: alphaMax= 40.
  real, parameter :: epsQ    = 1.e-14
  real, parameter :: epsN    = 1.e-3
  real, parameter :: epsZ    = 1.e-32

!  Q         mass mixing ratio
!  N         total concentration
!  Z         reflectivity
!  Cx        (pi/6)*RHOx
!  rho       air density
!  a         alpha (returned as solveAlpha)
!  g         function g(a)= [(6+a)(5+a)(4+a)]/[(3+a)(2+a)(1+a)],
!              where g = (Cx/(rho*Q))**2.*(Z*N)


  if (Q==0. .or. N==0. .or. Z==0. .or. Cx==0. .or. rho==0.) then
  ! For testing/debugging only; this module should never be called
  ! if the above condition is true.
    print*,'*** STOPPED in MODULE ### solveAlpha *** '
    print*,'*** : ',Q,N,Z,Cx*1.9099,rho
    stop
  endif

  IF (Q>epsQ .and. N>epsN .and. Z>epsZ ) THEN

     tmp1= Cx/(rho*Q)
     g   = tmp1*Z*tmp1*N    ! g = (Z*N)*[Cx / (rho*Q)]^2

 !Note: The above order avoids OVERFLOW, since tmp1*tmp1 is very large

!----------------------------------------------------------!
! !Solve alpha numerically: (brute-force; for testing only)
!      a= 0.
!      g2= 999.
!      do i=0,4000
!         a1= i*0.01
!         g1= (6.+a1)*(5.+a1)*(4.+a1)/((3.+a1)*(2.+a1)*(1.+a1))
!         if(abs(g-g1)<abs(g-g2)) then
!            a = a1
!            g2= g1
!         endif
!      enddo
!----------------------------------------------------------!

!Piecewise-polynomial approximation of g(a) to solve for a:
     if (g>=20.) then
       a= 0.
     else
       g2= g*g
       if (g<20.  .and.g>=13.31) a= 3.3638e-3*g2 - 1.7152e-1*g + 2.0857e+0
       if (g<13.31.and.g>=7.123) a= 1.5900e-2*g2 - 4.8202e-1*g + 4.0108e+0
       if (g<7.123.and.g>=4.200) a= 1.0730e-1*g2 - 1.7481e+0*g + 8.4246e+0
       if (g<4.200.and.g>=2.946) a= 5.9070e-1*g2 - 5.7918e+0*g + 1.6919e+1
       if (g<2.946.and.g>=1.793) a= 4.3966e+0*g2 - 2.6659e+1*g + 4.5477e+1
       if (g<1.793.and.g>=1.405) a= 4.7552e+1*g2 - 1.7958e+2*g + 1.8126e+2
       if (g<1.405.and.g>=1.230) a= 3.0889e+2*g2 - 9.0854e+2*g + 6.8995e+2
       if (g<1.230) a= alphaMax
     endif

     solveAlpha= max(0.,min(a,alphaMax))

  ELSE

     solveAlpha= 0.

  ENDIF

END FUNCTION solveAlpha
