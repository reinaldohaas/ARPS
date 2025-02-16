!########################################################################
!########################################################################
!#########                                                      #########
!#########              FUNCTIONS for raindrops                 #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

REAL FUNCTION rhfactorCG(qr,rho,Zerhf)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Same as SUBROUTINE rhfactor except for no truncation for raindrop
! size at 8 mm.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/29/2005
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL,INTENT(IN) :: qr, rho, Zerhf

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   rhfactorCG = Zerhf*((qr*rho) ** 1.77)

END FUNCTION rhfactorCG


REAL FUNCTION rvfactorCG(qr,rho,Zervf)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Same as SUBROUTINE rhfactor except for no truncation for raindrop
! size at 8 mm.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/29/2005
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL,INTENT(IN) :: qr, rho, Zervf

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

   rvfactorCG = Zervf*((qr*rho) ** 1.635)

END FUNCTION rvfactorCG

!########################################################################
!########################################################################
!#########                                                      #########
!#########                  FUNCTIONS for snow                  #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

REAL FUNCTION sfactor(Mfactor,qs,rho,sqpow,Zesf)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This user defined function calculate reflectivity factor
! from mixing ratio of snow. 
! For equation sets, see subroutine ZhhFromDualPol.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 12/3/2004
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL,INTENT(IN) :: Mfactor, qs, rho, sqpow, Zesf

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  sfactor = Zesf*Mfactor*(qs*rho)**sqpow

END FUNCTION sfactor


REAL FUNCTION sfactorPos(qs,rho,Zesposf)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Preexisting ARPS formular for wet snow. Formular is given in
! suboutine "reflec_ferrier".
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 2/3/2004
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL,INTENT(IN) :: qs, rho, Zesposf

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  sfactorPos = Zesposf*((qs*rho) ** 1.75)

END FUNCTION sfactorPos


REAL FUNCTION sfactorNeg(qs,rho,Zesnegf)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Preexisting ARPS formular for dry snow. Formular is given in
! suboutine "reflec_ferrier".
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 2/3/2004
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL,INTENT(IN) :: qs, rho, Zesnegf

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  sfactorNeg = Zesnegf*((qs*rho) ** 1.75)

END FUNCTION sfactorNeg

REAL FUNCTION shfactor(qs,rho,Zeshf)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This user defined function calculate partial reflectivity factor
! of snow for horizontal polarization. 
!
!                    4        2     
!           4 * lambda  * lapha  * No 
!   Zes = ------------------------------- * gamma(2*beta + 1)
!             4       2      (2*beta + 1)
!           pi  * |Kw|  * gam
!
!           where, alpha = 3.01e-5,   beta = 2.246
!
! Formula and constants are given by Dr. Guifu Zhang
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 10/29/2005
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL,INTENT(IN) :: qs, rho, Zeshf

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  shfactor = Zeshf*(qs*1.e6*rho)**1.669

END FUNCTION shfactor

REAL FUNCTION svfactor(qs,rho,Zesvf)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This user defined function calculate partial reflectivity factor
! of snow for vertical polarization. 
!
!                    4        2     
!           4 * lambda  * lapha  * No 
!   Zes = ------------------------------- * gamma(2*beta + 1)
!             4       2      (2*beta + 1)
!           pi  * |Kw|  * gam
!
!           where, alpha = 2.89e-5,   beta = 2.26
!
! Formula and constants are given by Dr. Guifu Zhang
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 10/29/2005
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL,INTENT(IN) :: qs, rho, Zesvf

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  svfactor = Zesvf*(qs*1.e6*rho)**1.678

END FUNCTION svfactor

REAL FUNCTION snow_alpha_a(fw)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This user defined function calculates alpha for Zhh
! for dry/melting snow.
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
! Declare local variables.
!-----------------------------------------------------------------------

  REAL,INTENT(IN) :: fw

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  snow_alpha_a = (0.194 + 7.094*fw + 2.135*fw**2. - 5.225*fw**3.)*10**-4.

END FUNCTION snow_alpha_a

REAL FUNCTION snow_alpha_b(fw)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This user defined function calculates alpha for Zvv
! for dry/melting snow.
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
! Declare local variables.
!-----------------------------------------------------------------------

  REAL,INTENT(IN) :: fw

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  snow_alpha_b = (0.191 + 6.916*fw - 2.841*fw**2. - 1.160*fw**3.)*10**-4.

END FUNCTION snow_alpha_b


!########################################################################
!########################################################################
!#########                                                      #########
!#########                  FUNCTIONS for hail                  #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

REAL FUNCTION hfactorPos(qh,rho,rqhpowf,Zehposf)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This user defined function calculate reflectivity factor
! from mixing ratio of melting hail.
! For equation sets, see subroutine ZhhFromDualPol.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 12/3/2004
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL,INTENT(IN) :: qh, rho, rqhpowf, Zehposf

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  hfactorPos = Zehposf*((qh*rho) ** rqhpowf)

END FUNCTION hfactorPos


REAL FUNCTION hfactorNeg(qh,rho,Zehnegf)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This user defined function calculate reflectivity factor
! from mixing ratio of dry hail.
! For equation sets, see subroutine ZhhFromDualPol.
!
! 
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 12/3/2004
! 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE
  
!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL,INTENT(IN) :: qh, rho, Zehnegf
  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  hfactorNeg = Zehnegf*((qh*rho) ** 1.75)

END FUNCTION hfactorNeg

REAL FUNCTION hail_alpha_a(fw)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This user defined function calculates alpha for Zhh
! for dry/melting hail.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/22/2007
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL,INTENT(IN) :: fw

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  hail_alpha_a = (0.191 + 2.39*fw - 12.57*fw**2. + 38.71*fw**3.   &
                  - 65.53*fw**4. + 56.16*fw**5. - 18.98*fw**6.)*10**-3.

END FUNCTION hail_alpha_a

REAL FUNCTION hail_alpha_b(fw)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This user defined function calculates alpha for Zvv
! for dry/melting hail.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/22/2007
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL,INTENT(IN) :: fw

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  hail_alpha_b = (0.165 + 1.72*fw - 9.92*fw**2. + 32.15*fw**3.   &
                  - 56.0*fw**4. + 48.83*fw**5. - 16.69*fw**6.)*10**-3.

END FUNCTION hail_alpha_b

SUBROUTINE partialRefIce(Cof,N0,Ai,Bi,Ci,alp_a,alp_b,rho_air, &
                         mxratio,rho_hydro,refIceHH,refIceVV)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This function calculates the partial reflectivity for each species
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/22/2007
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare variables.
!-----------------------------------------------------------------------
  REAL,INTENT(IN) :: Cof,N0,Ai,Bi,Ci,alp_a,alp_b,rho_air
  REAL,INTENT(IN) :: mxratio,rho_hydro

  REAL,INTENT(OUT) :: refIceHH,refIceVV

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  refIceHH = Cof*(Ai*alp_a**2+Bi*alp_b**2+2*Ci*alp_a*alp_b) &
             *(rho_air*mxratio/rho_hydro)**1.75/N0**0.75
  refIceVV = Cof*(Bi*alp_a**2+Ai*alp_b**2+2*Ci*alp_a*alp_b) &
             *(rho_air*mxratio/rho_hydro)**1.75/N0**0.75

END SUBROUTINE partialRefIce

SUBROUTINE partialRhoIce(Cof,N0,Ci,Di,alp_a,alp_b,rho_air,rho_0, &
                         mxratio,rho_hydro,refIceHV)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This function calculates the cross components, Z_hv, for each species
! for rho_hv calculation.
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 12/16/2007
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare variables.
!-----------------------------------------------------------------------
  REAL,INTENT(IN) :: Cof,N0,Ci,Di,alp_a,alp_b,rho_air
  REAL,INTENT(IN) :: mxratio,rho_hydro,rho_0

  REAL,INTENT(OUT) :: refIceHV

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  refIceHV = Cof*(Ci*alp_a**2+Ci*alp_b**2+2*Di*alp_a*alp_b*rho_0) &
             *(rho_air*mxratio/rho_hydro)**1.75/N0**0.75

END SUBROUTINE partialRhoIce

SUBROUTINE partialKdpIce(iceCof,Ck,alp_k,rho_air,mxratio,rho_hydro,  &
                         KdpIce)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This function calculates the partial Kdp for each species
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 1/22/2007
!
!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare variables.
!-----------------------------------------------------------------------
  REAL,INTENT(IN) :: iceCof,Ck,alp_k,rho_air
  REAL,INTENT(IN) :: mxratio,rho_hydro

  REAL,INTENT(OUT) :: KdpIce

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  KdpIce = iceCof*Ck*alp_k*(rho_air*mxratio/rho_hydro)

END SUBROUTINE partialKdpIce

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
!#########              FUNCTION  combine                       #########
!#########                                                      #########
!########################################################################
!########################################################################

REAL FUNCTION combine(rmax, rmin, temp, val1, val2)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Returns the reflectivity factor as a linear combination of wet and dry 
! formulars from 0 to 1 depending on the temperature.
! Variables:
!        rmax       higher temperature
!        rmin       lower temperature
!        val1       value from dry formular
!        val2       value from wet formular
!
!-----------------------------------------------------------------------
!
! AUTHOR:  Youngsun Jung, 12/6/2004
!
!-----------------------------------------------------------------------
! Include global variables only for dual-pol calculations
!-----------------------------------------------------------------------

  USE DUALPARA

!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL, INTENT(IN) :: rmax, rmin, val1, val2, temp
  REAL :: const

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  const = (rmax-temp)/(rmax-rmin)
  combine = const * val1 + (1-const)*val2

END FUNCTION combine

!########################################################################
!########################################################################
!#########                                                      #########
!#########              FUNCTION calculate_kdp                  #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

REAL FUNCTION calculate_kdp(rho,qr,qs,qh,drysnowopt)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! This subroutine estimates specific differential phase.
! For details, see "Kdp" in convertZ.f90.
!
!
! AUTHOR:  Youngsun Jung, 2/28/2007
!
!-----------------------------------------------------------------------
! Include global variables only for dual-pol calculations
!-----------------------------------------------------------------------

  USE DUALPARA

!-----------------------------------------------------------------------
! Force explicit declarations.
!-----------------------------------------------------------------------

  IMPLICIT NONE

!-----------------------------------------------------------------------
! Declare arguments.
!-----------------------------------------------------------------------

  REAL, INTENT(IN) :: rho ! Air density (kg m**-3)
  REAL, INTENT(IN) :: qr  ! Rain mixing ratio (kg kg**-1)
  REAL, INTENT(IN) :: qs  ! Snow mixing ratio (kg kg**-1)
  REAL, INTENT(IN) :: qh  ! Hail mixing ratio (kg kg**-1)

!-----------------------------------------------------------------------
! Declare local variables.
!-----------------------------------------------------------------------

  REAL :: rcomp,scomp,hcomp,sumcomp
  INTEGER :: drysnowopt

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  calculate_kdp = 0.0
  sumcomp = 0.0
  rcomp = 0.0
  scomp = 0.0

!-----------------------------------------------------------------------
! Check for bad air density value.
!-----------------------------------------------------------------------

        IF (rho > 0.0) THEN        ! outer if starts

!-----------------------------------------------------------------------
! Calculate kdp contribution from rain.
!-----------------------------------------------------------------------
          IF (drysnowopt < 4) THEN

            IF (qr > 0.0) THEN
              rcomp = constKdpr * ((qr*rho) ** 1.41)
            END IF

!-----------------------------------------------------------------------
! Calculate kdp contribution from snow.
!-----------------------------------------------------------------------

            IF (qs > 0.0) THEN
              scomp = constKdpsZ * ((qs*1.e6*rho) ** 1.27)
            ENDIF

!-----------------------------------------------------------------------
! Now add the contributions to kdp
!-----------------------------------------------------------------------

            sumcomp = rcomp + scomp
            calculate_kdp = sumcomp

          ELSE IF ((qs > 0.0 .or. qr > 0.0 .or. qh > 0.0) &
                    .and. drysnowopt == 4) THEN 
            calculate_kdp = rainIceKdp(qr,qs,qh,rho)
          END IF

        END IF                             ! outer if ends

  RETURN
END FUNCTION calculate_kdp
