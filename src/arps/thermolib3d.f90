!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   FUNCTION F_ES                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

FUNCTION f_es( p, t )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the saturation specific humidity using enhanced Teten's
!  formula.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  01/08/1998
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    p        Pressure (Pascal)
!    t        Temperature (K)
!
!  OUTPUT:
!
!    f_es     Saturation water vapor pressure (Pa)
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  REAL :: p         ! Pressure (Pascal)
  REAL :: t         ! Temperature (K)
  REAL :: f_es      ! Saturation water vapor pressure (Pa)
!
!-----------------------------------------------------------------------
!
!  Function f_es and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_esl, f_esi

!fpp$ expand (f_esl)
!fpp$ expand (f_esi)
!!dir$ inline always f_esl, f_esi
!*$*  inline routine (f_esl, f_esi)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF ( t >= 273.15 ) THEN      ! for water
    f_es = f_esl( p,t )
  ELSE                            ! for ice
    f_es = f_esi( p,t )
  END IF

  RETURN
END FUNCTION f_es

!
!-----------------------------------------------------------------------
!
!  Calculate the saturation water vapor over liquid water using
!  enhanced Teten's formula.
!
!-----------------------------------------------------------------------
!

FUNCTION f_esl( p, t )

  IMPLICIT NONE

  REAL :: p         ! Pressure (Pascal)
  REAL :: t         ! Temperature (K)
  REAL :: f_esl     ! Saturation water vapor pressure over liquid water

  REAL :: f

  INCLUDE 'globcst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  f = satfwa + satfwb * p
  f_esl = f * satewa * EXP( satewb*(t-273.15)/(t-satewc) )

  RETURN
END FUNCTION f_esl
!
!-----------------------------------------------------------------------
!
!  Calculate the saturation water vapor over ice using enhanced
!  Teten's formula.
!
!-----------------------------------------------------------------------
!

FUNCTION f_esi( p, t )

  IMPLICIT NONE

  REAL :: p         ! Pressure (Pascal)
  REAL :: t         ! Temperature (K)
  REAL :: f_esi     ! Saturation water vapor pressure over ice (Pa)

  REAL :: f

  INCLUDE 'globcst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  f = satfia + satfib * p
  f_esi = f * sateia * EXP( sateib*(t-273.15)/(t-sateic) )

  RETURN
END FUNCTION f_esi
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   FUNCTION F_DESDT                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION f_desdt( t )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate d(es)/dt/es using enhanced Teten's formula. See function
!  f_es.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  01/12/1998
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    t        Temperature (K)
!
!  OUTPUT:
!
!    f_desdt  d(es)/dt/es (1/K)
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  REAL :: t         ! Temperature (K)
  REAL :: f_desdt   ! d(es)/dt/es (1/K)
!
!-----------------------------------------------------------------------
!
!  Function f_desdtl and f_desdti and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_desdtl, f_desdti

!fpp$ expand (f_desdtl)
!fpp$ expand (f_desdti)
!!dir$ inline always f_desdtl, f_desdti
!*$*  inline routine (f_desdtl, f_desdti)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF ( t >= 273.15 ) THEN      ! for water
    f_desdt = f_desdtl(t)
  ELSE                         ! for ice
    f_desdt = f_desdti(t)
  END IF

  RETURN
  END FUNCTION f_desdt
!
!-----------------------------------------------------------------------
!
!  Calculate d(esl)/dt/esl using enhanced Teten's formula.
!  See function f_esl.
!
!-----------------------------------------------------------------------
!

  FUNCTION f_desdtl( t )

  IMPLICIT NONE

  REAL :: t         ! Temperature (K)
  REAL :: f_desdtl  ! d(esl)/dt/esl (1/K)

  INCLUDE 'globcst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  f_desdtl = satewb*(273.15-satewc)/(t-satewc)**2

  RETURN
  END FUNCTION f_desdtl
!
!-----------------------------------------------------------------------
!
!  Calculate d(esi)/dt/esi using enhanced Teten's formula.
!  See function f_esi.
!
!-----------------------------------------------------------------------
!

  FUNCTION f_desdti( t )

  IMPLICIT NONE

  REAL :: t         ! Temperature (K)
  REAL :: f_desdti  ! d(esi)/dt/esi (1/K)

  INCLUDE 'globcst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  f_desdti = sateib*(273.15-sateic)/(t-sateic)**2

  RETURN
  END FUNCTION f_desdti
!
!-----------------------------------------------------------------------
!
!  Calculate d(qvs)/dt/qvs using enhanced Teten's formula.
!  See function f_es.
!
!           Rd
!  rddrv = ----
!           Rv
!
!  qvs = rddrv * es / (p - ( 1.0 - rddrv ) * es )
!
!    1    d(qvs)           1 - rddrv             1    d(es)
!  ----- -------- = ( 1 + ----------- * qvs ) * ---- -------
!   qvs     dT               rddrv               es    dT
!
!-----------------------------------------------------------------------
!

  FUNCTION f_dqvsdt( p, t )

  IMPLICIT NONE

  REAL :: p         ! Pressure (Pa)
  REAL :: t         ! Temperature (K)
  REAL :: f_dqvsdt  ! d(esi)/dt/esi (1/K)

  REAL :: qvsat, desdt
!
!-----------------------------------------------------------------------
!
!  Function f_desdtl and f_desdti and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_desdt, f_qvsat

!fpp$ expand (f_desdt)
!fpp$ expand (f_qvsat)
!!dir$ inline always f_desdt, f_qvsat
!*$*  inline routine (f_desdt, f_qvsat)

  INCLUDE 'phycst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  qvsat = f_qvsat( p,t )
  desdt = f_desdt( t )
  f_dqvsdt = ( 1. + qvsat*(1.-rddrv)/rddrv ) * desdt

  RETURN
  END FUNCTION f_dqvsdt
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   FUNCTION F_QVSAT                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

FUNCTION f_qvsat( p, t )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the saturation specific humidity using enhanced Teten's
!  formula.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  01/08/1998
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    p        Pressure (Pascal)
!    t        Temperature (K)
!
!  OUTPUT:
!
!    f_qvsat  Saturation water vapor specific humidity (kg/kg).
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  REAL :: p         ! Pressure (Pascal)
  REAL :: t         ! Temperature (K)
  REAL :: f_qvsat   ! Saturation water vapor specific humidity (kg/kg)
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  Function f_es and inline directive for Cray PVP
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
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  f_qvsat = rddrv * f_es(p,t) / (p-(1.0-rddrv)*f_es(p,t))

  RETURN
END FUNCTION f_qvsat
!
!-----------------------------------------------------------------------
!
!  Calculate the saturation specific humidity over liquid water using
!  enhanced Teten's formula.
!
!-----------------------------------------------------------------------
!

  FUNCTION f_qvsatl( p, t )

  IMPLICIT NONE

  REAL :: p         ! Pressure (Pascal)
  REAL :: t         ! Temperature (K)
  REAL :: f_qvsatl  ! Saturation specific humidity over liquid
                    ! water (kg/kg)

  REAL :: fesl

  INCLUDE 'phycst.inc'

  REAL :: f_esl

!fpp$ expand (f_esl)
!!dir$ inline always f_esl
!*$*  inline routine (f_esl)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  fesl = f_esl( p, t )
  f_qvsatl = rddrv * fesl / (p-(1.0-rddrv)*fesl)

  RETURN
  END FUNCTION f_qvsatl
!
!-----------------------------------------------------------------------
!
!  Calculate the saturation specific humidity over ice using
!  enhanced Teten's formula.
!
!-----------------------------------------------------------------------
!

  FUNCTION f_qvsati( p, t )

  IMPLICIT NONE

  REAL :: p         ! Pressure (Pascal)
  REAL :: t         ! Temperature (K)
  REAL :: f_qvsati  ! Saturation specific humidity over ice (kg/kg)

  REAL :: fesi

  INCLUDE 'phycst.inc'

  REAL :: f_esi

!fpp$ expand (f_esi)
!!dir$ inline always f_esi
!*$*  inline routine (f_esi)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  fesi = f_esi( p, t )
  f_qvsati = rddrv * fesi / (p-(1.0-rddrv)*fesi)

  RETURN
  END FUNCTION f_qvsati
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   FUNCTION F_MRSAT                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION f_mrsat( p, t )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the saturation water vapor mixing ratio using enhanced
!  Teten's formula.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  01/08/1998
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    p        Pressure (Pascal)
!    t        Temperature (K)
!
!  OUTPUT:
!
!    f_mrsat  Saturation water vapor mixing ratio (kg/kg).
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  REAL :: p         ! Pressure (Pascal)
  REAL :: t         ! Temperature (K)
  REAL :: f_mrsat   ! Saturation water vapor mixing ratio (kg/kg)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  REAL :: fes
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  Function f_es and inline directive for Cray PVP
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
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  fes = f_es( p,t )
  f_mrsat = rddrv * fes / (p-fes)

  RETURN
  END FUNCTION f_mrsat
!
!-----------------------------------------------------------------------
!
!  Calculate the saturated water vapor mixing ratio over liquid water
!  using enhanced Teten's formula.
!
!-----------------------------------------------------------------------
!

  FUNCTION f_mrsatl( p, t )

  IMPLICIT NONE

  REAL :: p         ! Pressure (Pascal)
  REAL :: t         ! Temperature (K)
  REAL :: f_mrsatl  ! Saturation water vapor mixing ratio over liquid
                    ! water (kg/kg)

  REAL :: fesl

  INCLUDE 'phycst.inc'

  REAL :: f_esl

!fpp$ expand (f_esl)
!!dir$ inline always f_esl
!*$*  inline routine (f_esl)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  fesl = f_esl( p,t )
  f_mrsatl = rddrv * fesl / (p-fesl)

  RETURN
  END FUNCTION f_mrsatl
!
!-----------------------------------------------------------------------
!
!  Calculate the saturated water vapor mixing ratio over ice
!  using enhanced Teten's formula.
!
!-----------------------------------------------------------------------
!

  FUNCTION f_mrsati( p, t )

  IMPLICIT NONE

  REAL :: p         ! Pressure (Pascal)
  REAL :: t         ! Temperature (K)
  REAL :: f_mrsati  ! Saturation water vapor mixing ratio over ice
                    ! (kg/kg)

  REAL :: fesi

  INCLUDE 'phycst.inc'

  REAL :: f_esi

!fpp$ expand (f_esi)
!!dir$ inline always f_esi
!*$*  inline routine (f_esi)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  fesi = f_esi( p,t )
  f_mrsati = rddrv * fesi / (p-fesi)

  RETURN
  END FUNCTION f_mrsati
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   FUNCTION F_TDEW                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION f_tdew( p, t, qv )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the dew point temperature using reversed enhanced
!  Tetan's formula.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  01/09/1998
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    p        Pressure (Pascal)
!    t        Temperature (K)
!    qv       Specific humidity (kg/kg)
!
!  OUTPUT:
!
!    f_tdew   Dew point temperature (K)
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  REAL :: p         ! Pressure (Pascal)
  REAL :: t         ! Temperature (K)
  REAL :: qv        ! Specific humidity (kg/kg)
  REAL :: f_tdew    ! Dew point temperature (K)
!
!-----------------------------------------------------------------------
!
!  Include file
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
!
!-----------------------------------------------------------------------
!
!  Function f_tdewl and f_tdewi and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_tdewl, f_tdewi

!fpp$ expand (f_tdewl)
!fpp$ expand (f_tdewi)
!!dir$ inline always f_tdewl, f_tdewi
!*$*  inline routine (f_tdewl, f_tdewi)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF ( t >= 273.15 ) THEN         ! for water
    f_tdew = f_tdewl( p,qv )
    IF ( f_tdew < 273.15 ) THEN  ! using ice formula
      f_tdew = f_tdewi( p,qv )
    END IF
  ELSE                            ! for ice
    f_tdew = f_tdewi( p,qv )
  END IF

  RETURN
  END FUNCTION f_tdew
!
!-----------------------------------------------------------------------
!
!  Calculate the dew point temperature over liquid water using the
!  reversed enhanced Teten's formula.
!
!-----------------------------------------------------------------------
!

  FUNCTION f_tdewl( p, qv )

  IMPLICIT NONE

  REAL :: p         ! Pressure (Pascal)
  REAL :: qv        ! Specific humidity (kg/kg)
  REAL :: f_tdewl   ! Dew point temperature over liquid water (K)

  REAL :: qvs, esl, f, ln

  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  qvs = MAX( qv, 1.e-8 )

  f = satfwa + satfwb * p
  esl = p*qvs/(rddrv+(1-rddrv)*qvs)
  ln = ALOG( esl/(f*satewa) )

  f_tdewl = ( ln*satewc - 273.15*satewb ) / ( ln - satewb )

  RETURN
  END FUNCTION f_tdewl
!
!-----------------------------------------------------------------------
!
!  Calculate the dew point temperature over ice using the
!  reversed enhanced Teten's formula.
!
!-----------------------------------------------------------------------
!

  FUNCTION f_tdewi( p, qv )

  IMPLICIT NONE

  REAL :: p         ! Pressure (Pascal)
  REAL :: qv        ! Specific humidity (kg/kg)
  REAL :: f_tdewi   ! Dew point temperature over ice (K)

  REAL :: qvs, esi, f, ln

  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  qvs = MAX( qv, 1.e-8 )

  f = satfia + satfib * p
  esi = p*qvs/(rddrv+(1-rddrv)*qvs)
  ln = ALOG( esi/(f*sateia) )

  f_tdewi = ( ln*sateic - 273.15*sateib ) / ( ln - sateib )

  RETURN
  END FUNCTION f_tdewi
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE F_PT2PTE                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION f_pt2pte( p, pt, qv )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the equivalent potential temperature of an air
!  parcel that is either saturated or unsaturated.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  01/14/1998  Re-created from subroutine EQUIPT
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    p        Base state pressure (Pascals)
!    pt       Potential temperature (degrees Kelvin)
!    qv       Water vapor specific humidity (kg/kg)
!
!  OUTPUT:
!
!    f_pt2pte Equivalent potential temperature (K)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  REAL :: p         ! Pressure (Pascals)
  REAL :: pt        ! Potential temperature (degrees Kelvin)
  REAL :: qv        ! Water vapor specific humidity (kg/kg)

  REAL :: f_pt2pte  ! Equivalent potential temperature (K)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  REAL :: t       ! Temperature
  REAL :: tdew    ! Dew point temperature
  REAL :: tl      ! Temperature at the (lifting) condensation level
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  Function f_tdew and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_tdew

!fpp$ expand (f_tdew)
!!dir$ inline always f_tdew
!*$*  inline routine (f_tdew)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  t = pt * (p/p0)**rddcp

  ! DTD: slightly more accurate computation of T that takes into account
  ! water vapor

!  t = pt * (p/p0)**(0.2854*(1.-(0.28e-3)*1000.*(qv/(1.-qv))))

  tdew = f_tdew( p,t,qv )    ! Dew point temperture.

!  f_pt2pte = pt * EXP( lathv*qv/(cp*tdew) )

! DTD: changed to more accurate formula from Bolton (1980)

  ! First calculate the temperature at the condensation level (TL, eqn. 15 in Bolton)

  tl = 1./((1./(tdew - 56.))+LOG(t/tdew)/800.)+56.

  ! Then calculate equivalent potential temperature using Bolton's approximation
  ! (eqn. 43, Bolton 1980)

  ! Convert qv to mixing ratio

  qv = qv/(1.-qv)

  ! Now compute pte

  f_pt2pte = (t*(p0/p)**(0.2854*(1.-(0.28e-3)*1000.*qv)))*EXP((3.376/tl-0.00254)*1000.*qv*(1+(0.81e-3)*1000.*qv))



  RETURN
  END FUNCTION f_pt2pte
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETQVS                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getqvs( nx,ny,nz, ibgn,iend,jbgn,jend,kbgn,kend,             &
           p, t, qvsat )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the saturation specific humidity using Teten's formula.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  3/17/1991.
!
!  MODIFICATION HISTORY:
!
!  5/02/92 (M. Xue)
!  Added full documentation.
!
!  5/03/92 (M. Xue)
!  Further documentation.
!
!  2/10/93 (K. Droegemeier)
!  Cleaned up documentation.
!
!  9/10/94 (Weygandt & Y. Lu)
!  Cleaned up documentation.
!
!  12/11/1997 (Yuhe Liu)
!  Rewrote the subroutine calling function F_QVSAT to calculate the
!  saturated specific humidity for array.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Array dimension in x-direction
!    ny       Array dimension in y-direction
!    nz       Array dimension in z-direction
!    ibgn     Starting index in x-direction
!    jbgn     Starting index in y-direction
!    kbgn     Starting index in z-direction
!    iend     Last index in x-direction
!    jend     Last index in y-direction
!    kend     Last index in z-direction
!
!    p        Pressure (Pascal)
!    t        Temperature (K)
!
!  OUTPUT:
!
!    qvsat    Saturation water vapor specific humidity (kg/kg).
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz         ! Input array dimension

  INTEGER :: ibgn             ! Starting index in x-direction
  INTEGER :: jbgn             ! Starting index in y-direction
  INTEGER :: kbgn             ! Starting index in z-direction
  INTEGER :: iend             ! Last index in x-direction
  INTEGER :: jend             ! Last index in y-direction
  INTEGER :: kend             ! Last index in z-direction

  REAL :: p  (nx,ny,nz)       ! Pressure (Pascal)
  REAL :: t  (nx,ny,nz)       ! Temperature (K)
  REAL :: qvsat(nx,ny,nz)     ! Saturation water vapor specific
                              ! humidity (kg/kg)
!
!-----------------------------------------------------------------------
!
!  Include file
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
!
!-----------------------------------------------------------------------
!
!  Function f_qvsat and inline directives for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_qvsat

!fpp$ expand (f_qvsat)
!!dir$ inline always f_qvsat
!*$*  inline routine (f_qvsat)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=kbgn,kend
    DO j=jbgn,jend
      DO i=ibgn,iend
        qvsat(i,j,k) = f_qvsat( p(i,j,k), t(i,j,k) )
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE getqvs
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETMRS                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getmrs( nx,ny,nz, ibgn,iend,jbgn,jend,kbgn,kend,             &
           p, t, mrsat )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the saturation mixing ratio using Teten's formula.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  01/09/1998
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Array dimension in x-direction
!    ny       Array dimension in y-direction
!    nz       Array dimension in z-direction
!    ibgn     Starting index in x-direction
!    jbgn     Starting index in y-direction
!    kbgn     Starting index in z-direction
!    iend     Last index in x-direction
!    jend     Last index in y-direction
!    kend     Last index in z-direction
!
!    p        Pressure (Pascal)
!    t        Temperature (K)
!
!  OUTPUT:
!
!    mrsat    Saturation water vapor Mixing ratio (kg/kg).
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz         ! Input array dimension

  INTEGER :: ibgn             ! Starting index in x-direction
  INTEGER :: jbgn             ! Starting index in y-direction
  INTEGER :: kbgn             ! Starting index in z-direction
  INTEGER :: iend             ! Last index in x-direction
  INTEGER :: jend             ! Last index in y-direction
  INTEGER :: kend             ! Last index in z-direction

  REAL :: p  (nx,ny,nz)       ! Pressure (Pascal)
  REAL :: t  (nx,ny,nz)       ! Temperature (K)
  REAL :: mrsat(nx,ny,nz)     ! Saturated water vapor mixing ratio (kg/kg)
!
!-----------------------------------------------------------------------
!
!  Include file
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
!
!-----------------------------------------------------------------------
!
!  Function f_mrsat and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_mrsat

!fpp$ expand (f_mrsat)
!!dir$ inline always f_mrsat
!*$*  inline routine (f_mrsat)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=kbgn,kend
    DO j=jbgn,jend
      DO i=ibgn,iend
        mrsat(i,j,k) = f_mrsat( p(i,j,k), t(i,j,k) )
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE getmrs
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GETDEW                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE getdew( nx,ny,nz, ibgn,iend,jbgn,jend,kbgn,kend,             &
           p, t, qv, tdew )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the dew point temperature (K) using enhanced Teten's
!  formula.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  01/09/1998
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Array dimension in x-direction
!    ny       Array dimension in y-direction
!    nz       Array dimension in z-direction
!    ibgn     Starting index in x-direction
!    jbgn     Starting index in y-direction
!    kbgn     Starting index in z-direction
!    iend     Last index in x-direction
!    jend     Last index in y-direction
!    kend     Last index in z-direction
!
!    p        Pressure (Pascal)
!    t        Temperature (K)
!    qv       Specific humidity (kg/kg)
!
!  OUTPUT:
!
!    tdew     Dew point temperature (K)
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz         ! Input array dimension

  INTEGER :: ibgn             ! Starting index in x-direction
  INTEGER :: jbgn             ! Starting index in y-direction
  INTEGER :: kbgn             ! Starting index in z-direction
  INTEGER :: iend             ! Last index in x-direction
  INTEGER :: jend             ! Last index in y-direction
  INTEGER :: kend             ! Last index in z-direction

  REAL :: p  (nx,ny,nz)       ! Pressure (Pascal)
  REAL :: t  (nx,ny,nz)       ! Temperature (K)
  REAL :: qv (nx,ny,nz)       ! Specific humidity (kg/kg)
  REAL :: tdew(nx,ny,nz)      ! Dew point temperature (K)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
!
!-----------------------------------------------------------------------
!
!  Function f_mrsat and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_tdew

!fpp$ expand (f_tdew)
!!dir$ inline always f_tdew
!*$*  inline routine (f_tdew)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=kbgn,kend
    DO j=jbgn,jend
      DO i=ibgn,iend
        tdew(i,j,k) = f_tdew( p(i,j,k), t(i,j,k), qv(i,j,k) )
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE getdew
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE PT2PTE                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE pt2pte( nx,ny,nz, ibgn,iend,jbgn,jend,kbgn,kend,             &
           p, pt, qv, pte )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Calculate the dew point temperature (K) using enhanced Teten's
!  formula.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  01/09/1998
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    nx       Array dimension in x-direction
!    ny       Array dimension in y-direction
!    nz       Array dimension in z-direction
!    ibgn     Starting index in x-direction
!    jbgn     Starting index in y-direction
!    kbgn     Starting index in z-direction
!    iend     Last index in x-direction
!    jend     Last index in y-direction
!    kend     Last index in z-direction
!
!    p        Pressure (Pascal)
!    pt       Temperature (K)
!    qv       Specific humidity (kg/kg)
!
!  OUTPUT:
!
!    pte      Equivalent potential temperature (K)
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz         ! Input array dimension

  INTEGER :: ibgn             ! Starting index in x-direction
  INTEGER :: jbgn             ! Starting index in y-direction
  INTEGER :: kbgn             ! Starting index in z-direction
  INTEGER :: iend             ! Last index in x-direction
  INTEGER :: jend             ! Last index in y-direction
  INTEGER :: kend             ! Last index in z-direction

  REAL :: p  (nx,ny,nz)       ! Pressure (Pascal)
  REAL :: pt (nx,ny,nz)       ! Temperature (K)
  REAL :: qv (nx,ny,nz)       ! Specific humidity (kg/kg)
  REAL :: pte(nx,ny,nz)       ! Equivalent potential temperature (K)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
!
!-----------------------------------------------------------------------
!
!  Function f_pt2pte and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_pt2pte

!fpp$ expand (f_pt2pte)
!!dir$ inline always f_pt2pte
!*$*  inline routine (f_pt2pte)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  DO k=kbgn,kend
    DO j=jbgn,jend
      DO i=ibgn,iend
        pte(i,j,k) = f_pt2pte( p(i,j,k), pt(i,j,k), qv(i,j,k) )
      END DO
    END DO
  END DO

  RETURN
END SUBROUTINE pt2pte
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   FUNCTION F_PCCL                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION f_pccl(pm,p,t,td,mrbar,n)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  this function returns the pressure at the convective condensation
!  level given the appropriate sounding data.
!
!  the algorithm is decribed on p.17 of stipanuk, g.s.,1973:
!  "algorithms for generating a skew-t log p diagram and computing
!  selected meteorological quantities," atmospheric sciences labora-
!  tory, u.s. army electronics command, white sands missile range, new
!  mexico 88002.
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!    g.s. stipanuk     1973            original version.
!    reference stipanuk paper entitled:
!         "algorithms for generating a skew-t, log p
!         diagram and computing selected meteorological
!         quantities."
!         atmospheric sciences laboratory
!         u.s. army electronics command
!         white sands missile range, new mexico 88002
!         33 pages
!    baker, schlatter  17-may-1982
!
!
!  MODIFICATIONS:
!
!  01/30/1998 (Yuhe Liu)
!  Modified to use Kelvin and Pascal as units of temperature and
!  pressure
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    p       pressure (Pa). note that p(i).gt.p(i+1).
!    t       temperature (K)
!    td      dew point (K)
!    n       number of levels in the sounding and the dimension of
!            p, t and td
!    pm      pressure (Pa) at upper boundary of the layer for
!            computing the mean mixing ratio. p(1) is the lower
!            boundary.
!
!  OUTPUT:
!
!    mrbar   mean mixing ratio (kg/kg) in the layer bounded by
!            pressures p(1) at the bottom and pm at the top
!    f_pccl  pressure (Pa) at the convective condensation level
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: n

  REAL :: t(n),p(n),td(n)
  REAL :: pm

  REAL :: mrbar
  REAL :: f_pccl
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k,l

  REAL :: pc,tq,del,x,a
  REAL :: mrsat1,mrsat2
!
!-----------------------------------------------------------------------
!
!  Function f_mrsat and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_mrsat
  REAL :: f_tmr

!fpp$ expand (f_mrsat)
!fpp$ expand (f_tmr)
!!dir$ inline always f_mrsat, f_tmr
!*$*  inline routine (f_mrsat, f_tmr)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (pm /= p(1)) GO TO 5

  mrbar = f_mrsat( p(1),td(1) )
  pc= pm

  IF (ABS(t(1)-td(1)) < 0.05) GO TO 45
  GO TO 25

  5   CONTINUE

  mrbar = 0.
  k = 0

  10   CONTINUE

  k = k+1
  IF (p(k) > pm) GO TO 10
  k = k-1
  j = k-1
  IF (j < 1) GO TO 20

!   compute the average mixing ratio....alog = natural log

  DO i= 1,j
    l = i+1
    mrsat1 = f_mrsat( p(i),td(i) )
    mrsat2 = f_mrsat( p(l),td(l) )
    mrbar = (mrsat1+mrsat2)*ALOG(p(i)/p(l)) + mrbar
  END DO

  20   CONTINUE
  l= k+1

!   estimate the dew point at pressure pm.

  mrsat1 = f_mrsat( p(k),td(k) )
  tq= td(k)+(td(l)-td(k))*(ALOG(pm/p(k)))/(ALOG(p(l)/p(k)))
  mrsat2 = f_mrsat( pm,tq )
  mrbar=  mrbar+(mrsat1+mrsat2 )*ALOG(p(k)/pm)
  mrbar=  mrbar/(2.*ALOG(p(1)/pm))


!   find level at which the mixing ratio line  mrbar crosses the
!   environmental temperature profile.

  25   CONTINUE

  DO j= 1,n
    i = n-j+1
    IF (p(i) < 300.) CYCLE

!   f_tmr = temperature (celsius) at pressure p (mb) along a mixing
!        ratio line given by  mrbar (g/kg)

    x = f_tmr(mrbar,p(i))-t(i)
    IF (x <= 0.) GO TO 35
  END DO

  f_pccl = 0.0

  RETURN

!  set up bisection routine

  35   l = i
  i = i+1
  del = p(l)-p(i)
  pc = p(i)+.5*del
  a = (t(i)-t(l))/ALOG(p(l)/p(i))
  DO j = 1,10
    del = del/2.
    x = f_tmr(mrbar,pc)-t(l)-a*(ALOG(p(l)/pc))

!   the sign function replaces the sign of the first argument
!   with that of the second.

    pc = pc+SIGN(del,x)
  END DO

  45   f_pccl = pc

  RETURN
  END FUNCTION f_pccl
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   FUNCTION F_CT                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION f_ct(mrbar,pc,ps)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This function returns the convective temperature ct (K)
!  given the mean mixing ratio mrbar (kg/kg) in the surface layer,
!  the pressure pc (Pa) at the convective condensation level (ccl)
!  and the surface pressure ps (Pa).
!
!-----------------------------------------------------------------------
!
!  AUTHOR: G.S. Stipanuk     1973            original version.
!
!    reference stipanuk paper entitled:
!         "algorithms for generating a skew-t, log p
!         diagram and computing selected meteorological
!         quantities."
!         atmospheric sciences laboratory
!         u.s. army electronics command
!         white sands missile range, new mexico 88002
!         33 pages
!    baker, schlatter  17-may-1982
!
!  MODIFICATION HISTORY:
!
!  01/30/1998 (Yuhe Liu)
!  Modified for ARPS to use the units of Kelvin, Pascal and kg/kg for
!  temperature, pressure and mixing ratio respectively
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    mrbar    Mean mixing ratio (kg/kg)
!    pc       Pressure at CCL (Pa)
!    ps       Pressure at surface (Pa)
!
!  OUTPUT:
!
!    f_ct     Convective temperature (K)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  REAL :: mrbar
  REAL :: pc
  REAL :: ps

  REAL :: f_ct

  REAL :: tc
!
!-----------------------------------------------------------------------
!
!  Include file
!
!-----------------------------------------------------------------------
!
  INCLUDE 'phycst.inc'
!
!-----------------------------------------------------------------------
!
!  Function f_tmr and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_tmr

!fpp$ expand (f_tmr)
!!dir$ inline always f_tmr
!*$*  inline routine (f_tmr)

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  compute the temperature (K) at the ccl.

  tc= f_tmr(mrbar,pc)
!
!-----------------------------------------------------------------------
!
!  compute the potential temperature (K) at ccl, i.e., the dry
!  adiabat through the ccl.
!
!  ptc = tc*((100000./pc)**rddcp)
!
!  and then compute the surface temperature on the same dry adiabat.
!
!  f_ct = ptc*((ps/100000.)**rddcp)
!
!-----------------------------------------------------------------------
!
  f_ct = tc*((ps/pc)**rddcp)

  RETURN
  END FUNCTION f_ct
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   FUNCTION F_TMR                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION f_tmr(p,mr)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  This function returns the temperature (Kelvin) on a mixing
!  ratio line mr (kg/kg) at pressure p (Pa). The formula is given in
!  table 1 on page 7 of Stipanuk (1973).
!
!-----------------------------------------------------------------------
!
!  AUTHOR: G.S. Stipanuk, 1973            original version.
!
!    reference stipanuk paper entitled:
!         "algorithms for generating a skew-t, log p
!         diagram and computing selected meteorological
!         quantities."
!         atmospheric sciences laboratory
!         u.s. army electronics command
!         white sands missile range, new mexico 88002
!         33 pages
!    baker, schlatter  17-may-1982
!
!  MODIFICATIONS:
!
!  01/30/1998 (Yuhe Liu)
!  Modified to use the units of Kelvin, Pascal and kg/kg for
!  temperature, pressure and mixing ratio respectively.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    p        Pressure (Pa)
!    mr       Mixing ratio (kg/kg)
!
!  OUTPUT:
!
!    f_tmr    Temperature (K)
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  REAL :: p              ! Pressure (Pa)
  REAL :: mr             ! Mixing ratio (kg/kg)

  REAL :: f_tmr          ! Temperature (K)
!
!-----------------------------------------------------------------------
!
!  Local variables
!
!-----------------------------------------------------------------------
!
  REAL :: x, pmb, mrg

  REAL :: c1,c2,c3,c4,c5,c6

  DATA c1/.0498646455/
  DATA c2/2.4082965/
  DATA c3/7.07475/
  DATA c4/38.9114/
  DATA c5/.0915/
  DATA c6/1.2035/
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  pmb = p/100.
  mrg = mr*1000.

  x= ALOG10(mrg*pmb/(622.+mrg))
  f_tmr= 10.**(c1*x+c2)-c3+c4*((10.**(c5*x)-c6)**2.)

  RETURN
  END FUNCTION f_tmr
