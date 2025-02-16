!
!-----------------------------------------------------------------------
!
!  Include file 'sfcphycst.inc' for ARPS
!
!  This file initializes values of some physical constants
!  that remain unchanged throughout a model run
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Vince Wong and X. Song
!  09/01/1994
!
!  MODIFICATION HISTORY:
!  02/07/96 (V.Wong and X.Song)
!  Added a lower limiter, zolimit, for zo, and an upper limiter, z1limit,
!  for depth of the surface layer z1.
!
!  05/01/97 (V. Wong and X. Tan)
!  Added ztz0, ratio of thermal roughness length to momentum
!  roughness length.
!
!-----------------------------------------------------------------------
!
!
!
!-----------------------------------------------------------------------
!
!  Constants used in Businger's equations
!
!-----------------------------------------------------------------------
!
  REAL :: kv        ! Von Karman constant
  REAL :: kvwtr     ! Von Karman constant over the sea
  PARAMETER ( kv = 0.4 , kvwtr = 0.35)
!
!  Land case (See Byun's paper)
!
  REAL :: prantl0l  ! Turbulent Prandtl number for neutral stability,
                    ! defined as 0.74 in Byun's paper.
  PARAMETER ( prantl0l = 0.74 )

  REAL :: gammahl
  REAL :: gammaml
  REAL :: betahl
  REAL :: betaml
  PARAMETER (gammaml=15.0,                                              &
             gammahl=9.0,                                               &
             betaml=4.7,                                                &
             betahl=betaml/prantl0l)
!
!  Sea case (Given by V. Wong)
!
  REAL :: prantl0w  ! Turbulent Prandtl number for neutral stability,
                    ! defined as 1.
  PARAMETER ( prantl0w = 1.0 )

  REAL :: gammahw
  REAL :: gammamw
  REAL :: betahw
  REAL :: betamw
  PARAMETER (gammamw=16.0,                                              &
             gammahw=gammamw,                                           &
             betamw=7.0,                                                &
             betahw=betamw/prantl0w)
!
!-----------------------------------------------------------------------
!
!  Following constants are hot-wired when calling sfcphy3d.f
!
!-----------------------------------------------------------------------
!
  REAL :: oned3,oned9,oned54
  PARAMETER (oned3  = 1.0/3.0, oned9  = 1.0/9.0, oned54 = 1.0/54.0)

  REAL :: c1l,c2l,c3l,c4l,c5l
  PARAMETER (c1l= 1./(gammaml*gammaml),                                 &
             c2l=3.0*(gammahl/gammaml),                                 &
             c3l=-2.0*c1l/gammaml,                                      &
             c4l=9.0*(-gammahl/gammaml+3)/gammaml,                      &
             c5l=1.0/(3.0*gammaml))

  REAL :: c1w,c2w,c3w,c4w,c5w
  PARAMETER (c1w= 1./(gammamw*gammamw),                                 &
             c2w=3.0*(gammahw/gammamw),                                 &
             c3w=-2.0*c1w/gammamw,                                      &
             c4w=9.0*(-gammahw/gammamw+3)/gammamw,                      &
             c5w=1.0/(3.0*gammamw))
!
!-----------------------------------------------------------------------
!
!  Define a constant blimit=0.5; Note that it was defined to be 0.9
!  in Deardroff's paper.
!
!-----------------------------------------------------------------------
!
  REAL :: blimit
  PARAMETER(blimit=0.50)
!
!-----------------------------------------------------------------------
!
!  Define a minimal total wind speed for the calculation of surface
!  heat and moisture fluxes.
!
!-----------------------------------------------------------------------
!
  REAL :: vsfcmin
  PARAMETER ( vsfcmin = 0.1 )
!
!-----------------------------------------------------------------------
!
!  Define constants zolimit and z1limit, which is the lower and upper
!  limits for the roughness zo and the depth of the surface layer z1.
!
!-----------------------------------------------------------------------
!
  REAL :: zolimit,z1limit
  PARAMETER(zolimit=1.5E-5,z1limit=300.0)
!
!-----------------------------------------------------------------------
!
!  Define constant ztz0, zt = ztz0*z0
!
!-----------------------------------------------------------------------
!
  REAL :: ztz0
  PARAMETER(ztz0=1.0)
