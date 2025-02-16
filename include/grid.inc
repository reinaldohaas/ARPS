!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                       GRID.INC                       ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Include file 'grid.inc' for ARPS
!
!  This file contains grid location and map projection parameters.
!
!  These parameters are allocated in named common blocks and are
!  accessible to subroutines that include this file.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  10/01/1991
!
!  MODIFICATION HISTORY:
!
!  2000/04/28 (Gene Bassett)
!    Split off from globcst.inc.
!
!  02/15/02 (J. Brotzge, D. Weber)
!  Added soil grid stretching option
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!  Model grid parameters.
!
!-----------------------------------------------------------------------
!
  REAL :: dx                 ! Grid spacing in x-direction in
                             ! computational as well as in physical
                             ! space (m)
  REAL :: dy                 ! Grid spacing in y-direction in
                             ! computational as well as in physical
                             ! space (m)

  INTEGER :: strhopt         ! Vertical grid stretching option.
                             ! = 0, no stretching in vectical.
                             ! = 1, vertical stretching with function
                             !      dzp=a+b*k**3 being used to specify
                             !      dz as a function of k.
                             ! = 2, vertical stretching with function
                             !      dzp=c+a*tanh(b*k) being used to
                             !      specify dz as a function of k.

  REAL :: dz                 ! Averaged grid spacing in vertical
                             ! direction in transformed computataional
                             ! space (m).
  REAL :: dzmin              ! Minimun grid spacing in vertical
                             ! direction in physcal space (m).
  REAL :: zrefsfc            ! The reference height of the surface
                             ! (ground level) (m)

  REAL :: dlayer1            ! The depth of the lower layer with uniform
                             ! (dz=dzmin) vertical spacing (m)

  REAL :: dlayer2            ! The depth of the mid layer with stetched
                             ! vertical spacing (m)

  REAL :: strhtune           ! Tuning parameter for stretching option 2
                             ! A Value between 0.2 and 5.0 is appropriate.
                             ! A larger value gives a more linear stretching.

  INTEGER :: soilstrhopt     ! Vertical grid stretching option within soil.
                             ! = 0, no stretching in vectical.
                             ! = 1, vertical stretching with function
                             !      dzp=a+b*k**3 being used to specify
                             !      dz as a function of k.
                             ! = 2, vertical stretching with function
                             !      dzp=c+a*tanh(b*k) being used to
                             !      specify dz as a function of k.

  REAL :: dzsoil             ! Averaged grid spacing in vertical
                             ! direction in transformed computataional
                             ! space (m).

  REAL :: soildzmin          ! Minimun grid spacing in vertical
                             ! direction in physcal space (m).

  REAL :: zrefsoil           ! The reference height of the surface
                             ! (below ground level) (m)

  REAL :: soildlayer1        ! The depth of the lower layer with uniform
                             ! (dz=dzmin) vertical spacing (m)

  REAL :: soildlayer2        ! The depth of the mid layer with stetched
                             ! vertical spacing (m)

  REAL :: soilstrhtune       ! Tuning parameter for stretching option 2
                             ! A Value between 0.2 and 5.0 is appropriate.
                             ! A larger value gives a more linear stretching.

  REAL :: zflat              ! The height at which the grid levels
                             ! becomes flat in the terrain-following
                             ! coordinate transformation (m).

  REAL :: dxinv              ! Reciprocal of dx (i.e. 1/dx)
  REAL :: dyinv              ! Reciprocal of dy (i.e. 1/dy)
  REAL :: dzinv              ! Reciprocal of dz (i.e. 1/dz)
  REAL :: dzsoilinv          ! Reciprocal of dzsoil (i.e. 1/dzsoil)
  REAL :: dzsoilinv2         ! Square of the reciprocal of dzsoil

  COMMON /arpsc060/ strhopt, soilstrhopt 

!  COMMON /arpsc061/ dx, dy, dz, nzsoil, dzsoil, dzmin, soildzmin, &
!                    dxinv, dyinv, dzinv, dzsoilinv, dzsoilinv2

  COMMON /arpsc061/ dx, dy, dz, dzsoil, dzmin, soildzmin, &
                    dxinv, dyinv, dzinv, dzsoilinv, dzsoilinv2

  
  COMMON /arpsc062/ zrefsfc,zrefsoil,dlayer1,dlayer2,soildlayer1,     &
                    soildlayer2,strhtune,soilstrhtune, zflat

  REAL :: xl                 ! Length of the physical model domain in
                             ! x-direction (m). In uniform spacing case
                             ! ,xl = (nx-3)*dx.
  REAL :: yl                 ! Length of the physical model domain in
                             ! y-direction (m). In uniform spacing case
                             ! ,yl = (ny-3)*dy.
  REAL :: zh                 ! Length of the physical model domain in
                             ! z-direction (m). In uniform spacing case
                             ! , zh = (nz-3)*dz.

  COMMON /arpsc070/ xl, yl, zh

  REAL :: ctrlat             ! Center latitude of the model domain (deg. N)
  REAL :: ctrlon             ! Center longitude of the model domain (deg. E)
  INTEGER :: crdorgnopt      ! Option for origin of the model physical domain
  REAL :: xorig              ! x-coordinate origin of the model
                             ! physical domain
  REAL :: yorig              ! y-coordinate origin of the model
                             ! physical domain
  REAL :: zorig              ! z-coordinate origin of the model
                             ! physical domain

  COMMON /arpsc071/ ctrlat, ctrlon, crdorgnopt, xorig, yorig, zorig

  REAL :: xgrdorg            ! x-coordinate of model grid in
                             ! absolute space
  REAL :: ygrdorg            ! y-coordinate of model grid in
                             ! absolute space
  COMMON /arpsc075/ xgrdorg, ygrdorg
!
!-----------------------------------------------------------------------
!
!  Parameters for map projections
!
!-----------------------------------------------------------------------
!
  REAL :: eradius         ! Mean earth radius in meters
  PARAMETER ( eradius = 6371000. )
!
  INTEGER :: mapproj         ! Type of map projection in the model
                             ! grid:
                             ! modproj = 0  No projection.
                             !         = 1  Polar Stereographic
                             !              projection.
                             !         = 2  Mercator projection.
                             !         = 3  Lambert projection.

  REAL :: trulat1         ! 1st real true latitude of map projection
  REAL :: trulat2         ! 2nd real true latitude of map projection
  REAL :: trulon          ! Real true longitude of map projection
  REAL :: sclfct          ! Map scale factor.
                          ! Distance on map, between two latitudes
                          ! trulat1 and trulat2,
                          ! is = (Distance on earth)*sclfct
                          ! For ARPS model runs,
                          ! generally this is 1.0

  REAL :: swlatu          ! latitude  at southwest corner for u-grid
  REAL :: swlonu          ! longitude at southwest corner for u-grid
  REAL :: nelatu          ! latitude  at northeast corner for u-grid
  REAL :: nelonu          ! longitude at northeast corner for u-grid

  REAL :: swlatv          ! latitude  at southwest corner for v-grid
  REAL :: swlonv          ! longitude at southwest corner for v-grid
  REAL :: nelatv          ! latitude  at northeast corner for v-grid
  REAL :: nelonv          ! longitude at northeast corner for v-grid

  REAL :: swlats          ! latitude  at southwest corner for s-grid
  REAL :: swlons          ! longitude at southwest corner for s-grid
  REAL :: nelats          ! latitude  at northeast corner for s-grid
  REAL :: nelons          ! longitude at northeast corner for s-grid

  COMMON  /arps400/ mapproj
  COMMON  /arps401/ trulat1, trulat2, trulon, sclfct
  COMMON  /arps402/ swlatu, swlonu, nelatu, nelonu,                     &
                    swlatv, swlonv, nelatv, nelonv,                     &
                    swlats, swlons, nelats, nelons
