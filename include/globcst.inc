!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                     GLOBCST.INC                      ######
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
!  Include file 'globcst.inc' for ARPS
!
!  This file contains most of the model run control parameters
!  and parameters used by various physics packages.
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
!  5/03/1992 (M. Xue and K. Droegemeier)
!    Major restructuring. Some parameter names changed.
!
!  7/21/1992 (M. Xue)
!    Added parameters for energy statistics calculations, initial
!    random perturbation option, model run controls and additional
!    history dump format.
!
!  3/10/1993 (M. Xue)
!    Added parameters for advection options, and for terrain
!    specifications.
!
!  2/12/94 (Y. Liu)
!    Added parameters for the surface energy budget model.
!
!  9/10/94 (Weygandt & Y. Lu)
!    Cleaned up documentation.
!
!  1/28/95 (G. Bassett)
!    Added new parameter, buoyopt
!
!  02/07/1995 (Yuhe Liu)
!    Deleted the table data array, veg(14) and added an input
!    parameter, veg0.
!
!  2/2/96 (Donghai Wang & yuhe Liu)
!    Added parameters for the map projection factor.
!
!  03/26/96 (Yuhe Liu)
!    Added parameters for radiation computing.
!
!  4/2/96 (Donghai Wang, X. Song and M. Xue)
!    Added parameters for the vertically implicit mixing.
!
!  5/7/96 (Donghai Wang and M. Xue)
!     Added a new option for Rayleigh damping.
!
!  10/21/97 (Donghai Wang)
!     Added two new options, buoy2nd and rhofctopt.
!
!  03/13/2002 (Eric Kemp)
!     Added bmjcum_acct to time accounting flags.
!
!  02/11/02 (J. Brotzge, D. Weber)
!     Added an option for a choice of soil schemes.
!
!  03/19/02 (D. Weber, M. Xue & Xin Jin)
!     Added three new options, impfallopt, fallvalpha, fallvbeta
!
!  08/27/02 (D. Weber)
!     Added new option fallopt for using various fall velocity schemes
!     and update coefficients.
!
!-----------------------------------------------------------------------
!
  INTEGER :: mgrdmax         ! Maximum of mgrid
  PARAMETER ( mgrdmax = 100 )

  INTEGER :: mgrid           ! Grid identifier
  INTEGER :: nestgrd         ! Nested grid option
                             ! = 0, no grid nesting
                             ! = 1, grid nesting
  COMMON /arpsc001/ mgrid, nestgrd
!
!-----------------------------------------------------------------------
!
!  runname: A character string containing the pertinent information
!  for this run. The initial characters before a blank space, a comma
!  or otherwise the first 6 characters will be used to construct
!  history data and restart data dump file names.  This string is
!  written into the output data file.
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=80  ) :: runname  ! Name of this run
  INTEGER :: lfnkey                ! Length of string (no. of characters)
                                   ! to be used for naming files (file name
                                   ! keys).

  INTEGER :: nocmnt                ! Number of comment lines
  CHARACTER (LEN=80 ) :: cmnt(50)  ! String of comments on this model run

  COMMON /arpsc010/ runname, cmnt
  COMMON /arpsc011/ lfnkey, nocmnt
!
!-----------------------------------------------------------------------
!
!   Model initialization control parameters.
!
!-----------------------------------------------------------------------
!
  INTEGER :: runmod          ! Model configuration option.
                             ! = 1, model run in 3-D mode
                             ! = 2, model run in 2-D mode (x-z plane)
                             ! = 3, model run in 2-D mode (y-z plane)
                             ! = 4, model run in 1-D mode (z direction)

  INTEGER :: hxopt           ! Calculate prior observations option.
                             ! = 0, disable
                             ! = 1, enable
  INTEGER :: memid           ! The member''s number
  INTEGER :: hx_interval     ! Time interval for conducting hx calculation(unit:sec)

  COMMON /arpsc020/ runmod,hxopt,memid,hx_interval
!
!-----------------------------------------------------------------------
!
!  Model initialization control parameters.
!
!-----------------------------------------------------------------------
!
  INTEGER :: initopt         ! Model initialization option.
                             ! = 1, self initialization
                             !      (e.g. analytical field).
                             ! = 2, restart run, initialize from a
                             !      restart file.
                             ! = 3, initialize from an external data
                             !      set.

  INTEGER :: timeopt         ! Option to check the consistence between
                             ! the user specified time and data time
                             ! = 1, check, but continue with the
                             !      specified time (see initime and tstart)
                             ! = 2, check, but continue with data time
                             ! else, check and stop for inconsistence

  INTEGER :: inibasopt       ! Initialization option for base state
                             ! fields.
                             ! = 1, Initialize the base state using a
                             !      single sounding.
                             ! = 2, isentropic atmosphere
                             ! = 3, isothermal atmosphere
                             ! = 4, atmosphere with a constant static
                             !      stability.
                             ! = 5, an analytical thermodynamic
                             !      sounding profile after Weisman and
                             !      Klemp (1982, WMR)
                             ! For options rather than 1, additional
                             ! parameters are hardwired inside
                             ! subroutine INIBASE.

  INTEGER :: viniopt         ! Initialization option for base state
                             ! velocity when inibasopt.ne.1.
                             ! = 1, ubar= ubar0, vbar= vbar0
                             ! = 2, ubar0 and vbar0 are specified by user
                             !      inside code. By default, ubar0=vbar0=0.0.

  INTEGER :: soilinitopt     ! Initialization option for soil variables
                             ! = 1, integrate soil model using initial
                             !      initial atmospheric forcing

  INTEGER :: tsfcopt         ! Initialization option for skin temperature
                             ! = 0, no adjustment
                             ! = 1, adjustment by offset to air temp (k=2)
                             ! = 2, adjustment by iterating energy budget

  INTEGER :: pt0opt          ! Initial potential temperature
                             ! perturbation
                             ! option for initopt=1.
                             ! = 1, initial perturbation is
                             !      bubble-shaped,
                             ! = 2, initial perturbation has a random
                             !      distribution.
                             ! = 3, random initial perturbation symmetric
                             !      about x and y axis.
                             ! = 4, initial perturbation as given in
                             !      Skamarock and Klemp (1994)

  INTEGER :: inifmt          ! Format option for the history dump type
                             ! restart file
                             ! See hdmpfmt.
  CHARACTER (LEN=256) :: inifile ! Filename of t-dependent variable
                             ! initialization for initopt=3.
  CHARACTER (LEN=256) :: inigbf ! Filename of base&grid variable
                             ! initialization for initopt=3.
  INTEGER :: inifunt         ! Fortran I/O unit number for initial
                             ! data input.

  INTEGER :: buoyopt         ! Flag for turning off buoyancy terms
                             ! (buoyopt = 0 dissables buoyancy).

  INTEGER :: buoy2nd         ! Option for the second order terms in
                             ! the linerized buoyancy terms
                             ! = 1, including the 2nd-order terms
                             ! = 0, only the 1st-order terms

  INTEGER :: rhofctopt       ! Option for removing the density approximation
                             ! in the pressure gradient force(PGF) terms
                             ! = 1, removing the approximation
                             ! = 0, using the approximation

  INTEGER :: bsnesq          ! Bousinessq approximation?
                             ! = 1, yes
                             ! = 0, no

  COMMON /arpsc030/ initopt,inibasopt,viniopt,pt0opt,inifunt,           &
                    inifmt,soilinitopt,timeopt,tsfcopt
  COMMON /arpsc031/ inifile,inigbf,buoyopt,buoy2nd,rhofctopt,           &
                    bsnesq

  REAL :: soiltintv          ! Time interval in soil model initialization
                             ! integration

  INTEGER, PARAMETER :: maxpert=100
  REAL :: ptpert0(maxpert)            ! The magnitude of initial potential
                                      ! temperature perturbation (K)

  REAL :: pt0radx(maxpert)            ! The radius of bubble in x-direction.
  REAL :: pt0rady(maxpert)            ! The radius of bubble in y-direction.
  REAL :: pt0radz(maxpert)            ! The radius of bubble in z-direction.
  REAL :: pt0ctrx(maxpert)            ! The x coordinate of the bubble center.
  REAL :: pt0ctry(maxpert)            ! The y coordinate of the bubble center.
  REAL :: pt0ctrz(maxpert)            ! The z coordinate of the bubble center.

  REAL :: ubar0              ! User specified base state wind in x-dir.
  REAL :: vbar0              ! User specified base state wind in y-dir.

  COMMON /arpsc040/ ptpert0, pt0radx, pt0rady, pt0radz,                 &
                    pt0ctrx, pt0ctry, pt0ctrz, ubar0,vbar0,             &
                    soiltintv
!
!-----------------------------------------------------------------------
!
!  Coordinate transforamtion and terrain option parameters:
!
!-----------------------------------------------------------------------
!
  INTEGER :: crdtrns         ! Coordinate transformation option
                             ! = 0, no coordinate transforamtion at all
                             ! = 1, with coordinate transforamtion.

  INTEGER :: ternopt         ! Model terrain option.
                             ! = 0, no terrain, the ground is flat;
                             ! = 1, bell-shaped mountain;
                             ! = 2, terrain data read from terrain data
                             !      base (not implemented yet)
  INTEGER :: mntopt          ! Option for choosing idealized mountain type.
  REAL :: hmount             ! The mountain height (m)
  REAL :: mntwidx            ! The half-width of bell-shaped mountain in x-dir.
  REAL :: mntwidy            ! The half-width of bell-shaped mountain in y-dir.
  REAL :: mntctrx            ! The x coordinate of the bell-shaped
                             ! mountain center.
  REAL :: mntctry            ! The y coordinate of the bell-shaped
                             ! mountain center.
  CHARACTER (LEN=256) :: terndta  ! Name of the terrain data file

  COMMON /arpsc045/ crdtrns, ternopt, mntopt
  COMMON /arpsc046/ hmount, mntwidx, mntwidy, mntctrx, mntctry
  COMMON /arpsc047/ terndta
!
!-----------------------------------------------------------------------
!
!  Environmental sounding file name and Fortran I/O input unit.
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256  ) :: sndfile  ! Name of the input sounding file.
  INTEGER :: sfunit          ! Fortran I/O unit used by sounding input.

  COMMON /arpsc050/ sndfile
  COMMON /arpsc051/ sfunit
!
!-----------------------------------------------------------------------
!
!  Options for solution used in w equation
!
!  Acoustic/sound wave speed definition used in the pressure
!  equation.
!
!-----------------------------------------------------------------------
!
  INTEGER :: vimplct         ! Vertical implicit option for
                             ! the w and p equations.
                             ! = 1, implicit solution
                             ! = 0, explicit solution

  INTEGER :: ptsmlstp        ! Option to solve potential temperature
                             ! eqaution.
                             ! = 0, solve potential temperature eq.
                             !      outside small time steps
                             ! = 1, solve potential temperature eq.
                             !      within small time steps

  INTEGER :: peqopt          ! Option parameter of p-equation formulations
                             ! = 1 or 2 (1 is as in 4.0 User''s guide).


  INTEGER :: csopt           ! Sound wave speed option used in the
                             ! pressure equation.
                             ! = 1, csound = cp/cv * rd * T
                             ! = 2, csound = cp/cv * rd * T * csfactr
                             ! = 3, csound = specfied constant.
  REAL :: csfactr            ! A multiplication factor for the sound speed.
  REAL :: csound             ! User specified constant sound speed.

  COMMON /arpsc080/ vimplct, ptsmlstp,csopt, csfactr, csound,peqopt
!
!-----------------------------------------------------------------------
!
!  Model integration time step and integration control parameters.
!
!-----------------------------------------------------------------------
!
  REAL :: dtbig              ! Model integration time step.
  REAL :: dtsml              ! Small time step for the acoustic wave
                             ! integrations.

  REAL :: tstart             ! Time when the model is initialized.
  REAL :: tstop              ! Time when the integration is to be stoped.
  REAL :: curtim             ! Current model time.
  REAL :: tacoef             ! Time average weighting coefficient used
                             ! in the vertically implicit solver.

  INTEGER :: nstep           ! Current number of time steps taken.
  INTEGER :: nsteps          ! The total number of time steps in
                             ! this run.
  INTEGER :: nsmstp          ! The number of small acoustic wave time
                             ! integration steps per large time step.

  COMMON /arpsc090/ dtbig, tstart, tstop, curtim, dtsml, tacoef
  COMMON /arpsc100/ nstep, nsteps, nsmstp
!
!-----------------------------------------------------------------------
!
!  Parameters of the reference calendar day and real time
!
!-----------------------------------------------------------------------
!
  INTEGER :: year            ! Year of the reference calendar day
  INTEGER :: month           ! Month of the reference calendar day
  INTEGER :: day             ! Day of the reference calendar day
  INTEGER :: hour            ! Hours of the reference real time
  INTEGER :: minute          ! Minutes of the reference real time
  INTEGER :: second          ! Seconds of the reference real time

  INTEGER :: jday            ! Julian day starting from Jan. 1 of the
                             ! year

  COMMON /arpsc101/ year, month, day, hour, minute, second, jday
!
!-----------------------------------------------------------------------
!
!  Robert-Asselin time filter coefficient for leapfrog time
!  integration:
!
!-----------------------------------------------------------------------
!
  REAL :: flteps             ! Robert-Asselin time filter coefficient
                             ! (non-dimensional).

  COMMON /arpsc110/ flteps
!
!-----------------------------------------------------------------------
!
!  Option parameters for momentum and scalar advections.
!
!-----------------------------------------------------------------------
!
  INTEGER :: tintegopt       ! Option for time integration scheme
                             ! = 1, leapfrog scheme
                             ! = 2, Runge-Kutta 3rd-order (RK3) scheme

  INTEGER :: madvopt         ! Option parameter for momentum advection
                             ! = 1, 2nd order centerd;
                             ! = 2, 4th order centerd;
                             ! = 3, ...

  INTEGER :: sadvopt         ! Option parameter for scalar advection
                             ! = 1, 2nd order centerd;
                             ! = 2, 4th order centerd;
                             ! = 3, Flux-corrected transport (FCT)

  INTEGER :: fctorderopt     ! Option of the spactial order of accuracy of
                             ! the FCT advection (sadvopt=4) and MPDCD
                             ! advection schemes (sadvopt=5)
                             ! = 1  2nd order
                             ! = 2  4th order

  INTEGER :: fctadvptprt     ! Option for FCT advection for potential temperature.
                             ! Used only when sadvopt=4.
                             ! = 0, FCT scheme is applied to ptbar+ptprt.
                             ! = 1, FCT scheme is applied to ptprt only.
                             ! = 2, FCT scheme is applied to
                             !      ptbar+ptprt-min(ptbar+ptprt).

  COMMON /arpsc111/ tintegopt, madvopt, sadvopt, fctorderopt, fctadvptprt
!
!-----------------------------------------------------------------------
!
!  Subgrid-scale turbulent mixing parameters.
!
!-----------------------------------------------------------------------
!
  INTEGER :: tmixopt         ! Control parameter for turbulent mixing
                             ! options.
                             ! = 0, zero turbulent mixing
                             ! = 1, constant mixing coefficient,
                             !      km = tmixcst
                             ! = 2, Smagorinsky mixing coefficient.
                             ! = 3, Smagorinsky mixing coefficient
                             !      plus a constant coeffcient, tmixcst
  INTEGER :: trbisotp        ! Option for isotropic subgrid scale
                             ! turbulence.
                             ! = 0, the turbulence is assumed to be
                             !      anisotropic;
                             ! = 1, the turbulence is assumed to be
                             !      isotropic.
  INTEGER :: tkeopt          ! Option for 1.5 order TKE formulation
                             ! used by tmixopt=4
                             ! = 1, after Wyngaard;
                             ! = 2, after Deardroff;
                             ! = 3, after Sun and Chang.
  INTEGER :: trbvimp         ! Option parameter for the vertically
                             ! implicit mixing
                             ! = 0, vertical explicit
                             ! = 1, vertical implicit
  INTEGER :: tmixvert        ! Option parameter for the computing
                             ! only the vertical turb. mixing terms
                             ! = 0, all terms computed
                             ! = 1, only vertical terms computed.

  REAL    :: alfcoef         ! Time average weighting coefficient used
                             ! in the vertically implicit mixing

  REAL    :: prantl          ! Turbulent prantl number for Smagorinsky turb.
  REAL    :: tmixcst         ! Value of the constant mixing coefficient

  REAL    :: kmlimit         ! Constant used to set a limit of km
                             !

  COMMON /arpsc120/ tmixopt,trbisotp,tkeopt,trbvimp,tmixvert
  COMMON /arpsc121/ prantl, tmixcst, kmlimit,alfcoef
!
!-----------------------------------------------------------------------
!
!  Computational mixing parameters.
!
!-----------------------------------------------------------------------
!
  INTEGER :: cmix2nd         ! 2nd order computational mixing option.
                             ! = 0, mixing off
                             ! = 1, mixing on
  REAL :: cfcm2h             ! 2nd order horizontal computational mixing
                             ! coefficient caled by horizontal grid
                             ! spacing (1/s).
  REAL :: cfcm2v             ! 2nd order vertical computational mixing
                             ! coefficient scaled by vertical grid
                             ! spacing (1/s).
  REAL :: cfcmh2             ! Dimensional 2nd order horizontal
                             ! computational mixing coefficient (m**2/s)
  REAL :: cfcmv2             ! Dimensional 2nd order vertical
                             ! computational mixing coefficient (m**2/s)

  INTEGER :: cmix4th         ! 4th order computational mixing option
                             ! = 0, mixing off
                             ! = 1, mixing on
  REAL :: cfcm4h             ! 4th order horizontal computational mixing
                             ! coefficient scaled by horizontal grid
                             ! spacing (1/s).
  REAL :: cfcm4v             ! 4th order vertical computational mixing
                             ! coefficient scaled by vertical grid
                             ! spacing (1/s).
  REAL :: cfcmh4             ! Dimensional 4th order horizontal
                             ! computational mixing coefficient (m**2/s)
  REAL :: cfcmv4             ! Dimensional 4th order vertical
                             ! computational mixing coefficient (m**2/s)
  REAL :: scmixfctr          ! Reduction factor of the computational
                             ! mixing coefficients for scalars relative to
                             ! those of velocities.
  INTEGER :: cmix_opt, scmix_opt  !  options for monotonic computational
                             ! mixing.

  COMMON /arpsc130/ cmix2nd, cmix4th
  COMMON /arpsc131/ cfcm2h, cfcm2v, cfcm4h, cfcm4v
  COMMON /arpsc132/ cfcmh2, cfcmv2, cfcmh4, cfcmv4
  COMMON /arpsc133/ scmixfctr,cmix_opt, scmix_opt
!
!-----------------------------------------------------------------------
!
!  Upper level Rayleigh damping parameters.
!
!-----------------------------------------------------------------------
!
  INTEGER :: raydmp          ! Rayleigh damping option.
                             ! = 0, damping off.
                             ! = 1, damping the difference between the
                             !      total and base state fields.
                             ! = 2, damping the difference between the
                             !      total and external fields defined
                             !      in the EXBC file.
  REAL :: zbrdmp             ! Height of the bottom of Rayleigh damping
                             ! (sponge) layer (m).
  REAL :: cfrdmp             ! Rayleigh damping coefficient (1/second).
  INTEGER :: rayklow         ! Lowest level with Rayleigh damping

  COMMON /arpsc140/ zbrdmp, cfrdmp
  COMMON /arpsc141/ raydmp, rayklow
!
!-----------------------------------------------------------------------
!
!  Acoustic wave divergence damping parameters.
!
!-----------------------------------------------------------------------
!
  INTEGER :: divdmp          ! Acoustic wave divergence damping option.
                             ! = 0, divergence damping off
                             ! = 1, divergence damping on
  REAL :: cdvdmph            ! Divergence damping coefficient (m**2/s)
                             ! in horizontal
  REAL :: cdvdmpv            ! Divergence damping coefficient (m**2/s)
                             ! in vertical
  REAL :: divdmpndh          ! Non-dimensional divergence damping
                             ! coefficient in horizontal
  REAL :: divdmpndv          ! Non-dimensional divergence damping
                             ! coefficient in vertical

  COMMON /arpsc150/ divdmp
  COMMON /arpsc151/ cdvdmph,cdvdmpv, divdmpndh,divdmpndv

!-----------------------------------------------------------------------
!
!  Pollutant concentration parameters. (michi)
!
!-----------------------------------------------------------------------
!
!michi

  INTEGER :: ccin            ! Concentration option
                             ! = 0, no pollutant concentration
                             ! = 1, pollutant concentration (no input cc)
                             ! = 2, pollutant concentration

  INTEGER :: cpoint          ! Number of emitted place

  INTEGER :: icc(50)         ! Grid position of emitted place in x direction
  INTEGER :: jcc(50)         ! Grid position of emitted place in y direction
  INTEGER :: kcc(50)         ! Grid position of emitted place in z direction

  REAL    :: ccemit(50)         ! Emitted concentration
  INTEGER :: ccstart(50)        ! Emitted start time
  INTEGER :: ccend(50)          ! Emitted end time

  COMMON /arpsc155/  ccin,cpoint,icc,jcc,kcc,ccemit,ccstart,ccend

!michi

!
!-----------------------------------------------------------------------
!
!  Moisture/ice microphysics parameters.
!
!-----------------------------------------------------------------------
!
  INTEGER :: mphyopt         ! Microphysics option.
                             ! = 0, no microphysics.
                             ! = 1, Kessler warm rain microphysics.
                             ! = 2, Ice microphysics.

  INTEGER :: moist           ! Moist processes option.
                             ! = 0, dry run, all processes
                             !      related to moisture are turned off
                             ! = 1, moist processes are activated

  INTEGER :: nmphystp        ! Microphysics calls interval, in timesteps
                             ! = 1, called in each time step
                             ! > 1, called every nmphystp time steps.

  INTEGER :: graupel_ON      ! Option to suppress graupel initiation
                             ! for mphyopt = 8,9,10,11
                             ! = 0    Graupel initiation is suppressed.
                             ! = 1    Graupel is allowed.

  INTEGER :: hail_ON         ! Option to suppress hail initiation
                             ! for mphyopt = 8,9,10,11
                             ! = 0  Hail initiation is suppressed.
                             ! = 1  Hail initiation is allowed.


  INTEGER :: cnvctopt        ! Option for convective cumulus
                             ! parameterization
                             ! = 0, no convective parameterization
                             ! = 1, Kuo scheme
                             ! = 2, Kuo scheme and warm rain microphysics
                             ! = 3, Kain and Fristch
  INTEGER :: ice             ! Ice microphysics option.
                             ! = 0, ice processes are turned off
                             ! = 1, ice processes are activated.

  INTEGER :: impfallopt      ! option for vertically implicit scheme.
                             ! = 0, explicit scheme
                             ! = 1, implicit scheme

  REAL :: fallvalpha         ! weight coefficient for time level n+1.

  REAL :: fallvbeta          ! weight coefficient for time level n.

  INTEGER :: fallopt         ! option for fall velocity scheme/coeff.
                             ! = 1, exisiting Lin fomulation
                             ! = 2, Ferrier (1994) with updated fall
                             ! velocity coefficients from observations.

  INTEGER :: subsatopt       ! Sub-saturation option.
                             ! = 0, condensation occurs when RH>=100%
                             ! = 1, condensation occurs when RH>=rhsat,
                             !      a user specified value.
                             ! = 2, rhsat is defined as a linear function
                             !      of grid size instead.
                             !

  COMMON /arpsc160/ mphyopt, moist, nmphystp, graupel_ON, hail_ON,      &
                    cnvctopt, ice, subsatopt,                           &
                    impfallopt, fallvalpha, fallvbeta, fallopt

!-----------------------------------------------------------------------
!
!  Cumulus physics parameters.
!
!-----------------------------------------------------------------------
!
  REAL :: kffbfct            ! A feedback factor for K-F scheme.
  INTEGER :: kfsubsattrig    ! The sub-saturation trigger for for K-F scheme.
  REAL :: wcldbs             ! Vertical motion needed at cloud base for
                             ! convection.

  REAL :: confrq             ! Frequency of conv param. updates in sec.

  REAL :: qpfgfrq            ! Frequency of grid param. updates in sec.

  INTEGER :: idownd          ! Downdraft flag
                             ! = 0, no downdrafts.
                             ! = 1, simple downdraft model.

  COMMON /arpsc161/ kffbfct,kfsubsattrig,wcldbs,confrq,qpfgfrq,idownd

!
!-----------------------------------------------------------------------
!
!  Sub-saturation parameters.
!
!-----------------------------------------------------------------------
!
  REAL :: rhsat              ! Threshold of RH for condensation

  REAL :: rhsatmin           ! Minimum threshold of RH for condesation
  REAL :: dx_rhsatmin        ! Grid size corresponding to rhsatmin
  REAL :: dx_rhsat100        ! Grid size corresponding to RH =100%

  COMMON /arpsc162/ rhsat, rhsatmin, dx_rhsatmin, dx_rhsat100

!
!-----------------------------------------------------------------------
!
!  Microphysics constants.
!
!-----------------------------------------------------------------------
!
  REAL :: autotr             ! Cloud water mixing ratio threshold for
                             ! autoconversion from cloud to rain water
                             ! (kg/kg).
  PARAMETER ( autotr = 1.0E-3 )

  REAL :: autort             ! Autoconversion rate (1/seconds).
  PARAMETER ( autort = 1.0E-3 )

  REAL :: accrrt             ! Accretion rate (1/seconds).
  PARAMETER ( accrrt = 2.2 )

  REAL :: rho0               ! Reference density (kg/meter**3). )
  PARAMETER ( rho0 = 1.225 )

!
!-----------------------------------------------------------------------
!
!  Saturation specific humidity parameters used in enhanced Teten''s
!  formula. (See A. Buck, JAM 1981)
!
!-----------------------------------------------------------------------
!
  REAL :: satfwa, satfwb
  PARAMETER ( satfwa = 1.0007 )
  PARAMETER ( satfwb = 3.46E-8 )  ! for p in Pa

  REAL :: satewa, satewb, satewc
  PARAMETER ( satewa = 611.21 )   ! es in Pa
  PARAMETER ( satewb = 17.502 )
  PARAMETER ( satewc = 32.18 )

  REAL :: satfia, satfib
  PARAMETER ( satfia = 1.0003 )
  PARAMETER ( satfib = 4.18E-8 )  ! for p in Pa

  REAL :: sateia, sateib, sateic
  PARAMETER ( sateia = 611.15 )   ! es in Pa
  PARAMETER ( sateib = 22.452 )
  PARAMETER ( sateic = 0.6 )
!
!-----------------------------------------------------------------------
!
!  The x and y component of the domain translation speed
!  (ground-relative).
!
!-----------------------------------------------------------------------
!
  REAL :: umove              ! Model domain translation speed in x
                             ! direction.
  REAL :: vmove              ! Model domain translation speed in y
                             ! direction.
  COMMON /arpsc170/ umove, vmove
!
!-----------------------------------------------------------------------
!
!  Surface physics parameters.
!
!-----------------------------------------------------------------------
!
  INTEGER :: sfcphy          ! Surface physics option.
                             ! = 0, no surface physics parameterization
                             ! = 1, simple bulk aerodynamic formulation
  INTEGER :: landwtr         ! Flag to indicate if land and water will be
                             ! distinquished when the surface roughness
                             ! length is calculated.  If yes,the
                             ! vagetation data (in sfcdtfl) has to be
                             ! read in for all active options of sfcopt
                             ! = 0, no distinquishing.
                             ! = 1, surface roughness length is
                             !      calculated differently for land
                             !      and water.

  INTEGER :: cdhwtropt       ! Option to set constant value of cdh
                             ! (cdq) over water.

  INTEGER :: wtrexist        ! Flag to check if any grid point is
                             ! covered by water.


  REAL :: cdmlnd             ! Bulk aerodynamic momentum drag
                             ! coefficient.
  REAL :: cdmwtr             ! Bulk aerodynamic momentum drag
                             ! coefficient.
  REAL :: cdhlnd             ! Bulk aerodynamic heat flux drag
                             ! coefficient.
  REAL :: cdhwtr             ! Bulk aerodynamic heat flux drag
                             ! coefficient.
  REAL :: cdqlnd             ! Bulk aerodynamic moisture flux drag
                             ! coefficient.
  REAL :: cdqwtr             ! Bulk aerodynamic moisture flux drag
                             ! coefficient.
                             ! cdm, cdh and cdq are all non-dimensional
  INTEGER :: pbldopt         ! Option for determining PBL depth.
  REAL :: pbldpth0           ! User-specified PBL depth.
  REAL :: lsclpbl0           ! Length scale
  INTEGER :: sflxdis         ! Option for distributing fluxes in PBL
                             ! = 0, without flux distribution;
                             ! = 1, with flux distribution for heat,
                             !      moist, and momentum;
                             ! = 2, with flux distribution for heat, moist.
  INTEGER :: tqflxdis        ! Option for distributing heat and moisture
                             ! fluxes qudratically in a specified
                             ! depth dtqflxdis
                             ! = 0, no distribution;
                             ! = 1, with distribution
  REAL :: dtqflxdis          ! Depth of flux distribution for tqflxdis=1
                             ! 200 m recommended.

  INTEGER :: smthflx         ! Option to smooth surface fluxes
                             ! = 0, no smoothing
                             ! = 1, smoothing
  INTEGER :: numsmth         ! Number of smoothing passes
                             ! >=1 when smthflx = 1

  INTEGER :: sfcdiag         ! Flag of surface diagnotic calculation

  REAL :: ztime              ! Zulu time for any give day
                             ! Used to read Mesonet data files

  COMMON /arpsc180/ sfcphy,landwtr,cdhwtropt,sfcdiag,pbldopt,           &
                    wtrexist,sflxdis,tqflxdis,smthflx,numsmth,ztime
  COMMON /arpsc181/ cdmlnd,cdmwtr, cdhlnd,cdhwtr,                       &
                    cdqlnd,cdqwtr,pbldpth0,lsclpbl0,dtqflxdis

  INTEGER :: sfcdat          ! Surface data input flag.
  INTEGER :: soilinit        ! Soil model varialble initialization flag.
  INTEGER :: sfcunit         ! Fortran I/O unit number for surface data
                             ! input.
  REAL :: dtsfc              ! Time step for surface energy budget model
  INTEGER :: nsfcst          ! Number  of steps for surface energy budget model

  COMMON /arpsc182/ sfcdat, soilinit, sfcunit, nsfcst
  COMMON /arpsc183/ dtsfc

  INTEGER :: soilmodel_forced ! Option for forcing soil scheme
                             ! = 0, Non-forced mode
                             ! = 1, Forced mode (e.g.,Mesonet forcing)

  INTEGER :: soilmodel_option  ! Option for soil scheme
                             ! = 1, Two-layer force-restore method
                             ! = 2, Multi-layer 'OUSoil' scheme

  INTEGER :: nstyp           ! Number of soil types
  INTEGER :: styp            ! Soil type
  INTEGER :: vtyp            ! Vegetation type
  REAL :: lai0               ! Leaf Area Index
  REAL :: roufns0            ! Surface roughness
  REAL :: veg0               ! Vegetation fraction

  REAL :: tsoilint(100)      ! Soil temperature
  REAL :: qsoilint(100)      ! Soil moisture

  REAL :: ptslnd0            ! Initial ground-level soil potential
                             ! temperature over land
  REAL :: ptswtr0            ! Initial ground-level soil potential
                             ! temperature over water
  REAL :: tsoil0             ! Initial deep soil temperature
  REAL :: wetsfc0            ! Initial surface soil moisture
  REAL :: wetdp0             ! Initial deep soil moisture
  REAL :: wetcanp0           ! Initial canopy moisture
  REAL :: snowdpth0          ! Snow depth (m)

  REAL :: ttprt              ! Offset of tsfc from sfc air temperature
  REAL :: tbprt              ! Offset of tsoil from sfc air temperature
  REAL :: wgrat              ! Saturation rate of sfc soil moisture
  REAL :: w2rat              ! Saturation rate of deep soil moisture

  INTEGER :: prtsoilflx

  COMMON /arpsc184/ soilmodel_forced,soilmodel_option,nstyp,styp,vtyp,  &
                    prtsoilflx

  COMMON /arpsc185/ lai0,roufns0,veg0,tsoilint,qsoilint,                &
                    ptslnd0,ptswtr0,tsoil0,                             &
                    wetsfc0,wetdp0,wetcanp0,snowdpth0,                  &
                    ttprt,tbprt,wgrat,w2rat

  CHARACTER (LEN=256 ) :: sfcdtfl  ! File name of surface data.
  CHARACTER (LEN=256 ) :: sitemeso ! File name of site file name.
  CHARACTER (LEN=256 ) :: siteflux ! File name of site file name.
  CHARACTER (LEN=256 ) :: siternet ! File name of site file name.
  CHARACTER (LEN=256 ) :: sitesoil ! File name of site file name.
  CHARACTER (LEN=256 ) :: siteveg ! File name of site file name.

  CHARACTER (LEN=256 ) :: soilinfl  ! File name of surface initialization data.

  COMMON /arpsc186/ sfcdtfl, sitemeso, siteflux, siternet, sitesoil,    &
                    siteveg,  soilinfl
!
!-----------------------------------------------------------------------
!
!  Coriolis parameters
!
!-----------------------------------------------------------------------
!
  INTEGER :: coriopt         ! Coriolis term option.
                             ! = 0, no coriolis terms
                             ! = 1, coriolis terms involving horizontal wind
                             ! = 2, coriolis terms involving both
                             !      horizontal and vertical wind
                             ! = 3, as 1, but Coriolis parameters are altitude
                             !      dependent and effects of spatial gradient
                             !      of map factor and of earth curvature are
                             !      also included
                             ! = 4, as 2, but Coriolis parameters are altitude
                             !      dependent and effects of spatial gradient
                             !      of map factor and of earth curvature are
                             !      also included


  INTEGER :: earth_curvature ! Earth curvature term option
                             ! = 0, earth curvature terms are neglected
                             ! = 1, those terms are included when coriopt = 3 or 4

  INTEGER :: coriotrm        ! A flag indicating if total or perturbation
                             ! wind is used in the Coriolis terms
                             ! in the horizontal momemtum equations.
                             ! = 1, total wind,
                             ! = 2, perturbation wind.

  REAL :: latitud            ! Model center latitude (degrees).
  REAL :: longitud           ! Model center longitude (degrees).

  COMMON /arpsc190/ coriopt, coriotrm, earth_curvature
  COMMON /arpsc191/ latitud, longitud

!----------------------------------------------------------------------
!
! The amplitude of the annual cycle of the difference in deep and skin
! layer soil seasonal-mean temperatures.
!
!----------------------------------------------------------------------

  INTEGER :: tsoil_offset
  REAL :: tsoil_offset_amplitude

  COMMON /arpsc197/ tsoil_offset, tsoil_offset_amplitude

!
!-----------------------------------------------------------------------
!
!  The name of the directory into which output files will be written:
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=256  ) :: dirname  ! The name of output directory
  INTEGER :: ldirnam       ! The length of the directory name string

  COMMON /arpsc200/ ldirnam
  COMMON /arpsc201/ dirname
!
!-----------------------------------------------------------------------
!
!  Model output parameters.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nfmtprt         ! time steps    between formatted print out
  REAL    :: tfmtprt         ! time interval between formatted print out

  INTEGER :: nmaxmin         ! time steps    between max/min statistics calc.
  REAL    :: tmaxmin         ! time interval between max/min statistics calc.

  INTEGER :: nenergy         ! time steps    between energy statistics calc.
  REAL    :: tenergy         ! time interval between energy statistics calc.

  INTEGER :: nrstout         ! time steps    between restart data dumps
  REAL    :: trstout         ! time interval between restart data dumps

  COMMON /arpsc210/ nfmtprt, nmaxmin, nenergy, nrstout
  COMMON /arpsc220/ tfmtprt, tmaxmin, tenergy, trstout

  INTEGER :: totout          ! Flag to indicate that the history dump
                             ! contains total values instead of
                             ! perturbation
                             ! = 0, Perturbations
                             ! = 1, total values

  INTEGER :: grdout          ! Grid output option.
                             ! = 0, no grid output
                             ! = 1, grid output
  INTEGER :: basout          ! Base state field output option.
                             ! = 0, no base state fields output
                             ! = 1, base state fields output
  INTEGER :: varout          ! Mass & perturbation wind output option.
                             ! = 0, no mass or perturbation wind output
                             ! = 1, mass & perturbation wind output
  INTEGER :: mstout          ! Moist variable output option.
                             ! = 0, no moisture variables output
                             ! = 1, qv, qc, qr, qi, qs and qh output.
  INTEGER :: rainout         ! Rain variable output option.
                             ! = 0, no rain variables output
                             ! = 1, raing and rainc output.
  INTEGER :: prcout          ! Precipitation rates output option.
                             ! = 0, no rain variables output
                             ! = 1, prcrate output.
  INTEGER :: iceout          ! Ice variable ouptu option.
                             ! = 0, no ice variables output
                             ! = 1, qi, qs and qh output
  INTEGER :: tkeout          ! TKE output option.
                             ! = 0, no output
                             ! = 1, output
  INTEGER :: trbout          ! Turbulence field (km) output option.
                             ! = 0, km not output
                             ! = 1, km output
  INTEGER :: sfcout          ! Surface variable output option.
                             ! = 0, no surface variables output
                             ! = 1, tsfc,tsoil,wetsfc,wetdp,wetsfc output
  INTEGER :: snowout         ! Snow cover output option
                             ! = 0, no snow cover output
                             ! = 1, snow cover output
  INTEGER :: landout         ! Surface property array output option.
                             ! = 0, no surface property array output
                             ! = 1, soiltyp,vegtyp,lai,roufns,veg output
  INTEGER :: radout          ! Radiation arrays output option.
                             ! = 0, no radiation arrays output
                             ! = 1, radfrc, radsw, rnflx output
  INTEGER :: flxout          ! Surface fluxes output option.
                             ! = 0, no surface fluxes output
                             ! = 1, usflx, vsflx, ptsflx, qvsflx output

  INTEGER :: exbcdmp         ! Flag to dump the ARPS fields into an
                             ! additional file to make the ARPS
                             ! external boundary data set (if lbcopt=2).
                             ! = 0, no EXBC dumping;
                             ! = 1, EXBC dumping;
                             ! = 3, HDF4 format dumping;

  INTEGER :: exbchdfcompr    ! HDF4 compression option for exbcdmp = 3.


  INTEGER :: extdadmp        ! Flag to dump the fields that contains
                             ! external data array interpolated to the
                             ! current model time. When lbcopt.ne.2,
                             ! reset extdadmp to 0.
                             ! = 0, no dumping
                             ! = 1, dump.

  INTEGER :: sfcdmp          ! Flag to dump the ARPS surface data into
                             ! a surface data file.
                             ! = 0, no surface data file
                             ! = 1, dump in Fortran unformatted format
                             ! = 3, dump in HDF4 uncompresses format.

  INTEGER :: soildmp         ! Flag to dump the ARPS soil data into
                             ! a soil data file.
                             ! = 0, no soil data file
                             ! = 1, dump in Fortran unformatted format
                             ! = 3, dump in HDF4 uncompresses format.

  INTEGER :: terndmp         ! Flag to dump the ARPS terrain data into
                             ! a terrain data file.
                             ! = 0, no terrain data file
                             ! = 1, dump in Fortran unformatted format
                             ! = 3, dump in HDF4 uncompresses format.

  INTEGER :: exbcfmt         ! Format of EXBC file to be read in.
                             ! = 1, Fortran unformatted
                             ! = 3, HDF4 (uncompressed or compressed).

  INTEGER :: incrfmt         ! Format of nudging file to be read in.
                             ! = 1, Fortran unformatted
                             ! = 3, HDF4 (uncompressed or compressed).

  INTEGER :: sfcfmt          ! Format of surface data file to be read in.
                             ! = 1, Fortran unformatted
                             ! = 3, HDF4 (uncompressed).

  INTEGER :: soilfmt         ! Format of soil data file to be read in.
                             ! = 1, Fortran unformatted
                             ! = 3, HDF4 (uncompressed).

  INTEGER :: ternfmt         ! Format of terrain data file to be read in.
                             ! = 1, Fortran unformatted
                             ! = 3, HDF4 (uncompressed).

  INTEGER :: qcexout         ! Flag to dump qc into EXBC file
                             ! = 0, no EXBC dumping
                             ! = 1, EXBC dumping if exbcdmp \= 0

  INTEGER :: qrexout         ! Flag to dump qr into EXBC file
                             ! = 0, no EXBC dumping
                             ! = 1, EXBC dumping if exbcdmp \= 0

  INTEGER :: qiexout         ! Flag to dump qi into EXBC file
                             ! = 0, no EXBC dumping
                             ! = 1, EXBC dumping if exbcdmp \= 0

  INTEGER :: qsexout         ! Flag to dump qs into EXBC file
                             ! = 0, no EXBC dumping
                             ! = 1, EXBC dumping if exbcdmp \= 0

  INTEGER :: qhexout         ! Flag to dump qh into EXBC file
                             ! = 0, no EXBC dumping
                             ! = 1, EXBC dumping if exbcdmp \= 0

  INTEGER :: qgexout         ! Flag to dump qg into EXBC file
                             ! = 0, no EXBC dumping
                             ! = 1, EXBC dumping if exbcdmp \= 0

  INTEGER :: nqexout         ! Flag to dump concentration number of hydrometero
                             ! particals into EXBC file
                             ! = 0, no EXBC dumping
                             ! = 1, EXBC dumping if exbcdmp \= 0

  INTEGER :: zqexout         ! Flag to dump reflectivity of hydrometero
                             ! particals into EXBC file
                             ! = 0, no EXBC dumping
                             ! = 1, EXBC dumping if exbcdmp \= 0

  INTEGER :: grdbasfout      ! grid & base output file flag

  COMMON /arpsc230/ grdout, basout, varout, mstout, rainout,            &
                    prcout, iceout, tkeout, trbout, sfcout,             &
                    landout,radout, flxout, totout, snowout,            &
                    exbcdmp,exbchdfcompr,extdadmp,                      &
                    qcexout,qrexout,qiexout,qsexout,qhexout,            &
                    qgexout,nqexout,zqexout,                            &
                    sfcdmp,soildmp,terndmp,exbcfmt,incrfmt,sfcfmt,      &
                    soilfmt,ternfmt,grdbasfout

!
!-----------------------------------------------------------------------
!
!  History data dump parameters.
!
!-----------------------------------------------------------------------
!
  INTEGER :: hdmpopt         ! History data dump option.
                             ! = 1, dump with constant time interval
                             ! = 2, dump at the time specified by user

  INTEGER :: hdmpfmt         ! History data dump format option.
                             ! = 0, No data dump is produced
                             ! = 1, Unformatted binary data dump
                             ! = 2, Formatted ascii data dump
                             ! = 3, NCSA HDF4 format data dump
                             ! = 4, Compressed binary data dump
                             ! = 5, Dump for Savi3D visualization
                             !      package
                             ! = 6, binary allowing data point skip
                             ! = 7, Uncompressed NetCDF format dump.
                             ! = 8, Compressed NetCDF format dump.
                             ! = 9, GrADS data dump
                             ! = 10, GRIB data dump

  INTEGER :: grbpkbit        ! Bit length to pack floating data
                             ! in GRIB format

  INTEGER :: hdfcompr        ! HDF4 compression option

  REAL :: thisdmp            ! time interval between history data dumps
  INTEGER :: nhisdmp         ! time steps    between history data dumps

  REAL :: tstrtdmp           ! time to start history data dump
  INTEGER :: nstrtdmp        ! time steps from beginning to tstrtdmp

  INTEGER :: numhdmp         ! number of history dumps specified by user
  INTEGER :: hdmpmax         ! maximum number of history dumps allowed
                             ! in user specifying option
  PARAMETER ( hdmpmax = 100 )

  REAL :: hdmptim(hdmpmax)   ! times of model specified by user to
                             ! dump out history data
  INTEGER :: hdmpstp(hdmpmax)
                             ! steps from beginning of model specified
                             ! by user to dump out history data
                             ! hdmpstp = nint( hdmptim/dtbig )

  INTEGER :: nchdmp          ! Fortran I/O channel for history
                             ! data dump.
  CHARACTER (LEN=256) :: hdmpfn  ! History data dump file name.
  INTEGER :: ldmpf           ! Length of the history data dump file name

  INTEGER :: istager         ! Flag to dump vectors on vector grids.
                             ! = 0, no stagering, dump on scalar grid
                             ! = 1, stagering, dump vectors on vectors grids

  COMMON /arpsc240/ hdmpfn
  COMMON /arpsc241/ ldmpf,nchdmp
  COMMON /arpsc242/ grbpkbit,istager,hdfcompr
  COMMON /arpsc243/ hdmpopt,hdmpfmt,nhisdmp,nstrtdmp,numhdmp,hdmpstp
  COMMON /arpsc244/ thisdmp,tstrtdmp,hdmptim
!
!-----------------------------------------------------------------------
!
!  More output control parameters
!
!-----------------------------------------------------------------------
!


  INTEGER :: imgopt          ! HDF image dumping option, 0 or 1
  INTEGER :: nimgdmp         ! time steps    between HDF image dumps
  REAL    :: timgdmp         ! time interval between HDF image dumps

  INTEGER :: pltopt          ! Graphic plotting option, 0 or 1
  INTEGER :: nplots          ! time steps    between graphic plotting calls
  REAL    :: tplots          ! time interval between graphic plotting calls
  INTEGER :: filcmprs        ! Option for output file compression, 0 or 1.
  INTEGER :: readyfl         ! Option for using a ready file, 0 or 1.

  COMMON /arpsc305/ imgopt, pltopt
  COMMON /arpsc310/ nimgdmp, nplots
  COMMON /arpsc320/ timgdmp, tplots
  COMMON /arpsc321/ filcmprs, readyfl
!
!-----------------------------------------------------------------------
!
!  Parameters for grid translation.
!
!-----------------------------------------------------------------------
!
  INTEGER :: cltkopt         ! Option for cell tracking.
                             ! = 0, no cell-tracking,
                             ! = 1, cell-tracking is turned on.
  INTEGER :: grdtrns         ! Domain translation option.
                             ! = 0, no change in domain translation.
                             ! = 1, domain translation based on cell
                             !      tracking.
                             ! = 2, automatic domain translation via
                             !      optimal procedure.
  REAL :: chkdpth            ! Depth (m) over which the optimal domain
                             ! translation speed is calculated for the
                             ! automatic domain translation.
  REAL :: twindow            ! Time window (s) between update of grid
                             ! translation components for automatic
                             ! domain translation.
  INTEGER :: nceltrk         ! time steps    between cell tracking calls
  REAL :: tceltrk            ! time interval between cell tracking calls
  REAL :: tcrestr            ! time required to restore cell center
                             !               to domain center

  COMMON /arpsc330/ cltkopt,grdtrns,nceltrk
  COMMON /arpsc331/ chkdpth,twindow,tceltrk,tcrestr
!
!-----------------------------------------------------------------------
!
!  Debug information printing level.
!
!-----------------------------------------------------------------------
!
  INTEGER :: lvldbg          ! Level of debug information printing.
                             ! =0, no printing;
                             ! =1, Print variables in big t-step;
                             ! =2, Print forcings in big t-step;
                             ! =3, Print variables in small t-step;
                             ! =4, Print forcings in small t-step;
                             ! =5, Print individual forcing terms
                             !     and misc. info;

  COMMON /arpsc325/ lvldbg
!
!-----------------------------------------------------------------------
!
!  Model restart data dump parameters.
!
!-----------------------------------------------------------------------
!
  INTEGER :: restrt          ! Restart run option.
                             ! = 0, non-restart run;
                             ! = 1, restart run;

  INTEGER :: rstiunt         ! Fortran I/O unit for restart data read.
  CHARACTER (LEN=256   ) :: rstinf  ! File name of restart input data.

  INTEGER :: rstount         ! Fortran I/O unit for restart data dump.
  CHARACTER (LEN=256  ) :: rstoutf  ! File name of restart data dump.

  COMMON /arpsc260/ restrt, rstount, rstiunt
  COMMON /arpsc270/ rstoutf
  COMMON /arpsc271/ rstinf
!
!-----------------------------------------------------------------------
!
!  The Fortran I/O channel (unit) used for max/min statistics output.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nchmax          ! Fortran I/O channel for max/min
                             ! statistics output.
  COMMON /arpsc280/ nchmax
!
!-----------------------------------------------------------------------
!
!  The Fortran I/O channel (unit) used by energy statistics output
!
!-----------------------------------------------------------------------
!
  INTEGER :: ncheng          ! Fortran I/O channel for energy
                             ! statistics output.
  COMMON /arpsc290/ ncheng
!
!-----------------------------------------------------------------------
!
!  Parameters used by Savi3D.
!
!-----------------------------------------------------------------------
!
  INTEGER :: grafh           ! Grafhandle used by savi3D
  INTEGER :: gridh           ! Gridhandle used by savi3D
  COMMON /arpsc300/ grafh, gridh

  INTEGER :: dsindex
  INTEGER :: gridid
  COMMON /savicmn/dsindex, gridid
!
!-----------------------------------------------------------------------
!
!  Parameters for map projections
!
!-----------------------------------------------------------------------
!
  INTEGER :: mpfctopt        ! Option parameter for map factor
                             ! = 0, map factor set to 1
                             ! = 1, map factor calculated according to
                             !      mapproj

  INTEGER :: mptrmopt        ! Option parameter for map factor terms in
                             ! momentum advection
                             ! = 0, ignore the terms
                             ! = 1, include the terms

  INTEGER :: maptest         ! Option parameter for map factor testing
                             ! with a symmetric field.
                             ! = 0, no test
                             ! = 1, symmetric field applied
  COMMON  /arps450/ mpfctopt, mptrmopt, maptest
!
!-----------------------------------------------------------------------
!
!  Parameters for radiation
!
!-----------------------------------------------------------------------
!
  INTEGER :: radopt          ! Option to switch on/off radiation physics
                             ! = 0, No radiation physics;
                             ! = 1, Radiation physics;

  INTEGER :: radstgr         ! Option to switch on/off radiation
                             ! computing at staggering points
                             ! = 0, No staggering
                             ! = 1, staggering

  INTEGER :: rlwopt          ! Option for longwave schemes
                             ! = 0, for high = .false.
                             ! = 1, for high = .true.

  !augustin
  INTEGER :: radshade        ! Option for shade computing
                             ! = 0, no shading
                             ! = 1, shading.
                             ! = 2, idealized shading
                             ! used in simulations of 2D valleys
                             ! north south oriented the shade will
                             ! be computed only in the middle of the valley


  INTEGER :: raddiag         ! Option to dump variables calculated in
                             ! radiation transfer equations
                             ! = 0, no dump
                             ! = 1, dump to GrADS format

  INTEGER :: nradstp         ! Time steps to update radiation
                             ! calculation

  REAL :: dtrad              ! Time interval to update radiation
                             ! calculation

  REAL :: cldh2m             ! Height in meter that separates high cloud
                             ! and middle cloud
  REAL :: cldm2l             ! Height in meter that separates middle cloud
                             ! and low cloud

  INTEGER :: ict             ! Vertical indices of height cldh2m
  INTEGER :: icb             ! Vertical indices of height cldm2l

  PARAMETER ( cldh2m = 5500.0, cldm2l = 3000.0 )   ! in meter
  !augustin

  COMMON  /arps500/ radopt,radstgr, rlwopt,radshade,raddiag,nradstp,      &
                    ict,icb
  COMMON  /arps501/ dtrad

  CHARACTER (LEN=20) :: arpsversion  ! A string containing the version number of ARPS.
  COMMON /arps601/ arpsversion

  !
  !  Accumulated CPU time usage by various processes or terms in the
  !  equations.
  !

  INTEGER, PARAMETER :: max_acct = 22
  INTEGER            :: use_acct

  INTEGER, PARAMETER :: misc_acct=1,init_acct=2,output_acct=3,   &
       rad_acct=4,sfcphy_acct=5,soil_acct=6,tinteg_acct=7,       &
       tmix_acct=8,cmix_acct=9,raydmp_acct=10,advuvw_acct=11,    &
       advs_acct=12,coriol_acct=13,buoy_acct=14,tkesrc_acct=15,  &
       bc_acct=16,qpfgrd_acct=17,cum_acct=18,smlstp_acct=19,     &
       microph_acct=20,qfall_acct=21,mp_acct=22
  ! Note: can add other accounting sets for adas, etc.

  DOUBLE PRECISION :: acct_cpu_time(max_acct), acct_wall_time(max_acct)
  DOUBLE PRECISION :: wall_init, cpu_init, total_cpu, total_wall, &
                      wall_acct, cpu_acct, wall_inter, cpu_inter

  COMMON /arps701/ acct_cpu_time, acct_wall_time,                &
                   wall_init, cpu_init, total_cpu, total_wall,   &
                   wall_acct, cpu_acct, wall_inter, cpu_inter

  INTEGER :: current_acct, interrupt_acct

  COMMON /arps702/ use_acct, current_acct, interrupt_acct

  REAL :: ptground,pttrop, ttrop, htrop, qvmixed, rhmixed, mixtop, zshear
  COMMON /wk_sounding/ ptground, pttrop, ttrop, htrop, qvmixed,         &
                       rhmixed,  mixtop, zshear

!-----------------------------------------------------------------------
!
! Microphysics scalar index
!
!-----------------------------------------------------------------------

  INTEGER :: nscalar, nscalarq
  INTEGER :: P_QC, P_QR, P_QI, P_QS, P_QH, P_QG
  INTEGER :: P_NC, P_NR, P_NI, P_NS, P_NH, P_NG
  INTEGER ::       P_ZR, P_ZI, P_ZS, P_ZH, P_ZG
  INTEGER :: P_CC

  CHARACTER(LEN=40) :: qnames(20)  ! Note that Q, N, Z has special meaning (especially in tinteg)
  CHARACTER(LEN=40) :: qdescp(20)  ! so do not use them for newly added variables

  COMMON /arps_scalar/ nscalar, nscalarq, qnames, qdescp,               &
                       P_QC, P_QR, P_QI, P_QS, P_QH, P_QG,              &
                       P_NC, P_NR, P_NI, P_NS, P_NG, P_NH,              &
                             P_ZR, P_ZI, P_ZS, P_ZG, P_ZH,              &
                       P_CC
