!
!##################################################################
!##################################################################
!######                                                      ######
!######                     ASSIM.INC                        ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!---------------------------------------------------------------------------
!
!  Include file 'ASSIM.INC' for variational adjustment & velocity
!  insertion/dynamic retrieval with host model ARPS40.
!
!
!  Variational adjustment & velocity insertion/dynamic retrieval
!  parameters are allocated in a named common block. These parameters
!  are accessible to all subroutines that include this file.
!
!--------------------------------------------------------------------------
!
!
!  AUTHOR: Alan Shapiro and Steve Lazarus
!  9/2/1993.
!
!  MODIFICATION HISTORY:
!
!  01/03/96 (Limin Zhao)
!  Added control parameters for blending ADAS data with Doppler
!  retrieval(common block /blnd/).
!
!  11/15/96 (Limin Zhao)
!  Added control parameters for sensitivity tests(Common block /opts/).
!
!-------------------------------------------------------------------------
!
  INTEGER :: assimopt  ! Option to run ARPS data assimilation
                       ! instead of ARPS model

  REAL :: xshift  ! x-coordinate of radar relative to the model grid (0,0,0)
  REAL :: yshift  ! y-coordinate of radar relative to the model grid (0,0,0)
  REAL :: zshift  ! z-coordinate of radar relative to the model grid (0,0,0)

  INTEGER :: nvf      ! Maximum number of input data files (minimum of 3).

  INTEGER :: dtyp     ! Input data file type
                      ! = 0, input files are model history data files
                      ! = 1, input files are pre-processed Lincoln Lab data
                      ! = 2, input files are other (user specified)

  INTEGER :: ii       ! Counter indicating input data file name and time

  INTEGER :: recovopt ! A switch to turn on dynamic retrieval
                      ! = 0, no retrieval.
                      ! = 1, insert velocities and retrieve p, T.

  INTEGER :: varopt   ! A switch to turn on variational adjustment
                      ! = 0, no variational adjustment.
                      ! = 1, insert velocities and perform variational
                      !   adjustment.

  INTEGER :: insrtopt ! A switch to turn on the direct velocity insertion
                      ! = 0, no direct insertion.
                      ! = 1, direct insertion.


  INTEGER :: ldirnm   ! Length of the directory name where input data resides

  INTEGER :: irecov   ! Recovery flag
                      ! = 0, insertion/recovery is not performed this time step
                      ! = 1, insertion/recovery is performed this time step

  INTEGER :: ivar     ! Variational Adjustment flag
                      ! = 0, no variational adjustment is performed this time
                      !   step
                      ! = 1, variational adjustment is performed this time step

  INTEGER :: insrt    ! Data insertion flag
                      ! = 0, no direct insertion performed this time step
                      ! = 1, direct insertion performed this time step

  CHARACTER (LEN=255) :: assimnm       ! Name of recovery logfile

  CHARACTER (LEN=255) :: dirnam        ! Directory from which the input data is read

  CHARACTER (LEN=25) :: assimdat(100)  ! Names of the input data files

  CHARACTER (LEN=255 ) :: adjdat(100)  ! Names of the adjusted model data files

!     character*80 gbfile        ! Name of grid/base state file (model
                                 ! input only)

  CHARACTER (LEN=255) :: adasdat(100)  ! Name of ADAS files
!
!-----------------------------------------------------------------------
!
!  Control parameters for blending.
!
!-----------------------------------------------------------------------
!
  INTEGER :: ip_wgt,nadas

  REAL :: v_err, vfill_err, adas_err,spval
!
!-----------------------------------------------------------------------
!
!  Options for sensitivity tests.
!
!-----------------------------------------------------------------------
!
  INTEGER :: bc_opt,itest,ig
!
!-----------------------------------------------------------------------
!
  COMMON /assimrecvry01/ recovopt,varopt,insrtopt,                      &
                         irecov,ivar,insrt,                             &
                         dtyp,nvf,ii,ldirnm
  COMMON /assimrecvry02/ xshift,yshift,zshift
  COMMON /assimrecvry03/ assimnm, dirnam, assimdat, adjdat

  COMMON /assimblend01/ ip_wgt,nadas
  COMMON /assimblend02/ v_err, vfill_err, adas_err, spval
  COMMON /assimblend03/ adasdat

  COMMON /assimopts01/ assimopt, bc_opt, itest, ig

  REAL :: voltim1, voltim2, voltim3
  COMMON /assimvoltim02/ voltim1, voltim2, voltim3

  CHARACTER (LEN=4) :: radid
  REAL :: latrad,lonrad,elvrad

  COMMON /assimradloc02/ latrad,lonrad,elvrad
  COMMON /assimradloc03/ radid
