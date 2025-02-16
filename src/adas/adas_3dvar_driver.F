PROGRAM adas_3dvar
!
!##################################################################
!##################################################################
!######                                                      ######
!######                   PROGRAM ADAS_3DVAR                 ######
!######             ARPS Data Analysis System                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!  PURPOSE:
!
!  This program is a driver for original arps data analysis system
!  (ADAS) and three-dimensional variational data assimilation
!  (3DVAR) on the ARPS grid.
!
!  To run adas:
!    makearps adas
!    bin/adas < arps.input
!  To run arps3dvar
!    makearps arps3dvar
!    bin/arps3dvar < arps.input
!
!  The driver establishes sizing of variables
!  and adjustable parameters.
!  It calls a routine to set up the grid, calls the data readers,
!  a routine to analyse the data and another to write it out.
!
!  AUTHOR:
!
!  For ADAS:
!  Keith Brewster, CAPS, March, 1994
!  OLAPS sfc analysis Keith Brewster, CAPS, March, 1994
!
!  FOR 3DVAR:
!  06/10/2000 (Jidong Gao)
!            Basic version, shared with ADAS input & output
!
!  MODIFICATION HISTORY:
!  01/17/96  (KB) Named changed to ADAS, dimensions imported from
!            adas.dims and input variables from namelist
!            rather than parameter statements in this code.
!  11/27/96  (KB) Cleaned-up to remove unused arrays.
!  01/25/96  (KB) Added tem variables passed to anxiter to accomodate
!            new version of anxiter to speed-up execution.
!  06/15/97  (KB) Update for version 4.3.0, new arguments to call
!            to INITVAR.
!  07/16/97  (KB) Added nt to calling list of INITVAR and made nt
!            equal to one to save memory.
!  05/26/98  (KB) Modified logic in superadiabatic check. (kb)
!  08/31/98  (KB) Added code for nudging to official version (kb)
!  09/22/98  Jason Levit
!            Updated call in initgrdvar for arps4.4.0, including
!            change of dimension for qcumsrc.
!  12/07/98  Added background-based surface correlations to
!            account for cumulus-parameterization small-scale changes.
!  12/20/98  Changed moisture analysis variable to qv.
!            Error stats are still in RH, converted to qv based on
!            background value.
!  12/30/98  Added processing of qr from retrieval data.
!
!  2000/04/14 (Gene Bassett)
!            Merged the ADAS input file back into the ARPS input file.
!
!  2000/05/09 (Gene Bassett)
!            Converted to F90, creating allocation and adas main
!            subroutines.
!
!  2000/07/28 (Ming Xue)
!            Converted to F90 free format. Did away with adaasamin for
!            better compilation efficiency. Use ALLOCATABLE instead of
!            POINTER allocation to avoid double memory usage.
!
!  07/23/2001 (K. Brewster)
!            Added dynamic array allocation for the IAU increment arrays.
!            Made corresponding changes to the increment handling
!            routines to account for the array allocation.
!
!  08/07/2001 (K. Brewster)
!            Added a preliminary calculation of vertical velocity
!            before the cloud analyis (before cmpclddrv call).
!
!  11/01/2001 (K. Brewster)
!            Renamed radar_ref_3d array to ref_mos_3d since its a mosaic.
!            Removed potentially very large 4-d temporary arrays for
!            radar by reorganizing logic for radar mosaic creation.
!
!  04/10/2002 (K. Brewster)
!            Added code for soil surface temperature adjustment.
!            Includes allocation and deallocation of arrays temporarily
!            needed for radiation calculations.
!
!  05/16/2002 (K. Brewster)
!            Added BMJ arrays to initgrdvar call to conform with
!            latest updates to that subroutine.  Added deallocate after
!            the call to initgrdvar for some arrays not needed after
!            that call.
!
!  06/10/2002 (Eric Kemp)
!            Addition of new soil model variables.
!
!  04/09/2003 (K. Brewster)
!            Modified I/O of source errors to stop if I/O problem.
!            Correction of problem in super-obbing.
!
!  03/01/2005 (S. Leyton and D. Weber)
!            Add MPI capabilities.
!
!  12/01/2005 (K. W. Thomas)
!            Update MPI capabilities.
!
!  10/20/2005 (Fanyou Kong)
!            Correct a bug that the subroutine calls of bcu, bcv, lbcw,
!            and bcsclr have inconsistent arguement list to the
!            corresponding subroutines.
!
!  03/13/2007 (Kevin W. Thomas)
!            Update code as a consequence of MPI'ing ADAS.  Note that
!            ARPS3DVAR isn't MPI'd yet, even though some of the hooks
!            are now present.
!  11/05/2008 (K. W. Thomas)
!            MPI radial velocity.  Code cleanup.  Fix a few MPI bugs.
!
!  04/30/2010 (Keith Brewster)
!            Added deallocate of radar column arrays before cloud code
!            to save memory.
!
!  11/17/2011 (Keith Brewster)
!            Added diagnostic output file, adasstn, containing locations
!            of all single and upper air observations used.   Modified
!            the statistics reported in adasstat diagnostic file.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
#ifdef arps3dvar
#define nvar_anx_MODE nvar_anx_3DVAR
#define nvar_rad_MODE nvar_rad_3DVAR
  USE Vmodule_3dvar
#else
#define nvar_anx_MODE nvar_anx_ADAS
#define nvar_rad_MODE nvar_rad_ADAS
#endif

  USE module_arpsArray
  USE module_soilArray
  USE module_anaIncArray
  USE module_analysisArrays

  IMPLICIT NONE
  INCLUDE 'alloc.inc'
  INCLUDE 'adas.inc'
  INCLUDE 'nudging.inc'
  INCLUDE 'exbc.inc'
  INCLUDE 'radcst.inc'
  !INCLUDE '3dvarcputime.inc'
  INCLUDE 'mp.inc'     ! Message passing parameters

  REAL :: cputimearray(10)
  COMMON /cputime3dvar/ cputimearray

  REAL :: exbcbuf( 1 ) ! EXBC buffer array (unused)
!
!-----------------------------------------------------------------------
!
!  Analysis scratch space
!
!-----------------------------------------------------------------------
!
  INTEGER, ALLOCATABLE :: knt(:,:)
  REAL,    ALLOCATABLE :: rngsqi(:)
  REAL,    ALLOCATABLE :: wgtsum(:,:)
  REAL,    ALLOCATABLE :: zsum(:,:)
  REAL,    ALLOCATABLE :: sqsum(:,:)
!
!-----------------------------------------------------------------------
!
!  Additional grid variables
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: xs(:), xslg(:)
  REAL, ALLOCATABLE :: ys(:), yslg(:)
  REAL, ALLOCATABLE :: zs(:,:,:)
  REAL, ALLOCATABLE :: latgr(:,:)
  REAL, ALLOCATABLE :: longr(:,:)
!
!-----------------------------------------------------------------------
!
!  Temporary arrays.  Some are only used for MPI
!
!-----------------------------------------------------------------------
!
  REAL,    ALLOCATABLE :: tem1  (:,:,:)     ! Temporary work array.
  REAL,    ALLOCATABLE :: tem2  (:,:,:)     ! Temporary work array.
  REAL,    ALLOCATABLE :: tem3  (:,:,:)     ! Temporary work array.
  REAL,    ALLOCATABLE :: tem4  (:,:,:)     ! Temporary work array.
  REAL,    ALLOCATABLE :: tem5  (:,:,:)     ! Temporary work array.
  REAL,    ALLOCATABLE :: tem6  (:,:,:)     ! Temporary work array.
  REAL,    ALLOCATABLE :: tem7  (:,:,:)     ! Temporary work array.
  REAL,    ALLOCATABLE :: tem8  (:,:,:)     ! Temporary work array.
  REAL,    ALLOCATABLE :: tem9  (:,:,:)     ! Temporary work array.
  REAL,    ALLOCATABLE :: tem10 (:,:,:)     ! Temporary work array.
  REAL,    ALLOCATABLE :: tem11 (:,:,:)     ! Temporary work array.
  INTEGER, ALLOCATABLE :: item1(:)          ! Temporary work array.
  INTEGER, ALLOCATABLE :: ibegin(:)
  INTEGER, ALLOCATABLE :: iend(:)

!
! These arrays are continually cleaned up, so they get recycled frequently
! between different data types.
!

  INTEGER, ALLOCATABLE :: tems1di(:)
  INTEGER, ALLOCATABLE :: temr1di(:)

  REAL,    ALLOCATABLE :: tems1dr(:)
  REAL,    ALLOCATABLE :: temr1dr(:)

  INTEGER, ALLOCATABLE :: tems2di(:,:)
  INTEGER, ALLOCATABLE :: temr2di(:,:)

  REAL,    ALLOCATABLE :: tems2dr(:,:)
  REAL,    ALLOCATABLE :: temr2dr(:,:)

  INTEGER, ALLOCATABLE :: tems3ditmp(:,:,:)
  INTEGER, ALLOCATABLE :: temr3ditmp(:,:,:)

  REAL,    ALLOCATABLE :: tems3drtmp(:,:,:)
  REAL,    ALLOCATABLE :: temr3drtmp(:,:,:)

!
! Single level arrays.
!
  REAL, ALLOCATABLE :: tems2drsng(:,:)
  REAL, ALLOCATABLE :: temr2drsng(:,:)

!
! Multi-level arrays.
!
  REAL, ALLOCATABLE :: tems2drua(:,:)
  REAL, ALLOCATABLE :: temr2drua(:,:)

  REAL, ALLOCATABLE :: tems3drua(:,:,:)
  REAL, ALLOCATABLE :: temr3drua(:,:,:)

!
! Radar arrays.
!
  REAL, ALLOCATABLE :: tems3dirad(:,:,:)
  REAL, ALLOCATABLE :: temr3dirad(:,:,:)

  REAL, ALLOCATABLE :: tems3drrad(:,:,:)
  REAL, ALLOCATABLE :: temr3drrad(:,:,:)
!
!-----------------------------------------------------------------------
!
!  Radiation variables only needed for adjtsfc
!
!-----------------------------------------------------------------------
!
  REAL, ALLOCATABLE :: rsirbm(:,:), rsirdf(:,:), rsuvbm(:,:), rsuvdf(:,:)
  REAL, ALLOCATABLE :: cosz(:,:), cosss(:,:)
  REAL, ALLOCATABLE :: fdirir(:,:), fdifir(:,:), fdirpar(:,:), fdifpar(:,:)
  REAL, ALLOCATABLE :: radbuf(:), sh(:,:)           !added by DTD

  REAL, ALLOCATABLE :: temrad1 (:,:,:)   ! Temporary work array.
  REAL, ALLOCATABLE :: tem12 (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem13 (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem14 (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem15 (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem16 (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem17 (:,:,:)     ! Temporary work array. added by DTD
  REAL, ALLOCATABLE :: tem18 (:,:,:)     ! Temporary work array.

  REAL, ALLOCATABLE :: tem1soil (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem2soil (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem3soil (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem4soil (:,:,:)     ! Temporary work array.
  REAL, ALLOCATABLE :: tem5soil (:,:,:)     ! Temporary work array.
!
!-----------------------------------------------------------------------
!
!  Cloud analysis variables other than model's ones:
!
!-----------------------------------------------------------------------
!
  REAL,    ALLOCATABLE :: ref_mos_3d (:,:,:)
                     ! 3D remapped radar reflectivity field as input
                     ! may be modified by sate. VIS as output
  REAL,    ALLOCATABLE :: refmax_mos_3d (:,:,:)
  REAL,    ALLOCATABLE :: rhv_mos_3d (:,:,:)
  REAL,    ALLOCATABLE :: zdr_mos_3d (:,:,:)
  REAL,    ALLOCATABLE :: kdp_mos_3d (:,:,:)
  REAL,    ALLOCATABLE :: dist_mos_3d (:,:,:)
  REAL,    ALLOCATABLE :: w_cld (:,:,:)      ! vertical velocity (m/s) in clouds
  INTEGER, ALLOCATABLE :: cldpcp_type_3d(:,:,:)
                                       ! cloud/precip type field

  INTEGER, ALLOCATABLE :: icing_index_3d(:,:,:) ! icing severity index
  LOGICAL, ALLOCATABLE :: l_mask_pcptype(:,:) ! analyze precip type
!
!-----------------------------------------------------------------------
!
!  Analysis variables
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: PTR_U  = 1, PTR_V = 2, PTR_P = 3, PTR_PT = 4,   &
                        PTR_QV = 5, PTR_W = 6

  REAL, ALLOCATABLE :: anx(:,:,:,:)
!
!-----------------------------------------------------------------------
!
!  Cross-category correlations
!
!-----------------------------------------------------------------------
!
  REAL,    ALLOCATABLE :: xcor(:,:)
  INTEGER, ALLOCATABLE :: icatg(:,:)
!
!-----------------------------------------------------------------------
!
!  ARPS dimensions:
!
!  nx, ny, nz: Dimensions of analysis grid.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nx       ! Number of grid points in the x-direction
  INTEGER :: ny       ! Number of grid points in the y-direction
  INTEGER :: nz       ! Number of grid points in the z-direction
  INTEGER :: nzsoil   ! Number of soil levels

  INTEGER :: nxlg     ! Number of grid points in the x-direction large grid
  INTEGER :: nylg     ! Number of grid points in the y-direction large grid

  INTEGER :: nxndg    ! Number of x grid points for IAU
  INTEGER :: nyndg    ! Number of y grid points for IAU
  INTEGER :: nzndg    ! Number of z grid points for IAU
!
!-----------------------------------------------------------------------
!
!  Soil types.
!
!-----------------------------------------------------------------------
!
  INTEGER :: nstyps            ! Number of soil types
!
!-----------------------------------------------------------------------
!
!  Misc
!
!-----------------------------------------------------------------------
!
  INTEGER :: nt                ! Number of time levels of data
  INTEGER :: ncat              ! Number of correlation categories
  REAL    :: rlim              ! Largest possible radius of influence
  REAL :: rlimv                ! Same, but for radial velocity only

!-----------------------------------------------------------------------
!
!  Analysis parameters
!
!-----------------------------------------------------------------------
!
  REAL,    PARAMETER :: rhmin=0.05  ! rh safety net value to prevent neg qv
  INTEGER, PARAMETER :: iqspr=3     ! Use qobs of pstn to combine x,y,elev
  INTEGER, PARAMETER :: klim= 1     ! Min of one other sfc station for QC
!
!-----------------------------------------------------------------------
!
!  Indices of specific observation variables
!
!-----------------------------------------------------------------------
!
  INTEGER, PARAMETER :: ipres=3   ! Pressure
  INTEGER, PARAMETER :: iptmp=4   ! Potential Temperature
  INTEGER, PARAMETER :: iqv=5     ! Specific Humidity
!
!-----------------------------------------------------------------------
!
!  ARPS cloud cover analysis variables:
!
!-----------------------------------------------------------------------
!
! Remapping parameters of the remapped radar data grid
!
  INTEGER :: strhopt_rad                  ! streching option
  INTEGER :: mapproj_rad                  ! map projection indicator
  REAL :: dx_rad,dy_rad,dz_rad,dzmin_rad  ! grid spcngs
  REAL :: ctrlat_rad,ctrlon_rad           ! central lat and lon
  REAL :: tlat1_rad,tlat2_rad,tlon_rad    ! true lat and lon
  REAL :: scale_rad                       ! map scale factor
!
!-----------------------------------------------------------------------
!
!  Common block that stores remapping parameters for the radar
!  data file.
!
!-----------------------------------------------------------------------
!
  COMMON /remapfactrs_rad/ strhopt_rad,mapproj_rad
  COMMON /remapfactrs_rad2/ dx_rad,dy_rad,dz_rad,dzmin_rad,             &
           ctrlat_rad,ctrlon_rad,                                       &
           tlat1_rad,tlat2_rad,tlon_rad,scale_rad
!
!-----------------------------------------------------------------------
!
!  ARPS include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'phycst.inc'
  INCLUDE 'bndry.inc'
  INCLUDE 'grid.inc'        ! Grid parameters
  INCLUDE 'adjust.inc'
!
!-----------------------------------------------------------------------
!
!  Surface Station variables
!
!-----------------------------------------------------------------------
!
!  INTEGER :: isrcsng(mx_sng)
!  INTEGER :: icatsng(mx_sng)
!  REAL :: latsng(mx_sng,ntime)
!  REAL :: lonsng(mx_sng,ntime)
!  REAL :: hgtsng(mx_sng,ntime)
!  REAL :: xsng(mx_sng)
!  REAL :: ysng(mx_sng)
!  REAL :: trnsng(mx_sng)
!  INTEGER :: timesng(mx_sng,ntime)
!  CHARACTER (LEN=5) :: stnsng(mx_sng,ntime)
!
!-----------------------------------------------------------------------
!
!  Surface (single-level) read-in observation variables
!
!-----------------------------------------------------------------------
!
!  CHARACTER (LEN=8) :: wx(mx_sng,ntime)
!  CHARACTER (LEN=8) :: csrcsng(mx_sng,ntime)
!  CHARACTER (LEN=1) :: store_emv(mx_sng,5,ntime)
!  CHARACTER (LEN=4) :: store_amt(mx_sng,5,ntime)
!  INTEGER :: kloud(mx_sng,ntime),idp3(mx_sng,ntime)
!  REAL :: store_hgt(mx_sng,5,ntime)
!  REAL :: obrdsng(mx_sng,nvar_sng,ntime)
!
!  REAL :: obsng  (nvar_anx_MODE,mx_sng)
!  REAL :: odifsng(nvar_anx_MODE,mx_sng)
!  REAL :: oanxsng(nvar_anx_MODE,mx_sng)
!#ifdef arps3dvar
!  REAL :: corsng(mx_sng,nvar_anx_3dvar)
!#else
!  REAL :: corsng(nvar_anx_adas,mx_sng)
!#endif
!  REAL :: thesng(mx_sng)
!  INTEGER :: ival(mx_sng)
!
!-----------------------------------------------------------------------
!
!  Upper Air Station variables
!
!-----------------------------------------------------------------------
!
!  INTEGER :: isrcua(mx_ua)
!  REAL :: xua(mx_ua)
!  REAL :: yua(mx_ua)
!  REAL :: trnua(mx_ua)
!  REAL :: elevua(mx_ua)
!  REAL :: hgtua(nz_ua,mx_ua)
!  INTEGER :: nlevsua(mx_ua)
!  CHARACTER (LEN=5) :: stnua(mx_ua)
!
!-----------------------------------------------------------------------
!
!  Upper-air observation variables
!
!-----------------------------------------------------------------------
!
!  REAL :: obsua (nvar_anx_MODE,nz_ua,mx_ua)
!  REAL :: odifua(nvar_anx_MODE,nz_ua,mx_ua)
!  REAL :: oanxua(nvar_anx_MODE,nz_ua,mx_ua)
!#ifdef arps3dvar
!  REAL :: corua(nz_ua,mx_ua,nvar_anx_3dvar)
!#else
!  REAL :: corua(nvar_anx_adas,nz_ua,mx_ua)
!#endif
!  REAL :: theua(nz_ua,mx_ua)
!
!-----------------------------------------------------------------------
!
!  Radar site variables
!
!-----------------------------------------------------------------------
!
!  CHARACTER (LEN=5) :: stnrad(mx_rad)
!  INTEGER :: isrcrad(0:mx_rad)
!  REAL :: latrad(mx_rad),lonrad(mx_rad)
!  REAL :: elvrad(mx_rad)
!
!-----------------------------------------------------------------------
!
!  Radar observation variables
!
!-----------------------------------------------------------------------
!
!  INTEGER :: irad(mx_colrad)
!  INTEGER :: nlevrad(mx_colrad)
!  REAL :: latradc(mx_colrad)
!  REAL :: lonradc(mx_colrad)
!  REAL :: xradc(mx_colrad)
!  REAL :: yradc(mx_colrad)
!  REAL :: trnradc(mx_colrad)
!  REAL :: distrad(mx_colrad)
!  REAL :: uazmrad(mx_colrad)
!  REAL :: vazmrad(mx_colrad)
!  REAL :: hgtradc(nz_rdr,mx_colrad)
!  REAL :: theradc(nz_rdr,mx_colrad)
!  REAL :: obsrad(nvar_radin,nz_rdr,mx_colrad)
!
!  REAL :: odifrad(nvar_rad_MODE,nz_rdr,mx_colrad)
!  REAL :: oanxrad(nvar_rad_MODE,nz_rdr,mx_colrad)
!#ifdef arps3dvar
!  REAL :: corrad(nz_rdr,mx_colrad,nvar_rad_3dvar)
!#else
!  REAL ::  corrad(nvar_rad_adas,nz_rdr,mx_colrad)
!#endif
!  REAL ::    dsdr(nz_rdr,mx_colrad)
!  REAL ::    dhdr(nz_rdr,mx_colrad)
!
!-----------------------------------------------------------------------
!
!  Retrieval radar variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=5) :: stnret(mx_ret)
  INTEGER :: isrcret(0:mx_ret)
  REAL :: latret(mx_ret)
  REAL :: lonret(mx_ret)
  REAL :: elvret(mx_ret)
!
!-----------------------------------------------------------------------
!
!  Retrieval observation variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iret(mx_colret)
  INTEGER :: nlevret(mx_colret)
  REAL :: latretc(mx_colret)
  REAL :: lonretc(mx_colret)
  REAL :: xretc(mx_colret)
  REAL :: yretc(mx_colret)
  REAL :: trnretc(mx_colret)
  REAL :: hgtretc(nz_ret,mx_colret)
  REAL :: theretc(nz_ret,mx_colret)
  REAL :: qrret(nz_ret,mx_colret)

  REAL :: obsret (nvar_anx_MODE,nz_ret,mx_colret)
  REAL :: odifret(nvar_anx_MODE,nz_ret,mx_colret)
  REAL :: oanxret(nvar_anx_MODE,nz_ret,mx_colret)
#ifdef arps3dvar
  REAL :: corret(nz_ret,mx_colret,nvar_anx_3dvar)
#else
  REAL :: corret(nvar_anx_adas,nz_ret,mx_colret)
#endif
!
!-----------------------------------------------------------------------
!
!  Namen de variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=6) :: var_sfc(nvar_sng)
  DATA var_sfc/ 't     ','td    ','dd    ','ff    ',                    &
                'urot  ','vrot  ','pstn  ','pmsl  ',                    &
                'alt   ','ceil  ','low   ','cover ',                    &
                'solar ','vis'/
  CHARACTER (LEN=6) :: var_ua(nvar_ua)
  DATA var_ua / 'press ',                                               &
                'temp  ','dewpt ','dd    ','ff    '/
#ifdef arps3dvar
  CHARACTER (LEN=6) :: var_anx(nvar_anx_3dvar)
  DATA var_anx/ 'u     ','v     ',                                      &
                'press ','theta ','qv    ','w     '/
#else
  CHARACTER (LEN=6) :: var_anx(nvar_anx_adas)
  DATA var_anx/ 'u     ','v     ',                                      &
                'press ','theta ','qv    '/
#endif
!
!-----------------------------------------------------------------------
!
!  Source-dependent parameters
!  Qsrc is the expected square error it is used for
!  setting the QC threshold (qclim) and for relative
!  weighting of observations.
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=80) :: srcback
  REAL :: hqback(nz_tab)
  INTEGER :: nlvqback

  CHARACTER (LEN=80) :: srcsng_full(nsrc_sng)

  CHARACTER (LEN=80) :: srcua_full(nsrc_ua)
  INTEGER :: nlvuatab(nsrc_ua)
  REAL :: huaqsrc(nz_tab,nsrc_ua)

  CHARACTER (LEN=80) :: srcrad_full(nsrc_rad)
  REAL :: qsrcrad(nvar_radin,nsrc_rad)

  CHARACTER (LEN=80) :: srcret_full(nsrc_ret)
  INTEGER :: nlvrttab(nsrc_ret)
  REAL :: hrtqsrc(nz_tab,nsrc_ret)

  REAL :: qback  (nvar_anx_MODE,nz_tab)
  REAL :: qsrcsng(nvar_anx_MODE,nsrc_sng)
  REAL :: qsrcua (nvar_anx_MODE,nz_tab,nsrc_ua)
  REAL :: qsrcret(nvar_anx_MODE,nz_tab,nsrc_ret)
!
!-----------------------------------------------------------------------
!
!  Quality Control Parameters and Variables
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: rmiss = -99.

!  INTEGER :: qualrdsng(mx_sng,nvar_sng,ntime)
  REAL    :: qcmulsng(nsrc_sng)

!  INTEGER :: qualsng(nvar_anx_MODE,mx_sng)
!  REAL    :: qobsng(nvar_anx_MODE,mx_sng)
  REAL    :: qcthrsng(nvar_anx_MODE,nsrc_sng)
  REAL    :: barqclim(nvar_anx_MODE,nsrc_sng)

  REAL :: qcmulua(nsrc_ua)

!  INTEGER :: qualua(nvar_anx_MODE,nz_ua,mx_ua)
!  REAL :: qobsua (nvar_anx_MODE,nz_ua,mx_ua)
  REAL :: qcthrua(nvar_anx_MODE,nsrc_ua)
!
  REAL :: qcmulrad(nsrc_rad)
  REAL :: qcthrrad(nvar_radin,nsrc_rad)
!  INTEGER :: qualrad(nvar_rad_MODE,nz_rdr,mx_colrad)
!  REAL :: qobsrad(nvar_rad_MODE,nz_rdr,mx_colrad)
!
  REAL :: qcmulret(nsrc_ret)
  INTEGER :: qualret(nvar_anx_MODE,nz_ret,mx_colret)
  REAL :: qobsret(nvar_anx_MODE,nz_ret,mx_colret)
  REAL :: qcthrret(nvar_anx_MODE,nsrc_ret)
!
!-----------------------------------------------------------------------
!
!  Climin is the minimum possible value of each observed variable.
!  Climax is the maximum possible value of each observed variable.
!  Used for surface data only.
!
!-----------------------------------------------------------------------
!
  REAL :: climin(nvar_sng),climax(nvar_sng)
  DATA climin /  -50.,    -50.,    0.,    0.,                           &
                -100.,   -100.,  700.,  880.,                           &
                 880.,      0.,    0.,    0.,                           &
                   0.,      0./
  DATA climax /  125.,    125.,  360.,  100.,                           &
                 100.,    100., 1100., 1090.,                           &
                1090.,  30000.,30000.,    1.,                           &
                1500.,    150./
!
!-----------------------------------------------------------------------
!
!  Filenames and such
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=12)  :: suffix
  CHARACTER (LEN=256) :: froot
  CHARACTER (LEN=256) :: qclist
  CHARACTER (LEN=256) :: qcsave
  CHARACTER (LEN=256) :: stats        ! runname + '.lst'
  INTEGER :: istatus,istat_radar
  INTEGER :: nobsng,nobsrd(ntime)
  INTEGER :: i4timanx
  INTEGER :: lenfnm,maxsuf,lensuf,dotloc,mxroot,lenroot
  INTEGER :: iqclist,iqcsave,iwstat

  CHARACTER (LEN=256) :: status_file  ! runname + '.adasstat'
  CHARACTER (LEN=256) :: stnlist_file  ! runname + '.adasstn'

  INTEGER n_used_sng(nsrc_sng), n_used_ua(nsrc_sng)
  INTEGER jsta, klev, ngoodlev

!
!-----------------------------------------------------------------------
!
!  MPI bookkeeping
!
!-----------------------------------------------------------------------
!

  INTEGER :: indexsng(mx_sng)    ! indexes each data point according to
                                 ! which local proc domain it resides
  INTEGER :: oindexsng(mx_sng)   ! obs number mapping index
  LOGICAL :: usesng(mx_sng)       ! whether each obs should be use in 3dvar
  INTEGER, ALLOCATABLE :: ksng(:) ! Number of obs owned by each processor
  INTEGER :: ksngmax              ! Max "ksng" value.

  INTEGER :: indexua(mx_ua)      ! indexes each data point according to
                                 ! which local proc domain it resides
  INTEGER :: oindexua(mx_ua)     ! obs number mapping index
  LOGICAL :: useua(mx_ua)        ! whether each obs should be use in 3dvar
  INTEGER, ALLOCATABLE :: kua(:) ! Number of obs owned by each processor
  INTEGER :: kuamax              ! Max "kua" value.


  INTEGER :: indexrad(mx_colrad) ! indexes each data point according to
                                 ! which local proc domain it resides
  INTEGER :: oindexrad(mx_colrad)! obs number mapping index
  LOGICAL :: userad(mx_colrad)    ! whether each obs should be use in 3dvar
  INTEGER, ALLOCATABLE :: krad(:) ! Number of obs owned by each processor
  INTEGER :: kradmax             ! Max "krad" value.
  INTEGER :: ncolrad_mpi         ! Real total, as "ncolrad" is now local.

  INTEGER :: indexret(mx_colret) ! indexes each data point according to
                                 ! which local proc domain it resides
                                 ! WARNING:  used not yet coded!!!

  INTEGER :: iproc,iprocv        ! Radii of influence "i" direction
  INTEGER :: jproc,jprocv        ! Radii of influence "j" direction

  INTEGER, ALLOCATABLE :: mpi_map(:,:) ! Map of communication entries
  INTEGER :: nmap                ! Number of entries in "mpi_map".

  INTEGER, ALLOCATABLE :: mpi_mapv(:,:) ! Same for radial velocity
  INTEGER :: nmapv               ! Number of entries in "mpi_mapv".

  INTEGER :: iradvel             ! Reduce memory requirements for velocity.
                                 ! Never change the set value.  This is a
                                 ! code development tool only!

!
!-----------------------------------------------------------------------
!
!  Physical constants
!
!-----------------------------------------------------------------------
!
  REAL, PARAMETER :: kts2ms=0.514444, mb2pa=100., f2crat=(5./9.)
!
!-----------------------------------------------------------------------
!
!  ptlapse is minimum potential temperature lapse rate to use in
!  superadiabatic check.  Used to ensure initial static stability.
!
!-----------------------------------------------------------------------

  REAL, PARAMETER :: ptlapse = (1.0/4000.)    ! deg K / m
!
!-----------------------------------------------------------------------
!
!  Function f_qvsat and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_qvsat

!fpp$ expand (f_qvsat)
!!dir$ inline always f_qvsat
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=8  ) :: rdsource
  CHARACTER (LEN=8  ) :: srcstr
  CHARACTER (LEN=256) :: basdmpfn
  CHARACTER (LEN=80 ) :: header
  INTEGER :: i,j,k,nq,ios,ktab,isrc,jsrc,ivar,ifile
  INTEGER :: grdbas,lbasdmpf,rbufsz
  INTEGER :: nobsua,ncolrad,ncolret,itype,kntgood
  REAL :: latua,lonua

  REAL :: tgrid,qvmin,qvsat,qsrcmax,rngsq
  REAL :: temp,tvbar,qvprt,qtot,pconst,pr1,pr2

  INTEGER :: ixrad(mx_rad),jyrad(mx_rad)

  INTEGER, PARAMETER  :: retqropt=1

  INTEGER :: imstat

  INTEGER :: cntl_var_rh
            ! = 1: relative humindity is used in analysis
            ! = 0: specific humindity is used in analysis

  INTEGER :: nvar_anx, nvar_rad
  INTEGER :: ipass, ierr
  REAL    :: f_cputime
  REAL    :: cpuin, cpuout

  CHARACTER (LEN=8) :: tmp
  INTEGER :: lentmp

#ifdef arps3dvar
  CHARACTER(LEN=9) :: PROGNAME='ARPS3DVAR'
#else
  CHARACTER(LEN=4) :: PROGNAME='ADAS'
#endif

  CHARACTER(LEN=256) :: namelist_filename
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------

  cputimearray=0.0
  nt = 1    ! Number of time levels of data
  ncat = 4  ! Number of correlation categories
  iradvel = 1 ! NEVER CHANGE!!!  Leave as 1.   kwt

!-----------------------------------------------------------------------
!  Set grid mode to non-nesting and grid index, mgrid, to 1.
!-----------------------------------------------------------------------

  nestgrd=0
  mgrid=1

!-----------------------------------------------------------------------
!  read in ARPS namelists
!-----------------------------------------------------------------------

  namelist_filename = ' '
  CALL initpara(nx,ny,nz,nzsoil,nstyps, namelist_filename)

  IF (mp_opt > 0) THEN
    nxlg = (nx - 3) * nproc_x + 3
    nylg = (ny - 3) * nproc_y + 3
  ELSE
    nxlg = nx
    nylg = ny
  END IF

  IF (lbcopt == 2) THEN
    IF (myproc == 0) THEN
      WRITE (*,*) "INITADAS: resetting lbcopt to 1 & lateral bc's to 4"
    END IF
    lbcopt = 1
    ebc = 4
    wbc = 4
    nbc = 4
    sbc = 4
  END IF

  nvar_anx = nvar_anx_MODE
  nvar_rad = nvar_rad_MODE

!-----------------------------------------------------------------------
!  Read in adas namelist parameters
!-----------------------------------------------------------------------
!
  CALL initadas(namelist_filename)

#ifdef arps3dvar
  CALL init3dvar(cntl_var_rh, npass, namelist_filename)
#else
  cntl_var_rh = 0
#endif

!
! ARPS3DVAR doesn't use this interface.
!

#ifdef adas
!
  IF (mp_opt > 0) THEN
    CALL radinf(nvar_anx,rlim)
!
!   See what our processor radius of influence is in each direction.
!   Note we subtract four, as we need a distance.
!
    iproc = INT( rlim / ( (nx-4) * dx )) + 1
    jproc = INT( rlim / ( (ny-4) * dy )) + 1

    IF (myproc == 0) WRITE(6,'(1x,2(a,I2))')  &
       'Process radius of influence x & y ',iproc,' x ',jproc
!
!   Handle radial velocity separately.
!
    IF (raduvobs > 0) THEN
      CALL radinfv(nvar_anx,rlimv)
!
!   See what our processor radius of influence is in each direction.
!   Note we subtract four, as we need a distance.
!
      iprocv = INT( rlimv / ( (nx-4) * dx )) + 1
      jprocv = INT( rlimv / ( (ny-4) * dy )) + 1

      IF (myproc == 0) WRITE(6,'(1x,2(a,I2))')  &
         'Process velocity radius of influence x & y ',iprocv,' x ',jprocv
    ELSE
      iprocv = 1
      jprocv = 1
    END IF

  ELSE
#endif
!
!   We'll allocate a trivial array.
!
    iproc = 1
    jproc = 1

    iprocv = 1
    jprocv = 1
#ifdef adas
  ENDIF
#endif

  IF(incrdmp > 0 ) THEN
    nxndg=nx
    nyndg=ny
    nzndg=nz
  ELSE
    nxndg=1
    nyndg=1
    nzndg=1
  END IF
!
!-----------------------------------------------------------------------
!
! Allocate analysis arrays
!
!-----------------------------------------------------------------------

  CALL allocate_adas_3dvar_arrays(progname,mx_sng,nvar_sng,nvar_anx,ntime,     &
                                  mx_ua,nz_ua,mx_rad,                          &
                                  mx_colrad,nvar_rad,nvar_radin,nz_rdr,        &
                                  raduvobs,radrhobs,istatus)
  IF (istatus /= 0) THEN
    WRITE(*,'(1x,/,A,I2,A)') 'ERROR: allocate_adas_3dvar_arrays, status= ',    &
                             istatus, '.'
    CALL arpsstop('FATAL: Unable to allocate array. Program ends.',1)
  END IF

  CALL allocatearpsArray(nx,ny,nz,nzsoil)
  CALL allocatesoilArray(nx,ny,nz,nstyps,nzsoil)
  CALL allocateanaIncArray(nxndg,nyndg,nzndg)

  ALLOCATE(knt(nvar_anx,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:knt")
  ALLOCATE(rngsqi(nvar_anx),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:rngsqi")
  ALLOCATE(wgtsum(nvar_anx,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:wgtsum")
  ALLOCATE(zsum(nvar_anx,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:zsum")
  ALLOCATE(sqsum(nvar_anx,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:sqsum")

  ALLOCATE(xs(nx),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:xs")
  ALLOCATE(ys(ny),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:ys")

  ALLOCATE(xslg(nxlg),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:xslg")
  ALLOCATE(yslg(nylg),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:yslg")

  ALLOCATE(zs(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:zs")
  ALLOCATE(latgr(nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:latgr")
  ALLOCATE(longr(nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:longr")

  ALLOCATE(tem1soil(nx,ny,nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:tem1soil")
  ALLOCATE(tem2soil(nx,ny,nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:tem2soil")
  ALLOCATE(tem3soil(nx,ny,nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:tem3soil")
  ALLOCATE(tem4soil(nx,ny,nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:tem4soil")
  ALLOCATE(tem5soil(nx,ny,nzsoil),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:tem5soil")

  ALLOCATE(tem1(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:tem1")
  ALLOCATE(tem2(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:tem2")
  ALLOCATE(tem3(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:tem3")
  ALLOCATE(tem4(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:tem4")
  ALLOCATE(tem5(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:tem5")
  ALLOCATE(tem6(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:tem6")
  ALLOCATE(tem7(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:tem7")
  ALLOCATE(tem8(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:tem8")
  ALLOCATE(tem9(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:tem9")
  ALLOCATE(tem10(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:tem10")
  ALLOCATE(tem11(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:tem11")
  ALLOCATE(ibegin(ny),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:ibegin")
  ALLOCATE(iend(ny),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:iend")

  ALLOCATE(anx(nx,ny,nz,nvar_anx),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:anx")

  ALLOCATE(xcor(ncat,ncat),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:xcor")
  ALLOCATE(icatg(nx,ny),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:icatg")

  nmap = ( 2 * iproc + 1 ) * ( 2 * jproc + 1 )
  ALLOCATE(mpi_map(nmap,2),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:mpi:map")

  nmapv = ( 2 * iprocv + 1 ) * ( 2 * jprocv + 1 )
  ALLOCATE(mpi_mapv(nmap,2),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:mpi:mapv")

  IF (myproc == 0) WRITE (*,*) "Finished allocating arrays"

!-----------------------------------------------------------------------
! Initialize the allocated arrays to zero.
!-----------------------------------------------------------------------

  x=0.0
  y=0.0
  z=0.0

  hterain=0.0
  mapfct =0.0

  zp  =0.0
  zpsoil  =0.0
  j1  =0.0
  j2  =0.0
  j3  =0.0
  j3soil  =0.0
  aj3x=0.0
  aj3y=0.0
  aj3z=0.0
  j3inv=0.0
  j3soilinv=0.0

  u    =0.0
  v    =0.0
  w    =0.0
  wcont=0.0
  ptprt=0.0
  pprt =0.0

  qv   =0.0
  qscalar = 0.0
  tke  =0.0

  ubar =0.0
  vbar =0.0
  ptbar=0.0
  pbar =0.0
  rhostr=0.0
  qvbar=0.0

  udteb=0.0
  udtwb=0.0
  vdtnb=0.0
  vdtsb=0.0

  pdteb=0.0
  pdtwb=0.0
  pdtnb=0.0
  pdtsb=0.0

  trigs1=0.0
  trigs2=0.0
  ifax1=0.0
  ifax2=0.0

  vwork1=0.0
  vwork2=0.0

  wsave1=0.0
  wsave2=0.0

  ppi   =0.0
  csndsq=0.0

  radfrc=0.0
  radsw =0.0
  rnflx =0.0
  radswnet =0.0
  radlwin =0.0

  soiltyp=0.0
  stypfrct=0.0
  vegtyp =0.0
  lai    =0.0
  roufns =0.0
  veg    =0.0

  qvsfc  =0.0
  tsoil  =0.0
  qsoil  =0.0
  wetcanp=0.0
  snowdpth=0.0

  ptcumsrc=0.0
  qcumsrc =0.0
  w0avg   =0.0
  nca    =0.0
  kfraincv =0.0
  cldefi = 0.0
  xland = 0.0
  bmjraincv = 0.0
  raing  =0.0
  rainc  =0.0
  prcrate=0.0

  usflx =0.0
  vsflx =0.0
  ptsflx=0.0
  qvsflx=0.0

  knt=0.0
  rngsqi=0.0
  wgtsum=0.0
  zsum=0.0
  sqsum=0.0

  xs=0.0
  ys=0.0
  xslg = 0.0
  yslg = 0.0
  zs=0.0
  latgr=0.0
  longr=0.0

  tem1=0.0
  tem2=0.0
  tem3=0.0
  tem4=0.0
  tem5=0.0
  tem6=0.0
  tem7=0.0
  tem8=0.0
  tem9=0.0
  tem10=0.0
  tem11=0.0

  qsrcsng = 0.0

  ibegin=0
  iend=0

  anx=0.0

  xcor=0.0
  icatg=0.0

!
!-----------------------------------------------------------------------
!
!  Read in quality control info.  Most of the information will not be
!  broadcast, as processor 0 will handle pre-processing of obs before
!  the data gets sent to the other processors.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Set expected _squared_ background error for each variable
!  This will depend on the particular background field used,
!  season, age of forecast, etc.
!
!-----------------------------------------------------------------------
!

  IF (myproc == 0) THEN
     WRITE(6,'(/,a,a)') ' Reading ', TRIM(backerrfil)
     OPEN(4,FILE=trim(backerrfil),STATUS='old')
     READ(4,'(a80)') srcback
     READ(4,'(a80)') header
     WRITE(6,'(//,a,/a)') 'Background Std Error for',srcback
     WRITE(6,'(1x,a)')                                                  &
         '   k    hgt(m)  u(m/s)  v(m/s) pres(mb) temp(K)  RH(%)'
     DO ktab=1,nz_tab
       READ(4,*,iostat=istatus) hqback(ktab),                           &
                qback(PTR_P, ktab),                                     &
                qback(PTR_PT,ktab),                                     &
                qback(PTR_QV,ktab),                                     &
                qback(PTR_U, ktab),                                     &
                qback(PTR_V, ktab)
       IF( istatus /= 0 ) EXIT
       WRITE(6,'(1x,i4,f10.2,5f8.2)') ktab,hqback(ktab),                &
                               (qback(ivar,ktab),ivar=1,5)
       qback(PTR_P,ktab)=100.*qback(PTR_P,ktab)
#ifdef arps3dvar
       IF( cntl_var_rh == 1 ) THEN
         qback(PTR_QV,ktab)=0.01*qback(PTR_QV,ktab)
       ELSE
         qback(PTR_QV,ktab)=0.01
       ENDIF
       qback(PTR_W,ktab)=qback(PTR_U,ktab)      ! add error for w. Gao
#else
       qback(PTR_QV,ktab)=0.01*qback(PTR_QV,ktab)
#endif
     END DO
     nlvqback=ktab-1
     CLOSE(4)

     DO jsrc=1,nsrc_sng
        IF(srcsng(jsrc) /= 'NULL') THEN
           WRITE(6,'(/,a,a)') ' Reading ', TRIM(sngerrfil(jsrc))
           OPEN(4,FILE=trim(sngerrfil(jsrc)),STATUS='old')
           READ(4,'(a8)') rdsource
           IF(rdsource /= srcsng(jsrc)) THEN
              WRITE(6,'(a,i4,a,a,a,a,a)')                               &
              ' Mismatch of source names',jsrc,                         &
              ' read ',rdsource,' expected ',srcsng(jsrc)
              CALL arpsstop("mismatch",1)
           END IF
           READ(4,'(a80)') srcsng_full(jsrc)
           READ(4,*) qcmulsng(jsrc)
           READ(4,'(a80)') header
           READ(4,*) qsrcsng(PTR_P, jsrc), qsrcsng(PTR_PT,jsrc),        &
                     qsrcsng(PTR_QV,jsrc), qsrcsng(PTR_U, jsrc),        &
                     qsrcsng(PTR_V, jsrc)
           WRITE(6,'(//,a,a,/a)') 'Single-level std error for ',        &
                            srcsng(jsrc),srcsng_full(jsrc)
           WRITE(6,'(a,f5.1)') ' QC multiplier: ',qcmulsng(jsrc)
           WRITE(6,'(1x,a)')                                            &
               '   u (m/s)  v (m/s)  pres(mb) temp(K) RH(%)'
           WRITE(6,'(1x,5f8.2)') (qsrcsng(ivar,jsrc),ivar=1,5)
           qsrcsng(PTR_P,jsrc)=100.*qsrcsng(PTR_P,jsrc)
#ifdef arps3dvar
           IF( cntl_var_rh == 1 ) THEN
             qsrcsng(PTR_QV,jsrc)=0.01*qsrcsng(PTR_QV,jsrc)
           ELSE
             qsrcsng(PTR_QV,jsrc)=0.004
           ENDIF
           qsrcsng(PTR_W,jsrc)=qsrcsng(PTR_U,jsrc)
#else
           qsrcsng(PTR_QV,jsrc)=0.01*qsrcsng(PTR_QV,jsrc)
#endif
           CLOSE(4)
        ELSE
          DO k=1,npass
            iusesng(jsrc,k) = 0
          END DO
        END IF
     END DO

     DO jsrc=1,nsrc_ua
        IF(srcua(jsrc) /= 'NULL') THEN
           WRITE(6,'(/,a,a)') ' Reading ', TRIM(uaerrfil(jsrc))
           OPEN(4,FILE=trim(uaerrfil(jsrc)),STATUS='old')
           READ(4,'(a8)') rdsource
           IF(rdsource /= srcua(jsrc)) THEN
              WRITE(6,'(a,i4,a,a,a,a,a)')                               &
                 ' Mismatch of source names',jsrc,                      &
                 ' read ',rdsource,' expected ',srcua(jsrc)
              CALL arpsstop("mismatch",1)
           END IF
           READ(4,'(a80)') srcua_full(jsrc)
           READ(4,*) qcmulua(jsrc)
           READ(4,'(a80)') header
           WRITE(6,'(//,a,a,/a)') 'UA data std error for ',             &
                            srcua(jsrc),srcua_full(jsrc)
           WRITE(6,'(a,f5.1)') ' QC multiplier: ',qcmulua(jsrc)
           WRITE(6,'(1x,a)')                                            &
              '   k    hgt(m)  u(m/s)  v(m/s) pres(mb) temp(K)  RH(%)'
           DO ktab=1,nz_tab
              READ(4,*,iostat=istatus) huaqsrc(ktab,jsrc),              &
                              qsrcua(PTR_P, ktab,jsrc),                 &
                              qsrcua(PTR_PT,ktab,jsrc),                 &
                              qsrcua(PTR_QV,ktab,jsrc),                 &
                              qsrcua(PTR_U, ktab,jsrc),                 &
                              qsrcua(PTR_V, ktab,jsrc)
              IF( istatus /= 0 ) EXIT
              WRITE(6,'(1x,i4,f10.2,5f8.2)') ktab,huaqsrc(ktab,jsrc),   &
                 (qsrcua(ivar,ktab,jsrc),ivar=1,5)
              qsrcua(PTR_P,ktab,jsrc)=100.*qsrcua(PTR_P,ktab,jsrc)
#ifdef arps3dvar
              IF( cntl_var_rh == 1 ) THEN
                qsrcua(PTR_QV,ktab,jsrc)=0.01*qsrcua(PTR_QV,ktab,jsrc)
              ELSE
                qsrcua(PTR_QV,ktab,jsrc)=0.004
              ENDIF
              qsrcua(PTR_W,ktab,jsrc)=qsrcua(PTR_U,ktab,jsrc)    !for w GAO
#else
              qsrcua(PTR_QV,ktab,jsrc)=0.01*qsrcua(PTR_QV,ktab,jsrc)
#endif
           END DO
           nlvuatab(jsrc)=ktab-1
           CLOSE(4)
        ELSE
          DO k=1,npass
            iuseua(jsrc,k) = 0
          END DO
        END IF
     END DO

     DO jsrc=1,nsrc_rad
        IF(srcrad(jsrc) /= 'NULL') THEN
          WRITE(6,'(/,a,a)') ' Reading ', TRIM(raderrfil(jsrc))
          OPEN(4,FILE=trim(raderrfil(jsrc)),STATUS='old')
          READ(4,'(a8)') rdsource
          IF(rdsource /= srcrad(jsrc)) THEN
            WRITE(6,'(a,i4,a,a,a,a,a)')                                 &
               ' Mismatch of source names',jsrc,                        &
               ' read ',rdsource,' expected ',srcrad(jsrc)
            CALL arpsstop("mismatch",1)
          END IF
          READ(4,'(a80)') srcrad_full(jsrc)
          READ(4,*) qcmulrad(jsrc)
          READ(4,'(a80)') header
          READ(4,*) qsrcrad(1,jsrc),                                    &
                    qsrcrad(2,jsrc),                                    &
                    qsrcrad(3,jsrc)
          WRITE(6,'(//,a,a,/a)') 'Radar data std error for ',           &
                          srcrad(jsrc),srcrad_full(jsrc)
          WRITE(6,'(a,f5.1)') ' QC multiplier: ',qcmulrad(jsrc)
          WRITE(6,'(1x,a)')                                             &
             '    ref(dBz) Vrad(m/s) SpWid(m/s)'
          WRITE(6,'(1x,4f8.2)')                                         &
             (qsrcrad(ivar,jsrc),ivar=1,nvar_radin)
          CLOSE(4)
        ELSE
          DO k=1,npass
            iuserad(jsrc,k) = 0
          END DO
        END IF
     END DO

     DO jsrc=1,nsrc_ret
        IF(srcret(jsrc) /= 'NULL') THEN
           WRITE(6,'(/,a,a)') ' Reading ', TRIM(reterrfil(jsrc))
           OPEN(4,FILE=trim(reterrfil(jsrc)),STATUS='old')
           READ(4,'(a8)') rdsource
           IF(rdsource /= srcret(jsrc)) THEN
             WRITE(6,'(a,i4,a,a,a,a,a)')                                &
                ' Mismatch of source names',jsrc,                       &
                ' read ',rdsource,' expected ',srcret(jsrc)
             CALL arpsstop("mismatch",1)
           END IF
           READ(4,'(a80)') srcret_full(jsrc)
           READ(4,*) qcmulret(jsrc)
           READ(4,'(a80)') header
           WRITE(6,'(//,a,a,/a)') 'Retrieval std error for ',           &
                               srcret(jsrc),srcret_full(jsrc)
           WRITE(6,'(a,f5.1)') ' QC multiplier: ',qcmulrad(jsrc)
           WRITE(6,'(1x,a)')                                            &
               '   k    hgt(m)  u(m/s)  v(m/s) pres(mb) temp(K)  RH(%)'
           DO ktab=1,nz_tab
              READ(4,*,iostat=istatus) hrtqsrc(ktab,jsrc),              &
                           qsrcret(3,ktab,jsrc), qsrcret(4,ktab,jsrc),  &
                           qsrcret(5,ktab,jsrc), qsrcret(1,ktab,jsrc),  &
                           qsrcret(2,ktab,jsrc)
              IF( istatus /= 0 ) EXIT
              WRITE(6,'(1x,i4,f10.2,5f8.2)') ktab,hrtqsrc(ktab,jsrc),   &
                              (qsrcret(ivar,ktab,jsrc),ivar=1,5)
              qsrcret(3,ktab,jsrc)=100.*qsrcret(3,ktab,jsrc)
#ifdef arps3dvar
              IF( cntl_var_rh == 1 ) THEN
                qsrcret(5,ktab,jsrc)=0.01*qsrcret(5,ktab,jsrc)
              ELSE
                qsrcret(5,ktab,jsrc)=0.004
              ENDIF
#else
              qsrcret(5,ktab,jsrc)=0.01*qsrcret(5,ktab,jsrc)
#endif
           END DO
           nlvrttab(jsrc)=ktab-1
           CLOSE(4)
        ELSE
          DO k=1,npass
            iuseret(jsrc,k) = 0
          END DO
        END IF
     END DO
!
!-----------------------------------------------------------------------
!
!  Change standard error to standard error variance by squaring
!  Calculate quality control thresholds
!
!-----------------------------------------------------------------------
!
    DO ktab=1,nlvqback
      DO ivar=1,nvar_anx
        qback(ivar,ktab)=qback(ivar,ktab)*qback(ivar,ktab)
      END DO
    END DO
!
    DO isrc=1,nsrc_sng
      DO ivar=1,nvar_anx
        qsrcsng(ivar,isrc)=qsrcsng(ivar,isrc)*qsrcsng(ivar,isrc)
        qcthrsng(ivar,isrc)=qcmulsng(isrc)*                             &
                         SQRT(qback(ivar,1)+qsrcsng(ivar,isrc))
        barqclim(ivar,isrc)=qcmulsng(isrc)*SQRT(qsrcsng(ivar,isrc))
      END DO
      qcthrsng(iqv,isrc)=min(qcthrsng(iqv,isrc),0.75)
      barqclim(iqv,isrc)=min(barqclim(iqv,isrc),0.75)
    END DO
!
    DO isrc=1,nsrc_ua
      DO ivar=1,nvar_anx
        qsrcmax=0.
        DO ktab=1,nlvuatab(isrc)
          qsrcua(ivar,ktab,isrc) =                                      &
              qsrcua(ivar,ktab,isrc)*qsrcua(ivar,ktab,isrc)
          qsrcmax=AMAX1(qsrcmax,qsrcua(ivar,ktab,isrc))
        END DO
        qcthrua(ivar,isrc)=qcmulua(isrc)*                               &
                           SQRT(qback(ivar,1)+qsrcmax)
      END DO
      qcthrua(iqv,isrc)=min(qcthrua(iqv,isrc),0.75)
    END DO
!
    DO isrc=1,nsrc_rad
      DO ivar=1,nvar_radin
        qsrcrad(ivar,isrc) =                                            &
            qsrcrad(ivar,isrc)*qsrcrad(ivar,isrc)
        qcthrrad(ivar,isrc)=qcmulrad(isrc)*                             &
                            SQRT(qback(ivar,1)+qsrcrad(ivar,isrc))
      END DO
    END DO
!
    DO isrc=1,nsrc_ret
      DO ivar=1,nvar_anx
        qsrcmax=0.
        DO ktab=1,nlvrttab(isrc)
          qsrcret(ivar,ktab,isrc) =                                     &
              qsrcret(ivar,ktab,isrc)*qsrcret(ivar,ktab,isrc)
          qsrcmax=AMAX1(qsrcmax,qsrcret(ivar,ktab,isrc))
        END DO
        qcthrret(ivar,isrc)=qcmulret(isrc)*                             &
                            SQRT(qback(ivar,1)+qsrcmax)
      END DO
      qcthrret(iqv,isrc)=min(qcthrret(iqv,isrc),0.75)
    END DO
  END IF

!
!  Broadcast read variables to all procs.  Some are dead, so we keep them
!  that way.  Others are only used by processor zero initially.
!
!  Also broadcast variables that may have been updated.
!

  CALL mpupdater(hqback, nz_tab)
  CALL mpupdater(qback, (nvar_anx*nz_tab))
  CALL mpupdatei(nlvqback, 1)

  CALL mpupdatec(srcsng,   (8*nsrc_sng))
  CALL mpupdater(qcthrsng, (nvar_anx*nsrc_sng))

  CALL mpupdater(qcthrua, (nvar_anx*nsrc_ua))

  CALL mpupdater(qcthrrad, (nvar_radin*nsrc_rad))

  CALL mpupdater(qcthrret, (nvar_anx*nsrc_ret))

  CALL mpupdater(qsrcrad,nvar_radin*nsrc_rad)

  CALL mpupdatei(iusesng,(nsrc_sng+1)*mx_pass)
  CALL mpupdatei(iuseua,(nsrc_ua+1)*mx_pass)
  CALL mpupdatei(iuserad,(nsrc_rad+1)*mx_pass)
  CALL mpupdatei(iuseret,(nsrc_ret+1)*mx_pass)

  CALL make_mpi_map(mpi_map,nmap,iproc,jproc,nx,ny)
  CALL make_mpi_map(mpi_mapv,nmapv,iprocv,jprocv,nx,ny)

!
!-----------------------------------------------------------------------
!
!  Set cross-correlations between numerical categories.
!  Roughly 1=clear,2=some evidence of outflow,3=precip,4=conv precip
!  This could be read in as a table.
!
!-----------------------------------------------------------------------
!
  DO j=1,ncat
    DO i=1,ncat
      xcor(i,j)=1.0-(IABS(i-j))*0.25
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Initialize grids, forming first guess fields based on
!  the model initial option specified in the ARPS init file.
!
!-----------------------------------------------------------------------
!

  CALL initgrdvar(nx,ny,nz,nzsoil,1,nstyps,1,                           &
                  x,y,z,zp,zpsoil,hterain,mapfct,                       &
                  j1,j2,j3,j3soil,aj3x,aj3y,aj3z,j3inv,j3soilinv,       &
                  u,v,w,wcont,ptprt,pprt,qv,qscalar,tke,                &
                  udteb, udtwb, vdtnb, vdtsb,                           &
                  pdteb,pdtwb ,pdtnb ,pdtsb,                            &
                  trigs1,trigs2,ifax1,ifax2,                            &
                  wsave1,wsave2,vwork1,vwork2,                          &
                  ubar,vbar,tem1,ptbar,pbar,tem10,tem10,                &
                  rhostr,tem10,qvbar,ppi,csndsq,                        &
                  soiltyp,stypfrct,vegtyp,lai,roufns,veg,               &
                  tsoil,qsoil,wetcanp,snowdpth,qvsfc,                   &
                  ptcumsrc,qcumsrc,w0avg,nca,kfraincv,                  &
                  cldefi,xland,bmjraincv,                               &
                  raing,rainc,prcrate,exbcbuf,                          &
                  radfrc,radsw,rnflx,radswnet,radlwin,                  &
                  usflx,vsflx,ptsflx,qvsflx,                            &
                  tem1soil,tem2soil,tem3soil,tem4soil,tem5soil,         &
                  tem1,tem2,tem3,tem4,tem5,                             &
                  tem6,tem7,tem8,tem9)
  curtim = tstart   ! Added by G. Ge on 3/3/200

  DEALLOCATE( tem1soil, tem2soil, tem3soil, tem4soil, tem5soil )
!
!-----------------------------------------------------------------------
!
!  Deallocate some arrays needed only for initgrdvar
!
!-----------------------------------------------------------------------
!
  DEALLOCATE(ptcumsrc)
  DEALLOCATE(qcumsrc)
  DEALLOCATE(w0avg)
  DEALLOCATE(kfraincv)
  DEALLOCATE(nca)
  DEALLOCATE(cldefi)
  DEALLOCATE(xland)
!
!-----------------------------------------------------------------------
!
!  Set location of scalar fields.
!
!-----------------------------------------------------------------------
!
  DO i=1,nx-1
    xs(i)=0.5*(x(i)+x(i+1))
  END DO
  xs(nx)=2.0*xs(nx-1)-xs(nx-2)
  DO j=1,ny-1
    ys(j)=0.5*(y(j)+y(j+1))
  END DO
  ys(ny)=2.0*ys(ny-1)-ys(ny-2)
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        zs(i,j,k)=0.5*(zp(i,j,k)+zp(i,j,k+1))
      END DO
    END DO
  END DO
  DO i=1,nx-1
    DO j=1,ny-1
      zs(i,j,nz)=2.*zs(i,j,nz-1)-zs(i,j,nz-2)
    END DO
  END DO

  CALL xytoll(nx,ny,xs,ys,latgr,longr)

!
!-----------------------------------------------------------------------
!
! In the MPI world, some subroutines will need to know characteristics
! of the large grid.  To simply later subroutine calls (one call not two
! surrounded by IF/ELSE/ENDIF, everybody that always needs the large
! grid will use the large grid variables.  Large grid info is used by
! early decision making when determining what obs are needed.
!
!-----------------------------------------------------------------------
!

  IF (mp_opt > 0 ) THEN
     CALL mpimerge1dx(xs,nx,xslg)
     CALL mpimerge1dy(ys,ny,yslg)

     CALL mpupdater(xslg,nxlg)
     CALL mpupdater(yslg,nylg)
  ELSE
     xslg = xs
     yslg = ys
  END IF

  CALL ctim2abss( year,month,day,hour,minute,second, i4timanx)
!
!-----------------------------------------------------------------------
!
!  Identify the background field correlation category for
!  each 2-d point.   This is based on precip rate, cumulus
!  parameterization switch, and surface relative humidity.
!
!-----------------------------------------------------------------------
!
  CALL setcat(nx,ny,nz,nscalar,ccatopt,zs,                              &
              ptprt,pprt,qv,qscalar,                                    &
              ptbar,pbar,                                               &
              prcrate,icatg)
!
!-----------------------------------------------------------------------
!
!  Load "background fields" for analysis
!  Note, analysis is done at scalar points except for W with arps3dvar
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        anx(i,j,k,PTR_U )=0.5*(u(i,j,k)+u(i+1,j,k))
        anx(i,j,k,PTR_V )=0.5*(v(i,j,k)+v(i,j+1,k))
        anx(i,j,k,PTR_P )=pbar(i,j,k)+pprt(i,j,k)
        anx(i,j,k,PTR_PT)=ptbar(i,j,k)+ptprt(i,j,k)
        anx(i,j,k,PTR_QV)=qv(i,j,k)
      END DO
    END DO
  END DO

#ifdef arps3dvar

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
         anx(i,j,k,PTR_W)= 0.5*(w(i,j,k)+w(i,j,k+1))
      END DO
    END DO
  END DO

#endif
!
!-----------------------------------------------------------------------
!
!  convert specific humidity to relative humidity ( modified by Ming Hu)
!
!-----------------------------------------------------------------------
!
  IF( cntl_var_rh == 1 ) THEN
    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          tem10(i,j,k)= anx(i,j,k,PTR_PT)*(anx(i,j,k,PTR_P)/p0)**rddcp
        END DO
      END DO
    END DO
    CALL getqvs(nx,ny,nz, 1,nx-1,1,ny-1,1,nz-1, anx(1,1,1,PTR_P),tem10,qv)
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          anx(i,j,k,PTR_QV)= anx(i,j,k,PTR_QV)/qv(i,j,k)
          IF( anx(i,j,k,PTR_QV) > 1.0 ) anx(i,j,k,PTR_QV)=1.0
          IF( anx(i,j,k,PTR_QV) < 0.0 ) anx(i,j,k,PTR_QV)=0.0
        END DO
      END DO
    END DO
  ENDIF


  IF (mp_opt > 0) THEN
    DO k=1,nvar_anx
      CALL mpsendrecv2dew(anx(1,1,1,k),nx,ny,nz,ebc,wbc,0,tem1)
      CALL mpsendrecv2dns(anx(1,1,1,k),nx,ny,nz,nbc,sbc,0,tem1)
    END DO
  END IF

!-----------------------------------------------------------------------
!
!  Save background fields on staggered grid in the increment
!  arrays.
!
!  PROBABLY NEEDS MPI WORK!
!
!-----------------------------------------------------------------------
!
  IF(incrdmp > 0) THEN
    CALL incrsave(nx,ny,nz,nxndg,nyndg,nzndg,                           &
                  u,v,w,pprt,ptprt,qv,                                  &
                  qscalar,                                              &
                  uincr,vincr,wincr,pincr,ptincr,qvincr,                &
                  qcincr,qrincr,qiincr,qsincr,qhincr,                   &
                  istatus)
  END IF
!
!-----------------------------------------------------------------------
!
!  Filename creation
!
!-----------------------------------------------------------------------
!
  iwstat = 0
  IF (myproc == 0) THEN
    WRITE(stats,'(a,a)') runname(1:lfnkey),'.lst'
    CALL getunit(iwstat)
    PRINT *, 'Writing ', TRIM(stats)
    OPEN(iwstat,IOSTAT=ios,FILE=trim(stats),STATUS='unknown')
    IF(ios /= 0) iwstat=0
  END IF

  CALL mpupdatei(iwstat,1)

  maxsuf=LEN(suffix)
  mxroot=LEN(froot)
  nobsng=0

  iqclist = 0
  iqcsave = 0

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Read and preprocess the single level data, all done on processor 0.
!  Broadcast the data to the rest of the ADAS world.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  IF (myproc == 0) THEN
    DO ifile=1,nsngfil
      WRITE(6,'(a,i4,2a)') 'Processing file ', ifile, ' ',             &
           TRIM(sngfname(ifile))
      lenfnm=LEN_TRIM(sngfname(ifile))
      CALL exsufx(sngfname(ifile),lenfnm,suffix,maxsuf,dotloc,lensuf)

      IF(lensuf == 3. .AND. suffix(1:3) == 'lso') THEN

        CALL exfroot(sngfname(ifile),lenfnm,froot,mxroot,lenroot)
        WRITE(qclist,'(a,a)') froot(1:lenroot),'.sqc'
        WRITE(qcsave,'(a,a)') froot(1:lenroot),'.lsq'
!
!-----------------------------------------------------------------------
!
!  Open the files for listing QC info
!  To suppress listing set the unit numbers to zero
!
!-----------------------------------------------------------------------
!

        CALL getunit(iqclist)
        PRINT *, 'Opening qclist: ', TRIM(qclist)
        OPEN(iqclist,IOSTAT=ios,FILE=trim(qclist),STATUS='unknown')
        IF(ios /= 0) iqclist=0
        CALL getunit(iqcsave)
        PRINT *, 'Opening qclist: ', TRIM(qcsave)
        OPEN(iqcsave,IOSTAT=ios,FILE=trim(qcsave),STATUS='unknown')
        IF(ios /= 0) iqcsave=0

!
!-----------------------------------------------------------------------
!
!  Read surface data, QC and convert units.
!
!-----------------------------------------------------------------------
!

        rngsq=sfcqcrng*sfcqcrng
        IF (myproc == 0) PRINT *, 'Calling prepsfcobs'
        CALL prepsfcobs(ntime,mx_sng,                                     &
            nvar_sng,nvar_anx,nsrc_sng,nobsng,ipres,iptmp,iqv,            &
            sngfname(ifile),sngtmchk(ifile),blackfil,                     &
            var_sfc,var_anx,srcsng,qsrcsng,                               &
            rmiss,iqspr,sprdist,nobsrd,timesng,                           &
            stnsng,csrcsng,isrcsng,latsng,lonsng,hgtsng,xsng,ysng,        &
            nxlg,nylg,xslg,yslg,                                          &
            wx,kloud,idp3,store_emv,store_hgt,store_amt,                  &
            obrdsng,obsng,qualrdsng,qobsng,qualsng,                       &
            ival,climin,climax,                                           &
            rngsq,klim,wlim,qcthrsng,barqclim,                            &
            knt,wgtsum,zsum,sqsum,iqclist,iqcsave,istatus)

        IF(iqclist /= 0) THEN
          CLOSE(iqclist)
          CALL retunit(iqclist)
        END IF
        IF(iqcsave /= 0) THEN
          CLOSE(iqcsave)
          CALL retunit(iqcsave)
        END IF

      ELSE

!
!-----------------------------------------------------------------------
!
!  Read other single-level data and convert units.  Note that we are
!  using the global domain here, even if we're in MPI mode, as only
!  processor 0 is doing this.
!
!-----------------------------------------------------------------------
!
        IF (myproc == 0) PRINT *, 'Calling prepsglobs'
        CALL prepsglobs(mx_sng,ntime,srcsng,nsrc_sng,                   &
              sngfname(ifile),stnsng,latsng,lonsng,xsng,ysng,           &
              hgtsng,obsng,qobsng,qualsng,isrcsng,qsrcsng,              &
              csrcsng,nxlg,nylg,nz,nvar_anx,anx,xslg,yslg,zp,           &
              tem1,tem2,tem3,tem4,tem5,tem6,                            &
              rmiss,nobsng,istatus)
      END IF
    END DO
  END IF ! myproc == 0

!
!-----------------------------------------------------------------------
!
!  MPI issues.  Broadcast the number of data points, then the data.  If
!  there are no data values, we don't waste anyone's time.
!
!-----------------------------------------------------------------------
!
  CALL mpupdatei(nobsng,1)
  ALLOCATE(item1(nobsng),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:item1:nobsng")
  ALLOCATE(ksng(nprocs),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:ksng:nprocs")

  IF (mp_opt > 0) THEN
!   CALL mpupdatei(nobsng,1)
    IF (nobsng > 0) THEN
      CALL mpupdatei(timesng,mx_sng*ntime)
      CALL mpupdatec(stnsng,mx_sng*ntime*5)
      CALL mpupdatec(csrcsng,mx_sng*ntime*8)
      CALL mpupdatei(isrcsng,mx_sng)
      CALL mpupdater(latsng,mx_sng*ntime)
      CALL mpupdater(lonsng,mx_sng*ntime)
      CALL mpupdater(hgtsng,mx_sng*ntime)
      CALL mpupdater(xsng,mx_sng)
      CALL mpupdater(ysng,mx_sng)
      CALL mpupdatec(wx,mx_sng*ntime*8)
      CALL mpupdatei(kloud,mx_sng*ntime)
      CALL mpupdatei(idp3,mx_sng*ntime)
      CALL mpupdatec(store_emv,mx_sng*5*ntime)
      CALL mpupdater(store_hgt,mx_sng*5*ntime)
      CALL mpupdatec(store_amt,mx_sng*5*ntime*4)
      CALL mpupdater(obrdsng,mx_sng*nvar_sng*ntime)
      CALL mpupdater(obsng,nvar_anx*mx_sng)
      CALL mpupdater(qobsng,nvar_anx*mx_sng)
      CALL mpupdatei(qualsng,nvar_anx*mx_sng)
      CALL mpupdater(knt,nvar_anx*nz)
      CALL mpupdater(wgtsum,nvar_anx*nz)
      CALL mpupdater(zsum,nvar_anx*nz)

!
!-----------------------------------------------------------------------
!      We need to know which processor "owns" which obs, so they can be
!      consulted on an "as needed" basis.
!-----------------------------------------------------------------------
!

      CALL mpiprocess(nobsng,indexsng,usesng,nprocs,ksng,ksngmax,       &
                      isrcsng,item1,nx,ny,xsng,ysng,xs,ys,0)

      CALL data_filt(nobsng,isrcsng,indexsng,nmap,mpi_map)

!
!-----------------------------------------------------------------------
!
!  Although we're similar to PREPSFCGLOBS, and need to do one final task
!  (compute heights), we now use the local grid variables, not the
!  global.
!
!-----------------------------------------------------------------------
!
      CALL prepsglmdcrs_sm(mx_sng,ntime,latsng,lonsng,                  &
                           xsng,ysng,hgtsng,csrcsng,                    &
                           nx,ny,nz,nvar_anx,anx,xs,ys,zp,              &
                           tem1,tem2,tem3,tem4,tem5,tem6,               &
                           rmiss,nobsng,indexsng,istatus)

      !
      ! Collect and distribute the data.
      !
      ALLOCATE(tems1dr(ksngmax),STAT=istatus)
      CALL check_alloc_status(istatus, "mpi_1dr_collect:tmps")
      ALLOCATE(temr1dr(ksngmax),STAT=istatus)
      CALL check_alloc_status(istatus, "mpi_1dr_collect:tmpr")

      CALL mpi_1dr_collect(hgtsng,nobsng,indexsng,                       &
                 nprocs, ksng, ksngmax, mpi_map, nmap, tems1dr, temr1dr)

!
! Variables temps1dr/temr1dr aren't deallocated yet, as we'll reused then
! below for more single point processing.
!

    END IF  ! nobsng > 0

  ELSE
    CALL mpiprocess(nobsng,indexsng,usesng,nprocs,ksng,ksngmax,         &
                    isrcsng,item1,nx,ny,xsng,ysng,xs,ys,0)

  END IF  ! mp_opt > 0

  DEALLOCATE (item1)
!
!-----------------------------------------------------------------------
!
!  Calculate initial obs differences for single-level data
!
!-----------------------------------------------------------------------
!

  IF (myproc == 0 ) PRINT *, 'Calling grdtosng'
  CALL grdtosng(nx,ny,nxlg,nylg,nz,nz_tab,mx_sng,nvar_anx,nobsng,       &
                    xs,ys,xslg,yslg,zp,icatg,anx,qback,hqback,nlvqback, &
                    tem1,tem2,tem3,tem4,tem5,tem6,                      &
                    stnsng,isrcsng,icatsng,hgtsng,xsng,ysng,            &
                    obsng,qobsng,qualsng,                               &
                    odifsng,oanxsng,thesng,trnsng,indexsng)
  IF (mp_opt > 0) CALL mpbarrier

!
!-----------------------------------------------------------------------
!      Some data needs to be collected and distributed.
!-----------------------------------------------------------------------
!

  IF ( mp_opt > 0 .AND. nobsng > 0 ) THEN
    ALLOCATE(tems2drsng(nvar_anx,ksngmax),STAT=istatus)
    CALL check_alloc_status(istatus, "mpi_2dr_collect:tmpssng")
    ALLOCATE(temr2drsng(nvar_anx,ksngmax),STAT=istatus)
    CALL check_alloc_status(istatus, "mpi_2dr_collect:tmprsng")

!
!  1D real temporary arrays already allocated, but 1D ints aren't.
!

    ALLOCATE(tems1di(ksngmax),STAT=istatus)
    CALL check_alloc_status(istatus, "mpi_1di_collect:tmps")
    ALLOCATE(temr1di(ksngmax),STAT=istatus)
    CALL check_alloc_status(istatus, "mpi_1di_collect:tmpr")

    CALL mpi_2dr_collect( qobsng, nvar_anx, mx_sng, nobsng, indexsng,  &
      nprocs, ksng, ksngmax, mpi_map, nmap, tems2drsng, temr2drsng)
    CALL mpi_2dr_collect( odifsng, nvar_anx, mx_sng, nobsng, indexsng,  &
      nprocs, ksng, ksngmax, mpi_map, nmap, tems2drsng, temr2drsng)
    CALL mpi_2dr_collect( oanxsng, nvar_anx, mx_sng, nobsng, indexsng,  &
      nprocs, ksng, ksngmax, mpi_map, nmap, tems2drsng, temr2drsng)

    CALL mpi_1di_collect( icatsng, nobsng, indexsng,                    &
      nprocs, ksng, ksngmax, mpi_map, nmap, tems1di, temr1di)

    CALL mpi_1dr_collect( thesng, nobsng, indexsng,                     &
      nprocs, ksng, ksngmax, mpi_map, nmap, tems1dr, temr1dr)
    CALL mpi_1dr_collect( trnsng, nobsng, indexsng,                     &
      nprocs, ksng, ksngmax, mpi_map, nmap, tems1dr, temr1dr)

!
!   2D space can't be deallocated, as the arrays get passed to "anxiter".
!

    DEALLOCATE(tems1di)
    DEALLOCATE(temr1di)
    DEALLOCATE(tems1dr)
    DEALLOCATE(temr1dr)
  ELSE
    ALLOCATE(tems2drsng(1,1),STAT=istatus)
    CALL check_alloc_status(istatus, "mpi_2dr_collect:tmpssng")
    ALLOCATE(temr2drsng(1,1),STAT=istatus)
    CALL check_alloc_status(istatus, "mpi_2dr_collect:tmprsng")
  END IF

!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Read upper-air data, QC and convert units
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
  IF (myproc == 0) THEN
    PRINT *, 'Calling prepuaobs'
    CALL prepuaobs(nx,ny,nz,nvar_anx,                                  &
                   nz_ua,mx_ua,nz_tab,nsrc_ua,mx_ua_file,              &
                   anx,xs,ys,zp,tem1,tem2,tem3,tem4,tem5,tem6,tem7,    &
                   nuafil,uafname,srcua,                               &
                   stnua,elevua,xua,yua,hgtua,obsua,                   &
                   qsrcua,huaqsrc,nlvuatab,                            &
                   qobsua,qualua,isrcua,nlevsua,                       &
                   rmiss,nobsua,istatus)
  END IF

!
!-----------------------------------------------------------------------
!
!  MPI issues.  Broadcast the number of data points, then the data.  If
!  there are no data values, we don't waste anyone's time.
!
!-----------------------------------------------------------------------
!

  CALL mpupdatei(nobsua,1)

  ALLOCATE(item1(nobsua),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:item1:nobsua")
  ALLOCATE(kua(nprocs),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:kua:nprocs")

  IF (mp_opt > 0) THEN
    IF (nobsua > 0) THEN
      CALL mpupdatec(stnua,mx_ua*5)
      CALL mpupdater(elevua,mx_ua)
      CALL mpupdater(xua,mx_ua)
      CALL mpupdater(yua,mx_ua)
      CALL mpupdater(hgtua,nz_ua*mx_ua)
      CALL mpupdater(obsua,nvar_anx*nz_ua*mx_ua)
      CALL mpupdater(qobsua,nvar_anx*nz_ua*mx_ua)
      CALL mpupdatei(nlvuatab,nsrc_ua)
      CALL mpupdatei(qualua,nvar_anx*nz_ua*mx_ua)
      CALL mpupdatei(isrcua,mx_ua)
      CALL mpupdatei(nlevsua,mx_ua)

!
!-----------------------------------------------------------------------
!      We need to know which processor "owns" which obs, so they can be
!      consulted on an "as needed" basis.
!-----------------------------------------------------------------------
!

      CALL mpiprocess(nobsua,indexua,useua,nprocs,kua,kuamax,           &
          isrcua,item1,nx,ny,xua,yua,xs,ys,0)

      CALL data_filt(nobsua,isrcua,indexua,nmap,mpi_map)

    END IF
  ELSE
    CALL mpiprocess(nobsua,indexua,useua,nprocs,kua,kuamax,             &
                    isrcua,item1,nx,ny,xua,yua,xs,ys,0)

  END IF

  DEALLOCATE (item1)
!
!-----------------------------------------------------------------------
!
!  Calculate initial obs differences for upper-air data
!
!-----------------------------------------------------------------------
!
  IF (myproc == 0) PRINT *, 'Calling grdtoua'

  CALL grdtoua(nx,ny,nxlg,nylg,nz,nz_tab,nz_ua,mx_ua,nvar_anx,nobsua,   &
              xs,ys,xslg,yslg,zp,anx,qback,hqback,nlvqback,             &
              tem1,tem2,tem3,tem4,tem5,tem6,                            &
              stnua,isrcua,elevua,xua,yua,hgtua,                        &
              obsua,qobsua,qualua,nlevsua,                              &
              odifua,oanxua,theua,trnua,indexua)
  IF (mp_opt > 0) CALL mpbarrier
!
!-----------------------------------------------------------------------
!      Some data needs to be collected and distributed.
!-----------------------------------------------------------------------
!

  IF ( mp_opt > 0 .AND. nobsua > 0 ) THEN
    ALLOCATE(tems3drua(nvar_anx,nz_ua,kuamax),STAT=istatus)
    CALL check_alloc_status(istatus, "mpi_3dr_collect:tmpsua")
    ALLOCATE(temr3drua(nvar_anx,nz_ua,kuamax),STAT=istatus)
    CALL check_alloc_status(istatus, "mpi_3dr_collect:tmprua")

    ALLOCATE(tems1dr(kuamax),STAT=istatus)
    CALL check_alloc_status(istatus, "mpi_1dr_collect:tmps")
    ALLOCATE(temr1dr(kuamax),STAT=istatus)
    CALL check_alloc_status(istatus, "mpi_1dr_collect:tmpr")

    ALLOCATE(tems2drua(nz_ua,kuamax),STAT=istatus)
    CALL check_alloc_status(istatus, "mpi_2dr_collect:tmpsua")
    ALLOCATE(temr2drua(nz_ua,kuamax),STAT=istatus)
    CALL check_alloc_status(istatus, "mpi_2dr_collect:tmprua")

    ALLOCATE(tems3ditmp(nvar_anx,nz_ua,kuamax),STAT=istatus)
    CALL check_alloc_status(istatus, "mpi_3di_collect:tems3ditmp")
    ALLOCATE(temr3ditmp(nvar_anx,nz_ua,kuamax),STAT=istatus)
    CALL check_alloc_status(istatus, "mpi_3di_collect:temr3ditmp")

    CALL mpi_3dr_collect( qobsua,  nvar_anx, nz_ua, mx_ua, nobsua, indexua, &
      nprocs, kua, kuamax, mpi_map, nmap, tems3drua, temr3drua)
    CALL mpi_3di_collect( qualua,  nvar_anx, nz_ua, mx_ua, nobsua, indexua, &
      nprocs, kua, kuamax, mpi_map, nmap, tems3ditmp, temr3ditmp)
    CALL mpi_3dr_collect( odifua,  nvar_anx, nz_ua, mx_ua, nobsua, indexua, &
      nprocs, kua, kuamax, mpi_map, nmap, tems3drua, temr3drua)
    CALL mpi_3dr_collect( oanxua,  nvar_anx, nz_ua, mx_ua, nobsua, indexua, &
      nprocs, kua, kuamax, mpi_map, nmap, tems3drua, temr3drua)
    CALL mpi_1dr_collect( trnua, nobsua, indexua,                           &
       nprocs, kua, kuamax, mpi_map, nmap, tems1dr, temr1dr)
    CALL mpi_2dr_collect( theua, nz_ua, mx_ua, nobsua, indexua,  &
      nprocs, kua, kuamax, mpi_map, nmap, tems2drua, temr2drua)

    DEALLOCATE(tems1dr)
    DEALLOCATE(temr1dr)
    DEALLOCATE(tems2drua)
    DEALLOCATE(temr2drua)
  ELSE
    ALLOCATE(tems3drua(1,1,1),STAT=istatus)
    CALL check_alloc_status(istatus, "mpi_3dr_collect:tmpsua")
    ALLOCATE(temr3drua(1,1,1),STAT=istatus)
    CALL check_alloc_status(istatus, "mpi_3dr_collect:tmprua")
  END IF

!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Read radar data, unfold and convert units
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
  IF (myproc == 0) PRINT *, 'Calling prepradar'
!
!-----------------------------------------------------------------------
!      We need to know which processor "owns" which obs, so they can be
!      consulted on an "as needed" basis.  Unlike non-radar data, the
!      processing is different, and is done in "prepradar" since all the
!      bookkeeping we need is done there in temporary arrays.  Since the
!      "krad" array is small, and gets passed around a few times, we'll
!      always allocate it.
!-----------------------------------------------------------------------
!
  ALLOCATE(krad(nprocs),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:rad:nprocs")

  IF(cloudopt > 0) THEN
    print *, ' allocating full ref_mos_3d',nx,ny,nz
    ALLOCATE(ref_mos_3d(nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:ref_mos_3d")
    ALLOCATE(refmax_mos_3d(nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:refmax_mos_3d")
    ALLOCATE(rhv_mos_3d(nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:rhv_mos_3d")
    ALLOCATE(zdr_mos_3d(nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:zdr_mos_3d")
    ALLOCATE(kdp_mos_3d(nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:kdp_mos_3d")
    ALLOCATE(dist_mos_3d(nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:dist_mos_3d")
  ELSE
    print *, ' allocating placeholder ref_mos_3d'
    ALLOCATE(ref_mos_3d(1,1,1),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:ref_mos_3d")
    ALLOCATE(refmax_mos_3d(1,1,1),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:refmax_mos_3d")
    ALLOCATE(rhv_mos_3d(1,1,1),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:rhv_mos_3d")
    ALLOCATE(zdr_mos_3d(1,1,1),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:zdr_mos_3d")
    ALLOCATE(kdp_mos_3d(1,1,1),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:kdp_mos_3d")
    ALLOCATE(dist_mos_3d(1,1,1),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:dist_mos_3d")
  END IF
  ref_mos_3d=-99.
  refmax_mos_3d=-99.
  rhv_mos_3d=-99.
  zdr_mos_3d=-99.
  kdp_mos_3d=-99.
  dist_mos_3d=-99.

  CALL prepradar(nx,ny,nz,nz_tab,nvar_anx,nvar_radin,                   &
             nvar_rad,mx_rad,nsrc_rad,nz_rdr,mx_colrad,mx_pass,         &
             raduvobs,radrhobs,radistride,radkstride,cloudopt,          &
             iuserad,npass,refrh,rhradobs,                              &
             xs,ys,zs,hterain,latgr,longr,anx,qback,hqback,nlvqback,    &
             nradfil,radfname,srcrad,isrcrad,qsrcrad,qcthrrad,          &
             stnrad,latrad,lonrad,elvrad,                               &
             latradc,lonradc,xradc,yradc,irad,nlevrad,                  &
             distrad,uazmrad,vazmrad,hgtradc,theradc,trnradc,dsdr,dhdr, &
             obsrad,oanxrad,odifrad,qobsrad,qualrad,                    &
             ncolrad,ncolrad_mpi,indexrad,oindexrad,                    &
             refmax_mos_3d,ref_mos_3d,                                  &
             rhv_mos_3d,zdr_mos_3d,kdp_mos_3d,dist_mos_3d,              &
             iprocv,jprocv,iradvel,                                     &
             istatus,istat_radar,tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8)

  ALLOCATE(item1(ncolrad),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:item1:ncolrad")

  IF ( mp_opt > 0 ) THEN

    IF ( ncolrad_mpi > 0 ) THEN
!
!-----------------------------------------------------------------------
!      Some data needs to be collected and distributed.
!-----------------------------------------------------------------------
!

      IF ( iradvel == 0) THEN
        CALL mpiprocess(ncolrad,indexrad,userad,nprocs,krad,kradmax,      &
                        isrcrad,item1,nx,ny,xradc,yradc,xs,ys,1)
      ELSE  ! note that the array indexrad is redefined from above
        CALL mpiprocess_complex(ncolrad_mpi,ncolrad,indexrad,oindexrad,   &
                        userad, nprocs,krad,kradmax,isrcrad,item1,        &
                        nx,ny,xradc,yradc,xs,ys,1,istatus)
      ENDIF

      ! 1D integer arrays

      ALLOCATE(tems1di(kradmax),STAT=istatus)
      CALL check_alloc_status(istatus, "mpi_1di_collect:tmps")
      ALLOCATE(temr1di(kradmax),STAT=istatus)
      CALL check_alloc_status(istatus, "mpi_1di_collect:tmpr")

      CALL mpi_1di_collect( irad, ncolrad, indexrad,                      &
        nprocs, krad, kradmax, mpi_mapv, nmapv, tems1di, temr1di)
      CALL mpi_1di_collect( nlevrad, ncolrad, indexrad,                   &
        nprocs, krad, kradmax, mpi_mapv, nmapv, tems1di, temr1di)

      DEALLOCATE(tems1di)
      DEALLOCATE(temr1di)

      ! 1D real arrays

      ALLOCATE(tems1dr(kradmax),STAT=istatus)
      CALL check_alloc_status(istatus, "mpi_1dr_collect:tmps")
      ALLOCATE(temr1dr(kradmax),STAT=istatus)
      CALL check_alloc_status(istatus, "mpi_1dr_collect:tmpr")

      CALL mpi_1dr_collect( distrad, ncolrad, indexrad,                   &
        nprocs, krad, kradmax, mpi_mapv, nmapv, tems1dr, temr1dr)
      CALL mpi_1dr_collect( uazmrad, ncolrad, indexrad,                   &
        nprocs, krad, kradmax, mpi_mapv, nmapv, tems1dr, temr1dr)
      CALL mpi_1dr_collect( vazmrad, ncolrad, indexrad,                   &
        nprocs, krad, kradmax, mpi_mapv, nmapv, tems1dr, temr1dr)
      CALL mpi_1dr_collect( latradc, ncolrad, indexrad,                   &
        nprocs, krad, kradmax, mpi_mapv, nmapv, tems1dr, temr1dr)
      CALL mpi_1dr_collect( lonradc, ncolrad, indexrad,                   &
        nprocs, krad, kradmax, mpi_mapv, nmapv, tems1dr, temr1dr)
      CALL mpi_1dr_collect( xradc, ncolrad, indexrad,                     &
        nprocs, krad, kradmax, mpi_mapv, nmapv, tems1dr, temr1dr)
      CALL mpi_1dr_collect( yradc, ncolrad, indexrad,                     &
        nprocs, krad, kradmax, mpi_mapv, nmapv, tems1dr, temr1dr)
      CALL mpi_1dr_collect( trnradc, ncolrad, indexrad,                   &
        nprocs, krad, kradmax, mpi_mapv, nmapv, tems1dr, temr1dr)

      DEALLOCATE(tems1dr)
      DEALLOCATE(temr1dr)

      ! 2D real arrays

      ALLOCATE(tems2dr(nz_rdr,kradmax),STAT=istatus)
      CALL check_alloc_status(istatus, "mpi_2dr_collect:tmps")
      ALLOCATE(temr2dr(nz_rdr,kradmax),STAT=istatus)
      CALL check_alloc_status(istatus, "mpi_2dr_collect:tmrs")

      CALL mpi_2dr_collect( hgtradc, nz_rdr, mx_colrad, ncolrad,          &
        indexrad, nprocs, krad, kradmax, mpi_mapv, nmapv, tems2dr, temr2dr)

      CALL mpi_2dr_collect( theradc, nz_rdr, mx_colrad, ncolrad,          &
        indexrad, nprocs, krad, kradmax, mpi_mapv, nmapv, tems2dr, temr2dr)

      DEALLOCATE(tems2dr)
      DEALLOCATE(temr2dr)

      ! 3D real variables, two passes due to two different sizes.  One group
      ! is preserved for later use.

      ALLOCATE(tems3drtmp(nvar_radin,nz_rdr,kradmax),STAT=istatus)
      CALL check_alloc_status(istatus, "mpi_3dr_collect:tems3drtmp")
      ALLOCATE(temr3drtmp(nvar_radin,nz_rdr,kradmax),STAT=istatus)
      CALL check_alloc_status(istatus, "mpi_3dr_collect:temr3drtmp")

      CALL mpi_3dr_collect(obsrad, nvar_radin, nz_rdr, mx_colrad, ncolrad,  &
        indexrad, nprocs, krad, kradmax, mpi_mapv, nmapv, tems3drtmp, temr3drtmp)

      DEALLOCATE(tems3drtmp)
      DEALLOCATE(temr3drtmp)

      ALLOCATE(tems3drrad(nvar_anx,nz_rdr,kradmax),STAT=istatus)
      CALL check_alloc_status(istatus, "mpi_3dr_collect:tmpsrad")
      ALLOCATE(temr3drrad(nvar_anx,nz_rdr,kradmax),STAT=istatus)
      CALL check_alloc_status(istatus, "mpi_3dr_collect:tmprrad")

      CALL mpi_3dr_collect(oanxrad, nvar_rad, nz_rdr, mx_colrad, ncolrad,  &
        indexrad, nprocs, krad, kradmax, mpi_mapv, nmapv, tems3drrad, temr3drrad)
      CALL mpi_3dr_collect(odifrad, nvar_rad, nz_rdr, mx_colrad, ncolrad,  &
        indexrad, nprocs, krad, kradmax, mpi_mapv, nmapv, tems3drrad, temr3drrad)
      CALL mpi_3dr_collect(qobsrad, nvar_rad, nz_rdr, mx_colrad, ncolrad,  &
        indexrad, nprocs, krad, kradmax, mpi_mapv, nmapv, tems3drrad, temr3drrad)

      ALLOCATE(tems3dirad(nvar_rad,nz_rdr,kradmax),STAT=istatus)
      CALL check_alloc_status(istatus, "mpi_3di_collect:tmpsrad")
      ALLOCATE(temr3dirad(nvar_rad,nz_rdr,kradmax),STAT=istatus)
      CALL check_alloc_status(istatus, "mpi_3di_collect:tmprrad")

      CALL mpi_3di_collect(qualrad, nvar_rad, nz_rdr, mx_colrad, ncolrad,  &
        indexrad, nprocs, krad, kradmax, mpi_mapv, nmapv, tems3dirad, temr3dirad)

      CALL data_filt(ncolrad,irad,indexrad,nmapv,mpi_mapv)

    ELSE    ! no observations
      ALLOCATE(tems3drrad(1,1,1),STAT=istatus)
      CALL check_alloc_status(istatus, "mpi_3dr_collect:tmpsrad")
      ALLOCATE(temr3drrad(1,1,1),STAT=istatus)
      CALL check_alloc_status(istatus, "mpi_3dr_collect:tmprrad")

      ALLOCATE(tems3dirad(1,1,1),STAT=istatus)
      CALL check_alloc_status(istatus, "mpi_3di_collect:tmpsrad")
      ALLOCATE(temr3dirad(1,1,1),STAT=istatus)
      CALL check_alloc_status(istatus, "mpi_3di_collect:tmprrad")
    END IF  ! there are observations

  ELSE  ! set array indexrad & userad for 3DVAR
    CALL mpiprocess(ncolrad,indexrad,userad,nprocs,krad,kradmax,        &
                    isrcrad,item1,nx,ny,xradc,yradc,xs,ys,1)
    ALLOCATE(tems3drrad(1,1,1),STAT=istatus)
    CALL check_alloc_status(istatus, "mpi_3dr_collect:tmpsrad")
    ALLOCATE(temr3drrad(1,1,1),STAT=istatus)
    CALL check_alloc_status(istatus, "mpi_3dr_collect:tmprrad")

    ALLOCATE(tems3dirad(1,1,1),STAT=istatus)
    CALL check_alloc_status(istatus, "mpi_3di_collect:tmpsrad")
    ALLOCATE(temr3dirad(1,1,1),STAT=istatus)
    CALL check_alloc_status(istatus, "mpi_3di_collect:tmprrad")

  END IF

  DEALLOCATE(item1)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Read retrieval data.
!
!  Retrieval code has *NOT* been tested, so it *WILL* not be MPI ready!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
  IF (myproc == 0) PRINT *, 'Calling prepretr'
  CALL prepretr(nx,ny,nz,nvar_anx,                                      &
                nz_ret,mx_ret,mx_colret,nz_tab,nsrc_ret,                &
                nretfil,retfname,                                       &
                isrcret,srcret,nlvrttab,qsrcret,hrtqsrc,                &
                stnret,latret,lonret,elvret,                            &
                latretc,lonretc,xretc,yretc,iret,nlevret,               &
                hgtretc,obsret,qrret,qobsret,qualret,                   &
                rmiss,ncolret,tem1,istatus)

!
!-----------------------------------------------------------------------
!
!  Now get the retrieval data so that each processor only
!  contains the data it needs for its local domain plus the radius of
!  influence, where applicable.
!
!  Not guaranteed MPI ready.
!
!-----------------------------------------------------------------------
!

  IF (myproc == 0) PRINT *, 'Calling grdtoret'
  CALL grdtoret(nx,ny,nz,nz_tab,                                        &
                nz_ret,mx_ret,mx_colret,nvar_anx,ncolret,               &
                xs,ys,zp,anx,qback,hqback,nlvqback,                     &
                tem1,tem2,tem3,tem4,tem5,tem6,                          &
                stnret,iret,xretc,yretc,hgtretc,                        &
                obsret,qobsret,qualret,nlevret,                         &
                odifret,oanxret,theretc,trnretc)
  IF (mp_opt > 0) CALL mpbarrier
!
!-----------------------------------------------------------------------
!
! convert observation from q to RH  !Ming Hu
!-----------------------------------------------------------------------
!
  IF( cntl_var_rh == 1 ) THEN
    CALL qv2rhobs(nvar_anx_3dvar,mx_sng,nobsng,obsng,qualsng,           &
                  nz_ua,mx_ua,nobsua,obsua,qualua,nlevsua,              &
                  nz_ret,mx_colret,ncolret,obsret,qualret,nlevret)
    !write(*,*) 'RH=',obsng(5,804),obsng(5,1194)
  ENDIF
!
!-----------------------------------------------------------------------
!
!  Quality-control observation differences
!
!-----------------------------------------------------------------------
!

  IF (myproc == 0) PRINT *, 'Calling qcdiff'
  CALL qcdiff(nvar_anx,nvar_rad,nvar_radin,mx_sng,nsrc_sng,             &
              indexsng,indexua,indexrad,indexret,                       &
              usesng,useua,userad,                                      &
              nz_ua,mx_ua,nsrc_ua,                                      &
              nz_rdr,mx_rad,mx_colrad,nsrc_rad,                         &
              nz_ret,mx_ret,mx_colret,nsrc_ret,                         &
              nobsng,nobsua,ncolrad,ncolret,var_anx,                    &
              stnsng,isrcsng,hgtsng,obsng,oanxsng,odifsng,              &
              qcthrsng,qualsng,                                         &
              stnua,isrcua,hgtua,obsua,oanxua,odifua,                   &
              qcthrua,qualua,nlevsua,                                   &
              stnrad,irad,isrcrad,hgtradc,obsrad,odifrad,               &
              qcthrrad,qualrad,nlevrad,                                 &
              stnret,iret,isrcret,hgtretc,                              &
              obsret,oanxret,odifret,                                   &
              qcthrret,qualret,nlevret,                                 &
              wgtsum,zsum)

  IF ( mp_opt > 0 .AND. nobsng > 0) THEN
    CALL mpi_2dr_collect( qualsng, nvar_anx, mx_sng, nobsng, indexsng,  &
      nprocs, ksng, ksngmax, mpi_map, nmap, tems2drsng, temr2drsng)
  END IF
  IF ( mp_opt > 0 .AND. nobsua > 0 ) THEN
    CALL mpi_3di_collect( qualua,  nvar_anx, nz_ua, mx_ua, nobsua, indexua, &
      nprocs, kua, kuamax, mpi_map, nmap, tems3ditmp, temr3ditmp)
    DEALLOCATE(tems3ditmp)
    DEALLOCATE(temr3ditmp)
  END IF
!
!-----------------------------------------------------------------------
!
!  Analyze data
!
!-----------------------------------------------------------------------
!

#ifdef arps3dvar

  DEALLOCATE( tems3dirad,temr3dirad,tems3drrad,temr3drrad )  ! we donot use them for 3dvar
  DEALLOCATE( tems3drua, temr3drua )
  DEALLOCATE( tems2drsng,temr2drsng )

!-----------------------------------------------------------------------
!
!  allocate variable array for 3DVAR only
!
!-----------------------------------------------------------------------
!
  CALL allocateArray4Var(nx,ny,nz)

  ALLOCATE(tem12(nx,ny,nz), STAT = istatus)
  CALL check_alloc_status(istatus, "before minimization: tem12")
  ALLOCATE(tem13(nx,ny,nz), STAT = istatus)
  CALL check_alloc_status(istatus, "before minimization: tem13")
  ALLOCATE(tem14(nx,ny,nz), STAT = istatus)
  CALL check_alloc_status(istatus, "before minimization: tem14")
  ALLOCATE(tem15(nx,ny,nz), STAT = istatus)
  CALL check_alloc_status(istatus, "before minimization: tem15")
  ALLOCATE(tem16(nx,ny,nz), STAT = istatus)
  CALL check_alloc_status(istatus, "before minimization: tem16")
  ALLOCATE(tem17(nx,ny,nz), STAT = istatus)
  CALL check_alloc_status(istatus, "before minimization: tem17")
  ALLOCATE(tem18(nx,ny,nz), STAT = istatus)
  CALL check_alloc_status(istatus, "before minimization: tem18")
!
!-----------------------------------------------------------------------
!
!  Cycle in npass
!
!-----------------------------------------------------------------------
!
  DO ipass=1,npass
    IF (myproc == 0) WRITE(6,'(1x,a,I2)') 'Analysis pass is:',ipass

    cpuin=f_cputime()
    CALL FLUSH (06)

    IF (myproc == 0) WRITE(6,'(1x,a)') 'Calling minimization ...'
    CALL minimization(nx,ny,nz,                                         &
            nvar_anx_3dvar,nvar_radin,nvar_rad_3dvar,nz_ua,nz_rdr,nz_ret, &
            mapfct,j1,j2,j3,aj3x,aj3y,aj3z,j3inv,rhostr,                &
            mx_sng,mx_ua,mx_rad,mx_colrad,mx_ret,mx_colret,             &
            nsrc_sng,nsrc_ua,nsrc_rad,nsrc_ret,ncat,                    &
            mx_pass,ipass,iwstat,xs,ys,zs,x,y,z,zp,hterain,             &
            icatg,xcor,var_anx,                                         &
            indexsng,usesng,xsng,ysng,hgtsng,thesng,                    &
            odifsng,qobsng,qualsng,isrcsng,icatsng,nobsng,              &
            indexua,useua,xua,yua,hgtua,theua,                          &
            odifua,qobsua,qualua,isrcua,nlevsua,nobsua,                 &
            indexrad,userad,elvrad,xradc,yradc,                         &
            distrad,uazmrad,vazmrad,hgtradc,theradc,dsdr,dhdr,          &
            obsrad,odifrad,qobsrad,qualrad,                             &
            irad,isrcrad,nlevrad,ncolrad,                               &
            xretc,yretc,hgtretc,theretc,                                &
            odifret,qobsret,qualret,                                    &
            iret,isrcret,nlevret,ncolret,                               &
            srcsng,srcua,srcrad,srcret,                                 &
            ianxtyp,iusesng,iuseua,iuserad,iuseret,                     &
            xyrange,kpvrsq,wlim,zrange,zwlim,                           &
            thrng,rngsqi,knt,wgtsum,zsum,                               &
            corsng,corua,corrad,corret,                                 &
            oanxsng,oanxua,oanxrad,oanxret,                             &
            nlvqback,hqback,qback,                                      &
            nsngfil,nuafil,nradfil,raduvobs,nretfil,                    &
            anx,tem1,tem2,tem3,tem4,tem5,tem6,tem7,                     &
            tem8,tem9,tem10,tem11,tem12,tem13,tem14,tem15,tem16,tem17,  &
            tem18,istatus)

    cpuout=f_cputime()
    IF (myproc == 0) write(*,*) '=========== out minimization ==========='

    IF (ipass /= npass ) THEN

!
!-----------------------------------------------------------------------
!
!  Calculate initial obs differences
!
!-----------------------------------------------------------------------
!
      !
      !  for single-level data
      !
      CALL grdtosng(nx,ny,nxlg,nylg,nz,nz_tab,mx_sng,nvar_anx,nobsng,   &
                    xs,ys,xslg,yslg,zp,icatg,anx,qback,hqback,nlvqback, &
                    tem1,tem2,tem3,tem4,tem5,tem6,                      &
                    stnsng,isrcsng,icatsng,hgtsng,xsng,ysng,            &
                    obsng,qobsng,qualsng,                               &
                    odifsng,oanxsng,thesng,trnsng,indexsng)
      !
      !  for upper level data
      !
      CALL grdtoua(nx,ny,nxlg,nylg,nz,nz_tab,nz_ua,mx_ua,nvar_anx,nobsua,&
                   xs,ys,xslg,yslg,zp,anx,qback,hqback,nlvqback,        &
                   tem1,tem2,tem3,tem4,tem5,tem6,                       &
                   stnua,isrcua,elevua,xua,yua,hgtua,                   &
                   obsua,qobsua,qualua,nlevsua,                         &
                   odifua,oanxua,theua,trnua,indexua)
      !
      !  for radar data
      !
      CALL prcssrad(nx,ny,nz,nz_tab,nvar_anx,nvar_radin,                &
                  nvar_rad,mx_rad,nsrc_rad,nz_rdr,mx_colrad,            &
                  raduvobs,radrhobs,refrh,rhradobs,                     &
                  xs,ys,zs,hterain,anx,qback,hqback,nlvqback,           &
                  srcrad,qsrcrad,qcthrrad,                              &
                  isrcrad,stnrad,latrad,lonrad,elvrad,                  &
                  latradc,lonradc,xradc,yradc,irad,nlevrad,             &
                  hgtradc,obsrad,oanxrad,odifrad,                       &
                  distrad,uazmrad,vazmrad,                              &
                  theradc,trnradc,qobsrad,qualrad,dsdr,dhdr,            &
                  ncolrad,indexrad,istatus,tem1)
!
      IF (myproc == 0) PRINT *, 'Calling qcdiff'
      CALL qcdiff(nvar_anx,nvar_rad,nvar_radin,mx_sng,nsrc_sng,         &
                  indexsng,indexua,indexrad,indexret,                   &
                  usesng, useua,userad,                                 &
                  nz_ua,mx_ua,nsrc_ua,                                  &
                  nz_rdr,mx_rad,mx_colrad,nsrc_rad,                     &
                  nz_ret,mx_ret,mx_colret,nsrc_ret,                     &
                  nobsng,nobsua,ncolrad,ncolret,var_anx,                &
                  stnsng,isrcsng,hgtsng,obsng,oanxsng,odifsng,          &
                  qcthrsng,qualsng,                                     &
                  stnua,isrcua,hgtua,obsua,oanxua,odifua,               &
                  qcthrua,qualua,nlevsua,                               &
                  stnrad,irad,isrcrad,hgtradc,obsrad,odifrad,           &
                  qcthrrad,qualrad,nlevrad,                             &
                  stnret,iret,isrcret,hgtretc,                          &
                  obsret,oanxret,odifret,                               &
                  qcthrret,qualret,nlevret,                             &
                  wgtsum,zsum)
    END IF

    cpuout=f_cputime()
    IF (myproc == 0) THEN
       write(*,*) '=========== out qcdiff ==========='
       write(*,*) 'The cpu time used by ',ipass,' analysis is:',cpuout-cpuin
       write(*,*) ' The cost function time is:', cputimearray(1)
       write(*,*) ' The gradient part is:     ', cputimearray(2)
       write(*,*) ' The minimization part is: ', cputimearray(3)
    END IF
    cputimearray=0.0

  END DO ! ipass

  CALL deallocateArray4Var
  DEALLOCATE( tem12, tem13, tem14, tem15, tem16, tem17 )

  IF (myproc == 0) WRITE(6,*) ' The end of 3DVar analysis '

#else

!!
!!!!!!!!!!!!!!!!!!!!!! This is a branch for adas analysis !!!!!!!!!!!!!!
!!
  IF (myproc == 0) WRITE(6,*) 'Calling anxiter'
  CALL FLUSH (6)

!
! Notice that we are passing "nxlg/nylg" instead of "nx/ny".
!

  CALL anxiter(nx,ny,nz,                                                &
          nvar_anx,nvar_radin,nvar_rad,nz_ua,nz_rdr,nz_ret,             &
          mx_sng,mx_ua,mx_rad,mx_colrad,mx_ret,mx_colret,               &
          nsrc_sng,nsrc_ua,nsrc_rad,nsrc_ret,ncat,                      &
          mx_pass,npass,iwstat,xs,ys,zs,icatg,xcor,var_anx,             &
          mpi_map,nmap,mpi_mapv,nmapv,                                  &
          xsng,ysng,hgtsng,thesng,trnsng,                               &
          obsng,odifsng,qobsng,qualsng,isrcsng,icatsng,nobsng,          &
          indexsng,nprocs,ksng,ksngmax,                                 &
          indexua,kua,kuamax,                                           &
          indexrad,krad,kradmax,                                        &
          xua,yua,hgtua,theua,trnua,                                    &
          obsua,odifua,qobsua,qualua,isrcua,nlevsua,nobsua,             &
          elvrad,xradc,yradc,                                           &
          distrad,uazmrad,vazmrad,hgtradc,theradc,trnradc,              &
          obsrad,odifrad,qobsrad,qualrad,                               &
          irad,isrcrad,nlevrad,ncolrad,ncolrad_mpi,                     &
          xretc,yretc,hgtretc,theretc,trnretc,                          &
          obsret,odifret,qobsret,qualret,                               &
          iret,isrcret,nlevret,ncolret,                                 &
          srcsng,srcua,srcrad,srcret,                                   &
          ianxtyp,iusesng,iuseua,iuserad,iuseret,                       &
          xyrange,kpvrsq,wlim,zrange,zwlim,                             &
          thrng,trnropt,trnrcst,trnrng,rngsqi,knt,wgtsum,zsum,          &
          corsng,corua,corrad,corret,                                   &
          oanxsng,oanxua,oanxrad,oanxret,                               &
          anx,tem1,tem2,tem3,ibegin,iend,                               &
          tems2drsng,temr2drsng,tems3drua,temr3drua,                    &
          tems3dirad,temr3dirad,tems3drrad,temr3drrad,istatus)

  DEALLOCATE( tems3dirad,temr3dirad,tems3drrad,temr3drrad )
  DEALLOCATE( tems3drua, temr3drua )
  DEALLOCATE( tems2drsng,temr2drsng )

#endif
!
!-----------------------------------------------------------------------
!
!  convert relative humidity to specific humidity( modified by Ming Hu)
!
!-----------------------------------------------------------------------
!
  IF( cntl_var_rh == 1 ) THEN
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tem10(i,j,k)= anx(i,j,k,PTR_PT)*(anx(i,j,k,PTR_P)/p0)**rddcp
        END DO
      END DO
    END DO
    CALL getqvs(nx,ny,nz, 1,nx-1,1,ny-1,1,nz-1, anx(1,1,1,PTR_P),tem10,qv)
    DO k=1,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          IF( anx(i,j,k,PTR_QV) > 1.0 ) anx(i,j,k,PTR_QV)=1.0
          IF( anx(i,j,k,PTR_QV) < 0.0 ) anx(i,j,k,PTR_QV)=0.0
          anx(i,j,k,PTR_QV)= anx(i,j,k,PTR_QV)*qv(i,j,k)
        END DO
      END DO
    END DO
  ENDIF
!
!-----------------------------------------------------------------------
!
!	Write to ARPSControl status file
!
!-----------------------------------------------------------------------
! Only u (index 1) is checked. Eventually should check all
! Only 1 level of UA data is checked.
!
! 10   Good
! 100  Good
! 200  Good - superob
! -9   Outside domain
! -99  Initializing value
! -199 Rejected in quality test
!
!  > 0    : used++ : in file, in domain, good point
!  = -199 : rej++  : in file, in domain, bad point
!  = -9   :        : in file, outside domain
!  = 0    :        : not in file
!
! Good data has isrc* > 0 and qual* > 0.

  IF (myproc == 0) THEN
    WRITE (status_file, '(2a)') runname(1:lfnkey),'.adasstat'
    WRITE (6,'(2a)') 'Creating status file: ', TRIM(status_file)
    OPEN (UNIT=9,FILE=TRIM(status_file),form='formatted',              &
          status='unknown')

    WRITE (stnlist_file, '(2a)') runname(1:lfnkey),'.adasstn'
    WRITE (6,'(2a)') 'Creating station list file: ', TRIM(stnlist_file)
    OPEN (UNIT=11,FILE=TRIM(stnlist_file),form='formatted',            &
          status='unknown')
  END IF
  !CALL mpudatec(stnlist_file,256)
!
! MPI is slightly more complicated.  Values for stations not owned by the
! processor may say valid when then aren't.
!
  n_used_sng = 0
  DO jsta=1,nobsng
    IF (mp_opt > 0) THEN
      IF (indexsng(jsta) .NE. myproc ) qualsng(:,jsta) = -9
    END IF
    IF( isrcsng(jsta) > 0 ) THEN
      kntgood=0
      DO ivar=1,nvar_anx
        IF(qualsng(ivar,jsta) > 0) kntgood = kntgood + 1
      END DO
      IF ( kntgood > 0) THEN
        n_used_sng(isrcsng(jsta)) = n_used_sng(isrcsng(jsta)) + 1
        jsrc=isrcsng(jsta)
        itype=(10*isrcsng(jsta))+1
        IF (myproc == 0) THEN
          WRITE(11,'(a5,1x,a4,a8,a8,9x,f8.4,f9.4,i5,i5)')               &
             stnsng(jsta,1),' -- ',srcsng(jsrc),' Station',             &
             latsng(jsta,1),lonsng(jsta,1),NINT(hgtsng(jsta,1)),itype
        ELSE
          !WRITE(*,'(1x,3a,I5,a)')                                       &
          ! '*** Should have been reported to station list file: ',      &
          ! TRIM(stnlist_file),' by process ',myproc,' ***'
          !WRITE(*,'(1x,a5,1x,a4,a8,a8,9x,f8.4,f9.4,i5,i5)')             &
          !   stnsng(jsta,1),' -- ',srcsng(jsrc),' Station',             &
          !   latsng(jsta,1),lonsng(jsta,1),NINT(hgtsng(jsta,1)),itype
          !WRITE(*,'(1x,a,I5,a)')  '***  end of report from process ',myproc,' ***'
        END IF
      END IF
    END IF
  END DO
!
! Modified code to report counts by source.  KB 2011/10/27
! This way is more flexible for the increased number of sources
! now being used and avoids issues with counting in SAO and mesonet
! categories.
!
  DO jsrc=1,nsrc_sng
    IF (srcsng(jsrc) == 'NULL') CYCLE

    IF (mp_opt > 0) THEN
      CALL mpsumi(n_used_sng(jsrc),1)
    END IF

    IF (myproc == 0) THEN
      CALL cvtsrcstr(srcsng(jsrc),srcstr)
      WRITE (9,'(3a,i6,a)') '$',TRIM(srcstr),                          &
            '_n_used = ',n_used_sng(jsrc),';'
    END IF
  END DO
!
! Report upper air counts.
!
  DO jsta=1,nobsua
    IF (srcua(jsrc) == 'NULL') CYCLE

    jsrc = isrcua(jsta)
    IF (jsrc < 1 .OR. jsrc > nsrc_ua) THEN
      IF (myproc == 0) THEN
        WRITE (9,*) '# Unable to process: irsrcua=',                   &
         jsrc, ' ', TRIM(stnua(jsta))
      END IF
      CYCLE
    END IF

    n_used_ua(isrcua(jsta)) = n_used_ua(isrcua(jsta)) + 1
    IF (myproc == 0) THEN

      tmp = TRIM(srcua(jsrc))
      lentmp=LEN_TRIM(tmp)
      DO i=1,lentmp
        IF(tmp(i:i) == " ") THEN
          tmp(i:i) = "_"
        END IF
      END DO

!
! "tmp" values can *MUST* start with a letter.  Therefore "88D" is illegal.
! Put an "R" in front of the name.
!

      IF (tmp(1:3) == "88D") THEN
        tmp = "R"//tmp
      ENDIF

      WRITE (9,'(5a,i6,a)') '${',TRIM(tmp),'}{stn_levels}{"',      &
              TRIM(stnua(jsta)),'"} = ',nlevsua(jsta),';'
      IF (INDEX(srcua(jsrc),'RAOB') > 0 .OR.                           &
          INDEX(srcua(jsrc),'CLASS') > 0) THEN
        itype=601
        CALL xytoll(1,1,xua(jsta),yua(jsta),latua,lonua)
        WRITE(11,'(a5,1x,a4,a8,a8,9x,f8.4,f9.4,i5,i5)')                &
            stnua(jsta),' -- ',srcua(jsrc),' Station',                 &
            latua,lonua,NINT(hgtua(1,jsta)),itype
      ELSE IF (INDEX(srcua(jsrc),'PRO') > 0 .OR.                       &
               INDEX(srcua(jsrc),'PRF') > 0 .OR.                       &
               INDEX(srcua(jsrc),'SODAR') > 0 ) THEN
        itype=701
        CALL xytoll(1,1,xua(jsta),yua(jsta),latua,lonua)
        WRITE(11,'(a5,1x,a4,a8,a8,9x,f8.4,f9.4,i5,i5)')                &
            stnua(jsta),' -- ',srcua(jsrc),' Station',                 &
            latua,lonua,NINT(hgtua(1,jsta)),itype
      ELSE IF (INDEX(srcua(jsrc),'TOWR') > 0) THEN
        itype=801
        CALL xytoll(1,1,xua(jsta),yua(jsta),latua,lonua)
        WRITE(11,'(a5,1x,a4,a8,a8,9x,f8.4,f9.4,i5,i5)')                &
            stnua(jsta),' -- ',srcua(jsrc),' Station',                 &
            latua,lonua,NINT(hgtua(1,jsta)),itype
      ELSE IF (INDEX(srcua(jsrc),'VAD') > 0) THEN
        itype=901
        CALL xytoll(1,1,xua(jsta),yua(jsta),latua,lonua)
        WRITE(11,'(a5,1x,a4,a8,a8,9x,f8.4,f9.4,i5,i5)')                &
            stnua(jsta),' -- ',srcua(jsrc),' Station',                 &
            latua,lonua,NINT(hgtua(1,jsta)),itype
      ELSE
        itype=911
        CALL xytoll(1,1,xua(jsta),yua(jsta),latua,lonua)
        WRITE(11,'(a5,1x,a4,a8,a8,9x,f8.4,f9.4,i5,i5)')                &
            stnua(jsta),' -- ',srcua(jsrc),' Station',                 &
            latua,lonua,NINT(hgtua(1,jsta)),itype
      END IF
    END IF
  END DO

  IF( myproc == 0 ) THEN
    DO jsrc=1,nsrc_ua
      IF (srcua(jsrc) == 'NULL') CYCLE
      CALL cvtsrcstr(srcua(jsrc),srcstr)
      WRITE (9,'(3a,i6,a)') '$',TRIM(srcstr),'_n_used = ',             &
             n_used_ua(jsrc),';'
    END DO
    CLOSE (9)
    CLOSE (11)
  END IF

  CALL FLUSH (6)
!
!-----------------------------------------------------------------------
!
!  Use reflectivity to establish cloud water
!
!-----------------------------------------------------------------------
!
  IF( radcldopt > 0 ) CALL radmcro(nx,ny,nz,ref_mos_3d,                 &
                 radqvopt,radqcopt,radqropt,radptopt,                   &
                 refsat,rhradobs,                                       &
                 refcld,cldrad,ceilopt,ceilmin,dzfill,                  &
                 refrain,radsetrat,radreflim,radptgain,                 &
                 xs,ys,zs,zp,j3,                                        &
                 anx(1,1,1,PTR_P),anx(1,1,1,PTR_PT),anx(1,1,1,PTR_QV),  &
                 ptbar,qvbar,rhostr,                                    &
                 qscalar,tem1)
!
!-----------------------------------------------------------------------
!
!  Transfer analysed fields to ARPS variables for writing.
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        pprt(i,j,k) =anx(i,j,k,PTR_P )-pbar(i,j,k)
        ptprt(i,j,k)=anx(i,j,k,PTR_PT)-ptbar(i,j,k)
        tgrid =anx(i,j,k,PTR_PT) * (anx(i,j,k,PTR_P)/p0) ** rddcp
        qvsat=f_qvsat( anx(i,j,k,PTR_P),tgrid)
        qvmin=rhmin*qvsat
        qv(i,j,k)=MAX(qvmin,MIN(anx(i,j,k,PTR_QV),qvsat))
      END DO
    END DO
  END DO

  IF (mp_opt > 0) THEN
    CALL mpsendrecv2dew(pprt,nx,ny,nz,ebc,wbc,0,tem1)
    CALL mpsendrecv2dns(pprt,nx,ny,nz,nbc,sbc,0,tem1)

    CALL mpsendrecv2dew(ptprt,nx,ny,nz,ebc,wbc,0,tem1)
    CALL mpsendrecv2dns(ptprt,nx,ny,nz,nbc,sbc,0,tem1)

    CALL mpsendrecv2dew(qv,nx,ny,nz,ebc,wbc,0,tem1)
    CALL mpsendrecv2dns(qv,nx,ny,nz,nbc,sbc,0,tem1)
  END IF

  IF(spradopt == 1) THEN

!-----------------------------------------------------------------------
!
!  Check for superadiabatic layers in each column.
!  For spradopt=1,
!  Where superadiabatic levels are found, make the next level down
!  cooler, applying the minimum potential temperature lapse rate,
!  ptlapse.
!
!-----------------------------------------------------------------------
!
    DO j=1,ny-1
      DO i=1,nx-1
        ptprt(i,j,nz-1)=anx(i,j,nz-1,PTR_PT)-ptbar(i,j,nz-1)
        DO k=nz-2,2,-1
          anx(i,j,k,PTR_PT)=MIN(anx(i,j,k,PTR_PT),                      &
                (anx(i,j,k+1,PTR_PT)-ptlapse*(zs(i,j,k+1)-zs(i,j,k))))
          ptprt(i,j,k)=anx(i,j,k,PTR_PT)-ptbar(i,j,k)
        END DO
      END DO
    END DO

  ELSE IF(spradopt == 2) THEN
!
!-----------------------------------------------------------------------
!
!  For spradopt=2,
!  Where superadiabatic levels are found, make the next level up
!  warmer, applying the minimum potential temperature lapse rate,
!  ptlapse.
!
!-----------------------------------------------------------------------
!
    DO j=1,ny-1
      DO i=1,nx-1
        ptprt(i,j,2)=anx(i,j,2,PTR_PT)-ptbar(i,j,2)
        DO k=3,nz-1
          anx(i,j,k,PTR_PT)=MAX(anx(i,j,k,PTR_PT),                      &
                (anx(i,j,k-1,PTR_PT)+ptlapse*(zs(i,j,k)-zs(i,j,k-1))))
          ptprt(i,j,k)=anx(i,j,k,PTR_PT)-ptbar(i,j,k)
        END DO
      END DO
    END DO

  END IF

!-----------------------------------------------------------------------
!
! Transfer analysis wind fields
!
!-----------------------------------------------------------------------

  DO k=1,nz-1
    DO j=1,ny-1
      DO i=2,nx-1
        u(i,j,k)=0.5*(anx(i,j,k,PTR_U)+anx(i-1,j,k,PTR_U))
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO j=1,ny-1
      u(1,j,k)  = u(2,j,k)
      u(nx,j,k) = u(nx-1,j,k)
    END DO
  END DO

  IF (mp_opt > 0) THEN
    CALL mpsendrecv2dew(u,nx,ny,nz,ebc,wbc,1,tem1)
  END IF

  DO k=1,nz-1
    DO j=2,ny-1
      DO i=1,nx-1
        v(i,j,k)=0.5*(anx(i,j,k,PTR_V)+anx(i,j-1,k,PTR_V))
      END DO
    END DO
  END DO

  DO k=1,nz-1
    DO i=1,nx-1
      v(i,1,k)  = v(i,2,k)
      v(i,ny,k) = v(i,ny-1,k)
    END DO
  END DO

  IF (mp_opt > 0) THEN
    CALL mpsendrecv2dns(v,nx,ny,nz,nbc,sbc,2,tem1)
  END IF

#ifdef arps3dvar

   DO k=2,nz
     DO j=1,ny-1
       DO i=1,nx-1
         w(i,j,k)=0.5*(anx(i,j,k,PTR_W)+anx(i,j,k-1,PTR_W))
       END DO
     END DO
   END DO

   DO j = 1, ny-1
     DO i = 1, nx-1
       w(i,j,1)  = w(i,j,2)
       w(i,j,nz) = w(i,j,nz-1)
     END DO
   END DO

#else

!-----------------------------------------------------------------------
!
!  Adjust wind fields to minimize 3-d divergence
!  This is the first run-through for "w" in order to obtain
!  a good estimate of w for the cloud analysis.  Note here wndadj
!  is set to "2" for integration of 2-d divergence for "w" and
!  the cloud_w opt is set to 0.
!
!-----------------------------------------------------------------------
!
  !
  ! w_cld isn't allocated yet!
  !
  ALLOCATE(w_cld(nx,ny,nz),STAT=istatus)
  CALL check_alloc_status(istatus, "adas:w_cld")
  CALL adjuvw( nx,ny,nz,u,v,w,                                          &
               wcont,ptprt,ptbar,                                       &
               zp,j1,j2,j3,aj3z,mapfct,rhostr,w_cld,                    &
               2,obropt,obrzero,0,                                      &
               tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8)
  DEALLOCATE(w_cld)
!
!-----------------------------------------------------------------------
!
!  Enforce vertical w boundary conditions
!
!-----------------------------------------------------------------------
!
  CALL rhouvw(nx,ny,nz,rhostr,tem1,tem2,tem3)
  CALL vbcw(nx,ny,nz,w,wcont,tbc,bbc,u,v,                               &
              rhostr,tem1,tem2,tem3,                                    &
              j1,j2,j3)
#endif

!
!-----------------------------------------------------------------------
!
!  Free-up some memory.
!
!-----------------------------------------------------------------------
!
  DEALLOCATE(anx)
  IF (ALLOCATED(irad))    DEALLOCATE(irad)
  IF (ALLOCATED(nlevrad)) DEALLOCATE(nlevrad)
  IF (ALLOCATED(latradc)) DEALLOCATE(latradc,lonradc,xradc,yradc,trnradc,distrad)
  IF (ALLOCATED(uazmrad)) DEALLOCATE(uazmrad,vazmrad,dsdr,dhdr)
  IF (ALLOCATED(hgtradc)) DEALLOCATE(hgtradc,theradc,obsrad,odifrad,oanxrad)
  IF (ALLOCATED(corrad))  DEALLOCATE(corrad)
!
!-----------------------------------------------------------------------
!
!  Call cloud analysis.
!
!-----------------------------------------------------------------------
!
  IF( cloudopt > 0 ) THEN

    ALLOCATE(w_cld(nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:w_cld")
    w_cld=0.
    ALLOCATE(cldpcp_type_3d(nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:cldpcp_type_3d")
    cldpcp_type_3d=0
    ALLOCATE(icing_index_3d(nx,ny,nz),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:icing_index_3d")
    icing_index_3d=0
    ALLOCATE(l_mask_pcptype(nx,ny),STAT=istatus)
    CALL check_alloc_status(istatus, "adas:l_mask_pcptype")
    l_mask_pcptype=.false.

    CALL cmpclddrv(nx,ny,nz,i4timanx,                                   &
              xs,ys,zs,j3,hterain,latgr,longr,                          &
              nobsng,indexsng,stnsng,isrcsng,csrcsng,xsng,ysng,         &
              nxlg,nylg,xslg,yslg,                                      &
              timesng,latsng,lonsng,hgtsng,                             &
              kloud,store_amt,store_hgt,                                &
              stnrad,isrcrad,latrad,lonrad,elvrad,ixrad,jyrad,          &
              pprt,ptprt,qv,qscalar,w,                                  &
              pbar,ptbar,qvbar,rhostr,                                  &
              istat_radar,ref_mos_3d,w_cld,cldpcp_type_3d,              &
              icing_index_3d,l_mask_pcptype,                            &
              tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,tem9,tem10,tem11)

    DEALLOCATE (ref_mos_3d, cldpcp_type_3d, icing_index_3d, l_mask_pcptype )
    ! w_cld should also be deallocated here?
  END IF

  IF(retqropt > 0 .AND. ncolret > 0) CALL retmcro(nx,ny,nz,             &
                 mx_ret,nsrc_ret,nvar_anx,nz_ret,mx_colret,             &
                 srcret,isrcret,iret,nlevret,                           &
                 xretc,yretc,hgtretc,qrret,ncolret,                     &
                 dzfill,                                                &
                 xs,ys,zs,qscalar(:,:,:,P_QR))  ! This will be problematic
                 ! Since P_QR may be not defined, but it is safe now because
                 ! we do not do retrieval analysis at present.
!
!-----------------------------------------------------------------------
!
!  Create rhobar from rhostr
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem6(i,j,k)=rhostr(i,j,k)/j3(i,j,k)
      END DO
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Calculate and store the sound wave speed squared in csndsq.
!  Use original definition of sound speed.
!
!-----------------------------------------------------------------------
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        csndsq(i,j,k)= cpdcv*pbar(i,j,k)/tem6(i,j,k)
      END DO
    END DO
  END DO

  IF(hydradj == 1 .OR. hydradj == 2) THEN
    pconst=0.5*g*dz
!
!-----------------------------------------------------------------------
!
!  Create thermal buoyancy at each scalar point,
!  which is stored in tem2
!
!-----------------------------------------------------------------------
!
    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          qvprt=qv(i,j,k)-qvbar(i,j,k)
          qtot = 0.0
          DO nq=1,nscalarq
            qtot=qtot+qscalar(i,j,k,nq)
          END DO
          tem2(i,j,k)=j3(i,j,k)*tem6(i,j,k)*                            &
                    g*( (ptprt(i,j,k)/ptbar(i,j,k)) +                   &
                        (qvprt/(rddrv+qvbar(i,j,k))) -                  &
                        ((qvprt+qtot)/(1.+qvbar(i,j,k))) )
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Average thermal buoyancy to w points
!  which is stored in tem3
!
!-----------------------------------------------------------------------
!
    CALL avgsw(tem2,nx,ny,nz,1,nx,1,ny, tem3)

    IF(hydradj == 1) THEN

      DO i=1,nx
        DO j=1,ny
          tem1(i,j,1)=pprt(i,j,1)
          DO k=2,nz-2
            tem1(i,j,k)=                                                &
                  ( (1. - (pconst*j3(i,j,k-1)/csndsq(i,j,k-1)) )*       &
                      tem1(i,j,k-1) +                                   &
                      dz*tem3(i,j,k) ) /                                &
                      (1. + (pconst*j3(i,j,k)/csndsq(i,j,k)) )
            IF(j == 26 .AND. MOD(i,10) == 0) THEN
              IF(k == nz-2) PRINT *, '            Point i= ',i, '  j=26'
              PRINT *, ' k,pprt,tem1,diff =',                           &
                  k,pprt(i,j,k),tem1(i,j,k),                            &
                  (tem1(i,j,k)-pprt(i,j,k))
            END IF
            pprt(i,j,k)=tem1(i,j,k)
          END DO
          pprt(i,j,nz-1)=pprt(i,j,nz-2)
        END DO
      END DO

    ELSE IF(hydradj == 2) THEN

      DO i=1,nx
        DO j=1,ny
          tem1(i,j,nz-1)=pprt(i,j,nz-1)
          DO k=nz-2,2,-1
            tem1(i,j,k)=                                                &
                  ( (1.+ (pconst*j3(i,j,k+1)/csndsq(i,j,k+1)) )*        &
                      tem1(i,j,k+1) -                                   &
                      dz*tem3(i,j,k+1) ) /                              &
                      (1.- (pconst*j3(i,j,k  )/csndsq(i,j,k  )) )
            IF(j == 26 .AND. MOD(i,10) == 0) THEN
              IF(k == nz-2) PRINT *, '            Point i= ',i, '  j=26'
              PRINT *, ' k,pprt,tem1,diff =',                           &
                  k,pprt(i,j,k),tem1(i,j,k),                            &
                  (tem1(i,j,k)-pprt(i,j,k))
            END IF
            pprt(i,j,k)=tem1(i,j,k)
          END DO
          pprt(i,j,1)=pprt(i,j,2)
        END DO
      END DO
    END IF
  ELSE IF (hydradj == 3) THEN
!
!-----------------------------------------------------------------------
!
!  Calculate total pressure, store in tem1.
!  Calculate virtual temperature, store in tem2.
!
!-----------------------------------------------------------------------
!
    DO k=1,nz
      DO j=1,ny
        DO i=1,nx
          tem1(i,j,k) = pbar(i,j,k)+pprt(i,j,k)
          temp = (ptbar(i,j,k)+ptprt(i,j,k))*                           &
                 ((tem1(i,j,k)/p0)**rddcp)
          tem2(i,j,k) = temp*(1.0+rvdrd*qv(i,j,k))/                     &
                             (1.0+qv(i,j,k))
        END DO
      END DO
    END DO
!
!-----------------------------------------------------------------------
!
!  Integrate hydrostatic equation, begining at k=2
!
!-----------------------------------------------------------------------
!
    pconst=-g/rd
    DO k=3,nz-1
      DO j=1,ny-1
        DO i=1,nx-1
          tvbar=0.5*(tem2(i,j,k)+tem2(i,j,k-1))
          tem1(i,j,k)=tem1(i,j,k-1)*                                    &
                     EXP(pconst*(zs(i,j,k)-zs(i,j,k-1))/tvbar)
          pprt(i,j,k)=tem1(i,j,k)-pbar(i,j,k)
        END DO
      END DO
    END DO
    DO j=1,ny-1
      DO i=1,nx-1
        tvbar=tem2(i,j,2)
        tem1(i,j,1)=tem1(i,j,2)*                                        &
                   EXP(pconst*(zs(i,j,1)-zs(i,j,2))/tvbar)
        pprt(i,j,1)=tem1(i,j,1)-pbar(i,j,1)
      END DO
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Adjust wind fields to minimize 3-d divergence
!
!-----------------------------------------------------------------------
!
#ifdef adas
  CALL adjuvw( nx,ny,nz,u,v,w,                                          &
               wcont,ptprt,ptbar,                                       &
               zp,j1,j2,j3,aj3z,mapfct,rhostr,w_cld,                    &
               wndadj,obropt,obrzero,cldwopt,                           &
               tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8)
 ! the use of w_cld is very suspicious?, it may not be allocated.

!
!-----------------------------------------------------------------------
!
!  Enforce vertical w boundary conditions
!
!-----------------------------------------------------------------------
!
  CALL rhouvw(nx,ny,nz,rhostr,tem1,tem2,tem3)
  CALL vbcw(nx,ny,nz,w,wcont,tbc,bbc,u,v,                               &
              rhostr,tem1,tem2,tem3,                                    &
              j1,j2,j3)
#endif
!
!-----------------------------------------------------------------------
!
!  Enforce boundary conditions
!  Use zero gradient for lateral bc, and rigid for top,bottom.
!
!-----------------------------------------------------------------------
!
  ebc = 3
  wbc = 3
  nbc = 3
  sbc = 3
  tbc = 1
  bbc = 1
!
  DO k=1,nz
    DO j=1,ny
      pdteb(j,k)=0.
      pdtwb(j,k)=0.
    END DO
    DO i=1,nx
      pdtnb(i,k)=0.
      pdtsb(i,k)=0.
    END DO
  END DO
!
  CALL bcu(nx,ny,nz,0., u,                                              &
           pdteb,pdtwb,pdtnb,pdtsb,ebc,wbc,nbc,sbc,tbc,bbc,             &
           ebc,wbc,nbc,sbc)  ! should be global value. Use local value instead,
                             ! since ADAS is still not parallized.

  CALL bcv(nx,ny,nz,0., v,                                              &
           pdteb,pdtwb,pdtnb,pdtsb,ebc,wbc,nbc,sbc,tbc,bbc,             &
           ebc,wbc,nbc,sbc)  ! should be global value.

  CALL lbcw(nx,ny,nz,0.0,w,wcont,                                       &
            pdteb,pdtwb,pdtnb,pdtsb,ebc,wbc,nbc,sbc,                    &
            ebc,wbc,nbc,sbc)  ! should be global value.

  CALL bcsclr(nx,ny,nz,0.,ptprt,ptprt,ptprt,                            &
              pdteb,pdtwb,pdtnb,pdtsb,ebc,wbc,nbc,sbc,tbc,bbc,          &
              ebc,wbc,nbc,sbc)  ! should be global value.

  CALL bcsclr(nx,ny,nz,0., qv,qv,qv,                                    &
           pdteb,pdtwb,pdtnb,pdtsb,ebc,wbc,nbc,sbc,tbc,bbc,             &
           ebc,wbc,nbc,sbc)  ! should be global value.

  DO nq=1,nscalar
    CALL bcsclr(nx,ny,nz,0., qscalar(:,:,:,nq),qscalar(:,:,:,nq),       &
             qscalar(:,:,:,nq),                                         &
             pdteb,pdtwb,pdtnb,pdtsb,ebc,wbc,nbc,sbc,tbc,bbc,           &
             ebc,wbc,nbc,sbc)  ! should be global value.
  END DO
!
!-----------------------------------------------------------------------
!
!  Apply extrapolated gradient bc to pressure for geostrophic
!  type balance to winds along lateral boundaries.  Consistent
!  with horizontal winds constant.
!
!-----------------------------------------------------------------------
!
  DO k=2,nz-2
    DO j=1,ny-1
      pprt(1,j,k)=2.*pprt(2,j,k)-pprt(3,j,k)
    END DO
  END DO

  DO k=2,nz-2
    DO j=1,ny-1
      pprt(nx-1,j,k)=2.*pprt(nx-2,j,k)-pprt(nx-3,j,k)
    END DO
  END DO

  DO k=2,nz-2
    DO i=1,nx-1
      pprt(i,ny-1,k)=2.*pprt(i,ny-2,k)-pprt(i,ny-3,k)
    END DO
  END DO

  DO k=2,nz-2
    DO i=1,nx-1
      pprt(i,1,k)=2.*pprt(i,2,k)-pprt(i,3,k)
    END DO
  END DO

!
!-----------------------------------------------------------------------
!
!  For top and bottom bc, apply hydrostatic pressure equation to
!  total pressure, then subtract pbar.
!
!-----------------------------------------------------------------------
!
  pconst=-g/rd
  DO j=1,ny-1
    DO i=1,nx-1
      pr2 = pbar(i,j,2)+pprt(i,j,2)
      temp = (ptbar(i,j,2)+ptprt(i,j,2))*                               &
             ((pr2/p0)**rddcp)
      tvbar = temp*(1.0+rvdrd*qv(i,j,2))/                               &
                         (1.0+qv(i,j,2))
      pr1=pr2*EXP(pconst*(zs(i,j,1)-zs(i,j,2))/tvbar)
      pprt(i,j,1)=pr1-pbar(i,j,1)

      pr2 = pbar(i,j,nz-2)+pprt(i,j,nz-2)
      temp = (ptbar(i,j,nz-2)+ptprt(i,j,nz-2))*                         &
             ((pr2/p0)**rddcp)
      tvbar= temp*(1.0+rvdrd*qv(i,j,nz-2))/                             &
                         (1.0+qv(i,j,nz-2))
      pr1=pr2*EXP(pconst*(zs(i,j,nz-1)-zs(i,j,nz-2))/tvbar)
      pprt(i,j,nz-1)=pr1-pbar(i,j,nz-1)
    END DO
  END DO
!
!-----------------------------------------------------------------------
!
!  Adjust surface (skin) temperature based on radiation and
!  current air temperature.
!
!-----------------------------------------------------------------------
!
  IF(tsfcopt > 0 ) THEN

    IF (radopt == 2) THEN
      rbufsz = n2d_radiat*nx*ny + n3d_radiat*nx*ny*nz
    ELSE
      rbufsz = 1
    END IF

    ALLOCATE(rsirbm(nx,ny))
    rsirbm = 0.
    ALLOCATE(rsirdf(nx,ny))
    rsirdf = 0.
    ALLOCATE(rsuvbm(nx,ny))
    rsuvbm = 0.
    ALLOCATE(rsuvdf(nx,ny))
    rsuvdf = 0.

    ALLOCATE(cosz(nx,ny))
    cosz = 0.
    ALLOCATE(cosss(nx,ny))
    cosss = 0.


    ALLOCATE(fdirir(nx,ny))
    fdirir = 0.
    ALLOCATE(fdifir(nx,ny))
    fdifir = 0.
    ALLOCATE(fdirpar(nx,ny))
    fdirpar= 0.
    ALLOCATE(fdifpar(nx,ny))
    fdifpar= 0.

    ALLOCATE(temrad1(nx,ny,nz))
    temrad1 = 0.
    ALLOCATE(tem12(nx,ny,nz))
    tem12 = 0.
    ALLOCATE(tem13(nx,ny,nz))
    tem13 = 0.
    ALLOCATE(tem14(nx,ny,nz))
    tem14 = 0.
    ALLOCATE(tem15(nx,ny,nz))
    tem15 = 0.
    ALLOCATE(tem16(nx,ny,nz))
    tem16 = 0.
    ALLOCATE(tem17(nx,ny,nz))
    tem17 = 0
    ALLOCATE(radbuf(rbufsz))
    radbuf = 0.
    ALLOCATE(sh(nx,ny))
    sh = 0.

    CALL adjtsfc(nx,ny,nz,nscalar,rbufsz,                               &
               ptprt,pprt,qv,qscalar,                                   &
               ptbar,pbar,ppi,rhostr,                                   &
               x,y,z,zp,j3inv,                                          &
               soiltyp, tsoil(1,1,1,0),qsoil(1,1,1,0),snowdpth,         &
               radfrc, radsw, rnflx, radswnet, radlwin,                 &
               rsirbm,rsirdf,rsuvbm,rsuvdf,                             &
               cosz, cosss,                                             &
               fdirir,fdifir,fdirpar,fdifpar,                           &
               tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,                 &
               tem9,temrad1,tem11,tem12,tem13,tem14,tem15,tem16,        &
               radbuf, sh, tem17)

    ! DTD: placed call to fix_soil_nstp here, to propagate adjusted soil temperatures above
    ! to the other soil types, if nstyp > 1
    ! The call to adtsfc above only adjusts the soil temperatures for one soil type

    CALL fix_soil_nstyp(nx,ny,nzsoil,1,nstyp,tsoil, qsoil,wetcanp)

    DEALLOCATE(rsirbm)
    DEALLOCATE(rsirdf)
    DEALLOCATE(rsuvbm)
    DEALLOCATE(rsuvdf)

    DEALLOCATE(cosz)
    DEALLOCATE(cosss)

    DEALLOCATE(fdirir)
    DEALLOCATE(fdifir)
    DEALLOCATE(fdirpar)
    DEALLOCATE(fdifpar)

    DEALLOCATE(temrad1)
    DEALLOCATE(tem12)
    DEALLOCATE(tem13)
    DEALLOCATE(tem14)
    DEALLOCATE(tem15)
    DEALLOCATE(tem16)
    DEALLOCATE(tem17)
    DEALLOCATE(radbuf)
    DEALLOCATE(sh)

  END IF

!
!-----------------------------------------------------------------------
!
!  Calculate analysis increments and write to file, if desired.
!
!-----------------------------------------------------------------------
!
  IF(incrdmp > 0) THEN
    CALL incrcalc(nx,ny,nz,nxndg,nyndg,nzndg,                           &
                  u,v,w,pprt,ptprt,qv,                                  &
                  qscalar,                                              &
                  uincr,vincr,wincr,pincr,ptincr,qvincr,                &
                  qcincr,qrincr,qiincr,qsincr,qhincr,                   &
                  istatus)

    CALL incrdump(nxndg,nyndg,nzndg,incrdmp,incrhdfcompr,incdmpf,       &
                  uincr,vincr,wincr,pincr,ptincr,qvincr,                &
                  qcincr,qrincr,qiincr,qsincr,qhincr,                   &
                  uincdmp,vincdmp,wincdmp,                              &
                  pincdmp,ptincdmp,qvincdmp,                            &
                  qcincdmp,qrincdmp,qiincdmp,qsincdmp,qhincdmp,         &
                  istatus)

  END IF
!
!-----------------------------------------------------------------------
!
!  Data dump of the model grid and base state arrays:
!
!-----------------------------------------------------------------------
!
  IF (mp_opt > 0) THEN
    CALL mpsendrecv2dew(u,nx,ny,nz,ebc,wbc,1,tem1)
    CALL mpsendrecv2dns(u,nx,ny,nz,nbc,sbc,1,tem1)

    CALL mpsendrecv2dew(v,nx,ny,nz,ebc,wbc,2,tem1)
    CALL mpsendrecv2dns(v,nx,ny,nz,nbc,sbc,2,tem1)

    CALL mpsendrecv2dew(w,nx,ny,nz,ebc,wbc,3,tem1)
    CALL mpsendrecv2dns(w,nx,ny,nz,nbc,sbc,3,tem1)

    CALL mpsendrecv2dew(ptprt,nx,ny,nz,ebc,wbc,0,tem1)
    CALL mpsendrecv2dns(ptprt,nx,ny,nz,nbc,sbc,0,tem1)

    CALL mpsendrecv2dew(pprt,nx,ny,nz,ebc,wbc,0,tem1)
    CALL mpsendrecv2dns(pprt,nx,ny,nz,nbc,sbc,0,tem1)

    CALL mpsendrecv2dew(qv,nx,ny,nz,ebc,wbc,0,tem1)
    CALL mpsendrecv2dns(qv,nx,ny,nz,nbc,sbc,0,tem1)

    DO nq=1,nscalarq
      CALL mpsendrecv2dew(qscalar(:,:,:,nq),nx,ny,nz,ebc,wbc,0,tem1)
      CALL mpsendrecv2dns(qscalar(:,:,:,nq),nx,ny,nz,nbc,sbc,0,tem1)
    END DO

  END IF
!
!-----------------------------------------------------------------------
!
  DO k=1,nz
    DO j=1,ny
      DO i=1,nx
        tem1(i,j,k)=0.
      END DO
    END DO
  END DO
!
  DO k=1,nz-1
    DO j=1,ny-1
      DO i=1,nx-1
        tem2(i,j,k)=rhostr(i,j,k)/j3(i,j,k)
      END DO
    END DO
  END DO
!
  IF( hdmpfmt == 5 ) GO TO 700
  IF( hdmpfmt == 9 ) GO TO 700
!
!-----------------------------------------------------------------------
!
!  Find a unique name basdmpfn(1:lbasdmpf) for the grid and
!  base state array dump file
!
!-----------------------------------------------------------------------
!

  CALL gtbasfn(runname(1:lfnkey),dirname,ldirnam,hdmpfmt,               &
               mgrid,nestgrd,basdmpfn,lbasdmpf)

  IF (myproc == 0) WRITE(6,'(1x,a,a)')                                  &
       'Data dump of grid and base state arrays into file ',            &
        basdmpfn(1:lbasdmpf)

  CALL setcornerll( nx,ny,x,y )   ! set the corner lat/lon

  grdbas = 1                ! Dump out grid and base state arrays only

!    blocking inserted for ordering i/o for message passing
  DO i=0,nprocs-1,dumpstride
!  DO i=0,0,max_fopen

    IF(myproc >= i.AND.myproc <= i+dumpstride-1)THEN

      CALL dtadump(nx,ny,nz,nzsoil,nstyps,                              &
               hdmpfmt,nchdmp,basdmpfn(1:lbasdmpf),                     &
               grdbas,filcmprs,                                         &
               u,v,w,ptprt,pprt,qv,                                     &
               qscalar,                                                 &
               tem1,tem1,tem1,                                          &
               ubar,vbar,tem1,ptbar,pbar,tem2,qvbar,                    &
               x,y,z,zp,zpsoil,                                         &
               soiltyp,stypfrct,vegtyp,lai,roufns,veg,                  &
               tsoil,qsoil,wetcanp,snowdpth,                            &
               raing,rainc,prcrate,                                     &
               radfrc,radsw,rnflx,radswnet,radlwin,                     &
               usflx,vsflx,ptsflx,qvsflx,                               &
               tem3,tem4,tem5)
     END IF
     IF (mp_opt > 0) CALL mpbarrier
   END DO

!
!-----------------------------------------------------------------------
!
!  Find a unique name hdmpfn(1:ldmpf) for history dump data set
!  at time 'curtim'.
!
!-----------------------------------------------------------------------
!
  700 CONTINUE

  CALL gtdmpfn(runname(1:lfnkey),dirname,                               &
               ldirnam,curtim,hdmpfmt,                                  &
               mgrid,nestgrd, hdmpfn, ldmpf)

  IF (myproc == 0)                                                      &
    WRITE(6,'(1x,a,a)') 'History data dump in file ',hdmpfn(1:ldmpf)

  grdbas = 0                ! No base state or grid array is dumped.

!    blocking inserted for ordering i/o for message passing
  DO i=0,nprocs-1,dumpstride
!  DO i=0,0,max_fopen

    IF(myproc >= i.AND.myproc <= i+dumpstride-1)THEN

      CALL dtadump(nx,ny,nz,nzsoil,nstyps,                              &
               hdmpfmt,nchdmp,hdmpfn(1:ldmpf),                          &
               grdbas,filcmprs,                                         &
               u,v,w,ptprt,pprt,qv,                                     &
               qscalar,                                                 &
               tem1,tem1,tem1,                                          &
               ubar,vbar,tem1,ptbar,pbar,tem2,qvbar,                    &
               x,y,z,zp,zpsoil,                                         &
               soiltyp,stypfrct,vegtyp,lai,roufns,veg,                  &
               tsoil,qsoil,wetcanp,snowdpth,                            &
               raing,rainc,prcrate,                                     &
               radfrc,radsw,rnflx,radswnet,radlwin,                     &
               usflx,vsflx,ptsflx,qvsflx,                               &
               tem3,tem4,tem5)
    END IF
    IF (mp_opt > 0) CALL mpbarrier
  END DO

  CALL flush(6)
  IF (mp_opt > 0) THEN
    IF (myproc == 0) WRITE(6,'(/,5x,a)') '==== '//PROGNAME//'_MPI terminated normally ===='
    CALL mpexit(0)
  ELSE
    WRITE(6,'(/,5x,a)') '==== '//PROGNAME//' terminated normally ===='
    STOP
  END IF

END PROGRAM adas_3dvar

!
SUBROUTINE exsufx(fname,lenfnm,suffix,maxsuf,dotloc,lensuf)
  IMPLICIT NONE
  INTEGER :: lenfnm
  INTEGER :: maxsuf
  CHARACTER (LEN=lenfnm) :: fname
  CHARACTER (LEN=maxsuf) :: suffix
  INTEGER :: lensuf
  INTEGER :: dotloc
  LOGICAL, PARAMETER :: back = .TRUE.
!
  INTEGER :: i,j
!
  dotloc=0
  lensuf=0
  DO i=1,maxsuf
    suffix(i:i)=' '
  END DO
  dotloc=INDEX(fname,'.',back)
  IF(dotloc > 0) THEN
    lensuf=MIN(maxsuf,lenfnm-dotloc)
    suffix=fname((dotloc+1):lenfnm)
  END IF
  RETURN
END SUBROUTINE exsufx
!
SUBROUTINE exfroot(fname,lenfnm,froot,mxroot,lenroot)
  IMPLICIT NONE
  INTEGER :: lenfnm
  INTEGER :: mxroot
  CHARACTER (LEN=lenfnm) :: fname
  CHARACTER (LEN=mxroot) :: froot
  INTEGER :: lenroot

  INTEGER :: i
  INTEGER :: slashloc
  INTEGER :: dotloc
  LOGICAL, PARAMETER :: back = .TRUE.

  slashloc=INDEX(fname,'/',back)
  dotloc=INDEX(fname,'.',back)
  lenroot=(dotloc-slashloc)-1
  froot=fname((slashloc+1):(dotloc-1))
  RETURN
END SUBROUTINE exfroot
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  SUBROUTINE SETCAT                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE setcat(nx,ny,nz,nscalar,ccatopt,zs,                          &
                  ptprt,pprt,qv,qscalar,ptbar,pbar,prcrate,icatg)
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Assign categories to grid locations according to the background
!  state to account for decorrelation across areas with active
!  parameterized convection and those without.
!
!
!-----------------------------------------------------------------------
!
!
!  AUTHOR: Keith Brewster
!  8/7/98
!
!  MODIFICATION HISTORY:
!
!  2/2/06  Kevin Thomas
!  MPI the counts.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!  OUTPUT:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'mp.inc'
  INTEGER :: nx,ny,nz,nscalar
  INTEGER :: ccatopt
  REAL :: zs    (nx,ny,nz)  ! The physical height coordinate defined at
                            ! w-point of the staggered grid.
  REAL :: ptprt (nx,ny,nz)  ! Perturbation potential temperature (K)
  REAL :: pprt  (nx,ny,nz)  ! Perturbation pressure (Pascal)
  REAL :: qv    (nx,ny,nz)  ! Water vapor specific humidity (kg/kg)
  REAL :: qscalar(nx,ny,nz,nscalar)

  REAL :: ptbar (nx,ny,nz)     ! Base state potential temperature (K)
  REAL :: pbar  (nx,ny,nz)     ! Base state pressure (Pascal).
  REAL :: prcrate(nx,ny,4)     ! precipitation rate (kg/(m**2*s))
                               ! prcrate(1,1,1) = total precipitation rate
                               ! prcrate(1,1,2) = grid scale precip. rate
                               ! prcrate(1,1,3) = cumulus precip. rate
                               ! prcrate(1,1,4) = microphysics precip. rate
  INTEGER :: icatg(nx,ny)
!
  REAL :: prctrw,bltop,rhrkul,thevrkul
  PARAMETER (prctrw=1.0E-04,    & !  (kg/(m**2*s))
         bltop=1000.,   & !  m
         rhrkul=0.8,                                                    &
             thevrkul=9.0)  !  K**2
!
  INTEGER :: i,j
!  integer knt,k
  INTEGER :: nxlg,nylg
  INTEGER :: knt1,knt2a,knt2b,knt3,knt4
  INTEGER :: is_dup
  REAL :: pr,temp,qvsat,rh
!  real rhsum,rhmean
!  real the,thesum,thesq,themean,thevar
!
!-----------------------------------------------------------------------
!
!  Function f_qvsat and inline directive for Cray PVP
!
!-----------------------------------------------------------------------
!
  REAL :: f_qvsat

!fpp$ expand (f_qvsat)
!!dir$ inline always f_qvsat

!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
!
  INCLUDE 'phycst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
  knt1=0
  knt2a=0
  knt2b=0
  knt3=0
  knt4=0
!
!-----------------------------------------------------------------------
!
!  Initialize with default of no precip or precip influence
!
!  We jump thru some extra hoops to make sure the MPI and non-MPI counts
!  are the same.  If we don't, overlap points can be counted twice.
!
!-----------------------------------------------------------------------
!
  DO j=1,ny-1
    DO i=1,nx-1
      icatg(i,j)=1
    END DO
  END DO
!
  IF(ccatopt == 1) THEN
    DO j=1,ny-1
      is_dup = 0
      IF ((mp_opt > 0) .AND. (j >= ny-2) .AND. (loc_y .NE. nproc_y)) is_dup=1
      DO i=1,nx-1
        IF ((mp_opt > 0) .AND. (i >= nx-2) .AND. (loc_x .NE. nproc_x)) is_dup=1
!
!-----------------------------------------------------------------------
!
!  Is convective precip turned on or heavy rain occuring?
!
!-----------------------------------------------------------------------
!
        IF (prcrate(i,j,3) > 0. .OR. prcrate(i,j,1) > prctrw) THEN
          IF (is_dup == 0) knt4=knt4+1
          icatg(i,j)=4
!      print *, ' set i,j:',i,j,' cat 4',' prcrate(1):',
!    :             prcrate(i,j,1),' prcrate(3):',prcrate(i,j,3)
!
!-----------------------------------------------------------------------
!
!  Is it currently raining?
!
!-----------------------------------------------------------------------
!
        ELSE IF ( prcrate(i,j,1) > 0. ) THEN
          IF (is_dup == 0) knt3=knt3+1
          icatg(i,j)=3
!      print *, ' set i,j:',i,j,' cat 3',' prcrate(1):',
!    :             prcrate(i,j,1),' prcrate(3):',prcrate(i,j,3)
!
!-----------------------------------------------------------------------
!
!  Evidence of rain-cooled air?
!  Lapse rate close to a moist adiabatic lapse rate in the lowest 2-km,
!  and a high relative humidity in that layer.  Or a high relative
!  humidity at the first level above ground.
!
!-----------------------------------------------------------------------
!
        ELSE
!
!-----------------------------------------------------------------------
!
!  For the layer 0-2km above first grid level (k=2).
!  Calculate the mean relative humidity.
!  Calculate the standard deviation of saturated equivalent
!  potential temperature.
!
!-----------------------------------------------------------------------
!
!      rhsum=0.
!      thesum=0.
!      thesq=0.
!      knt=0
!      DO 40 k=2,nz-1
!        IF(k.gt.2 .and. (zs(i,j,k)-zs(i,j,2)).gt.bltop) GO TO 41
!        knt=knt+1
!        pr=pprt(i,j,k)+pbar(i,j,k)
!        temp=(ptprt(i,j,k)+ptbar(i,j,k)) * (pr/p0)**rddcp
!        qvsat=f_qvsat( pr,temp )
!        rh=min(1.,(qv(i,j,k)/qvsat))
!        the=F_PT2PTE( pr,(ptprt(i,j,k)+ptbar(i,j,k)),qvsat)
!        rhsum=rhsum+rh
!        thesum=thesum+the
!        thesq=thesq+(the*the)
!  40     CONTINUE
!  41     CONTINUE
!      IF(knt.gt.0) THEN
!        rhmean=rhsum/float(knt)
!        themean=thesum/float(knt)
!      ELSE
!        rhmean=0.
!        themean=0.
!      END IF
!      IF(knt.gt.1) THEN
!        thevar=(thesq-((thesum*thesum)/float(knt)))/float(knt-1)
!      ELSE
!        thevar=9999.
!      END IF
!
!      IF (rhmean.gt.rhrkul .and. thevar.lt.thevrkul) THEN
!        knt2a=knt2a+1
!        icatg(i,j)=2
!        print *, ' set i,j:',i,j,' cat 2',' prcrate(1):',
!    :             prcrate(i,j,1),' prcrate(3):',prcrate(i,j,3)
!      ELSE
!
          pr=pprt(i,j,2)+pbar(i,j,2)
          temp=(ptprt(i,j,2)+ptbar(i,j,2)) * (pr/p0)**rddcp
          qvsat=f_qvsat( pr,temp )
          rh=qv(i,j,2)/qvsat
          IF(rh > 0.8) THEN
            IF (is_dup == 0) knt2b=knt2b+1
            icatg(i,j)=2
!          print *, ' set i,j:',i,j,' cat 2*',' RH(2) = ',rh
          END IF

!      END IF
!      IF (mod(i,20).eq.0 .and. mod(j,10).eq.0)  THEN
!        print *, i,j,' rhmean=',rhmean,
!    :               ' thevar=',thevar,' icat=',icatg(i,j)
!        print *, ' '
!      END IF

        END IF
      END DO
    END DO
    IF (mp_opt > 0 ) THEN
      CALL mpsumi(knt2a,1)
      CALL mpsumi(knt2b,1)
      CALL mpsumi(knt3,1)
      CALL mpsumi(knt4,1)
      nxlg = (nx - 3) * nproc_x + 3
      nylg = (ny - 3) * nproc_y + 3
      knt1=((nxlg-1)*(nylg-1))-knt2a-knt2b-knt3-knt4
    ELSE
      knt1=((nx-1)*(ny-1))-knt2a-knt2b-knt3-knt4
    END IF

    IF (myproc == 0 ) THEN
      WRITE(6,'(1x,a)') 'Precip correlation category assignment counts'
      WRITE(6,'(5x,a,i8)') 'cat  1: ',knt1
      WRITE(6,'(5x,a,i8)') 'cat 2a: ',knt2a
      WRITE(6,'(5x,a,i8)') 'cat 2b: ',knt2b
      WRITE(6,'(5x,a,i8)') 'cat  3: ',knt3
      WRITE(6,'(5x,a,i8)') 'cat  4: ',knt4
    END IF
  END IF

  CALL iedgfill(icatg,1,nx,1,nx-1, 1,ny,1,ny-1,1,1,1,1)

  RETURN
END SUBROUTINE setcat
SUBROUTINE cvtsrcstr(srcname,srcstr)
!
! Sets name to upper case and replaces spaces with underscores.
!
  IMPLICIT NONE
  CHARACTER(LEN=8),INTENT(IN) :: srcname
  CHARACTER(LEN=8),INTENT(OUT) :: srcstr

  CHARACTER(LEN=1),PARAMETER :: underscore = '_'
  INTEGER, PARAMETER :: iaspace = 32
  INTEGER :: i,lensrc,jchar

  srcstr = srcname

!
! "srcstr" values can *MUST* start with a letter.  Therefore "88D" is illegal.
! Put an "R" in front of the name.
!

  IF (srcstr(1:3) == "88D") THEN
    srcstr = "R"//srcname
  ENDIF

  lensrc=LEN_TRIM(srcstr)

  DO i=1,lensrc
    jchar=IACHAR(srcstr(i:i))
    IF(jchar == iaspace ) THEN
      srcstr(i:i) = underscore
    ELSE IF( jchar > 96 .AND. jchar < 123) THEN
      srcstr(i:i) = ACHAR( (jchar-32) )
    END IF
  END DO
  RETURN
END SUBROUTINE cvtsrcstr
