!########################################################################
!########################################################################
!#########                                                      #########
!#########                 PROGRAM casa2arpsppi                 #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########          School of Meteorology / CASA / CAPS         #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

PROGRAM casa2arpsppi

!------------------------------------------------------------------------
!
! PURPOSE:
!
! Reads CASA radar data from NetCDF files, interpolates in time within
! a given assimilation window, then writes out results on ARPS columns
! and radar elevation levels (.gridtilt format) for use in arpsenkf.
!
!
! REQUIRED INPUT:
!
! A listfile generated using the scantier2a program developed by Keith
! Brewster
!
!------------------------------------------------------------------------
!
! AUTHOR:
!
! Nate Snook (written during 17 Dec. 2007 -- 30 Jan. 2008)
!
! 03/05/2010 (Y. Wang)
! Added tem_byte and changed interfaces of get_ncrad2ainfo and rd2atiltcdf
! to catch up with the latest package update.
!
! MODIFICATIONS:
!
! None as of yet.  Please note the abundant comments.  When modifying
! this code, please comment generously.  You'll thank yourself in six
! months when you open this code up again, and your successors won't be
! cursing your name long after you've moved on to other things.
!
!------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------
!
! Include files
!
!------------------------------------------------------------------------

  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
  INCLUDE 'netcdf.inc'

!------------------------------------------------------------------------
!
! Dimensions and parameters
!
!------------------------------------------------------------------------

  INTEGER, PARAMETER :: n_radials = 360   ! Max # of radials in an elevation.
                                          ! 360 is a natural choice, but could be something else.

  REAL, PARAMETER :: delta_azim = (n_radials / 360.0)
                                          ! What is the spacing between azimuths in degrees?

  REAL, PARAMETER :: missing  = -99999.   ! Assign this value for missing data.

  REAL, PARAMETER :: bad_data = -88888.   ! Assign this value for data found to
                                          ! be bad during processing.

  REAL, PARAMETER :: outside_window = -77777.
                                          ! Assign this value for data present,
                                          ! but outside the assimilation window.

  INTEGER, PARAMETER :: n_elevs = 8       ! The number of elevations.
                                          ! Could read this in, but I'm having allocation problems.

  INTEGER, PARAMETER :: n_times = 1600    ! Again, could read this in, but would be tricky.

  INTEGER, PARAMETER :: itime1970 = 315619200
                                          ! The offset time for 1970.  As in 'ncrad2arps'.

!------------------------------------------------------------------------
!
! Local variables
!
!------------------------------------------------------------------------

  CHARACTER(LEN=16) :: times(n_times)     ! The times of the netcdf files we will use.
                                          ! Dims: (n_times)

  INTEGER int_times(n_times)              ! An integer version of the CHARACTER variable time.

  INTEGER d_times(n_times)                ! The time (in seconds) between this data
                                          ! time and the previous one.

  CHARACTER(LEN=32) :: files(n_times)     ! The names of those netcdf files.
                                          ! Dims: (n_times)

  INTEGER :: elevs(n_elevs)               ! The elevation angles we have.
                                          ! Dims: (n_elevs)

  REAL    :: elevs_real(n_elevs)          ! A real (not integer) version of elevs.
                                          ! Dims: (n_elevs)

  LOGICAL :: data_status(n_radials,n_elevs,n_times)
                                          ! Is there data at a given time for a radial?
                                          ! Dims: (n_radials, n_elevs, n_times)

  INTEGER :: time_before(n_radials,n_elevs)! What earlier time was found for each radial?
                                           ! Dims: (n_radials, n_elevs)

  INTEGER :: time_after(n_radials,n_elevs) ! What later time was found for each radial?
                                           ! Dims: (n_radials, n_elevs)

  INTEGER :: id_before(n_radials,n_elevs)  ! What is the time index before?
                                           ! Dims: (n_radials, n_elevs)

  INTEGER :: id_after(n_radials,n_elevs)   ! What is the time index after?
                                           ! Dims: (n_radials, n_elevs)

  INTEGER :: elapsed_before(n_radials,n_elevs) ! Stores time elapsed between
                                               ! prev data and assim. time for each radial.

  INTEGER :: elapsed_after(n_radials,n_elevs)  ! Stores time elapsed between
                                               ! assim. time and next data for each radial.

  REAL, ALLOCATABLE :: z_3d_before(:,:,:)      ! 3D reflectivity volume for earlier time.
                                               ! Dims: (n_gates, n_radials, n_elevs)

  REAL, ALLOCATABLE :: z_3d_after(:,:,:)       ! 3D reflectivity volume for later time.
                                               ! Dims: (n_gates, n_radials, n_elevs)

  REAL, ALLOCATABLE :: z_3d_valid(:,:,:)       ! Interpolated 3D reflectivity volume.
                                               ! Dims: (n_gates, n_radials, n_elevs)

  REAL, ALLOCATABLE :: vr_3d_before(:,:,:)     ! 3D Vr volume for earlier time.
                                               ! Dims: (n_gates, n_radials, n_elevs)

  REAL, ALLOCATABLE :: vr_3d_after(:,:,:)      ! 3D Vr volume for later time.
                                               ! Dims: (n_gates, n_radials, n_elevs)

  REAL, ALLOCATABLE :: vr_3d_valid(:,:,:)      ! Interpolated 3D Vr volume.
                                               ! Dims: (n_gates, n_radials, n_elevs)

  REAL, ALLOCATABLE :: gridtilt_z(:,:,:)       ! Gridtilt reflectivity.
                                               ! Dims: (nx, ny, n_elevs)

  REAL, ALLOCATABLE :: gridtilt_vr(:,:,:)      ! Gridtilt radial velocity.
                                               ! Dims: (nx, ny, n_elevs)

  INTEGER :: first_radial             ! First radial that data is present for a given elevation.
  INTEGER :: first_radials(n_times)   ! An array to store the first_radial data.
  INTEGER :: last_radial              ! Last radial that data is present for a given elevation.
  INTEGER :: last_radials(n_times)    ! An array to store the last_radial data.
  REAL    :: elevation                ! The elevation associated with first_radial and last_radial.
  INTEGER :: elev_index               ! The array index associated with the given elevation angle.
  INTEGER :: sector_size              ! How big was the sector?
  INTEGER :: sector_sizes(n_times)              !An array to store the sector_size data.

  CHARACTER(LEN=80) :: radar_dir      ! directory containing our radar data.
                                      ! Length 80 important for use with NetCDF Reader.

  CHARACTER(LEN=4) :: radar_name      ! Name of the radar we're working with
                                      ! (e.g. KSAO, KCYR, KRSP, KLWE)

  INTEGER :: timeloopcount        !A counter to determine how many times we looped through the time loop.
  INTEGER :: intswap              !An integer for swapping variables

!Note -- the following four variables could be set in such a way as to be inconsistent.  AVOID THIS!
  INTEGER :: n_assim              !The number of assimilation cycles desired.
  INTEGER :: t_first_assim        !The time of the first data assimilation.  Format: hhmmss
  INTEGER :: t_last_assim         !The time of the last data assimilation.  Format: hhmmss
  INTEGER :: t_assim              !The assimilation time.  Format: hhmmss
  INTEGER :: dt_assim             !The interval between assimilation periods, IN SECONDS

  INTEGER :: cur_assim            !Time of the current assimilation cycle in hhmmss format
  INTEGER :: cur_assim_sec        !Time of current assimilation cycle IN SECONDS SINCE MIDNIGHT
  INTEGER :: elapsed_time         !Seconds elapsed since the index nearest the assimilation time.
  INTEGER :: index_reference      !Time index immediately before an assimilation cycle, for reference.
  INTEGER :: index_before         !Nearest time index before the current assimilation cycle where data exists.
  INTEGER :: index_after          !Nearest time index after the current assimilation cycle where data exists.
  INTEGER :: index_cdf            !What time index are we at while looping through NetCDF files?
  INTEGER :: begin_index          !What's the earliest index we could possibly need data from for this cycle?
  INTEGER :: end_index            !What's the latest index we could possibly need data from for this cycle?
  INTEGER :: cur_index            !What index are we at now?
  INTEGER :: elapsed_counter      !A counter to use in determining min_index and max_index.
                                  !   Remember:  Time indicies are defined by when radar files are present.

  LOGICAL :: found_flag           !A flag for use when searching for nearest times.
  LOGICAL :: forward_flag         !A flag to tell us which direction we're scanning.
  LOGICAL :: timeflag             !A flag for calculating time_arpsppi
  LOGICAL :: debugflag            !A flag for use in debugging.

  CHARACTER :: h1, m1, s1, h2, m2, s2           !For time conversion between CHARACTER and INTEGER
  INTEGER :: ih1, im1, is1, ih2, im2, is2       !For time conversion between INTEGER and CHARACTER
  CHARACTER :: ph1, pm1, ps1, ph2, pm2, ps2     !For finding dtime

!Note -- here we are using integers from 1 to 360 for our 3-D arrays of Z and Vr.  Any azimuthal...
!...resolution could be used by altering the casa2arpsppi code -- it doesn't have to be 1 degree.
  REAL, ALLOCATABLE :: refl_on_integer(:,:)     !Z on azimuths that are integers from 1 to 360.  Dims:(n_gates,true_nradials)
  REAL, ALLOCATABLE :: radv_on_integer(:,:)     !Vr on azimuths that are integers from 1 to 360.  Dims:(n_gates,true_nradials)
  REAL, ALLOCATABLE :: refl_integer_360(:,:)    !refl_on_integer(:,:) put onto a full 360 degree scan.
  REAL, ALLOCATABLE :: radv_integer_360(:,:)    !radv_on_integer(:,:) put onto a full 360 degree scan.
  REAL, ALLOCATABLE :: int_azims(:)             !The angles of the integer azimuths.  Dims:(true_nradials)
  INTEGER :: j_int_azim                         !For looping through the first spatial interpolation (to desired radials)
  INTEGER :: i_obs_azim                         !In the first spatial interpolation, what observed azimuth are we at?
  INTEGER :: k_cur_gate                         !For processing an individual radial in the first spatial interpolation.
  REAL :: dazim_left                            !for the linear interpolation
  REAL :: dazim_right                           !for the linear interpolation

  INTEGER :: scan_increment        !What increment do we use between radials when transferring to the full 360 scan?
  INTEGER :: start_radial          !Where do we begin for the transfer to 360 scan?
  INTEGER :: end_radial            !Where do we end for the transfer to 360 scan?

!------------------------------------------------------------------------
!
! ARPS subroutine variables
!
!------------------------------------------------------------------------

  INTEGER :: nx, ny, nz            !ARPS domain dimensions in the x, y, and z directions (number of gridpoints).
  INTEGER :: nzsoil                !Number of soil gridpoints (gridpoints below surface--in the -z direction)
  INTEGER :: nstypes               !Number of soil types.  Needed for calling 'initpara'.
  REAL, ALLOCATABLE :: x(:)        !x-coord. of the physical computational grid.  Defined at u-point on staggered grid.  Dims:(nx)
  REAL, ALLOCATABLE :: y(:)        !y-coord. of the physical computational grid.  Defined at v-point on staggered grid.  Dims:(ny)
  REAL, ALLOCATABLE :: z(:)        !z-coord. of the physical computational grid.  Defined at w-point on staggered grid.  Dims:(nz)
  REAL, ALLOCATABLE :: zp(:,:,:)   !Vertical coordinatel of gridpoints in physical space.  Dims:(nz)
  REAL, ALLOCATABLE :: zpsoil(:,:,:) !z-coord. of soil model computational grid.  Defined at center of soil layer.  Dims:(nz)
  REAL, ALLOCATABLE :: h_terrain(:,:)           !Height of the terrain in meters.  Dims:(nx,ny)
  REAL, ALLOCATABLE :: mapfct(:,:,:)            !Map factors at scalar, u and v points. For calling 'inigrid'.  Dims:(nx,ny,8)
  REAL, ALLOCATABLE :: j1(:,:,:)   !Coordinate transformation Jacobian.  For calling 'inigrid'.  Dims:(nx,ny,nz)
  REAL, ALLOCATABLE :: j2(:,:,:)   !Coordinate transformation Jacobian.  For calling 'inigrid'.  Dims:(nx,ny,nz)
  REAL, ALLOCATABLE :: j3(:,:,:)   !Coordinate transformation Jacobian.  For calling 'inigrid'.  Dims:(nx,ny,nz)
  REAL, ALLOCATABLE :: j3soil(:,:,:)            !Coordinate transformation Jacobian.  For calling 'inigrid'.  Dims:(nx,ny,nzsoil)
  REAL, ALLOCATABLE :: j3soilinv(:,:,:)         !Inverse of 'j3soil'.  For calling 'inigrid'.  Dims:(nx,ny,nzsoil)
  REAL, ALLOCATABLE :: tem1d1(:), tem1d2(:)     !Temporary 1D arrays.  For calling 'inigrid'.  Dims:(nz)
  REAL, ALLOCATABLE :: tem3d(:,:,:)             !Temporary 3D array.  For calling 'inigrid'.  Dims:(nx,ny,nz)
  REAL, ALLOCATABLE :: xs(:)       !???  A function of x-coord. needed for calling 'radcoord'.
  REAL, ALLOCATABLE :: ys(:)       !???  A function of y-coord. needed for calling 'radcoord'.
  REAL, ALLOCATABLE :: zps(:,:,:)  !???  A function of zp needed for calling 'radcoord'.
  REAL :: radar_x                  !x-coord. of the radar on the ARPS grid.  Calculated in 'radcoord' call.
  REAL :: radar_y                  !y-coord. of the radar on the ARPS grid.  Calculated in 'radcoord' call.

  INTEGER, ALLOCATABLE :: kntazm(:)             !Number of radials at each elevation angle. Needed for 'remap_arpsppi'  Dims:(n_elevs)
  INTEGER, ALLOCATABLE :: kntgat_ref(:,:)       !???-Z  Needed for 'rmpinit'  Dims:(n_radials,n_elevs)
  INTEGER, ALLOCATABLE :: kntgat_vel(:,:)       !???-Vr  Needed for 'rmpinit'  Dims:(n_radials,n_elevs)

  REAL :: range_min, range_max     !Minimum and maximum range, in meters, to process data from in 'remap_arpsppi'.
  INTEGER :: time_arpsppi          !Assimilation time to feed into 'abss2ctim'  All our data is considered valid at this time.
  REAL, ALLOCATABLE :: rngvol(:,:) !Range to gate in 3D volume for Z and Vr.  Used in 'wrttilts'  Dims:(n_gates,n_elevs)
  REAL, ALLOCATABLE :: azmvol(:,:) !Azimuth angle in 3D volume for Z and Vr.  Used in 'remap_arpsppi'  Dims:(n_radials,n_elevs)
  REAL, ALLOCATABLE :: elvvol(:,:) !Elevation angle in 3D volume for both Z and Vr.  Used in 'remap_arpsppi'  Dims:(n_radials,n_elevs)

  REAL, ALLOCATABLE :: gridtilthigh_ref(:,:,:)  !???  Gridtilt output from 'remap2dcts' for Z.  Dims:(nx,ny,n_elevs)
  REAL, ALLOCATABLE :: gridrange_ref(:,:,:)     !???  Gridtilt output from 'remap2dcts' for Z.  Dims:(nx,ny,n_elevs)
  REAL, ALLOCATABLE :: gridslr_ref(:,:)         !???  Gridtilt output from 'remap2dcts' for Z.  Dims:(nx,ny)
  REAL, ALLOCATABLE :: gridazm_ref(:,:)         !???  Gridtilt output from 'remap2dcts' for Z.  Dims:(nx,ny)
  REAL, ALLOCATABLE :: gridtiltval_ref(:,:,:)   !Gridtilt reflectivity output by 'remap2dcts'.  Dims:(nx,ny,n_elevs)
  REAL, ALLOCATABLE :: gridtilttime_ref(:)      !Valid time for each elevation for Z.           Dims:(n_elevs)
  REAL, ALLOCATABLE :: elevmean_ref(:)          !Mean reflectivity for each elevation, I guess. Dims:(n_elevs)

  REAL, ALLOCATABLE :: gridtilthigh_vel(:,:,:)  !???  Gridtilt output from 'remap2dcts' for Vr.    Dims:(nx,ny,n_elevs)
  REAL, ALLOCATABLE :: gridrange_vel(:,:,:)     !???  Gridtilt output from 'remap2dcts' for Vr.    Dims:(nx,ny,n_elevs)
  REAL, ALLOCATABLE :: gridslr_vel(:,:)         !???  Gridtilt output from 'remap2dcts' for Vr.    Dims:(nx,ny)
  REAL, ALLOCATABLE :: gridazm_vel(:,:)         !???  Gridtilt output from 'remap2dcts' for Vr.    Dims:(nx,ny)
  REAL, ALLOCATABLE :: gridtiltval_vel(:,:,:)   !Gridtilt radial velocity output by 'remap2dcts'.  Dims:(nx,ny,n_elevs)
  REAL, ALLOCATABLE :: gridtilttime_vel(:)      !Valid time for each elevation for Vr.             Dims:(n_elevs)
  REAL, ALLOCATABLE :: elevmean_vel(:)          !Mean radial velocity for each elevation, I guess. Dims:(n_elevs)

  INTEGER :: iyr, iday, imon, ihr, imin, isec   !For use in 'abss2ctim', also needed in 'wrttilts' and 'wrtgridtilt'
  CHARACTER (LEN = 100) :: rfname               !Name to use for the radar file.  Passed to 'wrttilts' and 'wrtgridtilt'
  CHARACTER (LEN = 100) :: rfnametilt           !Similar to rfname.  Used in 'wrttilts' and wrtgridtilt'

  REAL :: gatesp                                !Gate spacing, in meters.  Needed for 'remap_arpsppi'.

!------------------------------------------------------------------------
!
! NetCDF variables
!
!------------------------------------------------------------------------

  !Introduced in call to 'get_ncraddims'
  CHARACTER(LEN=256) :: filename   !The NetCDF file to use -- chosen from files(ntimes)
  INTEGER :: iradtype              !Radar type -- read in from NetCDF
  INTEGER :: n_gates               !The number of gates, read in from NetCDF
  INTEGER :: true_nradials         !The actual number of radials, read in from NetCDF
                                   !(earlier, assumed to be 360)
  INTEGER :: nvar                  !Number of variables present in NetCDF file.
                                   !Probably not of interest to us.
  INTEGER :: istatus               !Make sure the NetCDF reading routines don't fail.

  !Introduced in call to 'get_ncrad2ainfo'
  REAL(KIND=8), ALLOCATABLE :: tem_double(:)
                                   !For use in get_ncrad2ainfo call... not sure what it does.
  REAL,         ALLOCATABLE :: tem_real(:,:)

  CHARACTER(LEN=32) :: radar_name_cdf           !Length is important when calling NetCDF routines.
  REAL :: radar_lat                !Latitude of the radar.  Read in from NetCDF.
  REAL :: radar_lon                !Longitude of the radar.  Read in from NetCDF.
  REAL :: radar_altitude           !Altitude (elevation) of the radar.  Read in from NetCDF.
  INTEGER :: itimcdf               !The time that the NetCDF file is valid.  Needed for NetCDF calls.
  INTEGER :: timetmp               !A temporary variable for storing time.
  REAL :: frtime                   !???.  Some variable needed for calling get_ncrad2ainfo and rd2atiltcdf.
  REAL :: elv                      !Elevation angle?  Needed for calling get_ncrad2ainfo.
  CHARACTER(LEN=NF_MAX_NAME), ALLOCATABLE :: ncvarname(:)  !Variable names.  For use with get_ncrad2ainfo.

!Introduced in call to 'rd2atiltcdf'
  REAL :: vnyq                     !nyquist velocity of the radar.  Read in from NetCDF
  REAL :: rfirstg                  !Range to the first gate, in meters.  Needed for calling rd2atiltcdf.
  INTEGER :: irfirstg              !Integer variant of rfirstg.  Needed later in 'volbuild'.
  REAL :: rmisval                  !Missing value used in rd2atiltcdf.  Is DIFFERENT FROM casa2arpsppi variable 'missing'.
  REAL :: rngfval                  !???.  A value used in rd2atiltcdf.  Not sure what.  Useless elsewhere in casa2arpsppi.
  REAL :: initime                  !???.  Used in the call to rd2atiltcdf.  Utterly useless elsewhere in casa2arpsppi.
  REAL :: beamwidth                ! Beamwidth.  Read in from NetCDF.
  INTEGER, ALLOCATABLE :: gcfst(:)

  REAL, ALLOCATABLE :: istrgate(:) !???.  Some variable needed for calling rd2atiltcdf.  Dims:(true_nradials)
  REAL, ALLOCATABLE :: azim(:)     !Azimuth angle?  Needed for calling rd2atiltcdf.  Dims:(true_nradials)
  REAL, ALLOCATABLE :: gtspc(:)    !Gate spacing.  Needed for calling rd2atiltcdf.  Dims:(true_nradials)
  REAL, ALLOCATABLE :: attref(:,:) !Attenuated reflectivity?  Read in from NetCDF.  Dims:(n_gates, true_nradials)
  REAL, ALLOCATABLE :: refl(:,:)   !Observed/Corrected Reflectivity?  Read in from NetCDF.  Dims:(n_gates, true_nradials)
  REAL, ALLOCATABLE :: radv(:,:)   !Observed Radial Velocity.  Read in from NetCDF.  Dims:(n_gates, true_nradials)
  REAL, ALLOCATABLE :: rhohv(:,:)  !Cross-polarization correlation coefficient.  Read in from NetCDF.  Dims:(n_gates, true_nradials)
  REAL, ALLOCATABLE :: tem2d(:,:)  !temporary array for dual-pol variables. Not used right now. Dims:(n_gates, true_nradials)
  INTEGER, ALLOCATABLE :: gateqc(:,:)   ! Gate Quality Control Indicator  Read in from NetCDF.  Dims:(n_gates, true_nradials)

!------------------------------------------------------------------------
!
! Input variables
!
!------------------------------------------------------------------------

  CHARACTER(LEN=64) :: inputfile   !Name of the input file for casa2arpsppi
  CHARACTER(LEN=80) :: parsetmp, parsetmp2      !Temporary arrays for parsing the listfile
  INTEGER :: n_radars              !Number of radars to process data for.

!------------------------------------------------------------------------
!
! Loop counter variables
!
!------------------------------------------------------------------------

  INTEGER :: p_radar               !A counter variable for the outer DOWHILE loop.
  INTEGER :: q_time                !A counter variable for looping through times for a given radar.
  INTEGER :: r_elevation           !A counter for looping through elevations.
  INTEGER :: s_radial              !A counter for looping through radials.
  INTEGER :: t_gate                !A counter for looping through gates.
  INTEGER :: q_assim               !A counter to loop through assimilation cycles.
  INTEGER :: count,i,j,k           !Generic counters

!------------------------------------------------------------------------
!
! Namelist Input
!
!------------------------------------------------------------------------

  NAMELIST /casa2arps/ inputfile, t_assim, range_min, range_max

  CHARACTER(LEN=256) :: namelist_filename

!------------------------------------------------------------------------
!
! Output variables
!
!------------------------------------------------------------------------

  ! There are no output variables as of yet

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Executable code begins here!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!------------------------------------------------------------------------
! Specify constants and initialize arrays
!------------------------------------------------------------------------

  n_radars = 4
  p_radar  = 0   ! We're beginning with the first radar.  This is our counter.
  q_time   = 0   ! ...and we're beginning from the first time.
  elevs    = (/ 1.0, 2.0, 3.0, 5.0, 7.0, 9.0, 11.0, 14.0 /)
                 ! Could read this in, I suppose...

  first_radials = 0      ! Initialize the first and last radial arrays...
  last_radials  = 0      ! ...so they're clean for later use.
  timetmp       = 0      ! Same with 'timetmp'.
  t_assim       = 013000

!n_assim = 13
!t_first_assim = 013500 !In hhmmss format.
!t_last_assim = 023000  !hhmmss format.
!dt_assim = 300         !THIS ONE IN SECONDS!

!-----------------For use with a SINGLE FILE ONLY!---------------!
  n_assim = 1
  t_first_assim = t_assim
  t_last_assim = t_assim
  dt_assim = 300
!----------------------------------------------------------------!

  data_status = .FALSE.  !Initialize the data status for all times and radials to false.
  timeflag    = .FALSE.  !Initialize the time flag to false.
  debugflag   = .FALSE.  !Initialize the debugging flag to false.
  timeloopcount = 0

  time_before = missing  !Initialize time before...
  time_after  = missing  !...and after to the missing data value.
  id_before   = missing  !Same with the indices.
  id_after    = missing

  range_min = 500.0      !meters
  range_max = 45000.0    !meters

!------------------------------------------------------------------------
! Allocate variables.
!------------------------------------------------------------------------

  !ALLOCATE(times(n_times))
  !ALLOCATE(elevs(n_elevs))
  !ALLOCATE(data_status(n_radials,n_elevs,n_times))
  !ALLOCATE(files(n_times))

!------------------------------------------------------------------------
! Begin by opening and parsing the listfile.
!------------------------------------------------------------------------

  inputfile = '/home/nsnook/arps5.2.8/radar/listfile_KSAO.txt'
                                                            !Could read this in.
  READ(5,casa2arps)

  OPEN(42, FILE=TRIM(inputfile), STATUS="unknown", FORM="formatted")
                                                            !open the listfile

  !DO WHILE (i_radar .le. n_radars)                         !Get data for each radar.
  p_radar = p_radar + 1                !Increment our radar loop counter
  DO WHILE (parsetmp .ne. 'Radar')     !This do loop searches for the radar directory.
    READ(42,'(A5, 12x, A64)') parsetmp, parsetmp2
  END DO

  radar_dir  = TRIM(parsetmp2)         !Our radar directory should be contained here.
  parsetmp   = ADJUSTR(radar_dir)      !We'll use the other 'parsetmp' to get the radar name.
  radar_name = parsetmp((LEN(parsetmp) - 3):LEN(parsetmp))
                                       ! It should be the last 4 letters.

  PRINT*, 'Processing data for: ', radar_name          !Output the current radar name...
  PRINT'(1x,A17,A64)', 'Radar Directory: ', radar_dir  !...and directory to the terminal.

  DO WHILE (parsetmp .ne. 'BEGIN')   !Now we need to continue until we reach the beginning...
    READ(42,*) parsetmp              !...of the file listing.  This is where we begin reading...
  END DO                             !...in lines containing filenames, times, and sectors.

  DO WHILE (parsetmp .ne. ' Volumes for this radar')
                                     ! 'All' appears at the end of a radar block,
                                     ! before a new radar begins.
    timeloopcount = timeloopcount + 1       !Increment our time loop counter.

    READ(42,'(3x,A27,3x,F4.1,7x,I3,7x,I3,7x,I3)') parsetmp, elevation,  &
            first_radial, last_radial, sector_size  !Read in lines of data.

    IF(parsetmp .eq. 'IN NEW VOLUME') THEN  !If it's the beginning of a new volume...
      timeloopcount = timeloopcount - 1     !...then undo the counter increment
                                            !so this read isn't counted...

      CYCLE                                 !...and skip the "BEGIN NEW VOLUME" line
                                            !and read more data.
    END IF

    IF(parsetmp .eq. ' Volumes for this radar') THEN
                              !If we've reached the end of data for a radar...
      CYCLE                   !...then break out of the loop and begin the next radar.
    END IF

    IF(elevation < 1.0) THEN  ! Correct for one extra space
                              ! due to formatting of the list file.
      elevation = elevation*10.0
    END IF

    if (timeloopcount > n_times) exit
    first_radials(timeloopcount) = first_radial  !Store the first and last radial data...
    last_radials(timeloopcount) = last_radial    !...for later use.
    sector_sizes(timeloopcount) = sector_size    !We'll need the sector size data later too.

    files(timeloopcount) = parsetmp          !First entry contains the filename.
    parsetmp = ADJUSTL(files(timeloopcount)) !Make sure the filename is shifted
                                             !left for our next step...
    times(timeloopcount) = parsetmp(15:20)   !...now the date and time should be
                                             !contained here.

!------------------------------------------------------------------------
! A CHARACTER time doesn't do us much good, so convert it to INTEGER.
!------------------------------------------------------------------------

    IF (timeloopcount .ne. 1) THEN
       ph1 = h1                   !Set the previous time's hours, minutes, and seconds...
       ph2 = h2                   !...before we update them for the current iteration.
       pm1 = m1
       pm2 = m2                   !We'll need these to calculated d_times...
       ps1 = s1                   !...which tells us how far apart our data are in time.
       ps2 = s2
    END IF

    h1 = parsetmp(15:15)          !tens place for hours
    h2 = parsetmp(16:16)          !ones place for hours
    m1 = parsetmp(17:17)          !tens place for minutes
    m2 = parsetmp(18:18)          !ones place for minutes
    s1 = parsetmp(19:19)          !tens place for seconds
    s2 = parsetmp(20:20)          !ones place for seconds

    int_times(timeloopcount) = ((100000*(ichar(h1)-48))                 &
                             + (10000*(ichar(h2)-48))                   &
                             + (1000*(ichar(m1)-48))                    &
                             + (100*(ichar(m2)-48))                     &
                             + (10*(ichar(s1)-48))                      &
                             + (ichar(s2)-48))
             ! In case you're wondering, to get the corresponding integer,
             ! you must take (ichar(x) - 48).

    IF (timeloopcount .eq. 1) THEN   !Unless this is the first iteration
                                     !through the loop...
       d_times(timeloopcount) = 0
    ELSE                             !...calculate how many seconds have
                                     !elapsed since the last data time.
       d_times(timeloopcount) = 36000*(ichar(h1)-ichar(ph1))            &
                              +  3600*(ichar(h2)-ichar(ph2))            &
                              +   600*(ichar(m1)-ichar(pm1))            &
                              +    60*(ichar(m2)-ichar(pm2))            &
                              +    10*(ichar(s1)-ichar(ps1))            &
                              +       (ichar(s2)-ichar(ps2))
             ! Here ichar is used in a difference, so the (-48) is unnecessary.
    END IF

!--------------------------Debugging Print Block------------------------
!      PRINT'(A17,I4,6x,A9,A15)', 'Processing File #',timeloopcount,    &
!                                 'At time: ', times(timeloopcount)
!      PRINT*, 'Integer time: ', int_times(timeloopcount)
!      PRINT*, 'Filename = ', files(timeloopcount)
!      PRINT*, 'Sector size = ', sector_size
!
!      IF ((first_radial - last_radial .eq. sector_size) .or.           &
!          (last_radial - first_radial .eq. sector_size)) THEN
!         PRINT*, 'Size matches sector -- did not cross zero.'
!      ELSE
!         PRINT*, 'Crossed zero during this scan'
!      END IF
!
!      IF (sector_size .eq. 359) THEN
!         PRINT*, 'Full 360 scan'
!      END IF
!
!      PRINT '(A45,3x,I3,3x,I3,3x,F4.1)',                               &
!            'First radial, last radial, and elevation are: ',          &
!            first_radial, last_radial, elevation
!      PRINT*, '-------------------------------------------------------'
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
! The listfile has been parsed.  Now use that data to fill the
! data_status array with information about which radials are present and
! which are missing, where TRUE is present, and FALSE is absent.
!------------------------------------------------------------------------

    IF (first_radial > last_radial) THEN    !We want the first radial to be the lower number...
       intswap = first_radial               !...the direction of the sector is unimportant...
       first_radial = last_radial           !...we only care where data is present.
       last_radial = intswap
    END IF

    SELECT CASE(INT(elevation))    !Determine which index our elevation corresponds to.
      CASE(1)
        elev_index = 1
      CASE(2)
        elev_index = 2
      CASE(3)
        elev_index = 3
      CASE(5)
        elev_index = 4
      CASE(7)
        elev_index = 5
      CASE(9)
        elev_index = 6
      CASE(11)
        elev_index = 7
      CASE(14)
        elev_index = 8
      CASE DEFAULT
        STOP 'Elevation index doex not match given values! -- Please check the listfile!'
    END SELECT

    IF (abs(last_radial - first_radial - sector_size) <3 )  THEN
             ! This will only be true when we have not crossed zero.
      DO s_radial=first_radial,last_radial,1
             ! If we did not cross zero, data is present
             ! between first_radial and last_radial.
        data_status(s_radial,elev_index,timeloopcount) = .TRUE.
             ! There is data here!
      END DO
    ELSE
      data_status(:,elev_index,timeloopcount) = .TRUE.
             ! If not, then data is absent between first_radial and last_radial...
!     PRINT*, "--",timeloopcount,"--"
!     PRINT*, data_status(:,elev_index,timeloopcount)
!      write(0,*) "********",timeloopcount, elev_index,first_radial,last_radial
      IF (first_radial < 1) CYCLE

      DO s_radial=first_radial,last_radial,1            !...so we need to run the loop in reverse.
         data_status(s_radial,elev_index,timeloopcount) = .FALSE.
                              !There is no data here (except for a full 360)
      END DO
      data_status(first_radial,elev_index,timeloopcount) = .TRUE.
                              !Data at the first_radial and the last radial...
      data_status(last_radial,elev_index,timeloopcount) = .TRUE.
                              !...are present, but the loop marked them absent.
    END IF

!----------------------Another print block...-----------------------------!
!      IF (timeloopcount .eq. 192) THEN
!         PRINT*, "Elevation and Time indices:", elev_index, timeloopcount
!         PRINT*, "First and last radials: ", first_radial, last_radial, "Sector Size: ", sector_size
!         PRINT*, "Data Status:"
!         DO count = 1, n_radials, 1
!            PRINT*, "Status is: ", data_status(count,elev_index,timeloopcount), "for ", count
!         END DO
!      END IF
!      PRINT*, data_status((first_radial - 1),elev_index,timeloopcount),
!              data_status(first_radial,elev_index,timeloopcount),
!              data_status((first_radial + 1), elev_index,timeloopcount),
!     ' --- ', data_status((last_radial - 1),elev_index,timeloopcount),
!              data_status(last_radial, elev_index, timeloopcount),
!              data_status((last_radial + 1), elev_index, timeloopcount)
!      PRINT*, '-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
!-------------------------------------------------------------------------!

  END DO

  cur_assim = t_first_assim        !Start with the first assimilation cycle

  DO q_assim=1,n_assim             !For each assimilation time...

!------------------------------------------------------------------------
! Find the time index corresponding to our current assimilation cycle.
!------------------------------------------------------------------------

    count = 1                               !...begin at the start of the times array...
    DO WHILE (int_times(count) < cur_assim) !...and loop through it until
                                            ! we find the index of the assimilation cycle.
      count = count + 1
      IF (count > n_times) THEN
        STOP "You broke it... assimilation time not within time range of input data."
      END IF
    END DO
    index_reference = count - 1             !The index we find is the one immediately
                                            !after our assimilation time.
    index_before = count - 1                !We'll start our search for earlier data
                                            !from the previous time...
    index_after = count                     !...and our search for later data from this time.

!-------------------------------Another print block!-----------------------!
    PRINT*, "Assimilation time = ", cur_assim
    PRINT*, "Nearest data time before =", int_times(index_before)
    PRINT*, "Nearest data time after =", int_times(index_after)
    PRINT*, "Reference time index is:", index_reference
    PRINT*, "--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--=--"
!--------------------------------------------------------------------------!

!------------------------------------------------------------------------
! Now, for each radial, search for the nearest time before and after
! where data is present.
!------------------------------------------------------------------------

    DO r_elevation = 1,n_elevs             !For each elevation...
      DO s_radial = 1,n_radials            !...and each radial, search for
                                           !the nearest earlier existing data.
        found_flag = .FALSE.               !Reset the "found" flag.
        index_before = index_reference     !Set the initial guess for nearest previous data time.
        elapsed_time = 0
        DO WHILE ( .NOT. found_flag )
           IF ( data_status(s_radial,r_elevation,index_before) ) THEN
             time_before(s_radial,r_elevation) = int_times(index_before)
                                           !This is the one we're looking for!
             id_before(s_radial,r_elevation) = index_before
             found_flag = .TRUE.           !Note that we've found it.
           ELSE
             index_before = index_before - 1
             elapsed_time = elapsed_time + d_times(index_before + 1)  !Update the time elapsed.
             IF (index_before < 1) THEN
                found_flag = .TRUE.       !If we're at the beginning of the array,
                                          !then there is no earlier data.
             END IF
           END IF
        END DO
        elapsed_before(s_radial,r_elevation) = elapsed_time
                                          !Store the total elapsed time... it'll be useful later.
      END DO
    END DO

    DO r_elevation = 1,n_elevs            !Now, pull the same trick as above...
       DO s_radial = 1,n_radials             !...but searching for the nearest data time after.
          found_flag = .FALSE.               !Reset the "found" flag.
          index_after = index_reference + 1  !initial guess for nearest following data time.
          elapsed_time = 0
          DO WHILE ( .NOT. found_flag )
            IF ( data_status(s_radial,r_elevation,index_after) ) THEN
                time_after(s_radial,r_elevation) = int_times(index_after)
                                          !This is the one we're looking for!
                id_after(s_radial,r_elevation) = index_after
                found_flag = .TRUE.
             ELSE
                index_after = index_after + 1
                elapsed_time = elapsed_time + d_times(index_after)
                                          !Update the elapsed time.
                IF (index_after > n_times) THEN
                   found_flag = .TRUE.    !If we're at the end of the array,
                                          !then there is no later data.
                END IF
             END IF
          END DO
          elapsed_after(s_radial,r_elevation) = elapsed_time
                                          !Store the total elapsed time... it'll be useful later
       END DO
    END DO

!------------------------------------------------------------------------
! Now, mark locations where data is not close to the assimilation window.
! Follow this procedure:
!
!      1)  If data is present within 1/2 cycle one one side but not on
!          the other, then search out to 1 full cycle on the side where
!          data was not found.  If you find it, use it.  If not, mark the
!          radial as 'outside window'.  Keep the data on the good side.
!
!      2)  If data is not present within 1/2 cycle on either side, mark
!          the radial as 'outside window' both before and after.
!------------------------------------------------------------------------

    DO r_elevation = 1,n_elevs
      DO s_radial = 1,n_radials
        IF (elapsed_after(s_radial,r_elevation) > dt_assim) THEN
                    !If the next data is more than a full cycle after, then...
          time_after(s_radial,r_elevation) = outside_window
                    !...mark this radial as having data after outside window...
          id_after(s_radial,r_elevation) = outside_window
                    !...for both time and time index.
          IF (elapsed_before(s_radial,r_elevation) > (dt_assim/2)) THEN
            time_before(s_radial,r_elevation) = outside_window
                    !But only mark data before outside window if it is...
            id_before(s_radial,r_elevation) = outside_window
          END IF    !...more than a half cycle away from assimilation time.
        END IF
        IF (elapsed_before(s_radial,r_elevation) > dt_assim) THEN
                    !Same as above but for a full cycle before.
           time_before(s_radial,r_elevation) = outside_window
           id_before(s_radial,r_elevation) = outside_window
           IF (elapsed_after(s_radial,r_elevation) > (dt_assim/2)) THEN
              time_after(s_radial,r_elevation) = outside_window
              id_after(s_radial,r_elevation) = outside_window
           END IF
         END IF
         IF ((elapsed_after(s_radial,r_elevation) > (dt_assim/2)) .and. &
             (elapsed_before(s_radial,r_elevation) > (dt_assim/2))) THEN
            time_before(s_radial, r_elevation) = outside_window
                    !And if data is more than 1/2 cycle away on both sides...
            id_before(s_radial,r_elevation) = outside_window
            time_after(s_radial, r_elevation) = outside_window
                    !...then mark this radial as having data outside window.
            id_after(s_radial,r_elevation) = outside_window
         END IF

!--------------------------------Another one of those debugging print blocks-----------------!
!            PRINT*, "For elevation #", r_elevation, "and radial #", s_radial
!
!            PRINT*, "Nearest data before: ", time_before(s_radial,r_elevation)
!            PRINT*, "Seconds from prev. data:", elapsed_before(s_radial,r_elevation)
!
!            PRINT*, "Nearest data after: ", time_after(s_radial,r_elevation)
!            PRINT*, "Seconds until next data:", elapsed_after(s_radial,r_elevation)
!
!            IF (elapsed_after(s_radial,r_elevation) > (dt_assim/2) .and. elapsed_before(s_radial,r_elevation) > (dt_assim/2)) THEN
!               PRINT*, "DATA REJECTED:  Greater than (dt_assim/2) on both sides"
!            END IF
!            IF (elapsed_after(s_radial,r_elevation) > dt_assim) THEN
!               IF (elapsed_before(s_radial,r_elevation) > (dt_assim/2)) THEN
!                  PRINT*, "DATA REJECTED:  Nearest data after > dt_assim away, and no nearby data before."
!               ELSE
!                  PRINT*, "NO INTERPOLATION:  Data after > dt_assim away, but good data before -- will use."
!               END IF
!            END IF
!            IF (elapsed_before(s_radial,r_elevation) > dt_assim) THEN
!               IF (elapsed_after(s_radial,r_elevation) > (dt_assim/2)) THEN
!                  PRINT*, "DATA REJECTED:  Nearest data before > dt_assim away, and no nearby data after."
!               ELSE
!                  PRINT*, "NO INTERPOLATION:  Data before > dt_assim away, but good data after -- will use."
!               END IF
!            END IF
!            PRINT*, "------------------------------------------------"
!---------------------------------------------------------------------------------------------!

      END DO
    END DO

!------------------------------------------------------------------------
! Before we start looping through times, we need to allocate some arrays.
! We assume the radar is operating with a constant number of range gates.
!------------------------------------------------------------------------
    PRINT*, "Did we make it here?"

    filename = TRIM(files(1))             !Take a peek at the first radar file

    CALL get_ncraddims(filename,radar_dir,iradtype,n_gates,true_nradials,  &
                       nvar,istatus)      !Call this to find out how many gates there are.

    ALLOCATE(z_3d_before(n_gates,n_radials,n_elevs))   !Start by allocating the arrays...
    ALLOCATE(z_3d_after(n_gates,n_radials,n_elevs))    !...to store the data for the final...
    ALLOCATE(z_3d_valid(n_gates,n_radials,n_elevs))    !...linear temporal interpolation.
    ALLOCATE(vr_3d_before(n_gates,n_radials,n_elevs))
    ALLOCATE(vr_3d_after(n_gates,n_radials,n_elevs))
    ALLOCATE(vr_3d_valid(n_gates,n_radials,n_elevs))

    z_3d_before = missing        !Set all the 3D arrays to "missing data" by default.
    z_3d_after = missing         !We'll fill them up later with good data where it is available.
    vr_3d_before = missing
    vr_3d_after = missing
    z_3d_valid = missing
    vr_3d_valid = missing

!------------------------------------------------------------------------
! Now we need to read in data from the NetCDF files.  Oh, joy.  Start by
! determining what window viable data might be in.
!------------------------------------------------------------------------

!What is the first index we're interested in?  Loop through time to find out!
    PRINT*, "Commencing search for data within window."

    cur_index = index_reference             !Start from the reference index.
    elapsed_counter = 0                     !Initialize our elpsed time counter
    DO WHILE (elapsed_counter <= dt_assim)  !Data could possibly be needed anytime
                                            !within one cycle length.
       elapsed_counter = elapsed_counter + d_times(cur_index - 1)
                                            !Start accumulating elapsed time.
       cur_index = cur_index - 1            !And continue doing so backwards
                                            !until you hit the start of the window.
    END DO
    begin_index = cur_index - 1             !We won't search for data more than
                                            !one file outside the window.

    PRINT*, "Earliest index is:", begin_index, "Time to reference is: ", elapsed_counter

! Now repeat this process, going forwards to find the last index we're interested in.
    cur_index = index_reference             !Start from the reference index.
    elapsed_counter = 0                     !Re-initialize our elapsed time counter
    DO WHILE (elapsed_counter <= dt_assim)  !Data could possibly be needed anytime
                                            !within one cycle length.
       elapsed_counter = elapsed_counter + d_times(cur_index)
                                            !Start accumulating elapsed time.
       cur_index = cur_index + 1            !And continue doing so forwards
                                            !until you hit the end of the window.
     END DO
     end_index = cur_index + 1              !We won't search for data more than
                                            !one file outside the window.

     PRINT*, "Latest index is:", end_index, "Time to reference is: ", elapsed_counter

!    PRINT*, "Reference index (immediately before assimilation): ", index_reference

write(0,*) 'loop index_cdf = ',begin_index, end_index

     DO index_cdf = begin_index, end_index, 1   !Loop through all files where there may be data.

        filename = TRIM(files(index_cdf))    !This is the file to use with the NetCDF calls

!------------------------------------------------------------------------
! Call "get_ncraddims", which will tell us the true number of gates and
! radials in the given CDF file.
!------------------------------------------------------------------------

        CALL get_ncraddims(filename,radar_dir,iradtype,n_gates,true_nradials,  &
                           nvar,istatus)
        PRINT*, "Radar File: ", filename
!         PRINT*, "Located in directory: ", radar_dir
!         PRINT*, "iradtype = ", iradtype
!         PRINT*, "n_gates, true_nradials = ", n_gates, true_nradials
!         PRINT*, "nvar = ", nvar, "istatus = ", istatus

!------------------------------------------------------------------------
! Call "get_ncrad2ainfo", which reads some important data about the CASA
! radar location (i.e. lat, lon, elevation), as well as some NetCDF
! variables we may or may not need to know about.
!------------------------------------------------------------------------

        ALLOCATE(tem_double(true_nradials))       !Allocate our temporary variable
        ALLOCATE(ncvarname(nvar))                 !Allocate our variable names array
                                                  !to the number of variables present.
        radar_name_cdf = radar_name               !Set our NetCDF radar name to
                                                  !match the correct radar name.

        CALL get_ncrad2ainfo(filename,radar_dir,iradtype,true_nradials,nvar,tem_double,   &
                             radar_name_cdf,radar_lat,radar_lon,radar_altitude,           &
                             itimcdf,frtime,elv,ncvarname,istatus)

        IF ( timeflag ) THEN      !If we just found index_reference
                                            !last time through the loop...
          time_arpsppi = ((itimcdf/2) + (timetmp/2)) + itime1970
                 !...then calculate time_arpsppi from the nearest two files...
!         time_arpsppi = itimcdf                 !THIS ONE FOR DEBUGGING!!!  ~NS 1/10/08
          timeflag = .FALSE.                !...and unset the time flag so we
                                            !don't trigger this code again.
          PRINT*, "Time before:  ", timetmp, "Time after:", itimcdf
          PRINT*, "time_arpsppi: ", time_arpsppi
        END IF

        IF (index_cdf .eq. index_reference) THEN  !If we are at the reference index, then...
          timetmp = itimcdf     !...store the current file time for getting time_arpsppi...
          timeflag = .TRUE.     !...and note this by setting the time flag.
        END IF

!         PRINT*, "Radar lat, lon: ", radar_lat, radar_lon
!         PRINT*, "Radar elevation (altitude): ", radar_altitude
!         PRINT*, "itimcdf (it's back!): ", itimcdf
        PRINT*, "Elevation Angle of this file: ", elv
!         PRINT*, "Variables present: ", ncvarname

       DEALLOCATE(ncvarname)                     !Deallocate this... we don't need it now.
       DEALLOCATE(tem_double)                    !Deallocate this... we don't need it now.

!------------------------------------------------------------------------
! Call "rd2atiltcdf", which reads in, among other things, reflectivity
! and radial velocity from our given file.  This is the important stuff
! right here.
!------------------------------------------------------------------------

       PRINT*, "Allocating arrays for 'rd2atiltcdf'"

       ALLOCATE(gcfst(true_nradials),         STAT = istatus)
       CALL check_alloc_status(istatus,'casa2arps:gcfst')
       ALLOCATE(istrgate(true_nradials),      STAT = istatus)
       CALL check_alloc_status(istatus,'casa2arps:istrgate')
       ALLOCATE(rhohv(n_gates,true_nradials), STAT = istatus)
       CALL check_alloc_status(istatus,'casa2arps:rhohv')
       ALLOCATE(tem2d(n_gates,true_nradials), STAT = istatus)
       CALL check_alloc_status(istatus,'casa2arps:tem2d')
       ALLOCATE(radv(n_gates,true_nradials),  STAT = istatus)
       CALL check_alloc_status(istatus,'casa2arps:radv')
       ALLOCATE(refl(n_gates,true_nradials),  STAT = istatus)
       CALL check_alloc_status(istatus,'casa2arps:refl')
       ALLOCATE(attref(n_gates,true_nradials),STAT = istatus)
       CALL check_alloc_status(istatus,'casa2arps:attref')
       ALLOCATE(gtspc(true_nradials),         STAT = istatus)
       CALL check_alloc_status(istatus,'casa2arps:gtspc')
       ALLOCATE(azim(true_nradials),          STAT = istatus)
       CALL check_alloc_status(istatus,'casa2arps:azim')
       ALLOCATE(gateqc(n_gates,true_nradials),   STAT = istatus)
       CALL check_alloc_status(istatus,'casa2arps:gateqc')
       ALLOCATE(tem_real(n_gates,true_nradials), STAT = istatus)
       CALL check_alloc_status(istatus,'casa2arps:tem_real')
       ALLOCATE(tem_double(true_nradials),    STAT = istatus)
       CALL check_alloc_status(istatus,'casa2arps:tem_double')

       CALL rd2atiltcdf(true_nradials,n_gates,iradtype,filename,radar_dir, &
                        rmisval,rngfval,itimcdf,frtime,initime,         &
                        vnyq,rfirstg,beamwidth,istrgate,                &
                        gcfst,azim,gtspc,attref,refl,radv,rhohv,        &
                        tem2d,tem2d,gateqc,tem_double,tem_real)

       !PRINT*, "Reflectivity and Radial Velocity data have been read in."
       !PRINT*, "Nyquist velocity is: ", vnyq
       !PRINT*, "Azimuth angles (radials) present are: ", azim
       !PRINT*, "For Elevation: ", elv, " and azimuth: ", azim(1)
       !PRINT*, "Radar beamwidth is: ", beamwidth(1)
       !PRINT*, "Reflectivity along the radial is: ", refl(:,1)
       !PRINT*, "Radial velocity along the radial is: ", radv(:,1)

       DEALLOCATE(tem_double)      ! Deallocate immediately those variables that we...
       DEALLOCATE(tem_real)
       DEALLOCATE(gcfst)
       DEALLOCATE(istrgate)        ! ...will not be working with at all.

!--------------------DEBUGGING SPOT PLACEMENT SECTION---------------------!
! This section was added to trace a rotational error in the code.  If the
! code works correctly, this will place a line of strong reflectivity
! directed DUE EAST of the radar.
!         DO j= 1,sector_sizes(index_cdf), 1
!            IF((azim(j) > 88) .and. (azim(j) < 92)) THEN
!               refl(:,j) = 65.1
!            END IF
!         END DO
!-------------ROTATIONAL ERROR IS LOCATED AFTER THIS POINT----------------!
!-------------------------------------------------------------------------!

!------------------------------------------------------------------------
! Okay, the data has been read in, but the radials are not uniformly
! spaced.  This is a problem for our time interpolation algorithm.
! We'll need to do linear interpolation to get data on evenly-spaced
! radials that we can store in our 3D array.  Depending on our choice of
! max_radials (default:360), we can either oversample or smooth our data.
!------------------------------------------------------------------------

       ALLOCATE(refl_on_integer(n_gates,sector_sizes(index_cdf)))  !Allocate arrays to store the interpolated...
       ALLOCATE(radv_on_integer(n_gates,sector_sizes(index_cdf)))  !...Z and Vr data.
       ALLOCATE(int_azims(sector_sizes(index_cdf))) !Also, an array so we know what azimuths we're interpolating to.

       int_azims(1) = first_radials(index_cdf)      !Begin from the first radial.
       j_int_azim = 1                                     !We'll loop through all the radials we desire data on.
       i_obs_azim = 1                                     !Start with the first observed azimuth, and go from there.

       write(0,*) 'loop index_cdf  = ',index_cdf
!------------------------------------------------------------------------
! We need to find out if the radar was scanning forwards (azimuth angle
! increasing) or backwards (azimuth angle decreasing), and orient our
! interpolation loop accordingly.  There are four cases to consider:
!     1) Scanning forwards, not crossing zero.
!     2) Scanning forwards, crossing zero.
!     3) Scanning backwards, not crossing zero.
!     4) Scanning backwards, crossing zero.
! Each of these cases must be treated differently.
!------------------------------------------------------------------------

         IF(azim(7) < azim(8)) THEN                         !If this is the case, we're scanning forward.
            forward_flag = .TRUE.
            IF(azim(1) > real(first_radials(index_cdf))) THEN !If the first observed radial is already after our first desired radial...
               refl_on_integer(:,1) = refl(:,1)                     !...then assume that the first observed radial is valid...
               radv_on_integer(:,1) = radv(:,1)                     !...for both reflectivity and radial velocity.
            END IF
            DO j_int_azim = 2,sector_sizes(index_cdf),1  !Loop through and fill out our array of desired azimuths.
               int_azims(j_int_azim) = first_radials(index_cdf) + (j_int_azim - 1)
               IF (int_azims(j_int_azim) > 360) THEN
                  int_azims(j_int_azim) = int_azims(j_int_azim) - 360   !Can't have azimuth values above 360... we must have crossed zero.
               END IF
            END DO

            j_int_azim = 1  !Now reset j_int_azim so we can use it again in the next loop

!-----------------------------------Debug PRINT Block #50010------------------------------!
!            PRINT*, "Desired radials are: ", int_azims
!            PRINT*, "Sector size, j_int_azim: ", sector_sizes(index_cdf), j_int_azim
!-----------------------------------------------------------------------------------------!

            DO WHILE (j_int_azim < sector_sizes(index_cdf))  !loop forward until we reach the end of our sector.
write(0,*) 'loop j_int_azim = ',j_int_azim
               !If the below statement is true, then we've found two radials to linearly interpolate between.
               IF((azim(i_obs_azim) < int_azims(j_int_azim)) .and. (azim(i_obs_azim + 1) > int_azims(j_int_azim))) THEN
!               PRINT*, "Found one!  Observed azims: ", azim(i_obs_azim), azim(i_obs_azim + 1), "Desired azim: ", int_azims(j_int_azim)
                  dazim_left = int_azims(j_int_azim) - azim(i_obs_azim)     !Calculate how far away the observation is on the left...
                  dazim_right = azim(i_obs_azim + 1) - int_azims(j_int_azim)    !...and on the right for linear interpolation.
                  DO k_cur_gate = 1,n_gates,1          !For each gate in the radial, perform the linear interpolation.
                     IF ((refl(k_cur_gate,i_obs_azim) .ne. rmisval) .and. (refl(k_cur_gate,i_obs_azim + 1) .ne. rmisval)) THEN
                        refl_on_integer(k_cur_gate,j_int_azim) = ((dazim_left*refl(k_cur_gate,i_obs_azim + 1)   &
                                     + dazim_right*refl(k_cur_gate,i_obs_azim)) / (dazim_left + dazim_right))
                        IF (abs((refl(k_cur_gate,i_obs_azim) - refl(k_cur_gate,i_obs_azim + 1))) > 100) THEN
                           refl_on_integer(k_cur_gate,j_int_azim) = bad_data  !If the dBZ values aren't within 100dBZ of one another, flag it.
                        END IF
                     ELSE IF (refl(k_cur_gate,i_obs_azim) .eq. rmisval) THEN  !missing data on the left, use right side.
                        refl_on_integer(k_cur_gate,j_int_azim) = refl(k_cur_gate,i_obs_azim + 1)
                     ELSE                                                     !missing data on the right, use left side.
                        refl_on_integer(k_cur_gate,j_int_azim) = refl(k_cur_gate,i_obs_azim)
                     END IF

                     IF ((radv(k_cur_gate,i_obs_azim) .ne. rmisval) .and. (radv(k_cur_gate,i_obs_azim + 1) .ne. rmisval)) THEN
                        radv_on_integer(k_cur_gate,j_int_azim) = ((dazim_left*radv(k_cur_gate,i_obs_azim + 1)  &
                                     + dazim_right*radv(k_cur_gate,i_obs_azim)) / (dazim_left + dazim_right))
                        IF (abs((radv(k_cur_gate,i_obs_azim) - radv(k_cur_gate,i_obs_azim + 1))) > 100) THEN
                           radv_on_integer(k_cur_gate,j_int_azim) = bad_data  !If the Vr values aren't within 100 of one another, flag it.
                        END IF
                     ELSE IF (radv(k_cur_gate,i_obs_azim) .eq. rmisval) THEN  !missing data on the left, use right side.
                        radv_on_integer(k_cur_gate,j_int_azim) = radv(k_cur_gate,i_obs_azim + 1)
                     ELSE                                                     !missing data on the right, use left side.
                        radv_on_integer(k_cur_gate,j_int_azim) = radv(k_cur_gate,i_obs_azim)
                     END IF
                  END DO
                  j_int_azim = j_int_azim + 1            !Increment j_int_azim and interpolate for the next one.
               ELSE IF ((azim(i_obs_azim) < int_azims(j_int_azim)) .and. (azim(i_obs_azim + 1) < int_azims(j_int_azim))) THEN
                  i_obs_azim = i_obs_azim + 1            !Increment i_obs_azim and try it again.
!                  PRINT*, "Nope, not this one... trying azims: ", azim(i_obs_azim), azim(i_obs_azim + 1)
               ELSE   !This means we must be starting *NOT* between two observed radials.  We need to fix this...
                  refl_on_integer(:,j_int_azim) = refl(:,i_obs_azim)  !...by using the nearest data for reflectivity...
                  radv_on_integer(:,j_int_azim) = radv(:,i_obs_azim)  !...and radial velocity...
                  j_int_azim = j_int_azim + 1 !...and incrementing j_int_azim so that we should be fine now.
!                  PRINT*, "Starting *not* between two observed radials -- assigning nearest values and incrementing desired radial."
               END IF
               !Now, add code to deal with a possible zero crossing.
               IF ((int_azims(j_int_azim) > 340) .and. (azim(i_obs_azim + 1) < 20)) THEN  !If this is true, we're about to cross zero.
!                  PRINT*, "Crossed zero!  Azims: ", azim(i_obs_azim), azim(i_obs_azim + 1), "Desired azim: ", int_azims(j_int_azim)
                  dazim_left = int_azims(j_int_azim) - azim(i_obs_azim)     !Calculate how far away the observation is on the left...
                  dazim_right = azim(i_obs_azim + 1) - (int_azims(j_int_azim) + 360)   !...and on the right for linear interpolation.
                  DO k_cur_gate = 1,n_gates,1          !For each gate in the radial, perform the linear interpolation.
                     IF ((refl(k_cur_gate,i_obs_azim) .ne. rmisval) .and. (refl(k_cur_gate,i_obs_azim + 1) .ne. rmisval)) THEN
                        refl_on_integer(k_cur_gate,j_int_azim) = ((dazim_left*refl(k_cur_gate,i_obs_azim + 1)  &
                                     + dazim_right*refl(k_cur_gate,i_obs_azim)) / (dazim_left + dazim_right))
                        IF (abs((refl(k_cur_gate,i_obs_azim) - refl(k_cur_gate,i_obs_azim + 1))) > 100) THEN
                           refl_on_integer(k_cur_gate,j_int_azim) = bad_data  !If the dBZ values aren't within 100dBZ of one another, flag it.
                        END IF
                     ELSE IF (refl(k_cur_gate,i_obs_azim) .eq. rmisval) THEN  !missing data on the left, use right side.
                        refl_on_integer(k_cur_gate,j_int_azim) = refl(k_cur_gate,i_obs_azim + 1)
                     ELSE                                                     !missing data on the right, use left side.
                        refl_on_integer(k_cur_gate,j_int_azim) = refl(k_cur_gate,i_obs_azim)
                     END IF

                     IF ((radv(k_cur_gate,i_obs_azim) .ne. rmisval) .and. (radv(k_cur_gate,i_obs_azim + 1) .ne. rmisval)) THEN
                        radv_on_integer(k_cur_gate,j_int_azim) = ((dazim_left*radv(k_cur_gate,i_obs_azim + 1)  &
                                     + dazim_right*radv(k_cur_gate,i_obs_azim)) / (dazim_left + dazim_right))
                        IF (abs((radv(k_cur_gate,i_obs_azim) - radv(k_cur_gate,i_obs_azim + 1))) > 100) THEN
                           radv_on_integer(k_cur_gate,j_int_azim) = bad_data  !If the Vr values aren't within 100 of one another, flag it.
                        END IF
                     ELSE IF (radv(k_cur_gate,i_obs_azim) .eq. rmisval) THEN  !missing data on the left, use right side.
                        radv_on_integer(k_cur_gate,j_int_azim) = radv(k_cur_gate,i_obs_azim + 1)
                     ELSE                                                     !missing data on the right, use left side.
                        radv_on_integer(k_cur_gate,j_int_azim) = radv(k_cur_gate,i_obs_azim)
                     END IF
                  END DO
                  j_int_azim = j_int_azim + 1    !Increment both the observed and desired...
                  i_obs_azim = i_obs_azim + 1    !...radial so that we're past zero.
                  IF (azim(i_obs_azim) > int_azims(j_int_azim)) THEN      !If this put us so that we're not between two observed radials...
                     refl_on_integer(:,j_int_azim) = refl(:,i_obs_azim)  !...then use the nearest data for reflectivity...
                     radv_on_integer(:,j_int_azim) = radv(:,i_obs_azim)  !...and radial velocity...
                     j_int_azim = j_int_azim + 1 !...and increment j_int_azim so that we should be fine now.
                  END IF
               END IF
            END DO

!------------------------------------------------------------------------
! This code block is similar to above, but handles scanning backwards.
!------------------------------------------------------------------------

         ELSE                                               !If we're not scanning forward, we're scanning backwards.
            forward_flag = .FALSE.                          !So set the "scanning forwards" flag to false.
            IF(azim(1) < real(first_radials(index_cdf))) THEN !If the first observed radial is already before our first desired radial...
               refl_on_integer(:,1) = refl(:,1)                     !...then assume that the first observed radial is valid.
               radv_on_integer(:,1) = radv(:,1)                     !...for both reflectivity and radial velocity.
            END IF
            DO j_int_azim = 2,sector_sizes(index_cdf),1  !Loop through and fill out our array of desired azimuths.
               int_azims(j_int_azim) = first_radials(index_cdf) - (j_int_azim - 1)
               IF (int_azims(j_int_azim) < 1) THEN
                  int_azims(j_int_azim) = int_azims(j_int_azim) + 360   !Can't have azimuth values below 1... we must have crossed zero.
               END IF
            END DO

            j_int_azim = 1  !Now reset j_int_azim so we can use it again in the next loop

!-----------------------------Debug PRINT Block #73072-----------------------------------!
!            PRINT*, "Desired radials are: ", int_azims
!            PRINT*, "Sector size, j_int_azim: ", sector_sizes(index_cdf), j_int_azim
!----------------------------------------------------------------------------------------!

            DO WHILE (j_int_azim < sector_sizes(index_cdf))
                                   !loop forward until we reach the end of our sector.
                                   !If the below statement is true, then we've found
                                   !two radials to linearly interpolate between.
write(0,*) 'j_int_azim2 = ',j_int_azim
!               i_obs_azim = MIN(i_obs_azim,true_nradials-1)
               IF((azim(i_obs_azim) > int_azims(j_int_azim)) .and.      &
                  (azim(i_obs_azim + 1) < int_azims(j_int_azim))) THEN
!                  PRINT*, "Found one!  Observed azims: ", azim(i_obs_azim), &
!                          azim(i_obs_azim + 1), "Desired azim: ", int_azims(j_int_azim)
                  dazim_left = int_azims(j_int_azim) - azim(i_obs_azim + 1)
                                   !Calculate how far away the observation is on the left...
                  dazim_right = azim(i_obs_azim) - int_azims(j_int_azim)
                                   !...and on the right for linear interpolation.
                  DO k_cur_gate = 1,n_gates,1
                                   !For each gate in the radial, perform the linear interpolation.
                     IF ((refl(k_cur_gate,i_obs_azim + 1) .ne. rmisval) .and. &
                         (refl(k_cur_gate,i_obs_azim) .ne. rmisval)) THEN
                        refl_on_integer(k_cur_gate,j_int_azim) =        &
                               ((dazim_left*refl(k_cur_gate,i_obs_azim) &
                         + dazim_right*refl(k_cur_gate,i_obs_azim + 1)) &
                               / (dazim_left + dazim_right))
                        IF (abs((refl(k_cur_gate,i_obs_azim) -          &
                                 refl(k_cur_gate,i_obs_azim + 1))) > 100) THEN
                           refl_on_integer(k_cur_gate,j_int_azim) = bad_data
                             !If the dBZ values aren't within 100dBZ of one another, flag it.
                        END IF
                     ELSE IF (refl(k_cur_gate,i_obs_azim + 1) .eq. rmisval) THEN  !missing data on the left, use right side.
                        refl_on_integer(k_cur_gate,j_int_azim) = refl(k_cur_gate,i_obs_azim)
                     ELSE                                                     !missing data on the right, use left side.
                        refl_on_integer(k_cur_gate,j_int_azim) = refl(k_cur_gate,i_obs_azim + 1)
                     END IF

                     IF ((radv(k_cur_gate,i_obs_azim + 1) .ne. rmisval) .and. (radv(k_cur_gate,i_obs_azim) .ne. rmisval)) THEN
                        radv_on_integer(k_cur_gate,j_int_azim) = ((dazim_left*radv(k_cur_gate,i_obs_azim)  &
                             + dazim_right*radv(k_cur_gate,i_obs_azim + 1)) / (dazim_left + dazim_right))
                        IF (abs((radv(k_cur_gate,i_obs_azim) - radv(k_cur_gate,i_obs_azim + 1))) > 100) THEN
                           radv_on_integer(k_cur_gate,j_int_azim) = bad_data  !If the Vr values aren't within 100 of one another, flag it.
                        END IF
                     ELSE IF (radv(k_cur_gate,i_obs_azim + 1) .eq. rmisval) THEN  !missing data on the left, use right side.
                        radv_on_integer(k_cur_gate,j_int_azim) = radv(k_cur_gate,i_obs_azim)
                     ELSE                                                     !missing data on the right, use left side.
                        radv_on_integer(k_cur_gate,j_int_azim) = radv(k_cur_gate,i_obs_azim + 1)
                     END IF
                  END DO
                  j_int_azim = j_int_azim + 1            !Increment j_int_azim and interpolate for the next one.
               ELSE IF ((azim(i_obs_azim) > int_azims(j_int_azim)) .and. (azim(i_obs_azim + 1) > int_azims(j_int_azim))) THEN
                  i_obs_azim = i_obs_azim + 1            !Increment i_obs_azim and try it again.
!                  PRINT*, "Nope, not this one... trying azims: ", azim(i_obs_azim), azim(i_obs_azim + 1)
               ELSE IF (int_azims(j_int_azim) .eq. 360) THEN  !We don't want to trigger the starting criteria (below) for a zero crossing.
!                  PRINT*, "Approaching zero crossing..."
                  i_obs_azim = i_obs_azim + 1
!                  PRINT*, "Nope, not this one... trying azims: ", azim(i_obs_azim), azim(i_obs_azim + 1)
               ELSE   !This means we must be starting *NOT* between two observed radials.  We need to fix this...
                  refl_on_integer(:,j_int_azim) = refl(:,i_obs_azim)  !...by using the nearest data for reflectivity...
                  radv_on_integer(:,j_int_azim) = radv(:,i_obs_azim)  !...and radial velocity...
                  j_int_azim = j_int_azim + 1 !...and incrementing j_int_azim so that we should be fine now.
!                  PRINT*, "Starting *not* between two observed radials -- assigning nearest values and incrementing desired radial."
               END IF
               !Now, add code to deal with a possible zero crossing.
               !If the following is true, we're about to cross zero.
               IF ((int_azims(j_int_azim) > 340) .and. (azim(i_obs_azim) < 20) .and. (azim(i_obs_azim + 1) > 340)) THEN
!                  PRINT*, "Crossed zero!  Azims: ", azim(i_obs_azim), azim(i_obs_azim + 1), "Desired azim: ", int_azims(j_int_azim)
                  dazim_left = int_azims(j_int_azim) - azim(i_obs_azim + 1) !Calculate how far away the observation is on the left...
                  dazim_right = azim(i_obs_azim) - (int_azims(j_int_azim) + 360)    !...and on the right for linear interpolation.
                  DO k_cur_gate = 1,n_gates,1          !For each gate in the radial, perform the linear interpolation.
                     IF ((refl(k_cur_gate,i_obs_azim + 1) .ne. rmisval) .and. (refl(k_cur_gate,i_obs_azim) .ne. rmisval)) THEN
                        refl_on_integer(k_cur_gate,j_int_azim) = ((dazim_left*refl(k_cur_gate,i_obs_azim)    &
                             + dazim_right*refl(k_cur_gate,i_obs_azim + 1)) / (dazim_left + dazim_right))
                        IF (abs((refl(k_cur_gate,i_obs_azim) - refl(k_cur_gate,i_obs_azim + 1))) > 100) THEN
                           refl_on_integer(k_cur_gate,j_int_azim) = bad_data  !If the dBZ values aren't within 100dBZ of one another, flag it.
                        END IF
                     ELSE IF (refl(k_cur_gate,i_obs_azim + 1) .eq. rmisval) THEN  !missing data on the left, use right side.
                        refl_on_integer(k_cur_gate,j_int_azim) = refl(k_cur_gate,i_obs_azim)
                     ELSE                                                     !missing data on the right, use left side.
                        refl_on_integer(k_cur_gate,j_int_azim) = refl(k_cur_gate,i_obs_azim + 1)
                     END IF

                     IF ((radv(k_cur_gate,i_obs_azim + 1) .ne. rmisval) .and. (radv(k_cur_gate,i_obs_azim) .ne. rmisval)) THEN
                        radv_on_integer(k_cur_gate,j_int_azim) = ((dazim_left*radv(k_cur_gate,i_obs_azim)    &
                             + dazim_right*radv(k_cur_gate,i_obs_azim + 1)) / (dazim_left + dazim_right))
                        IF (abs((radv(k_cur_gate,i_obs_azim) - radv(k_cur_gate,i_obs_azim + 1))) > 100) THEN
                           radv_on_integer(k_cur_gate,j_int_azim) = bad_data  !If the Vr values aren't within 100 of one another, flag it.
                        END IF
                     ELSE IF (radv(k_cur_gate,i_obs_azim + 1) .eq. rmisval) THEN  !missing data on the left, use right side.
                        radv_on_integer(k_cur_gate,j_int_azim) = radv(k_cur_gate,i_obs_azim)
                     ELSE                                                     !missing data on the right, use left side.
                        radv_on_integer(k_cur_gate,j_int_azim) = radv(k_cur_gate,i_obs_azim + 1)
                     END IF
                  END DO
                  j_int_azim = j_int_azim + 1    !Increment both the observed and desired...
                  i_obs_azim = i_obs_azim + 1    !...radial so that we're past zero.
                  IF (azim(i_obs_azim) < int_azims(j_int_azim)) THEN      !If this put us so that we're not between two observed radials...
                     refl_on_integer(:,j_int_azim) = refl(:,i_obs_azim)  !...then use the nearest data for reflectivity...
                     radv_on_integer(:,j_int_azim) = radv(:,i_obs_azim)  !...and radial velocity...
                     j_int_azim = j_int_azim + 1 !...and increment j_int_azim so that we should be fine now.
                  END IF
               END IF
            END DO
         END IF

!--------------------DEBUGGING SPOT PLACEMENT SECTION---------------------!
!      DO j = 89,91,1
!         refl_on_integer(:,j) = 65.1
!      END DO
!-------------ROTATIONAL ERROR IS LOCATED BEFORE THIS POINT---------------!
!-------------------------------------------------------------------------!

!--------------------------------Yet another debugging print block-------------------------!
!         PRINT*, "For third desired radial, ", int_azims(3), "Reflectivity data are:"
!         PRINT*, refl_on_integer(:,3)
!         PRINT*, "Radial Velocity data are:"
!         PRINT*, radv_on_integer(:,3)
!
!         PRINT*, "For comparison, third observed radial, ", azim(3), "Reflectivity:"
!         PRINT*, refl(:,3)
!         PRINT*, "Observed Radial Velocity:"
!         PRINT*, radv(:,3)
!------------------------------------------------------------------------------------------!


!------------------------------------------------------------------------
! There; now the data is on uniform radials.  Now, our next question is:
! What part of this data (if any) do we need to store in either the
! "before" or "after" 3D array?  Begin by taking the uniform radial data
! and placing it on a full 360 degree elevation scan for easy access.
!------------------------------------------------------------------------

         ALLOCATE(refl_integer_360(n_gates,n_radials)) !Allocate for the full 360 scans...
         ALLOCATE(radv_integer_360(n_gates,n_radials)) !...in both Z and Vr.

         start_radial = int_azims(1)                           !Set our first guess for the first...
         end_radial = int_azims(sector_sizes((index_cdf)))     !...and last radials present for looping purposes.

         refl_integer_360 = missing                !Initialize the entire 360 scan arrays...
         radv_integer_360 = missing                !...using the 'missing data' value

         IF ( forward_flag ) THEN        !If we are scanning forward, then...
            scan_increment = delta_azim            !...increment forward through radials.
            IF (start_radial > end_radial) THEN    !If we crossed zero...
               end_radial = end_radial + 360       !...adjust the end radial to make a loop possible.
            END IF
         ELSE                                      !Otherwise we're scanning backward...
            scan_increment = 0 - delta_azim        !...so increment backward through radials.
            IF (start_radial < end_radial) THEN    !If we crossed zero...
               end_radial = end_radial - 360       !...adjust the end radial to make a loop possible.
            END IF

         END IF

         PRINT*, "Filling 360 array with the data..."
!         DO count = start_radial, end_radial, scan_increment
         count = start_radial                      !We want to make sure to start at the right location when assigning data.
         DO k = 1, sector_sizes(index_cdf), NINT(delta_azim)
            IF (count < 1) THEN                    !We might need to add 360 if we crossed zero going backwards.
               refl_integer_360(:,(count+360)) = refl_on_integer(:,(k))   !Now scan through the interpolated data...
               radv_integer_360(:,(count+360)) = radv_on_integer(:,(k))   !... and place it into the 360 degree arrays.
            ELSE IF (count > 360) THEN             !We might need to subtract 360 if we crossed zero going forwards.
               refl_integer_360(:,(count-360)) = refl_on_integer(:,(k))   !Now scan through the interpolated data...
               radv_integer_360(:,(count-360)) = radv_on_integer(:,(k))   !... and place it into the 360 degree arrays.
            ELSE                                   !Or we might not have crossed zero at all.
               refl_integer_360(:,count) = refl_on_integer(:,k)   !Now scan through the interpolated data...
               radv_integer_360(:,count) = radv_on_integer(:,k)   !... and place it into the 360 degree arrays.
            END IF
            count = count + scan_increment         !Move on to the next radial
         END DO

!------------------------------------------------------------------------
! Now, check with the time arrays to determine where our current
! elevation scan is needed in either the before or after 3D array of
! Z and Vr data.
!------------------------------------------------------------------------

         SELECT CASE(INT(elv))        !Determine which index our elevation corresponds to.
         CASE(1)
            elev_index = 1
         CASE(2)
            elev_index = 2
         CASE(3)
            elev_index = 3
         CASE(5)
            elev_index = 4
         CASE(7)
            elev_index = 5
         CASE(9)
            elev_index = 6
         CASE(11)
            elev_index = 7
         CASE(14)
            elev_index = 8
         CASE DEFAULT                 !If it doesn't match any elevation angle... we've got a problem.
            STOP 'Elevation index doex not match given values! -- Please check the listfile!'
         END SELECT

         DO count = 1,n_radials,1     !For each radial, loop through and see if this data is needed.
            IF (id_before(count,elev_index) .eq. index_cdf) THEN   !If the data before is needed, then...
               z_3d_before(:,count,elev_index) = refl_integer_360(:,count)      !...plop it into the 3D array.
               vr_3d_before(:,count,elev_index) = radv_integer_360(:,count)
               PRINT*, "***Data before filled at time: ", int_times(index_cdf), "(Index is: ", index_cdf, ")"
               PRINT*, "***Radial: ", count, "Elevation: ", elv
            END IF
            IF (id_after(count,elev_index) .eq. index_cdf) THEN   !If the data after is needed, then...
               z_3d_after(:,count,elev_index) = refl_integer_360(:,count)      !...plop it into the 3D array.
               vr_3d_after(:,count,elev_index) = radv_integer_360(:,count)
               PRINT*, "***Data after filled at time: ", int_times(index_cdf), "(Index is: ", index_cdf, ")"
               PRINT*, "***Radial: ", count, "Elevation: ", elv
            END IF
         END DO

         gatesp = gtspc(1)  !Set gate spacing to the value we read in from the NetCDF file.

         DEALLOCATE(azim)             !Deallocate the other seven NetCDF variables so that...
         DEALLOCATE(gtspc)            !...we can reallocate them in the next trip...
                                      !...through the loop for the next radar datafile.
         DEALLOCATE(attref)
         DEALLOCATE(refl)
         DEALLOCATE(radv)
         DEALLOCATE(rhohv)
         DEALLOCATE(tem2d)

         DEALLOCATE(refl_on_integer)  !Deallocate local arrays will be redefined for the next file.
         DEALLOCATE(radv_on_integer)
         DEALLOCATE(int_azims)
         DEALLOCATE(refl_integer_360)
         DEALLOCATE(radv_integer_360)


     END DO   !Finished looping through all times when data could be present for this assimilation cycle.
!------------------------------------------------------------------------
! The 3D arrays have been filled!  Now, perform linear interpolation in
! time to get a single 3D volume valid at the assimilation time for both
! Z and Vr.
!------------------------------------------------------------------------

      DO k = 1,n_elevs,1         !We'll need to perform this interpolation in a...
         DO j = 1,n_radials,1    !...triply-nested loop over elevations, radials, and gates...
            DO i = 1,n_gates,1   !...to get our 3D volume valid at the assimilation time.

               IF ((z_3d_before(i,j,k) > -66666.) .and. (z_3d_after(i,j,k) > -66666.)) THEN
                  z_3d_valid(i,j,k) = ((elapsed_before(j,k)*z_3d_after(i,j,k)                &
                         + elapsed_after(j,k)*z_3d_before(i,j,k)) / (elapsed_before(j,k) + elapsed_after(j,k)))
                  IF (abs((z_3d_before(i,j,k) - z_3d_after(i,j,k))) > 100) THEN
                     z_3d_valid(i,j,k) = bad_data  !If the dBZ values aren't within 100dBZ of one another, flag it.
                  END IF
               ELSE IF ((z_3d_before(i,j,k) < -66666.) .and. (z_3d_after(i,j,k) > -66666.)) THEN  !missing data before, use data after.
                  z_3d_valid(i,j,k) = z_3d_after(i,j,k)
               ELSE IF ((z_3d_after(i,j,k) < -66666.) .and. (z_3d_before(i,j,k) > -66666.)) THEN  !missing data after, use data before.
                  z_3d_valid(i,j,k) = z_3d_before(i,j,k)
               END IF

               IF ((vr_3d_before(i,j,k) > -66666.) .and. (vr_3d_after(i,j,k) > -66666.)) THEN
                  vr_3d_valid(i,j,k) = ((elapsed_before(j,k)*vr_3d_after(i,j,k)              &
                        + elapsed_after(j,k)*vr_3d_before(i,j,k)) / (elapsed_before(j,k) + elapsed_after(j,k)))
                  IF (abs((vr_3d_before(i,j,k) - vr_3d_after(i,j,k))) > 100) THEN
                     vr_3d_valid(i,j,k) = bad_data  !If the Vr values aren't within 100 of one another, flag it.
                  END IF
               ELSE IF ((vr_3d_before(i,j,k) < -66666.) .and. (vr_3d_after(i,j,k) > -66666.)) THEN  !missing data before, use data after.
                  vr_3d_valid(i,j,k) = vr_3d_after(i,j,k)
               ELSE IF ((vr_3d_after(i,j,k) < -66666.) .and. (vr_3d_before(i,j,k) > -66666.)) THEN  !missing data after, use data before.
                  vr_3d_valid(i,j,k) = vr_3d_before(i,j,k)
               END IF

            END DO
         END DO
      END DO

!--------------------DEBUGGING SPOT PLACEMENT SECTION---------------------!
!      DO j = 89,91,1
!         z_3d_valid(:,j,:) = 65.1
!      END DO
!-------------------------------------------------------------------------!


!------------------------------------------------------------------------
! Linear interpolation in time is complete and we have 3D Z and Vr
! arrays valid at our analysis time.  Now we need to determine where the
! radar is, where that falls on our ARPS grid, and how to put our 3D
! arrays onto the ARPS grid in a .gridtilt (arpsppi) format.
!------------------------------------------------------------------------


!---------------------------Debugging print block!----------------------------!
!      DO k = 1, n_elevs, 1
!         DO j = 3, n_radials, 3
!               PRINT'(a1,i3,a1,i3,a1,i3,a1,i3,a2,3x,10f10.2)', '(',n_gates-30,'-',n_gates-20,',',j,',',k,'):', z_3d_valid(n_gates-30:n_gates-20,j,k)
!         END DO
!      END DO
!      STOP
!-----------------------------------------------------------------------------!


      PRINT*, "Calling 'initpara'."

      namelist_filename = ' '
      CALL initpara(nx, ny, nz, nzsoil, nstypes, namelist_filename)  !Call 'initpara' to read in ARPS domain dimensions.

      PRINT*, "Back from 'initpara'."

      ALLOCATE(x(nx))              !Allocate ARPS variables needed...
      ALLOCATE(y(ny))              !...for calling 'inigrid'.  'inigrid'...
      ALLOCATE(z(nz))              !...will give us the information we need...
      ALLOCATE(zp(nx,ny,nz))             !...to call radcoord and locate the radar...
      ALLOCATE(zpsoil(nx,ny,nzsoil))     !...on the ARPS model grid.
      ALLOCATE(h_terrain(nx,ny))
      ALLOCATE(mapfct(nx,ny,8))
      ALLOCATE(j1(nx,ny,nz))
      ALLOCATE(j2(nx,ny,nz))
      ALLOCATE(j3(nx,ny,nz))
      ALLOCATE(j3soil(nx,ny,nzsoil))
      ALLOCATE(j3soilinv(nx,ny,nzsoil))
      ALLOCATE(tem1d1(nz))
      ALLOCATE(tem1d2(nz))
      ALLOCATE(tem3d(nx,ny,nz))

      PRINT*, "Calling 'inigrd'."

      CALL inigrd(nx, ny, nz, nzsoil, x, y, z, zp, zpsoil,          &
                  h_terrain, mapfct, j1, j2, j3, j3soil,           &
                  j3soilinv, tem1d1, tem1d2, tem3d)

      PRINT*, "Back from 'inigrd'."

      PRINT*, "nx, ny, nz:", nx, ny, nz
!      IF (debugflag .eq. .true.) THEN
!         STOP
!      END IF
!      debugflag = .true.

!------------------------------------------------------------------------
! Now, call 'radcoord' to determine where exactly the CASA radar is on
! the ARPS grid.
!------------------------------------------------------------------------

      ALLOCATE(xs(nx))             !Three more ARPS variables are needed for this call.  Allocate them.
      ALLOCATE(ys(ny))
      ALLOCATE(zps(nx,ny,nz))

      PRINT*, "Calling 'radcoord'."

      CALL radcoord(nx, ny, nz, x, y, z, zp, xs, ys, zps,            &
                    radar_lat, radar_lon, radar_x, radar_y)

      PRINT*, "Back from 'radcoord'."

!------------------------------------------------------------------------
! To put the radar data onto the ARPS grid, we'll need to call the
! subroutine 'remap_arpsppi'.  It puts the data onto the x-y ARPS grid
! but leaves data on radar elevation levels, which is what we want.
!------------------------------------------------------------------------

      ALLOCATE(kntazm(n_elevs))           !Allocate 'remap_arpsppi' variables.
      ALLOCATE(kntgat_ref(n_radials,n_elevs))
      ALLOCATE(kntgat_vel(n_radials,n_elevs))

      ALLOCATE(rngvol(n_gates,n_elevs))
      ALLOCATE(azmvol(n_radials,n_elevs))
      ALLOCATE(elvvol(n_radials,n_elevs))

      ALLOCATE(gridtilthigh_ref(nx,ny,n_elevs))
      ALLOCATE(gridrange_ref(nx,ny,n_elevs))
      ALLOCATE(gridslr_ref(nx,ny))
      ALLOCATE(gridazm_ref(nx,ny))
      ALLOCATE(gridtiltval_ref(nx,ny,n_elevs))
      ALLOCATE(gridtilttime_ref(n_elevs))
      ALLOCATE(elevmean_ref(n_elevs))

      ALLOCATE(gridtilthigh_vel(nx,ny,n_elevs))
      ALLOCATE(gridrange_vel(nx,ny,n_elevs))
      ALLOCATE(gridslr_vel(nx,ny))
      ALLOCATE(gridazm_vel(nx,ny))
      ALLOCATE(gridtiltval_vel(nx,ny,n_elevs))
      ALLOCATE(gridtilttime_vel(n_elevs))
      ALLOCATE(elevmean_vel(n_elevs))

      PRINT*, "Allocations successful."

      irfirstg = NINT(rfirstg)  !Get the integer equivalent of the range to the first gate.

      PRINT*, "Setting up 'vol' arrays."

      DO k=1,n_elevs,1
         kntazm(k) = n_radials  ! Set up the 'kntazm' array -- we want to
                                ! process 360 degrees of data at each elevation level.
         DO j=1,n_radials,NINT(delta_azim)
            azmvol(j,k)=j       ! Set up the 'azmvol' array -- it contains
                                ! the azimuth angles for each radial and elevation.
            elvvol(j,k)=k       ! Set up the 'elvvol' array -- it contains
                                ! the elevation angles for each radial and elevation.
         END DO
         DO i=1,n_gates,1
            rngvol(i,k)=((i-1)*gatesp) + irfirstg ! Set up the 'rngvol' array --
                                                  ! contains range to each gate in each elevation.
         END DO
      END DO

!------------------------------------------------------------------------
! Finally, we are ready to remap the data into a gridtilt format.  Call
! remap_arpsppi to do this.
!------------------------------------------------------------------------

      PRINT*, "Calling 'remap_arpsppi' for reflectivity."

!First, call 'remap_arpsppi' for reflectivity.
      CALL remap_arpsppi(n_gates,n_radials,n_elevs,nx,ny,                         &
                         kntazm,n_elevs,radar_lat,radar_lon,radar_x,              &
                         radar_y,radar_altitude,delta_azim,range_min,range_max,   &
                         gatesp,azmvol,elvvol,z_3d_valid,time_arpsppi,xs,ys,      &
                         gridtilthigh_ref,gridrange_ref,gridslr_ref,gridazm_ref,  &
                         gridtiltval_ref,gridtilttime_ref,elevmean_ref)

      PRINT*, "Back -- now calling 'remap_arpsppi' for velocity."

!Now, call it again for radial velocity.
      CALL remap_arpsppi(n_gates,n_radials,n_elevs,nx,ny,                         &
                         kntazm,n_elevs,radar_lat,radar_lon,radar_x,              &
                         radar_y,radar_altitude,delta_azim,range_min,range_max,   &
                         gatesp,azmvol,elvvol,vr_3d_valid,time_arpsppi,xs,ys,     &
                         gridtilthigh_vel,gridrange_vel,gridslr_vel,gridazm_ref,  &
                         gridtiltval_vel,gridtilttime_vel,elevmean_vel)

!-------------------------DEBUG PRINT BLOCK, YET AGAIN!------------------------!
!      DO k = 1, n_elevs, 1
!         PRINT*, 'For elevation: ', k
!         PRINT*, gridtiltval_ref(:,:,k)
!         PRINT*, '-------------------------------------------------------'
!      END DO
!      STOP
!------------------------------------------------------------------------------!

!------------------------------------------------------------------------
! We now have the arrays that 'wrttilts' and 'wrtgridtilt' need to
! produce '.tilts' and '.gridtilt' files -- the latter of which are the
! ones needed by ARPS EnKF for data assimilation.
!------------------------------------------------------------------------

      PRINT*, "Calling 'abss2ctim'."

!Call 'abss2ctim' to get a unique time for use in the output subroutines.
      CALL abss2ctim(time_arpsppi, iyr, imon, iday, ihr, imin, isec)

!      PRINT*, "Year: ", iyr, "Month:", imon
!      PRINT*, "Day: ", iday, "Hour:", ihr
!      PRINT*, "Min: ", imin, "Sec:", isec

      PRINT*, "Calling 'wrttilts'."

      WRITE(rfname,'(a4,a1,i4.4,2(i2.2),a1,3(i2.2))')                   &
            radar_name,'.', iyr,imon,iday,'.',ihr,imin,isec
      !Set a unique rfname for these output files.

!      PRINT*, "rfname = ", rfname

      !Call 'wrttilts' to write out the interpolated radar data on each radar elevation ('.tilts' files).
      CALL wrttilts(rfname,n_gates,n_gates,n_radials,n_elevs,           &
                    rngvol,azmvol,elvvol,vr_3d_valid,                   &
                    rngvol,azmvol,elvvol,z_3d_valid,                    &
                    radar_altitude,radar_lat,radar_lon,                 &
                    iyr,imon,iday,ihr,imin,isec )

      PRINT*, "Calling 'wrtgridtilt'."
!      PRINT*, "time_arpsppi:", time_arpsppi

      !Call 'wrtgridtilt' to write out the radar data on model gridpoints and radar elevation levels ('.gridtilt' files).
      ! NOTE: need to check for new interface
      CALL wrtgridtilt(rfname,n_elevs,nx,ny,n_elevs,radar_name,                        &
                       radar_lat,radar_lon,radar_x,radar_y,radar_altitude,             &
                       delta_azim,range_min,range_max,                                 &
                       gridtilthigh_ref,gridrange_ref,gridslr_ref,gridazm_ref,         &
                       gridtiltval_ref,gridtilttime_ref,gridtilttime_vel,elevmean_ref, &
                       gridtilthigh_vel,gridrange_vel,gridslr_vel,gridazm_vel,         &
                       gridtiltval_vel,elevmean_vel)

!Call 'write_arpsppi' to write out the radar data directly into the format needed by ppiplt and arpsenkf
!without needing to go through 'tintp'.

      elevs_real = real(elevs)     !Map the elevations from INTEGER to REAL for use in 'write_arpsppi'

      PRINT*, "Calling 'write_arpsppi'"

      CALL write_arpsppi(time_arpsppi,iyr,imon,iday,ihr,imin,isec,cur_assim, &
                         n_elevs,nx,ny,radar_name,radar_lat,radar_lon,       &
                         radar_x,radar_y,radar_altitude,delta_azim,          &
                         range_min,range_max,elevs_real,                     &
                         gridtilthigh_vel,gridrange_vel,                     &
                         gridtiltval_vel,gridtiltval_ref)

!------------------------------------------------------------------------
! This completes our data processing and output for this assimilation
! cycle.  Now deallocate all variables that need to be re-used next
! cycle to avoid allocation errors.
!------------------------------------------------------------------------

    DEALLOCATE(kntazm)           !Deallocate 'remap_arpsppi' variables.
    DEALLOCATE(kntgat_ref)
    DEALLOCATE(kntgat_vel)

    DEALLOCATE(rngvol)
    DEALLOCATE(azmvol)
    DEALLOCATE(elvvol)

    DEALLOCATE(gridtilthigh_ref)
    DEALLOCATE(gridrange_ref)
    DEALLOCATE(gridslr_ref)
    DEALLOCATE(gridazm_ref)
    DEALLOCATE(gridtiltval_ref)
    DEALLOCATE(gridtilttime_ref)
    DEALLOCATE(elevmean_ref)

    DEALLOCATE(gridtilthigh_vel)
    DEALLOCATE(gridrange_vel)
    DEALLOCATE(gridslr_vel)
    DEALLOCATE(gridazm_vel)
    DEALLOCATE(gridtiltval_vel)
    DEALLOCATE(gridtilttime_vel)
    DEALLOCATE(elevmean_vel)

    DEALLOCATE(x)                !Deallocate the ARPS variables.
    DEALLOCATE(y)
    DEALLOCATE(z)
    DEALLOCATE(zp)
    DEALLOCATE(zpsoil)
    DEALLOCATE(h_terrain)
    DEALLOCATE(mapfct)
    DEALLOCATE(j1)
    DEALLOCATE(j2)
    DEALLOCATE(j3)
    DEALLOCATE(j3soil)
    DEALLOCATE(j3soilinv)
    DEALLOCATE(tem1d1)
    DEALLOCATE(tem1d2)
    DEALLOCATE(tem3d)
    DEALLOCATE(xs)
    DEALLOCATE(ys)
    DEALLOCATE(zps)

    DEALLOCATE(z_3d_before)      !Deallocate our local 3D arrays for this cycle.
    DEALLOCATE(z_3d_after)
    DEALLOCATE(z_3d_valid)
    DEALLOCATE(vr_3d_before)
    DEALLOCATE(vr_3d_after)
    DEALLOCATE(vr_3d_valid)

!------------------------------------------------------------------------
! Now, the painful process of updating the current assimilation time.
! We're using an hhmmss format (since that's how it was read in from the
! filename, and because we'll be working extensively with filenames later
! on), so we need to convert between seconds and hhmmss.
!------------------------------------------------------------------------

    PRINT*, "+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+"
    PRINT*, "Current hhmmss assimilation time:", cur_assim

    CALL hhmmss_to_sec(cur_assim,cur_assim_sec)   !convert the hhmmss time to seconds since midnight.

    PRINT*, "This is equivalent to ",cur_assim_sec, " seconds since midnight."
    cur_assim_sec = cur_assim_sec + dt_assim      !increment the current time to the next cycle.

    CALL sec_to_hhmmss(cur_assim_sec,cur_assim)   !Convert back to hhmmss for use with filenames.
    PRINT*, "Incrementing assimilation time by ", dt_assim, " seconds"
    PRINT*, "New hhmmss assimilation time:", cur_assim
    PRINT*, "+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+"

  END DO    !End the loop over multiple assimilation times.

  !END DO    !End the loop over multiple radars.

  CLOSE(42)

END PROGRAM casa2arpsppi


!#######################################################################
!
!           SUBROUTINE hhmmss_decompose_integer
!                  Author: Nate Snook
!#######################################################################
!
! A subroutine to take an hhmmss integer time and decompose it into
! individual integer components. Currently not implemented in the
! casa2arpsppi program, but could be useful at a later time.
!-----------------------------------------------------------------------

SUBROUTINE hhmmss_decompose_integer(thetime,hour10,hour1,min10,min1,sec10,sec1)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: thetime      !The integer time to be decomposed
  INTEGER, INTENT(OUT) :: hour10      !Hours tens place
  INTEGER, INTENT(OUT) :: hour1       !Hours ones place
  INTEGER, INTENT(OUT) :: min10       !Minutes tens place
  INTEGER, INTENT(OUT) :: min1        !Minutes ones place
  INTEGER, INTENT(OUT) :: sec10       !Seconds tens place
  INTEGER, INTENT(OUT) :: sec1        !Seconds ones place

  sec1 = MOD(thetime, 10)
  sec10 = (MOD(thetime, 100) - sec1)/10
  min1 = (MOD(thetime,1000) - (10*sec10) - sec1)/100
  min10 = (MOD(thetime,10000) - (100*min1) - (10*sec10) - sec1)/1000
  hour1 = (MOD(thetime,100000) - (1000*min10) - (100*min1) - (10*sec10) - sec1)/10000
  hour10 = (thetime - (10000*hour1) - (1000*min10) - (100*min1) - (10*sec10) - sec1)/100000

  RETURN
END SUBROUTINE hhmmss_decompose_integer

!#######################################################################
!               SUBROUTINE hhmmss_to_sec
!                  Author: Nate Snook
!#######################################################################
!
! A subroutine to take an hhmmss integer time and convert it into
! seconds since midnight.
!-----------------------------------------------------------------------

SUBROUTINE hhmmss_to_sec(time_hhmmss,time_sec)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: time_hhmmss   !The time given in hhmmss
  INTEGER, INTENT(OUT) :: time_sec     !The equivalent in seconds since midnight (hhmmss = 000000)

  INTEGER :: hour10
  INTEGER :: hour1
  INTEGER :: min10
  INTEGER :: min1
  INTEGER :: sec10
  INTEGER :: sec1

  CALL hhmmss_decompose_integer(time_hhmmss,hour10,hour1,min10,min1,sec10,sec1)  !First decompose the time.

  time_sec = sec1 + (10*sec10) + (60*min1) + (600*min10) + (3600*hour1) + (36000*hour10)

  RETURN
END SUBROUTINE hhmmss_to_sec


!#######################################################################
!                    SUBROUTINE sec_to_hhmmss
!                       Author: Nate Snook
!#######################################################################
!
! A subroutine to take an integer time in seconds since midnight and
! convert it to hhmmss format.
!-----------------------------------------------------------------------

SUBROUTINE sec_to_hhmmss(time_sec,time_hhmmss)
  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: time_hhmmss   !The time given in hhmmss
  INTEGER, INTENT(IN) :: time_sec       !The equivalent in seconds since midnight (hhmmss = 000000)

  INTEGER :: hour10                     !These six variables are the same as in the...
  INTEGER :: hour1                      !hhmmss_decompose_integer subroutine.
  INTEGER :: min10
  INTEGER :: min1
  INTEGER :: sec10
  INTEGER :: sec1

  !Now, find the individual elements of the hhmmss time.
  sec1 = MOD(time_sec,10)
  sec10 = (MOD(time_sec,60) - sec1)/10
  min1 = (MOD(time_sec,600) - (10*sec10) - sec1)/60
  min10 = (MOD(time_sec,3600) - (60*min1) - (10*sec10) - sec1)/600
  hour1 = (MOD(time_sec,36000) - (600*min10) - (60*min1) - (10*sec10) - sec1)/3600
  hour10 = (time_sec - (3600*hour1) - (600*min10) - (60*min1) - (10*sec10) - sec1)/36000

  !And recompose them into an hhmmss integer.
  time_hhmmss = (100000*hour10) + (10000*hour1) + (1000*min10) + (100*min1) + (10*sec10) + sec1

  RETURN
END SUBROUTINE sec_to_hhmmss
