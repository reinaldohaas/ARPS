  INCLUDE "netcdf.inc"
!--------------------------------------------------------------
  INTEGER :: me_version_number
  PARAMETER(me_version_number = 1)

!    a standard function "MECheckVersion" determines whether or
!    not a file is a MERAF file.  It returns either ME_NOT_CDF,
!    ME_NOT_MEKA, or the Version number of the MERAF File

  INTEGER :: me_version_check_error,me_not_cdf,me_not_meka
  PARAMETER(me_version_check_error = -1,me_not_cdf= -2)
  PARAMETER(me_not_meka = -3)

!----------------------------------------------------------------
  INTEGER :: me_null_component
  PARAMETER(me_null_component = -1)
!----------------------------------------------------------------
  INTEGER :: me_verbose,me_quiet
  PARAMETER(me_verbose = 1)
  PARAMETER(me_quiet = 0)
!----------------------------------------------------------------

!  ME_READ and ME_WRITE are the two flags currently allowed when opening
! a dataset or scheme to indicate how the files are to be opened...
! there is currently no support for multiple writers on the same file.
  INTEGER :: me_read,me_write
  PARAMETER(me_read = ncnowrit)
  PARAMETER(me_write = ncwrite)

!----------------------------------------------------------------

  INTEGER :: true, false
  PARAMETER(true = 1)
  PARAMETER(false = 0)


!----------------------------------------------------------------

! ME_ALL can be passed by user to any of the Variable GetWrite functions
!    in the "To" place to have the data written to the whole dimension
!    (from the "From" point)
  INTEGER :: me_all
  PARAMETER(me_all = -8)

!----------------------------------------------------------------

! MERAF_HEADER is used for the GetHistory(long Frame,char R_History)
!    in the place of a frame number to indicate the header is
!    to be read and not the datafiles...
!
  INTEGER :: meraf_header
  PARAMETER(meraf_header = -2)
!----------------------------------------------------------------

! ConvertTimeU asks the user for a UTC offset this parameter
!    uses the local timezone for the UTC offset...?

  INTEGER :: melocaltimezone
  PARAMETER(melocaltimezone = 999)

!------------------------------------------------------------------
!    MAX_NAME_LENGTH is used everywhere for string sizes
  INTEGER :: max_path_length,max_name_length
  PARAMETER(max_path_length = 256)
  PARAMETER(max_name_length = 256)
!----------------------------------------------------------------

!  MAX_COORDSYSTEM_LENGTH, MAX_UNITS_LENGTH may be any int and are only
!    used to determine max lengths for attributes, which do not
!    affect any other lengths...

!------------------------------------------------------------------
  INTEGER :: max_time_length,max_dataclass_length
  INTEGER :: max_coordsystem_length,max_units_length
  INTEGER :: max_date_length,max_comment_length
  INTEGER :: max_filename_length,max_err_length,max_type_length
  INTEGER :: max_errorstring_length
  PARAMETER(max_time_length = 100)
  PARAMETER(max_dataclass_length = 20)
  PARAMETER(max_coordsystem_length  = 20)
  PARAMETER(max_units_length = 20)
  PARAMETER(max_date_length = 100)
  PARAMETER(max_comment_length = 256)
  PARAMETER(max_filename_length = 266)
  PARAMETER(max_err_length = 40)
  PARAMETER(max_type_length = 16)
  PARAMETER(max_errorstring_length = 80)
!----------------------------------------------------------------

  INTEGER :: next_frame
  PARAMETER(next_frame = 0)

!----------------------------------------------------------------
!----------------------------------------------------------------
!  ME_NO_MAXMIN,ME_DATA_END_MAXMIN,and ME_FRAME_END_MAXMIN
!    are used by the ConfigureMaxMin calls to determine when
!    or if the MaxMin records are created and written...
  INTEGER :: me_invalid_maxmin_control
  INTEGER :: me_no_maxmin
  INTEGER :: me_data_end_maxmin
  INTEGER :: me_frame_end_maxmin
  PARAMETER(me_invalid_maxmin_control = 0)
  PARAMETER(me_no_maxmin = 1)
  PARAMETER(me_data_end_maxmin = 2)
  PARAMETER(me_frame_end_maxmin = 4)
!----------------------------------------------------------------

! Scheme types

!----------------------------------------------------------------

! UpdateFlags are used in the Update Call to determine when an update
!    is necessary
!
  INTEGER :: update_all,update_new_frame,update_time
  PARAMETER(update_all = 0)
  PARAMETER(update_new_frame = 1)
  PARAMETER(update_time = 2)

!--------------------------------------------------------------
! Variable Types .. METype... not nc_types!
  INTEGER :: mescalar,mevector,meglobal,meclock,mepntatt
  INTEGER :: megeneric
  PARAMETER(mescalar = 1)
  PARAMETER(mevector = 2)
  PARAMETER(meglobal = 3)
  PARAMETER(meclock = 4)
  PARAMETER(mepntatt = 5)
  PARAMETER(megeneric = 6)

! File Control MEMode, MELIBERAL will cause the open files to stay
!    open and will eat up RAM, MECONSERVE closes out files after
!    the frames they hold are no longer being read..
  INTEGER :: meliberal,meconserve
  PARAMETER(meliberal = 0)
  PARAMETER(meconserve = 2)

! Coordinate Systems
  INTEGER :: mecartesian,mespherical,meregcartesian
  PARAMETER(mecartesian = 1)
  PARAMETER(mespherical = 2)
  PARAMETER(meregcartesian = 4)

! Coordinate Units
  INTEGER :: memeters,mekilometers,melightyears,meotherunits
  PARAMETER(memeters = 1)
  PARAMETER(mekilometers = 2)
  PARAMETER(melightyears = 3)
  PARAMETER(meotherunits = 4)

!-------------------------------------
! When GetXYZUnits is called the two options for the LocFlag
! Are:
  INTEGER :: meoffsetunits,mebaseunits
  PARAMETER(meoffsetunits = 0)
  PARAMETER(mebaseunits  = 1)
!-------------------------------------
! Location Configurations
  INTEGER :: me_non_time_var,me_time_var
  PARAMETER(me_non_time_var = 0)
  PARAMETER(me_time_var = 1)

! ModelRunning...
  INTEGER :: me_model_running,me_model_stopped
  PARAMETER(me_model_running = 1)
  PARAMETER(me_model_stopped = 0)

! Record Locking Flags
  INTEGER :: me_rec_lock,me_no_rec_lock
  PARAMETER(me_rec_lock = 5)
  PARAMETER(me_no_rec_lock = 0)

  INTEGER :: meerr,mesuccess,meiderr
  PARAMETER(meerr = -1)
  PARAMETER(mesuccess = 1)
  PARAMETER(meiderr = -1)

!--------------------------------------------------------------
!----- Default DataSet settings...
  INTEGER :: minvarname,defltmxpaths,defltmxschemes
  INTEGER :: defltmxvars,defltmxvects,defltmxvarname
  INTEGER :: defltmxpathname
  PARAMETER (minvarname = 35)
  PARAMETER (defltmxpaths = 16)
  PARAMETER (defltmxschemes = 32)
  PARAMETER (defltmxvars = 100)
  PARAMETER (defltmxvects = 30)
  PARAMETER (defltmxvarname = 64)
  PARAMETER (defltmxpathname = 128)
!--------------------------------------------------------------

