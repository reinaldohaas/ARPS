!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE open_ncd_file               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE open_ncd_file(filename,ncidout)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Open a WRF file and return NetCDF file handler.
!    NOTE: it is required to call close_wrf_file explicitly to close
!          the opened file in your calling program.
!
!------------------------------------------------------------------

  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN)  :: filename
  INTEGER,          INTENT(OUT) :: ncidout

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------
  INTEGER :: istatus
  LOGICAL :: fexists

  INCLUDE 'netcdf.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Begining of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  INQUIRE(FILE = filename, EXIST = fexists)
  IF (fexists) THEN
    istatus = NF_OPEN(TRIM(filename),NF_NOWRITE,ncidout)
    CALL nf_handle_error(istatus,'NF_OPEN in open_wrf_file')
  ELSE
    WRITE(6,'(2a)') 'File not found: ', filename
    STOP 'open_wrf_file'
  ENDIF

  RETURN

END SUBROUTINE open_ncd_file
!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE close_ncd_file               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE close_ncd_file(ncid)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!     Close the WRF file which is opened using open_wrf_file.
!
!------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ncid

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------
!
  INTEGER :: istatus

  INCLUDE 'netcdf.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  istatus = NF_CLOSE(ncid)
  CALL nf_handle_error(istatus,'NF_CLOSE in close_wrf_file')

  RETURN
END SUBROUTINE close_ncd_file
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE get_ncd_frames_per_outfile     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_ncd_frames_per_outfile(filename,numframes,istatus)

!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Get the size of unlimitted dimension in the file
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!  Yunheng Wang (03/26/2007)
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------

  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN)  :: filename
  INTEGER,          INTENT(OUT) :: numframes
  INTEGER,          INTENT(OUT) :: istatus

!------------------------------------------------------------------
!
!  Misc. local variables
!
!------------------------------------------------------------------

  INTEGER :: ncid
  INTEGER :: dimid
  INTEGER :: timelen

  INCLUDE 'netcdf.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus   = 0
  numframes = 0

  CALL open_ncd_file(TRIM(filename),ncid)

  istatus = NF_INQ_DIMID(ncid,'Time',dimid)
  CALL nf_handle_error(istatus,'NF_INQ_DIMID in get_ncd_frames_per_outfile')
  istatus = NF_INQ_DIMLEN(ncid,dimid,timelen)
  CALL nf_handle_error(istatus,'NF_INQ_DIMLEN in get_ncd_frames_per_outfile')
  IF( timelen < 1) THEN
    WRITE(6,'(1x,3a,I2)') 'ERROR: The unlimited dimension in the file ',  &
                          TRIM(filename),' is bad, timelen = ',timelen, '.'
    istatus = -1
  ELSE
    numframes = timelen
  END IF

  CALL close_ncd_file(ncid)

  RETURN
END SUBROUTINE get_ncd_frames_per_outfile
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE get_ncd_next_times             ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_ncd_next_time(ncid,itime,timestr,istatus)
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Read the the Date String in the WRF outputs at specified time
!
!-----------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: ncid     ! NetCDF file handler
  INTEGER, INTENT(IN)  :: itime    ! Time dimension value
                                   ! this is the unlimited dimension
  CHARACTER(LEN=*), INTENT(OUT) :: timestr
  INTEGER,          INTENT(OUT) :: istatus

!------------------------------------------------------------------
!
!  Misc. local variables
!
!------------------------------------------------------------------
  INTEGER :: dimid, varid
  INTEGER :: timelen, strlen, dateStrLen

  INCLUDE 'netcdf.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = NF_INQ_DIMID(ncid,'Time',dimid)
  CALL nf_handle_error(istatus,'NF_INQ_DIMID in get_wrf_Times')
  istatus = NF_INQ_DIMLEN(ncid,dimid,timelen)
  CALL nf_handle_error(istatus,'NF_INQ_DIMLEN in get_wrf_Times')
  IF(itime > timelen .OR. itime < 1) THEN
    WRITE(6,'(a,I2,a,I2)') ' The unlimited dimension in the file is ',  &
                           timelen, ' itime is ',itime
    STOP 'Time_out_of_bound'
  END IF

  istatus = NF_INQ_DIMID(ncid,'DateStrLen',dimid)
  CALL nf_handle_error(istatus,'NF_INQ_DIMID in get_wrf_Times')
  istatus = NF_INQ_DIMLEN(ncid,dimid,dateStrLen)
  CALL nf_handle_error(istatus,'NF_INQ_DIMLEN in get_wrf_Times')

  strlen = LEN(timestr)
  IF(strlen < dateStrLen) THEN
    WRITE(6,'(a)') ' Date string is not long enough to hold Times.'
    STOP 'timestr_too_short'
  END IF

  istatus = NF_INQ_VARID(ncid, 'Times', varid)
  CALL nf_handle_error(istatus,'NF_INQ_VARID in get_wrf_Times')
  istatus = NF_GET_VARA_TEXT(ncid,varid,(/1,itime/),                    &
                             (/dateStrLen,1/),timestr)
  CALL nf_handle_error(istatus,'NF_GET_VARA_TEXT in get_wrf_Times')

  RETURN
END SUBROUTINE get_ncd_next_time
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE get_wrf_att                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_ncd_dom_ti_integer(ncid,element,val,ireturn)
!-----------------------------------------------------------------------
!
! PURPOSE
!
!   Retieve WRF grib information from the NetCDF file which are stored
!   as Global attributes.
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,      INTENT(IN)  :: ncid
  CHARACTER(*), INTENT(IN)  :: element
  INTEGER,      INTENT(OUT) :: val
  INTEGER,      INTENT(OUT) :: ireturn

  INCLUDE 'netcdf.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ireturn = NF_GET_ATT_INT(ncid,NF_GLOBAL,element,val)
  CALL nf_handle_error(ireturn,'NF_GET_ATT_INT ('//TRIM(element)//') in getwrfd')

  RETURN
END SUBROUTINE get_ncd_dom_ti_integer

SUBROUTINE get_ncd_dom_ti_real(ncid,element,val,ireturn)
!-----------------------------------------------------------------------
!
! PURPOSE
!
!   Retieve WRF grib information from the NetCDF file which are stored
!   as Global attributes.
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,      INTENT(IN)  :: ncid
  CHARACTER(*), INTENT(IN)  :: element
  REAL,         INTENT(OUT) :: val
  INTEGER,      INTENT(OUT) :: ireturn

  INCLUDE 'netcdf.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ireturn = NF_GET_ATT_REAL(ncid,NF_GLOBAL,element,val)
  CALL nf_handle_error(ireturn,'NF_GET_ATT_REAL in getwrfd')

  RETURN
END SUBROUTINE get_ncd_dom_ti_real

SUBROUTINE get_ncd_dom_ti_char(ncid,element,val,ireturn)
!-----------------------------------------------------------------------
!
! PURPOSE
!
!   Retieve WRF grib information from the NetCDF file which are stored
!   as Global attributes.
!
!-----------------------------------------------------------------------
  IMPLICIT NONE

  INTEGER,      INTENT(IN)  :: ncid
  CHARACTER(*), INTENT(IN)  :: element
  CHARACTER(*), INTENT(OUT) :: val
  INTEGER,      INTENT(OUT) :: ireturn

  INCLUDE 'netcdf.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ireturn = NF_GET_ATT_TEXT(ncid,NF_GLOBAL,element,val)
  CALL nf_handle_error(ireturn,'NF_GET_ATT_TEXT in get_ncd_dom_ti_char.')

  RETURN
END SUBROUTINE get_ncd_dom_ti_char
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE get_ncd_1d                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_ncd_1d(ncid,itime,varname,nx,var1d,istatus)
!-----------------------------------------------------------------------
!
! PURPOSE:
!
!   Read in a 1D array from the WRF NetCDF file.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER,          INTENT(IN)  :: ncid
  INTEGER,          INTENT(IN)  :: itime
  CHARACTER(LEN=*), INTENT(IN)  :: varname
  INTEGER,          INTENT(IN)  :: nx
  REAL,             INTENT(OUT) :: var1d(nx)
  INTEGER,          INTENT(OUT) :: istatus
!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------

  INCLUDE 'netcdf.inc'

  INTEGER                    :: varid
  CHARACTER(LEN=NF_MAX_NAME) :: namein
  INTEGER                    :: vartype, ndims,natts,dimlen
  INTEGER                    :: dimids(NF_MAX_VAR_DIMS)


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = NF_INQ_VARID(ncid,varname,varid)
  CALL nf_handle_error(istatus,'NF_INQ_VARID ('//TRIM(varname)//') in get_wrf_1d')

  istatus = NF_INQ_VAR(ncid,varid,namein,vartype,ndims,dimids,natts)
  CALL nf_handle_error(istatus,'NF_INQ_VAR ('//TRIM(varname)//') in get_wrf_1d')

  IF(vartype /= NF_FLOAT) THEN
    WRITE(6,'(3a)') 'Variable ',varname, ' is not REAL.'
    STOP 'WRONG_VAR_TYPE'
  END IF

  IF(ndims /= 2) THEN
    WRITE(6,'(3a)') 'Variable ', varname, ' is not a 1D array.'
    STOP 'WRONG_VAR_DIMENSIONS'
  END IF

  istatus = NF_INQ_DIMLEN(ncid,dimids(1),dimlen)
  CALL nf_handle_error(istatus,'NF_INQ_DIMLEN ('//TRIM(varname)//') in get_wrf_1d')
  IF(dimlen /= nx) THEN
    WRITE(6,'(3a,I3,a,I3)') 'The dimension of variable ', varname,    &
                    ' is ',dimlen, ' and it should be ',nx
    STOP 'WRONG_DIM_length'
  END IF

  istatus = NF_INQ_DIMLEN(ncid,dimids(2),dimlen)
  CALL nf_handle_error(istatus,'NF_INQ_DIMLEN ('//TRIM(varname)//') in get_wrf_1d')
  IF(dimlen < itime) THEN
    WRITE(6,'(a,I3,a,I3)') 'The total records number is ', dimlen,     &
                    ' however, the required time level is ',itime
    STOP 'itime_tool_large'
  END IF

  istatus = NF_GET_VARA_REAL(ncid,varid,(/1,itime/),(/nx,1/),var1d)
  CALL nf_handle_error(istatus,'NF_GET_VARA_REAL ('//TRIM(varname)//') in get_wrf_1d')

  RETURN
END SUBROUTINE get_ncd_1d
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE get_ncd_scalar                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_ncd_scalar(ncid,itime,varname,var,istatus)
!-----------------------------------------------------------------------
!
! PURPOSE:
!
!   Read in a scalar from the WRF NetCDF file.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER,          INTENT(IN)  :: ncid
  INTEGER,          INTENT(IN)  :: itime
  CHARACTER(LEN=*), INTENT(IN)  :: varname
  REAL,             INTENT(OUT) :: var
  INTEGER,          INTENT(OUT) :: istatus
!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------

  INCLUDE 'netcdf.inc'

  INTEGER                    :: varid
  CHARACTER(LEN=NF_MAX_NAME) :: namein
  INTEGER                    :: vartype, ndims,natts,dimlen
  INTEGER                    :: dimids(NF_MAX_VAR_DIMS)

  !REAL :: varin(1)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = NF_INQ_VARID(ncid,varname,varid)
  CALL nf_handle_error(istatus,'NF_INQ_VARID ('//TRIM(varname)//') in get_ncd_scalar')

  istatus = NF_INQ_VAR(ncid,varid,namein,vartype,ndims,dimids,natts)
  CALL nf_handle_error(istatus,'NF_INQ_VAR ('//TRIM(varname)//') in get_ncd_scalar')

  IF(vartype /= NF_FLOAT) THEN
    WRITE(6,'(3a)') 'Variable ',varname, ' is not REAL.'
    STOP 'WRONG_VAR_TYPE'
  END IF

  IF(ndims /= 1) THEN
    WRITE(6,'(3a)') 'Variable ', varname, ' is not a scalar.'
    STOP 'WRONG_VAR_DIMENSIONS'
  END IF

  istatus = NF_INQ_DIMLEN(ncid,dimids(1),dimlen)
  CALL nf_handle_error(istatus,'NF_INQ_DIMLEN in get_wrf_1d')
  IF(dimlen < itime) THEN
    WRITE(6,'(a,I3,a,I3)') 'The total records number is ', dimlen,     &
                    ' however, the required time level is ',itime
    STOP 'itime_tool_large'
  END IF

  istatus = NF_GET_VARA_REAL(ncid,varid,(/itime/),(/1/),var)
  !istatus = NF_GET_VAR_REAL(ncid,varid,varin)
  CALL nf_handle_error(istatus,'NF_GET_VAR_REAL in get_wrf_scalar')

  RETURN
END SUBROUTINE get_ncd_scalar
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE get_ncd_2d                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_ncd_2d(ncid,itime,varname,nx,ny,var2d,istatus)
!-----------------------------------------------------------------------
!
! PURPOSE:
!
!   Read in a 2D array from the WRF NetCDF file.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER,          INTENT(IN)  :: ncid
  INTEGER,          INTENT(IN)  :: itime
  CHARACTER(LEN=*), INTENT(IN)  :: varname
  INTEGER,          INTENT(IN)  :: nx
  INTEGER,          INTENT(IN)  :: ny
  REAL,             INTENT(OUT) :: var2d(nx,ny)
  INTEGER,          INTENT(OUT) :: istatus
!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------

  INCLUDE 'netcdf.inc'

  INTEGER                    :: varid
  CHARACTER(LEN=NF_MAX_NAME) :: namein
  INTEGER                    :: vartype, ndims,natts,dimlen
  INTEGER                    :: dimids(NF_MAX_VAR_DIMS)


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = NF_INQ_VARID(ncid,varname,varid)
  CALL nf_handle_error(istatus,'NF_INQ_VARID in get_wrf_2d with '//TRIM(varname))

  istatus = NF_INQ_VAR(ncid,varid,namein,vartype,ndims,dimids,natts)
  CALL nf_handle_error(istatus,'NF_INQ_VAR in get_wrf_2d')

  IF(vartype /= NF_FLOAT) THEN
    WRITE(6,'(3a)') 'Variable ',varname, ' is not REAL.'
    STOP 'WRONG_VAR_TYPE'
  END IF

  IF(ndims /= 3) THEN
    WRITE(6,'(3a)') 'Variable ', varname, ' is not a 2D array.'
    STOP 'WRONG_VAR_DIMENSIONS'
  END IF

  istatus = NF_INQ_DIMLEN(ncid,dimids(1),dimlen)
  CALL nf_handle_error(istatus,'NF_INQ_DIMLEN in get_wrf_2d')
  IF(dimlen /= nx) THEN
    WRITE(6,'(3a,I3,a,I3)') 'First dimension of variable ', varname,    &
                    ' is ',dimlen, ' and it should be ',nx
    STOP 'WRONG_DIM_length'
  END IF

  istatus = NF_INQ_DIMLEN(ncid,dimids(2),dimlen)
  CALL nf_handle_error(istatus,'NF_INQ_DIMLEN in get_wrf_2d')
  IF(dimlen /= ny) THEN
    WRITE(6,'(3a,I3,a,I3)') 'Second dimension of variable ', varname,   &
                    ' is ',dimlen, ' and it should be ',ny
    STOP 'WRONG_DIM_length'
  END IF

  istatus = NF_INQ_DIMLEN(ncid,dimids(3),dimlen)
  CALL nf_handle_error(istatus,'NF_INQ_DIMLEN in get_wrf_2d')
  IF(dimlen < itime) THEN
    WRITE(6,'(a,I3,a,I3)') 'The total records number is ', dimlen,     &
                    ' however, the required time level is ',itime
    STOP 'itime_tool_large'
  END IF

  istatus = NF_GET_VARA_REAL(ncid,varid,(/1,1,itime/),(/nx,ny,1/),var2d)
  CALL nf_handle_error(istatus,'NF_GET_VARA_REAL in get_wrf_2d')

  RETURN
END SUBROUTINE get_ncd_2d
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE get_ncd_2di                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_ncd_2di(ncid,itime,varname,nx,ny,var2d,istatus)
!-----------------------------------------------------------------------
!
! PURPOSE:
!
!   Read in a 2D INTEGER array from the WRF NetCDF file.
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER,          INTENT(IN)  :: ncid
  INTEGER,          INTENT(IN)  :: itime
  CHARACTER(LEN=*), INTENT(IN)  :: varname
  INTEGER,          INTENT(IN)  :: nx
  INTEGER,          INTENT(IN)  :: ny
  INTEGER,          INTENT(OUT) :: var2d(nx,ny)
  INTEGER,          INTENT(OUT) :: istatus
!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------

  INCLUDE 'netcdf.inc'

  INTEGER                    :: varid
  CHARACTER(LEN=NF_MAX_NAME) :: namein
  INTEGER                    :: vartype, ndims,natts,dimlen
  INTEGER                    :: dimids(NF_MAX_VAR_DIMS)


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = NF_INQ_VARID(ncid,varname,varid)
  CALL nf_handle_error(istatus,'NF_INQ_VARID in get_wrf_2di')

  istatus = NF_INQ_VAR(ncid,varid,namein,vartype,ndims,dimids,natts)
  CALL nf_handle_error(istatus,'NF_INQ_VAR in get_wrf_2di')

  IF(vartype /= NF_INT) THEN
    WRITE(6,'(3a)') 'Variable ',varname, ' is not INTEGER.'
    STOP 'WRONG_VAR_TYPE'
  END IF

  IF(ndims /= 3) THEN
    WRITE(6,'(3a)') 'Variable ', varname, ' is not a 2D array.'
    STOP 'WRONG_VAR_DIMENSIONS'
  END IF

  istatus = NF_INQ_DIMLEN(ncid,dimids(1),dimlen)
  CALL nf_handle_error(istatus,'NF_INQ_DIMLEN in get_wrf_2di')
  IF(dimlen /= nx) THEN
    WRITE(6,'(3a,I3,a,I3)') 'First dimension of variable ', varname,    &
                    ' is ',dimlen, ' and it should be ',nx
    STOP 'WRONG_DIM_length'
  END IF

  istatus = NF_INQ_DIMLEN(ncid,dimids(2),dimlen)
  CALL nf_handle_error(istatus,'NF_INQ_DIMLEN in get_wrf_2di')
  IF(dimlen /= ny) THEN
    WRITE(6,'(3a,I3,a,I3)') 'Second dimension of variable ', varname,   &
                    ' is ',dimlen, ' and it should be ',ny
    STOP 'WRONG_DIM_length'
  END IF

  istatus = NF_INQ_DIMLEN(ncid,dimids(3),dimlen)
  CALL nf_handle_error(istatus,'NF_INQ_DIMLEN in get_wrf_2di')
  IF(dimlen < itime) THEN
    WRITE(6,'(a,I3,a,I3)') 'The total records number is ', dimlen,     &
                    ' however, the required time level is ',itime
    STOP 'itime_tool_large'
  END IF

  istatus = NF_GET_VARA_INT(ncid,varid,(/1,1,itime/),(/nx,ny,1/),var2d)
  CALL nf_handle_error(istatus,'NF_GET_VARA_INT in get_wrf_2di')

  RETURN
END SUBROUTINE get_ncd_2di
!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE get_ncd_3d                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_ncd_3d(ncid,itime,varname,nx,ny,nz,var3d,istatus)
!-----------------------------------------------------------------------
!
!  PURPOSE:
!    Read in a 3D array from the WRF NetCDF file
!
!-----------------------------------------------------------------------
  IMPLICIT NONE
  INTEGER,          INTENT(IN)  :: ncid
  INTEGER,          INTENT(IN)  :: itime
  CHARACTER(LEN=*), INTENT(IN)  :: varname
  INTEGER,          INTENT(IN)  :: nx
  INTEGER,          INTENT(IN)  :: ny
  INTEGER,          INTENT(IN)  :: nz
  REAL,             INTENT(OUT) :: var3d(nx,ny,nz)
  INTEGER,          INTENT(OUT) :: istatus

!-----------------------------------------------------------------------
!
! Misc. Local variables
!
!-----------------------------------------------------------------------
  INCLUDE 'netcdf.inc'

  INTEGER                    :: varid
  CHARACTER(LEN=NF_MAX_NAME) :: namein
  INTEGER                    :: vartype, ndims,natts,dimlen
  INTEGER                    :: dimids(NF_MAX_VAR_DIMS)

  INTEGER, PARAMETER         :: VAR_NOTEXIST = -1

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = NF_INQ_VARID(ncid,varname,varid)
  IF(istatus == NF_ENOTVAR) THEN
     WRITE(6,'(3a)') ' WARNING: variable ',varname,' does not exist.'
     var3d(:,:,:) = -999.0
     istatus = VAR_NOTEXIST
     RETURN
  END IF
  CALL nf_handle_error(istatus,'NF_INQ_VARID in get_wrf_3d')

  istatus = NF_INQ_VAR(ncid,varid,namein,vartype,ndims,dimids,natts)
  CALL nf_handle_error(istatus,'NF_INQ_VAR in get_wrf_3d')

  IF(vartype /= NF_FLOAT) THEN
    WRITE(6,'(3a)') 'Variable ',varname, ' is not REAL.'
    STOP 'WRONG_VAR_TYPE'
  END IF

  IF(ndims /= 4) THEN
    WRITE(6,'(3a)') 'Variable ', varname, ' is not a 3D array.'
    STOP 'WRONG_VAR_DIMENSIONS'
  END IF

  istatus = NF_INQ_DIMLEN(ncid,dimids(1),dimlen)
  CALL nf_handle_error(istatus,'NF_INQ_DIMLEN in get_wrf_3d')
  IF(dimlen /= nx) THEN
    WRITE(6,'(3a,I3,a,I3)') 'First dimension of variable ', varname,    &
                    ' is ',dimlen, ' and it should be ',nx
    STOP 'WRONG_DIM_length'
  END IF

  istatus = NF_INQ_DIMLEN(ncid,dimids(2),dimlen)
  CALL nf_handle_error(istatus,'NF_INQ_DIMLEN in get_wrf_3d')
  IF(dimlen /= ny) THEN
    WRITE(6,'(3a,I3,a,I3)') 'Second dimension of variable ', varname,   &
                    ' is ',dimlen, ' and it should be ',ny
    STOP 'WRONG_DIM_length'
  END IF

  istatus = NF_INQ_DIMLEN(ncid,dimids(3),dimlen)
  CALL nf_handle_error(istatus,'NF_INQ_DIMLEN in get_wrf_3d')
  IF(dimlen /= nz) THEN
    WRITE(6,'(3a,I3,a,I3)') 'Third dimension of variable ', varname,   &
                    ' is ',dimlen, ' and it should be ',nz
    STOP 'WRONG_DIM_length'
  END IF

  istatus = NF_INQ_DIMLEN(ncid,dimids(4),dimlen)
  CALL nf_handle_error(istatus,'NF_INQ_DIMLEN in get_wrf_3d')
  IF(dimlen < itime) THEN
    WRITE(6,'(a,I3,a,I3)') 'The total records number is ', dimlen,     &
                    ' however, the required time level is ',itime
    STOP 'itime_tool_large'
  END IF

  istatus = NF_GET_VARA_REAL(ncid,varid,(/1,1,1,itime/),               &
                             (/nx,ny,nz,1/),var3d)
  IF(istatus /= NF_NOERR .OR. istatus /= NF_EEXIST) THEN
    CALL nf_handle_error(istatus,'NF_GET_VARA_REAL in get_wrf_3d')
  END IF

  RETURN
END SUBROUTINE get_ncd_3d
!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE get_ncd_dimension            ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_ncd_dimension(ncid,nx_ext,ny_ext,nz_ext,nzsoil_ext,istatus)
!
!------------------------------------------------------------------------
!
! PURPOSE:
!
!   Read dimension parameters and the first time string from
!   WRF output file.
!
!   Please note this subroutine will open a file to read and then close
!   it. So it does not leave any opened handler for NetCDF file.
!
!------------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER,          INTENT(IN)  :: ncid
  INTEGER,          INTENT(OUT) :: nx_ext
  INTEGER,          INTENT(OUT) :: ny_ext
  INTEGER,          INTENT(OUT) :: nz_ext
  INTEGER,          INTENT(OUT) :: nzsoil_ext
  INTEGER,          INTENT(OUT) :: istatus

!------------------------------------------------------------------------
!
!  Misc. Local variables
!
!------------------------------------------------------------------------

  INTEGER :: dimid

  INCLUDE 'netcdf.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Begining of executable code
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
! Get WRF dimensions
!
!-----------------------------------------------------------------------

  istatus = NF_INQ_DIMID(ncid,'west_east_stag',dimid)
  CALL nf_handle_error(istatus,'NF_INQ_DIMID in get_ncd_dimension')
  istatus = NF_INQ_DIMLEN(ncid,dimid,nx_ext)
  CALL nf_handle_error(istatus,'NF_INQ_DIMLEN in get_ncd_dimension')

  istatus = NF_INQ_DIMID(ncid,'south_north_stag',dimid)
  CALL nf_handle_error(istatus,'NF_INQ_DIMID in get_ncd_dimension')
  istatus = NF_INQ_DIMLEN(ncid,dimid,ny_ext)
  CALL nf_handle_error(istatus,'NF_INQ_DIMLEN in get_ncd_dimension')

  istatus = NF_INQ_DIMID(ncid,'bottom_top_stag',dimid)
  CALL nf_handle_error(istatus,'NF_INQ_DIMID in get_ncd_dimension')
  istatus = NF_INQ_DIMLEN(ncid,dimid,nz_ext)
  CALL nf_handle_error(istatus,'NF_INQ_DIMLEN in get_ncd_dimension')

  istatus = NF_INQ_DIMID(ncid,'soil_layers_stag',dimid)
  CALL nf_handle_error(istatus,'NF_INQ_DIMID in get_ncd_dimension')
  istatus = NF_INQ_DIMLEN(ncid,dimid,nzsoil_ext)
  CALL nf_handle_error(istatus,'NF_INQ_DIMLEN in get_ncd_dimension')

  RETURN
END SUBROUTINE get_ncd_dimension
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE nf_handle_error           ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE nf_handle_error(ierr,sub_name)

!------------------------------------------------------------------------
!
! PURPOSE:
!   Write error message to the standard output if ierr contains an error
!
!------------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER,          INTENT(IN) :: ierr
  CHARACTER(LEN=*), INTENT(IN) :: sub_name
  CHARACTER(LEN=80) :: errmsg

  INCLUDE 'netcdf.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF(ierr /= NF_NOERR) THEN
    errmsg = NF_STRERROR(ierr)
    WRITE(6,*) 'NetCDF error: ',errmsg
    WRITE(6,*) 'Program stopped while calling ', sub_name
    STOP
  END IF

  RETURN
END SUBROUTINE nf_handle_error
