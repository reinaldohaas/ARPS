!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE open_ncd_for_write             ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE open_ncd_for_write(filename,LargeFile,ifile,wrfversion,      &
                  nx,ny,nz,nzsoil,bdywidth,nout,istatus)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!   Prepare NetCDF output file to be ready for write. After the call
!   the file should still in DEFINE mode.
!
!------------------------------------------------------------------

  USE wrf_metadata           ! use wrfversion from the module

  IMPLICIT NONE

  CHARACTER(LEN=*), INTENT(IN)  :: filename
  LOGICAL,          INTENT(IN)  :: LargeFile
  INTEGER,          INTENT(IN)  :: ifile
  REAL,             INTENT(IN)  :: wrfversion
  INTEGER,          INTENT(IN)  :: nx
  INTEGER,          INTENT(IN)  :: ny
  INTEGER,          INTENT(IN)  :: nz
  INTEGER,          INTENT(IN)  :: nzsoil
  INTEGER,          INTENT(IN)  :: bdywidth
  INTEGER,          INTENT(OUT) :: nout
  INTEGER,          INTENT(OUT) :: istatus

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INCLUDE 'netcdf.inc'

  INTEGER :: lenstr
  INTEGER :: dimunlim_id, dimwe_stagid, dimwe_unstagid
  INTEGER :: dimsn_stagid, dimsn_unstagid
  INTEGER :: dimbt_stagid, dimbt_unstagid
  INTEGER :: dimdatestr_id, dimscalar_id
  INTEGER :: dimbdywd_id, dimland_id, dimsoil_id, dimnzsoil_id
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  lenstr = LEN_TRIM(filename)
  
  IF (LargeFile) THEN
    istatus = NF_CREATE(TRIM(filename),IOR(NF_CLOBBER,NF_64BIT_OFFSET),  &
                        nout)                                ! CDF2
  ELSE
    istatus = NF_CREATE(TRIM(filename),NF_CLOBBER,nout)      ! CDF1
  END IF

  CALL nf_handle_error(istatus,'NF_CREATE')
  
  ! define dimensions
  istatus = NF_DEF_DIM(nout,'Time',NF_UNLIMITED,dimunlim_id)
  CALL nf_handle_error(istatus,'NF_DEF_DIM')
  istatus = NF_DEF_DIM(nout,'DateStrLen',19,dimdatestr_id)
  CALL nf_handle_error(istatus,'NF_DEF_DIM')
  
  istatus = NF_DEF_DIM(nout,'west_east_stag',nx,dimwe_stagid)
  CALL nf_handle_error(istatus,'NF_DEF_DIM')
  istatus = NF_DEF_DIM(nout,'west_east',nx-1,dimwe_unstagid)
  CALL nf_handle_error(istatus,'NF_DEF_DIM')
  istatus = NF_DEF_DIM(nout,'south_north_stag',ny,dimsn_stagid)
  CALL nf_handle_error(istatus,'NF_DEF_DIM')
  istatus = NF_DEF_DIM(nout,'south_north',ny-1,dimsn_unstagid)
  CALL nf_handle_error(istatus,'NF_DEF_DIM')

  IF (ifile == 1 .OR. ifile == 2) THEN 
    istatus = NF_DEF_DIM(nout,'bottom_top_stag',nz,dimbt_stagid)
    CALL nf_handle_error(istatus,'NF_DEF_DIM')
    istatus = NF_DEF_DIM(nout,'bottom_top',nz-1,dimbt_unstagid)
    CALL nf_handle_error(istatus,'NF_DEF_DIM')
  END IF
  
  IF (ifile == 1) THEN        ! WRF input file
    IF (wrfversion < 3.0) THEN
    istatus = NF_DEF_DIM(nout,'DIM0008',bdywidth,dimbdywd_id)
    CALL nf_handle_error(istatus,'NF_DEF_DIM')
    ELSE
      istatus = NF_DEF_DIM(nout,'DIM0009',bdywidth,dimbdywd_id)
      CALL nf_handle_error(istatus,'NF_DEF_DIM')
    END IF

    IF (wrfversion >= 3.1) THEN
      istatus = NF_DEF_DIM(nout,'land_cat_stag',LanduseCategories,dimbdywd_id)
      CALL nf_handle_error(istatus,'NF_DEF_DIM')

      istatus = NF_DEF_DIM(nout,'soil_cat_stag',SoilCategories,dimbdywd_id)
      CALL nf_handle_error(istatus,'NF_DEF_DIM')
    END IF
 
!    istatus = NF_DEF_DIM(nout,'ext_scalar',1,dimscalar_id)
!    CALL nf_handle_error(istatus,'NF_DEF_DIM')
  
    istatus = NF_DEF_DIM(nout,'soil_layers_stag',nzsoil,dimnzsoil_id)
    CALL nf_handle_error(istatus,'NF_DEF_DIM')
  
  ELSE IF (ifile == 2) THEN   ! WRF bdy file
    istatus = NF_DEF_DIM(nout,'bdy_width',bdywidth,dimbdywd_id)
    CALL nf_handle_error(istatus,'NF_DEF_DIM')
  ELSE IF (ifile == 3) THEN   ! WRF static file
    istatus = NF_DEF_DIM(nout,'land_cat',LanduseCategories,dimbdywd_id)
    CALL nf_handle_error(istatus,'NF_DEF_DIM')
    istatus = NF_DEF_DIM(nout,'soil_cat',SoilCategories,dimbdywd_id)
    CALL nf_handle_error(istatus,'NF_DEF_DIM')
    istatus = NF_DEF_DIM(nout,'month',12,dimbdywd_id)
    CALL nf_handle_error(istatus,'NF_DEF_DIM')
  ELSE
    ! Do nothing, should not reach here
  END IF

  IF (ifile == 3) THEN
    CALL define_wrf_static(nout,nx,ny,LanduseCategories,SoilCategories,istatus)
  ELSE

    CALL define_wrf_variables_V2(nout,ifile,wrfversion,nx,ny,nz,nzsoil, &
                                 bdywidth,istatus)
  END IF

  RETURN
END SUBROUTINE open_ncd_for_write
!
!##################################################################
!##################################################################
!######                                                      ######
!######           SUBROUTINE close_ncd_for_write             ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE close_ncd_for_write(nch,istatus)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!   Close the output file.
!
!------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: nch
  INTEGER, INTENT(OUT) :: istatus

  INCLUDE 'netcdf.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  istatus = NF_CLOSE(nch)
  CALL nf_handle_error(istatus,'NF_CLOSE')

  RETURN
END SUBROUTINE close_ncd_for_write

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

  IMPLICIT NONE
  INTEGER,          INTENT(IN) :: ierr
  CHARACTER(LEN=*), INTENT(IN) :: sub_name

  CHARACTER(LEN=80) :: errmsg

  INCLUDE 'netcdf.inc'

  IF(ierr /= NF_NOERR) THEN
    errmsg = NF_STRERROR(ierr)
    WRITE(6,'(/2a)') 'NetCDF error: ',errmsg
    CALL arpsstop('Program stopped while calling '//sub_name,1)
  END IF

  RETURN
END SUBROUTINE nf_handle_error
!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE write_global_attribute         ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE enter_ncd_define(ncid,ireturn)

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: ncid
  INTEGER, INTENT(OUT) :: ireturn

  INCLUDE 'netcdf.inc'

  ireturn = nf_redef(ncid)

  RETURN
END SUBROUTINE enter_ncd_define

SUBROUTINE exit_ncd_define(ncid,ireturn)

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: ncid
  INTEGER, INTENT(OUT) :: ireturn

  INCLUDE 'netcdf.inc'

  ireturn = nf_enddef(ncid)

  RETURN
END SUBROUTINE exit_ncd_define

SUBROUTINE put_ncd_dom_ti_char(ncid,attname,attstr,ireturn)

  IMPLICIT NONE

  INTEGER,      INTENT(IN)  :: ncid
  CHARACTER(*), INTENT(IN)  :: attname
  CHARACTER(*), INTENT(IN)  :: attstr
  INTEGER,      INTENT(OUT) :: ireturn

  INCLUDE 'netcdf.inc'

  INTEGER :: lenstr

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Begin of executable  code ... ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  lenstr = LEN_TRIM(attstr)

  ireturn = NF_PUT_ATT_TEXT(ncid,NF_GLOBAL,TRIM(attname),lenstr,attstr)
  CALL nf_handle_error(ireturn,'NF_PUT_ATT_TEXT in put_ncd_dom_ti_char.')

  RETURN
END SUBROUTINE put_ncd_dom_ti_char

SUBROUTINE put_ncd_dom_ti_integer(ncid,attname,attval,ireturn)

  IMPLICIT NONE

  INTEGER,      INTENT(IN)  :: ncid
  CHARACTER(*), INTENT(IN)  :: attname
  INTEGER,      INTENT(IN)  :: attval
  INTEGER,      INTENT(OUT) :: ireturn

  INCLUDE 'netcdf.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Begin of executable  code ... ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ireturn = NF_PUT_ATT_INT(ncid,NF_GLOBAL,TRIM(attname),NF_INT,1,attval)

  RETURN
END SUBROUTINE put_ncd_dom_ti_integer
  
SUBROUTINE put_ncd_dom_ti_real(ncid,attname,attval,attsiz,ireturn)

  IMPLICIT NONE

  INTEGER,      INTENT(IN)  :: ncid
  CHARACTER(*), INTENT(IN)  :: attname
  INTEGER,      INTENT(IN)  :: attsiz
  REAL,         INTENT(IN)  :: attval(attsiz)
  INTEGER,      INTENT(OUT) :: ireturn

  INCLUDE 'netcdf.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Begin of executable  code ... ...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ireturn = NF_PUT_ATT_REAL(ncid,NF_GLOBAL,TRIM(attname),NF_FLOAT,      &
                            attsiz,attval)

  RETURN
END SUBROUTINE put_ncd_dom_ti_real
  
SUBROUTINE write_var_meta(ncid,varid,varmeta)

  USE wrf_metadata

  IMPLICIT NONE
  INTEGER,                INTENT(IN) :: ncid
  TYPE(wrf_var_metadata), INTENT(IN) :: varmeta
  INTEGER,                INTENT(IN) :: varid

  INCLUDE 'netcdf.inc'

  INTEGER  :: istatus, lenstr

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code ......
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  istatus = NF_PUT_ATT_INT(ncid,varid,'FieldType',NF_INT,1,varmeta%fieldType) ! REAL
  istatus = NF_PUT_ATT_TEXT(ncid,varid,'MemoryOrder',3,varmeta%memoryOrder)
!  istatus = NF_PUT_ATT_REAL(nout,varid,'_FillValue',NF_FLOAT,1,rmissing) ! REAL

  lenstr = LEN_TRIM(varmeta%description)
  istatus = NF_PUT_ATT_TEXT(ncid,varid,'description',lenstr,varmeta%description)

  lenstr = LEN_TRIM(varmeta%units)
  IF(lenstr < 1) lenstr = 1
  istatus = NF_PUT_ATT_TEXT(ncid,varid,'units',lenstr,varmeta%units)

  lenstr = LEN_TRIM(varmeta%stagger)
  IF(lenstr < 1) lenstr = 1
  istatus = NF_PUT_ATT_TEXT(ncid,varid,'stagger',lenstr,varmeta%stagger)

  lenstr = LEN_TRIM(varmeta%coordinates)
  IF(lenstr > 2) istatus = NF_PUT_ATT_TEXT(ncid,varid,'coordinates',    &
                                           lenstr,varmeta%coordinates)

  RETURN
END SUBROUTINE write_var_meta

!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE write_times_str                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE put_ncd_dom_td_char(ncid,varname,DateStr,ntime,istatus)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write a Date String for NetCDF WRF data file
!    Must call after the calling of wrf_define_variables
!
!------------------------------------------------------------------

  IMPLICIT NONE

  INTEGER,           INTENT(IN)  :: ncid
  CHARACTER(LEN=*),  INTENT(IN)  :: varname
  CHARACTER(LEN=19), INTENT(IN)  :: DateStr
  INTEGER,           INTENT(IN)  :: ntime
  INTEGER,           INTENT(OUT) :: istatus

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INTEGER :: varid
  INTEGER :: dimids(2), dimlens(2)
  INTEGER :: lenstr
  INTEGER :: nt

  INCLUDE 'netcdf.inc'
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  istatus = NF_INQ_VARID(ncid,varname,varid)
  CALL nf_handle_error(istatus,'NF_INQ_VARID in write_times_str.')

  lenstr = LEN_TRIM(DateStr)

  nt = ntime
  IF(nt < 1) nt = 1

  istatus = NF_INQ_VARDIMID(ncid,varid,dimids)

  istatus = NF_INQ_DIMLEN(ncid,dimids(1),dimlens(1))
  istatus = NF_INQ_DIMLEN(ncid,dimids(2),dimlens(2))   ! unlimit dimension

!  IF (nt > dimlens(2) ) THEN
!    WRITE(6,*) 'Length of Unlimited record does not match with the data time'
!    STOP
!  END IF
  IF(lenstr /= dimlens(1)) THEN
    WRITE(6,*) 'Length of DateStr does not match with the dimension size'
    STOP
  END IF

  istatus = NF_PUT_VARA_TEXT(ncid,varid,(/1,nt/),(/dimlens(1),1/),DateStr)
  
  RETURN
END SUBROUTINE put_ncd_dom_td_char
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE write_ncd_real              ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE write_ncd_real(nout,varname,var,istatus)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write a REAL scalar to the output file.
!
!------------------------------------------------------------------
  USE wrf_metadata
  IMPLICIT NONE
  
  INTEGER, INTENT(IN)          :: nout   ! output channel,  NetCDF id
  REAL,    INTENT(IN)          :: var
  CHARACTER(LEN=*), INTENT(IN) :: varname
  INTEGER, INTENT(OUT)         :: istatus

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INTEGER :: varid
  INTEGER :: dim_id
  INTEGER :: dimlen

  INCLUDE 'netcdf.inc'
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

! get variable id and dimension length

  istatus = NF_INQ_VARID(nout,varname,varid)
  CALL nf_handle_error(istatus,'NF_INQ_VARID in write1d.')

  istatus = NF_INQ_VARDIMID(nout,varid,dim_id)

  istatus = NF_INQ_DIMLEN(nout,dim_id,dimlen)   ! unlimit dimension

! Write data

  istatus = NF_PUT_VARA_REAL(nout,varid,(/1/),(/1/),var)

  CALL nf_handle_error(istatus,'NF_PUT_VARA_REAL in write_ncd_real.')

  RETURN
END SUBROUTINE write_ncd_real
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE write_ncd_int                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE write_ncd_int(nout,varname,var,istatus)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write an integer scalar to the output file.
!
!------------------------------------------------------------------
  USE wrf_metadata
  IMPLICIT NONE
  
  INTEGER,          INTENT(IN) :: nout   ! output channel,  NetCDF file id
  INTEGER,          INTENT(IN) :: var
  CHARACTER(LEN=*), INTENT(IN) :: varname
  INTEGER,          INTENT(OUT):: istatus

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INTEGER :: varid
  INTEGER :: dim_id
  INTEGER :: dimlen

  INCLUDE 'netcdf.inc'
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! get variable id and dimension length

  istatus = NF_INQ_VARID(nout,varname,varid)
  CALL nf_handle_error(istatus,'NF_INQ_VARID in write_ncd_int.')

  istatus = NF_INQ_VARDIMID(nout,varid,dim_id)

  istatus = NF_INQ_DIMLEN(nout,dim_id,dimlen)   ! unlimit dimension

! Write data

  istatus = NF_PUT_VARA_INT(nout,varid,(/1/),(/1/),var)
  CALL nf_handle_error(istatus,'NF_PUT_VARA_INT in write_ncd_int.')

  RETURN
END SUBROUTINE write_ncd_int
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE write_ncd_1d                ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE write_ncd_1d(nout,varname,var1d,ndim,istatus)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write 1D vector to the output file.
!
!------------------------------------------------------------------
  USE wrf_metadata
  IMPLICIT NONE
  
  INTEGER, INTENT(IN)          :: nout   ! output channel,  NetCDF id
  INTEGER, INTENT(IN)          :: ndim
  REAL,    INTENT(IN)          :: var1d(ndim)
  CHARACTER(LEN=*), INTENT(IN) :: varname
  INTEGER, INTENT(OUT)         :: istatus

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INTEGER :: varid
  INTEGER :: dim_ids(2)
  INTEGER :: dimlens(2)

  INCLUDE 'netcdf.inc'
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  WRITE(6,FMT='(2a)',ADVANCE='NO') '  Writing 1D variable ', varname

! get variable id and dimension length

  istatus = NF_INQ_VARID(nout,varname,varid)
  CALL nf_handle_error(istatus,'NF_INQ_VARID in write1d.')

  istatus = NF_INQ_VARDIMID(nout,varid,dim_ids)

  istatus = NF_INQ_DIMLEN(nout,dim_ids(1),dimlens(1))
  istatus = NF_INQ_DIMLEN(nout,dim_ids(2),dimlens(2))   ! unlimit dimension

  IF(dimlens(1) /= ndim) THEN
    WRITE(6,'(/2(a,i3))') 'Mismatched dimension size. In file = ',       &
                      dimlens(1), ' Passed in = ', ndim
    STOP
  END IF

! Write data

  istatus = NF_PUT_VARA_REAL(nout,varid,(/1,1/),(/ndim,1/),var1d)

  CALL nf_handle_error(istatus,'NF_PUT_VARA_REAL in write1d')

!  WRITE(6,'(a)') '         ...  DONE.'

  RETURN
END SUBROUTINE write_ncd_1d
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE write_ncd_1di                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE write_ncd_1di(nout,varname,var1di,ndim,istatus)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write 1D integer vector to the output file.
!
!------------------------------------------------------------------
  USE wrf_metadata
  IMPLICIT NONE
  
  INTEGER,          INTENT(IN) :: nout   ! output channel,  NetCDF file id
  INTEGER,          INTENT(IN) :: ndim
  INTEGER,          INTENT(IN) :: var1di(ndim)
  CHARACTER(LEN=*), INTENT(IN) :: varname
  INTEGER,          INTENT(OUT):: istatus

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INTEGER :: varid
  INTEGER :: dim_ids(2)
  INTEGER :: dimlens(2)

  INCLUDE 'netcdf.inc'
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  WRITE(6,FMT='(2a)',ADVANCE='NO') '  Writing 1D integer variable ', varname

! get variable id and dimension length

  istatus = NF_INQ_VARID(nout,varname,varid)
  CALL nf_handle_error(istatus,'NF_INQ_VARID in write1d.')

  istatus = NF_INQ_VARDIMID(nout,varid,dim_ids)

  istatus = NF_INQ_DIMLEN(nout,dim_ids(1),dimlens(1))
  istatus = NF_INQ_DIMLEN(nout,dim_ids(2),dimlens(2))   ! unlimit dimension

  IF(dimlens(1) /= ndim) THEN
    WRITE(6,'(/a)') 'Mismatched dimension size.'
    STOP
  END IF

! Write data

  istatus = NF_PUT_VARA_INT(nout,varid,(/1,1/),(/ndim,1/),var1di)
  CALL nf_handle_error(istatus,'NF_PUT_VARA_INT in write1di.')

!  WRITE(6,'(a)') '         ...  DONE.'

  RETURN
END SUBROUTINE write_ncd_1di
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE write_ncd_2d                  ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE write_ncd_2d(nout,varname,var2d,ndimx,ndimy,istatus)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write 2D array to the output file.
!
!------------------------------------------------------------------
  USE wrf_metadata
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: nout
  INTEGER, INTENT(IN) :: ndimx,ndimy
  REAL,    INTENT(IN) :: var2d(ndimx,ndimy)
  CHARACTER(LEN=*), INTENT(IN) :: varname
  INTEGER, INTENT(OUT) :: istatus

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INTEGER :: varid
  INTEGER :: dim_ids(3)
  INTEGER :: dimlens(3)

  INCLUDE 'netcdf.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  WRITE(6,FMT='(2a)',ADVANCE='NO') '  Writing 2D variable ', varname

! get variable id and dimension length

  istatus = NF_INQ_VARID(nout,varname,varid)
  CALL nf_handle_error(istatus,'NF_INQ_VARID in write2d.')

  istatus = NF_INQ_VARDIMID(nout,varid,dim_ids)

  istatus = NF_INQ_DIMLEN(nout,dim_ids(1),dimlens(1))
  istatus = NF_INQ_DIMLEN(nout,dim_ids(2),dimlens(2))
  istatus = NF_INQ_DIMLEN(nout,dim_ids(3),dimlens(3))   ! unlimit dimension

  IF(dimlens(1) /= ndimx) THEN
    WRITE(6,'(/a)') 'Mismatched dimension size in X direction.'
    WRITE(6,*) 'Input X dimension = ',ndimx, ' Data X dimension =',dimlens(1)
    STOP
  END IF

  IF(dimlens(2) /= ndimy) THEN
    WRITE(6,'(/a)') 'Mismatched dimension size in Y direction.'
    WRITE(6,*) 'Input Y dimension = ',ndimx, ' Data Y dimension =',dimlens(2)
    STOP
  END IF

! Write data

  istatus = NF_PUT_VARA_REAL(nout,varid,(/1,1,1/),                     &
                              (/ndimx,ndimy,1/),var2d)
  CALL nf_handle_error(istatus,'NF_PUT_VARA_REAL in write2d')

!  WRITE(6,'(a)') '         ...  DONE.'

  RETURN
END SUBROUTINE write_ncd_2d
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE write_ncd_2di                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE write_ncd_2di(nout,varname,var2d,ndimx,ndimy,istatus)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write 2D array to the output file.
!
!------------------------------------------------------------------
  USE wrf_metadata
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: nout
  INTEGER, INTENT(IN) :: ndimx,ndimy
  INTEGER, INTENT(IN) :: var2d(ndimx,ndimy)
  CHARACTER(LEN=*), INTENT(IN) :: varname
  INTEGER, INTENT(OUT):: istatus

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INTEGER :: varid
  INTEGER :: dim_ids(3)
  INTEGER :: dimlens(3)

  INCLUDE 'netcdf.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  WRITE(6,FMT='(2a)',ADVANCE='NO') '  Writing 2D integer variable ', varname

! get variable id and dimension length

  istatus = NF_INQ_VARID(nout,varname,varid)
  CALL nf_handle_error(istatus,'NF_INQ_VARID in write2di.')

  istatus = NF_INQ_VARDIMID(nout,varid,dim_ids)

  istatus = NF_INQ_DIMLEN(nout,dim_ids(1),dimlens(1))
  istatus = NF_INQ_DIMLEN(nout,dim_ids(2),dimlens(2))
  istatus = NF_INQ_DIMLEN(nout,dim_ids(3),dimlens(3))   ! unlimit dimension

  IF(dimlens(1) /= ndimx) THEN
    WRITE(6,'(/a)') 'Mismatched dimension size in X direction.'
    WRITE(6,*) 'Input X dimension = ',ndimx, ' Data X dimension =',dimlens(1)
    STOP
  END IF

  IF(dimlens(2) /= ndimy) THEN
    WRITE(6,'(/a)') 'Mismatched dimension size in Y direction.'
    STOP
  END IF

! Write data

  istatus = NF_PUT_VARA_INT(nout,varid,(/1,1,1/),                     &
                              (/ndimx,ndimy,1/),var2d)
  CALL nf_handle_error(istatus,'NF_PUT_VARA_INT in write2di.')

!  WRITE(6,'(a)') '         ...  DONE.'

  RETURN
END SUBROUTINE write_ncd_2di
!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE write_ncd_3d                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE write_ncd_3d(nout,varname,var3d,ndimx,ndimy,ndimz,istatus)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write 3D array to the output file.
!
!------------------------------------------------------------------
  USE wrf_metadata

  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: nout   
  INTEGER, INTENT(IN) :: ndimx,ndimy,ndimz
  REAL,    INTENT(IN) :: var3d(ndimx,ndimy,ndimz)
  INTEGER, INTENT(OUT):: istatus

  CHARACTER(LEN=*), INTENT(IN) :: varname

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INTEGER :: varid
  INTEGER :: dim_ids(4)
  INTEGER :: dimlens(4)

  INCLUDE 'netcdf.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  WRITE(6,FMT='(2a)',ADVANCE='NO') '  Writing 3D variable ', varname

! get variable id and dimension length

  istatus = NF_INQ_VARID(nout,varname,varid)
  CALL nf_handle_error(istatus,'NF_INQ_VARID in write3d.')

  istatus = NF_INQ_VARDIMID(nout,varid,dim_ids)

  istatus = NF_INQ_DIMLEN(nout,dim_ids(1),dimlens(1))
  istatus = NF_INQ_DIMLEN(nout,dim_ids(2),dimlens(2))
  istatus = NF_INQ_DIMLEN(nout,dim_ids(3),dimlens(3))
  istatus = NF_INQ_DIMLEN(nout,dim_ids(4),dimlens(4))   ! unlimit dimension

  IF(dimlens(1) /= ndimx) THEN
    WRITE(6,'(/a)') 'Mismatched dimension size in X direction.'
    STOP
  END IF

  IF(dimlens(2) /= ndimy) THEN
    WRITE(6,'(/a)') 'Mismatched dimension size in Y direction.'
    STOP
  END IF

  IF(dimlens(3) /= ndimz) THEN
    WRITE(6,'(/a)') 'Mismatched dimension size in the 3rd dimension.'
    STOP
  END IF

! Write data

  istatus = NF_PUT_VARA_REAL(nout,varid,(/1,1,1,1/),                     &
                              (/ndimx,ndimy,ndimz,1/),var3d)
  CALL nf_handle_error(istatus,'NF_PUT_VARA_REAL in write3d')

!  WRITE(6,'(a)') '         ...  DONE.'

  RETURN
END SUBROUTINE write_ncd_3d
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE write_ncd_bdy                 ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE write_ncd_bdy(nout,ndimx,ndimy,ndimz,ndimbdy,ndimt,          &
                    varname,bdys,bdyn,bdyw,bdye,istatus)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write the 4 lateral boudnary arrays
!
!------------------------------------------------------------------
  USE wrf_metadata

  IMPLICIT NONE
  
  INTEGER, INTENT(IN)          :: nout   
  INTEGER, INTENT(IN)          :: ndimx,ndimy,ndimz,ndimbdy,ndimt
  CHARACTER(LEN=*), INTENT(IN) :: varname
  REAL,    INTENT(IN)          :: bdys(ndimx,ndimz,ndimbdy)
  REAL,    INTENT(IN)          :: bdyn(ndimx,ndimz,ndimbdy)
  REAL,    INTENT(IN)          :: bdyw(ndimy,ndimz,ndimbdy)
  REAL,    INTENT(IN)          :: bdye(ndimy,ndimz,ndimbdy)
  INTEGER, INTENT(OUT)         :: istatus

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INTEGER :: varid(4)
  INTEGER :: dimlens(4)
  INTEGER :: dim_ids(4)

  INTEGER :: n
  CHARACTER(LEN=30) :: vname

  INCLUDE 'netcdf.inc'

!                                            West  East  South North
  CHARACTER(LEN=2), PARAMETER :: appd(4) = (/'XS', 'XE', 'YS', 'YE'/)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  vname(1:30) = ' '

  DO n = 1, 4

    vname = varname//appd(n)
    WRITE(6,*) ' Writing lateral boundary arrays - ', vname

    ! get variable id and dimension length

    istatus = NF_INQ_VARID(nout,TRIM(vname),varid(n))
    CALL nf_handle_error(istatus,'NF_INQ_VARID in writebdy.')

    istatus = NF_INQ_VARDIMID(nout,varid,dim_ids)

    istatus = NF_INQ_DIMLEN(nout,dim_ids(1),dimlens(1))
    istatus = NF_INQ_DIMLEN(nout,dim_ids(2),dimlens(2))
    istatus = NF_INQ_DIMLEN(nout,dim_ids(3),dimlens(3))
    istatus = NF_INQ_DIMLEN(nout,dim_ids(4),dimlens(4))   ! unlimit dimension

    IF(dimlens(2) /= ndimz) THEN
      WRITE(6,*) 'Mismatched dimension size in Z direction.'
      STOP
    END IF

    IF(dimlens(3) /= ndimbdy) THEN
      WRITE(6,*) 'Mismatched dimension size in bdy_width.'
      STOP
    END IF

!    IF(dimlens(4) < ndimt) THEN
!      WRITE(6,*) 'Mismatched dimension size in bdy_width.'
!      STOP
!    END IF
  END DO


! Write data

  istatus = NF_PUT_VARA_REAL(nout,varid(1),(/1,1,1,ndimt/),               &
                             (/ndimy,ndimz,ndimbdy,1/),bdyw)    ! WEST
  CALL nf_handle_error(istatus,'NF_PUT_VARA_REAL in writebdy')

  istatus = NF_PUT_VARA_REAL(nout,varid(2),(/1,1,1,ndimt/),               &
                             (/ndimy,ndimz,ndimbdy,1/),bdye)    ! EAST
  CALL nf_handle_error(istatus,'NF_PUT_VARA_REAL in writebdy')

  istatus = NF_PUT_VARA_REAL(nout,varid(3),(/1,1,1,ndimt/),               &
                             (/ndimx,ndimz,ndimbdy,1/),bdys)    ! SOUTH
  CALL nf_handle_error(istatus,'NF_PUT_VARA_REAL in writebdy')

  istatus = NF_PUT_VARA_REAL(nout,varid(4),(/1,1,1,ndimt/),               &
                             (/ndimx,ndimz,ndimbdy,1/),bdyn)    ! NORTH
  CALL nf_handle_error(istatus,'NF_PUT_VARA_REAL in writebdy')

  RETURN
END SUBROUTINE write_ncd_bdy
!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE write_ncd_bdy2d               ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE write_ncd_bdy2d(nout,ndimx,ndimy,ndimbdy,ndimt,          &
                      varname,bdys,bdyn,bdyw,bdye,istatus)
!
!------------------------------------------------------------------
!
!  PURPOSE:
!
!    Write the 4 lateral boudnary arrays
!
!------------------------------------------------------------------
  USE wrf_metadata

  IMPLICIT NONE
  
  INTEGER, INTENT(IN)          :: nout   
  INTEGER, INTENT(IN)          :: ndimx,ndimy,ndimbdy,ndimt
  CHARACTER(LEN=*), INTENT(IN) :: varname
  REAL,    INTENT(IN)          :: bdys(ndimx,ndimbdy)
  REAL,    INTENT(IN)          :: bdyn(ndimx,ndimbdy)
  REAL,    INTENT(IN)          :: bdyw(ndimy,ndimbdy)
  REAL,    INTENT(IN)          :: bdye(ndimy,ndimbdy)
  INTEGER, INTENT(OUT)         :: istatus

!------------------------------------------------------------------
!
!  Misc. local variable
!
!------------------------------------------------------------------

  INTEGER :: varid(4)
  INTEGER :: dimlens(3)
  INTEGER :: dim_ids(3)

  INTEGER :: n
  CHARACTER(LEN=30) :: vname

  INCLUDE 'netcdf.inc'

!                                            West  East  South North
  CHARACTER(LEN=2), PARAMETER :: appd(4) = (/'XS', 'XE', 'YS', 'YE'/)
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  vname(1:30) = ' '

  DO n = 1, 4

    vname = varname//appd(n)
    WRITE(6,*) ' Writing lateral boundary arrays - ', vname

    ! get variable id and dimension length

    istatus = NF_INQ_VARID(nout,TRIM(vname),varid(n))
    CALL nf_handle_error(istatus,'NF_INQ_VARID in writebdy2d.')

    istatus = NF_INQ_VARDIMID(nout,varid,dim_ids)

    istatus = NF_INQ_DIMLEN(nout,dim_ids(1),dimlens(1))
    istatus = NF_INQ_DIMLEN(nout,dim_ids(2),dimlens(2))
    istatus = NF_INQ_DIMLEN(nout,dim_ids(3),dimlens(3))   ! unlimit dimension

    IF(dimlens(1) /= ndimx .AND. dimlens(1) /= ndimy) THEN
      WRITE(6,*) 'Mismatched dimension size in X/Y direction.'
      STOP
    END IF

    IF(dimlens(2) /= ndimbdy) THEN
      WRITE(6,*) 'Mismatched dimension size in bdy_width.'
      STOP
    END IF

  END DO


! Write data

  istatus = NF_PUT_VARA_REAL(nout,varid(1),(/1,1,ndimt/),               &
                                   (/ndimy,ndimbdy,1/),bdyw)    ! WEST
  CALL nf_handle_error(istatus,'NF_PUT_VARA_REAL in writebdy2d')

  istatus = NF_PUT_VARA_REAL(nout,varid(2),(/1,1,ndimt/),               &
                                   (/ndimy,ndimbdy,1/),bdye)    ! EAST
  CALL nf_handle_error(istatus,'NF_PUT_VARA_REAL in writebdy')

  istatus = NF_PUT_VARA_REAL(nout,varid(3),(/1,1,ndimt/),               &
                                   (/ndimx,ndimbdy,1/),bdys)    ! SOUTH
  CALL nf_handle_error(istatus,'NF_PUT_VARA_REAL in writebdy')

  istatus = NF_PUT_VARA_REAL(nout,varid(4),(/1,1,ndimt/),               &
                                   (/ndimx,ndimbdy,1/),bdyn)    ! NORTH
  CALL nf_handle_error(istatus,'NF_PUT_VARA_REAL in writebdy')

  RETURN
END SUBROUTINE write_ncd_bdy2d

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
  CALL nf_handle_error(istatus,'NF_INQ_VARID in get_wrf_2d')

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
  IF(istatus /= NF_NOERR .AND. istatus /= NF_EEXIST) THEN
    PRINT *, 'nx,ny,nz,itime = ',nx,ny,nz,itime
    CALL nf_handle_error(istatus,'NF_GET_VARA_REAL in get_wrf_3d for '//TRIM(varname))
  END IF

  RETURN
END SUBROUTINE get_ncd_3d
