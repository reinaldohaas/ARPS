
!########################################################################
!########################################################################
!######                                                            ######
!######              SUBROUTINE RD_VERIF_STATS_DIFF                ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################
  
SUBROUTINE rd_verif_stats_diff (filename, if_diff, diffile, date, &
                     amodel,agrid,adx,anx,any, fmodel,fgrid,fdx, &
                     vdx, vdy, vnx, vny, vmapproj, &
                     vtrulat1, vtrulat2, vtrulon, vsclfct, vctrlat, &
                     vctrlon, &
                     corner_lat, corner_lon,  &
                     nlevel, ntime, nvar_p, nvar_sfc, &
                     timesec, pressure, missing, mastercounter, &
                     counter_p, bias_p, rms_p, diff_p, &
                     counter_sfc, bias_sfc, rms_sfc, diff_sfc, &
                     varid_p, varname_p, varunit_p, &
                     varid_sfc, varname_sfc, varunit_sfc)

!----------------------------------------------------------------------- 
!
! PURPOSE:
!
! Reads in an HDF file containing verification statistics produced by 
! program verifgrid.  Based on subroutine wrt_verif_stats_diff.
!
! AUTHOR:  Eric Kemp, May 2000
!
!----------------------------------------------------------------------- 
!
! Variable declarations
!
!-----------------------------------------------------------------------

  USE hdf_constants

  IMPLICIT NONE

  INTEGER, INTENT(INOUT) :: date(6), vnx,vny, anx,any, vmapproj, &
      mastercounter,if_diff
  INTEGER, INTENT(INOUT) :: nlevel, ntime, nvar_p, nvar_sfc
  INTEGER, INTENT(OUT) :: timesec(ntime), pressure(nlevel)
  REAL, INTENT(OUT) :: vdx,vdy, adx, fdx, corner_lat(2,2),corner_lon(2,2),&
                     vtrulat1, vtrulat2, vtrulon, vsclfct, vctrlat,&
                     vctrlon, missing
  REAL, INTENT(OUT), DIMENSION(nlevel,ntime,nvar_p) :: &
      counter_p, bias_p, rms_p
  REAL, INTENT(OUT) :: diff_p(anx,any,nlevel,ntime,nvar_p), &
                    diff_sfc(anx,any,ntime,nvar_p)   
  REAL, INTENT(OUT), DIMENSION(ntime,nvar_sfc) :: &
      counter_sfc, bias_sfc, rms_sfc
  CHARACTER(LEN=*), INTENT(OUT) :: filename, diffile, &
      amodel, agrid, fmodel, fgrid
  CHARACTER(LEN=*), INTENT(OUT), DIMENSION(nvar_p) :: &
      varid_p, varname_p, varunit_p
  CHARACTER(LEN=*), INTENT(OUT), DIMENSION(nvar_sfc) :: &
      varid_sfc, varname_sfc, varunit_sfc

!-----------------------------------------------------------------------
!
! HDF variables and parameters.
!
!-----------------------------------------------------------------------
      
  INTEGER,PARAMETER :: maxrank = 5
  INTEGER, DIMENSION(maxrank) :: dims, start, edges, stride

  INTEGER :: sd_id,sds_id,sd_index

!-----------------------------------------------------------------------
!
! Miscellaneous variables
!
!-----------------------------------------------------------------------
  
  CHARACTER :: string*132
 
  INTEGER :: status, i
 
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!     
! Beginning of executable code... 
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
! Get the header information.
! 
!-----------------------------------------------------------------------

  CALL rd_verif_header (sd_id, missing, filename, date, &
                     amodel,agrid,adx, fmodel,fgrid,fdx, &
                     vdx, vdy, vnx, vny, vmapproj, &
                     vtrulat1, vtrulat2, vtrulon, vsclfct, vctrlat,&
                     vctrlon, &
                     corner_lat, corner_lon,  &
                     nlevel, ntime, nvar_p, nvar_sfc, &
                     timesec, pressure, mastercounter, &
                     counter_p, counter_sfc, &
                     varid_p, varname_p, varunit_p, &
                     varid_sfc, varname_sfc, varunit_sfc)

  dims(1) = nlevel
  dims(2) = ntime
  dims(3) = nvar_p
  sd_index = sfn2index(sd_id,'counter_p')      
  CALL read_sds (sd_id, sd_index, counter_p, 3, dims)
  sd_index = sfn2index(sd_id,'bias_p')
  CALL read_sds (sd_id, sd_index, bias_p, 3, dims)
  sd_index = sfn2index(sd_id,'rms_p')
  CALL read_sds (sd_id, sd_index, rms_p, 3, dims)

  dims(1) = ntime
  dims(2) = nvar_sfc
  sd_index = sfn2index(sd_id,'counter_sfc')
  CALL read_sds (sd_id, sd_index, counter_sfc, 2, dims)
  sd_index = sfn2index(sd_id,'bias_sfc')
  CALL read_sds (sd_id, sd_index, bias_sfc, 2, dims)
  sd_index = sfn2index(sd_id,'rms_sfc')
  CALL read_sds (sd_id, sd_index, rms_sfc, 2, dims)

  status = sfend(sd_id)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfend')

!-----------------------------------------------------------------------
!
! Read in the difference fields
!
!-----------------------------------------------------------------------

  IF (if_diff == 0) THEN
    PRINT *, 'rd_verif_stats_diff: Successful completion'
    RETURN
  END IF

  ! Note that diff_p and diff_sfc have dimensions anx,any
  PRINT *, 'rd_verif_stats_diff: Differences to file ', TRIM(diffile)
  CALL rd_verif_header ( sd_id, missing, diffile, date, &
                     amodel,agrid,adx, fmodel,fgrid,fdx, &
                     adx, adx, anx, any, vmapproj, &
                     vtrulat1, vtrulat2, vtrulon, vsclfct, vctrlat, &
                     vctrlon, &
                     corner_lat, corner_lon,  &
                     nlevel, ntime, nvar_p, nvar_sfc, &
                     timesec, pressure, mastercounter, &
                     counter_p, counter_sfc, &
                     varid_p, varname_p, varunit_p, &
                     varid_sfc, varname_sfc, varunit_sfc)
  PRINT *, 'rd_verif_stats_diff: Done call wrt_verif_header'

  dims(1) = anx
  dims(2) = any
  dims(3) = nlevel
  dims(4) = ntime
  dims(5) = nvar_p
  sd_index = sfn2index(sd_id,'diff_p')
  CALL read_sds (sd_id, sd_index, diff_p, 5, dims)

  dims(1) = anx
  dims(2) = any
  dims(3) = ntime
  dims(4) = nvar_sfc
  sd_index = sfn2index(sd_id,'diff_sfc')
  CALL read_sds (sd_id, sd_index, diff_sfc, 4, dims)

  status = sfend(sd_id)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfend')

  PRINT *, 'rd_verif_stats_diff: Successful completion'
  WRITE(6,*)'-----------------------------------------------------------'

END SUBROUTINE rd_verif_stats_diff

!########################################################################
!########################################################################
!######                                                            ######
!######                SUBROUTINE RD_VERIF_HEADER                  ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################

SUBROUTINE rd_verif_header (sd_id, missing, filename, date, &
                     amodel,agrid,adx, fmodel,fgrid,fdx, &
                     vdx, vdy, vnx, vny, vmapproj, &
                     vtrulat1, vtrulat2, vtrulon, vsclfct, vctrlat,&
                     vctrlon, &
                     corner_lat, corner_lon,  &
                     nlevel, ntime, nvar_p, nvar_sfc, &
                     timesec, pressure, mastercounter, &
                     counter_p, counter_sfc, &
                     varid_p, varname_p, varunit_p, &
                     varid_sfc, varname_sfc, varunit_sfc)

!----------------------------------------------------------------------- 
!
! PURPOSE:
!
! Reads header for verification statistics file in HDF format.  Based
! on rd_verif_header.
!     
! AUTHOR:  Eric Kemp, May 2000
!
!-----------------------------------------------------------------------
!
! Variable declarations
!
!-----------------------------------------------------------------------

  USE hdf_constants
                     
  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: sd_id
  INTEGER, INTENT(OUT) :: date(6), vnx,vny, vmapproj, mastercounter
  INTEGER, INTENT(IN) :: nlevel, ntime, nvar_p, nvar_sfc
  INTEGER, INTENT(OUT) :: timesec(ntime), pressure(nlevel)
  REAL, INTENT(OUT) :: vdx,vdy, adx, fdx, corner_lat(2,2), &
                     corner_lon(2,2),&
                     vtrulat1, vtrulat2, vtrulon, vsclfct, vctrlat,&
                     vctrlon, &
                     missing
  REAL, INTENT(OUT), DIMENSION(nlevel,ntime,nvar_p) :: counter_p
  REAL, INTENT(OUT), DIMENSION(ntime,nvar_sfc) :: counter_sfc
  CHARACTER(LEN=*), INTENT(OUT) :: filename, amodel, agrid, fmodel, &
                                   fgrid
  CHARACTER(LEN=*), INTENT(OUT), DIMENSION(nvar_p) :: &
      varid_p, varname_p, varunit_p
  CHARACTER(LEN=*), INTENT(IN), DIMENSION(nvar_sfc) :: &
      varid_sfc, varname_sfc, varunit_sfc

!-----------------------------------------------------------------------
!
! Miscellaneous variables
!
!-----------------------------------------------------------------------
  
  INTEGER :: sds_id, status, i, dims(5)
  CHARACTER :: string*132
  INTEGER :: attr_index
  CHARACTER*10 :: attr_name
  INTEGER :: data_type,n_values
  INTEGER :: sd_index

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!                    
! Open HDF file.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'-----------------------------------------------------------'
  PRINT *, 'wrt_verif_header: Creating ', TRIM(filename)
  sd_id = sfstart(filename, DFACC_READ)
  IF (sd_id == FAIL) CALL hdf_fail (1, sd_id, sd_id, 'sfstart')

!-----------------------------------------------------------------------
!
! Read file attributes describing the verification data.
!
!-----------------------------------------------------------------------
 
  attr_index =sffattr(sd_id,'forecast_source')
  status = sfgainfo(sd_id, attr_index, attr_name, data_type, n_values)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfgainfo')
  status = sfrattr(sd_id, attr_index, fmodel)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfrattr')
  WRITE(6,*) 'fmodel = ',TRIM(fmodel)

  attr_index =sffattr(sd_id,'forecast_grid')
  status = sfgainfo(sd_id, attr_index, attr_name, data_type, n_values)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfgainfo')
  status = sfrattr(sd_id, attr_index, fgrid)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfrattr')
  WRITE(6,*) 'fgrid = ',TRIM(fgrid)

  attr_index =sffattr(sd_id,'analysis_source')
  status = sfgainfo(sd_id, attr_index, attr_name, data_type, n_values)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfgainfo')
  status = sfrattr(sd_id, attr_index, amodel)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfrattr')
  WRITE(6,*) 'amodel = ',TRIM(amodel)

  attr_index =sffattr(sd_id,'analysis_grid')
  status = sfgainfo(sd_id, attr_index, attr_name, data_type, n_values)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfgainfo')
  status = sfrattr(sd_id, attr_index, agrid)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfrattr')
  WRITE(6,*) 'agrid = ',TRIM(agrid)

  attr_index =sffattr(sd_id,'start_date')
  status = sfgainfo(sd_id, attr_index, attr_name, data_type, n_values)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfgainfo')
  status = sfrattr(sd_id, attr_index, date)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfrattr')
  WRITE(6,*) 'Start Date = ',date

  attr_index =sffattr(sd_id,'vdx')
  status = sfgainfo(sd_id, attr_index, attr_name, data_type, n_values)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfgainfo')
  status = sfrattr(sd_id, attr_index, vdx)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfrattr')
  WRITE(6,*) 'vdx = ',vdx

  attr_index =sffattr(sd_id,'vdy')
  status = sfgainfo(sd_id, attr_index, attr_name, data_type, n_values)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfgainfo')
  status = sfrattr(sd_id, attr_index, vdy)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfrattr')
  WRITE(6,*) 'vdy = ',vdy

  attr_index =sffattr(sd_id,'vnx')
  status = sfgainfo(sd_id, attr_index, attr_name, data_type, n_values)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfgainfo')
  status = sfrattr(sd_id, attr_index, vnx)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfrattr')
  WRITE(6,*) 'vnx = ',vnx

  attr_index =sffattr(sd_id,'vny')
  status = sfgainfo(sd_id, attr_index, attr_name, data_type, n_values)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfgainfo')
  status = sfrattr(sd_id, attr_index, vny)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfrattr')
  WRITE(6,*) 'vny = ',vny

  attr_index =sffattr(sd_id,'vmapproj')
  status = sfgainfo(sd_id, attr_index, attr_name, data_type, n_values)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfgainfo')
  status = sfrattr(sd_id, attr_index, vmapproj)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfrattr')
  WRITE(6,*) 'vmapproj = ',vmapproj

  attr_index =sffattr(sd_id,'vtrulat1')
  status = sfgainfo(sd_id, attr_index, attr_name, data_type, n_values)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfgainfo')
  status = sfrattr(sd_id, attr_index, vtrulat1)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfrattr')
  WRITE(6,*) 'vtrulat1 = ',vtrulat1

  attr_index =sffattr(sd_id,'vtrulat2')
  status = sfgainfo(sd_id, attr_index, attr_name, data_type, n_values)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfgainfo')
  status = sfrattr(sd_id, attr_index, vtrulat2)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfrattr')
  WRITE(6,*) 'vtrulat2 = ',vtrulat2

  attr_index =sffattr(sd_id,'vtrulon')
  status = sfgainfo(sd_id, attr_index, attr_name, data_type, n_values)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfgainfo')
  status = sfrattr(sd_id, attr_index, vtrulon)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfrattr')
  WRITE(6,*) 'vtrulon = ',vtrulon

  attr_index =sffattr(sd_id,'vctrlat')
  status = sfgainfo(sd_id, attr_index, attr_name, data_type, n_values)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfgainfo')
  status = sfrattr(sd_id, attr_index, vctrlat)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfrattr')
  WRITE(6,*) 'vctrlat = ',vctrlat

  attr_index =sffattr(sd_id,'vctrlon')
  status = sfgainfo(sd_id, attr_index, attr_name, data_type, n_values)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfgainfo')
  status = sfrattr(sd_id, attr_index, vctrlon)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfrattr')
  WRITE(6,*) 'vctrlon = ',vctrlon

  attr_index =sffattr(sd_id,'adx')
  status = sfgainfo(sd_id, attr_index, attr_name, data_type, n_values)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfgainfo')
  status = sfrattr(sd_id, attr_index, adx)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfrattr')
  WRITE(6,*) 'adx = ',adx

  attr_index =sffattr(sd_id,'fdx')
  status = sfgainfo(sd_id, attr_index, attr_name, data_type, n_values)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfgainfo')
  status = sfrattr(sd_id, attr_index, fdx)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfrattr')
  WRITE(6,*) 'fdx = ',fdx

  attr_index =sffattr(sd_id,'mastercounter')
  status = sfgainfo(sd_id, attr_index, attr_name, data_type, n_values)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfgainfo')
  status = sfrattr(sd_id, attr_index, mastercounter)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfrattr')
  WRITE(6,*) 'mastercounter = ',mastercounter

!-----------------------------------------------------------------------
!
! Corner lat/lon, time
!
!-----------------------------------------------------------------------

  sd_index = sfn2index(sd_id,'corner_lat')
  dims(:) = 2
  CALL read_sds (sd_id,sd_index,corner_lat,2,dims) 
  
  sd_index = sfn2index(sd_id,'corner_lon')
  dims(:) = 2
  CALL read_sds (sd_id,sd_index,corner_lon,2,dims)

  sd_index = sfn2index(sd_id,'pressure')
  CALL read_sds (sd_id,sd_index,pressure,1,nlevel)

  sd_index = sfn2index(sd_id,'time')
  CALL read_sds (sd_id,sd_index,timesec,1,ntime)

!-----------------------------------------------------------------------
!
! Get data descriptions
!
!-----------------------------------------------------------------------

  sd_index = sfn2index(sd_id,'id_p')
  dims(2) = nvar_p
  dims(1) = LEN(varid_p(1))
  CALL read_sds_char (sd_id, sd_index, varid_p, 2, dims)
  DO i=1,nvar_p
    PRINT *, i, TRIM(varid_p(i))
  END DO

  sd_index = sfn2index(sd_id,'units_p')      
  dims(1) = LEN(varunit_p(1))
  CALL read_sds_char (sd_id, sd_index, varunit_p, 2, dims) 
  DO i=1,nvar_p
    PRINT *, i, TRIM(varunit_p(i))
  END DO

  sd_index = sfn2index(sd_id,'names_p')      
  dims(1) = LEN(varname_p(1))
  CALL read_sds_char (sd_id, sd_index, varname_p, 2, dims) 
  DO i=1,nvar_p
    PRINT *, i, TRIM(varname_p(i))
  END DO

  sd_index = sfn2index(sd_id,'id_sfc')
  dims(2) = nvar_sfc
  dims(1) = LEN(varid_sfc(1))
  CALL read_sds_char (sd_id, sd_index, varid_sfc, 2, dims)
  DO i=1,nvar_sfc
    PRINT *, i, TRIM(varid_sfc(i))
  END DO

  sd_index = sfn2index(sd_id,'units_sfc')
  dims(1) = LEN(varunit_sfc(1))
  CALL read_sds_char (sd_id, sd_index, varunit_sfc, 2, dims)
  DO i=1,nvar_sfc
    PRINT *, i, TRIM(varunit_sfc(i))
  END DO

  sd_index = sfn2index(sd_id,'names_sfc')
  dims(1) = LEN(varname_sfc(1))
  CALL read_sds_char (sd_id, sd_index, varname_sfc, 2, dims)
  DO i=1,nvar_sfc
    PRINT *, i, TRIM(varname_sfc(i))
  END DO

  PRINT *, 'rd_verif_header: Done'
  WRITE(6,*)'-----------------------------------------------------------'

END SUBROUTINE rd_verif_header

