!########################################################################
!########################################################################
!######                                                            ######
!######                SUBROUTINE RD_VERIF_DIMS                    ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################

SUBROUTINE rd_verif_dims(nx,ny,nlevel,ntime,nvar_p,nvar_sfc,filename)

!-----------------------------------------------------------------------
! 
! PURPOSE:
! 
! Retrieves dimensions of verification data in HDF files.
!
! AUTHOR:  Eric Kemp, March 2000
!
!----------------------------------------------------------------------- 

  USE hdf_constants

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! Variable declarations
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(INOUT) :: nx,ny,nlevel,ntime,nvar_p,nvar_sfc
  CHARACTER*(*), INTENT(IN) :: filename
 
!-----------------------------------------------------------------------
!
! HDF parameters and variables
!
!-----------------------------------------------------------------------

  INTEGER,PARAMETER :: maxrank = 5
  INTEGER, DIMENSION(maxrank) :: dims, start, edges, stride
  INTEGER :: sd_id,sds_id
 
!-----------------------------------------------------------------------
!
! Other variables
!
!-----------------------------------------------------------------------
  
  INTEGER :: i
  INTEGER,PARAMETER :: missing = -9999
 
  INTEGER, PARAMETER :: max_metadata = 6
  CHARACTER*10, PARAMETER :: metadata(max_metadata) =       &
     (/'nx        ','ny        ','nlevel    ','ntime     ', &
       'nvar_p    ','nvar_sfc  '/)
  CHARACTER*10 :: attr_name
  INTEGER :: attr_index,status,n_file_attrs,data_type,n_datasets,n_values
  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
! Open HDF file.
!
!-----------------------------------------------------------------------
 
  WRITE(6,*)'-----------------------------------------------------------' 
  PRINT *, 'rd_verif_dims: Opening ', TRIM(filename)
  sd_id = sfstart(filename, DFACC_READ)
  IF (sd_id == FAIL) CALL hdf_fail (1, sd_id, sd_id, 'sfstart')

!-----------------------------------------------------------------------
!
! Determine the number of data sets in the file and the number of
! file attributes.
!
!-----------------------------------------------------------------------

  status = sffinfo(sd_id, n_datasets, n_file_attrs)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sffinfo')
  PRINT *, 'rd_verif_dims: There are', n_datasets, 'datasets and', &
           n_file_attrs, 'file attributes'
  
!-----------------------------------------------------------------------
!
! Read file attributes describing the model data.
!
!-----------------------------------------------------------------------
 
  nx = missing
  ny = missing
  nlevel = missing
  nvar_p = missing
  nvar_sfc = missing
           
  DO i = 1,max_metadata

    PRINT *, 'rd_verif_dims: Getting attribute ', i, TRIM(metadata(i))
    attr_index = sffattr(sd_id,metadata(i))
    status = sfgainfo(sd_id, attr_index, attr_name, data_type, n_values)
    IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfgainfo')
 
    SELECT CASE (TRIM(attr_name))
      CASE ("nx")
        status = sfrattr(sd_id, attr_index, nx)
      CASE ("ny")
        status = sfrattr(sd_id, attr_index, ny)
      CASE ("nlevel")
        status = sfrattr(sd_id, attr_index, nlevel)
      CASE ("ntime")
        status = sfrattr(sd_id, attr_index, ntime)
      CASE ("nvar_p")   
        status = sfrattr(sd_id, attr_index, nvar_p)
      CASE ("nvar_sfc")
        status = sfrattr(sd_id, attr_index, nvar_sfc)
    END SELECT
    IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfrattr')
  END DO
    
  IF (nx == missing .OR. ny == missing .OR. nlevel == missing .OR. &
      ntime == missing .OR. nvar_p == missing .OR. nvar_sfc == missing) &
  THEN
    WRITE(6,*)'RD_VERIF_DIMS:  Missing dimension data.  Aborting...'
    WRITE(6,*)'nx = ',nx,' ny = ',ny,' nlevel = ',nlevel
    WRITE(6,*)'ntime = ',ntime,' nvar_p = ',nvar_p,' nvar_sfc = ',nvar_sfc
    STOP
  END IF

  WRITE(6,*)'rd_verif_dims: Done'
  WRITE(6,*)'-----------------------------------------------------------' 
      
END SUBROUTINE rd_verif_dims

!########################################################################
!########################################################################
!######                                                            ######
!######                 SUBROUTINE WRTQPFSTATS                     ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################
  
SUBROUTINE WRTQPFSTATS(nx,ny,ntime,nthresh,nstat,missing,               &
                       filename,model,grid,                             &
                       date,timesec,mapproj,scale,trulat1,trulat2,      &
                       trulon,dx,dy,scswlat,scswlon,                    &
                       mapswlat,mapswlon,mapnwlat,mapnwlon,             &
                       mapselat,mapselon,mapnelat,mapnelon,             &
                       thresh,threshunit,                               &
                       ctstatsid,ctstatsname,ctstats,                   &
                       CTA,CTB,CTC,CTD)

!----------------------------------------------------------------------- 
!
! PURPOSE:
!
! Write contingency table statistics to HDF file.
!
! AUTHOR:  Eric Kemp, March 2000
!
! MODIFICATION HISTORY:
!
! Eric Kemp, 3/31/2000
! Added lat/lon coordinates of corner points for NCL mapping.
!
!----------------------------------------------------------------------- 
  
  USE hdf_constants
 
  IMPLICIT NONE
 
!----------------------------------------------------------------------- 
!
! Variable declarations
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nx,ny,nstat,ntime,nthresh
  REAL, INTENT(IN) :: missing
  CHARACTER(LEN=*), INTENT(IN) :: filename
  CHARACTER(LEN=*), INTENT(IN) :: model,grid
  INTEGER, INTENT(IN) :: date(6),timesec(ntime),mapproj
  REAL, INTENT(IN) :: scale,trulat1,trulat2,trulon,dx,dy
  REAL, INTENT(IN) :: mapswlat,mapswlon,mapnwlat,mapnwlon, &
                      mapselat,mapselon,mapnelat,mapnelon
  REAL, INTENT(IN) :: scswlat,scswlon
  REAL, INTENT(IN), DIMENSION(nthresh) :: thresh
  CHARACTER(LEN=*), INTENT(IN) :: threshunit
  CHARACTER(LEN=*),INTENT(IN),DIMENSION(nstat) :: ctstatsid, &
       ctstatsname
  REAL, INTENT(IN) :: ctstats(ntime,nthresh,nstat)
  REAL, INTENT(IN), DIMENSION(ntime,nthresh) :: CTA,CTB,CTC,CTD

!-----------------------------------------------------------------------
!
! HDF parameters and variables
!
!-----------------------------------------------------------------------
  
  INTEGER,PARAMETER :: maxrank = 5
  INTEGER, DIMENSION(maxrank) :: dims, start, edges, stride
  INTEGER :: sd_id,sds_id
  
!-----------------------------------------------------------------------
!
! Other variables
!
!-----------------------------------------------------------------------
  
  CHARACTER*256 :: string
  INTEGER :: i, j, status

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! 
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
!-----------------------------------------------------------------------
!
! Open HDF file.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'-----------------------------------------------------------'
  PRINT *, 'wrtqpfstats: Opening ', TRIM(filename)
  sd_id = sfstart(filename, DFACC_CREATE)
  IF (sd_id == FAIL) CALL hdf_fail (1, sd_id, sd_id, 'sfstart')

!-----------------------------------------------------------------------
! 
! Add scalar attributes describing the model data.
!
!-----------------------------------------------------------------------
  
  PRINT *, 'wrtqpfstats: Adding scalar attributes'
  status = sfsnatt(sd_id, 'start_date', DFNT_INT32, 6, date)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfscatt(sd_id, 'model', DFNT_CHAR8,LEN_TRIM(model), model)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfscatt')
  status = sfscatt(sd_id, 'grid', DFNT_CHAR8,LEN_TRIM(grid), grid)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfscatt')
  status = sfsnatt(sd_id, 'dx', DFNT_FLOAT32, 1, dx)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'dy', DFNT_FLOAT32, 1, dy)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'nx', DFNT_INT32, 1, nx)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'ny', DFNT_INT32, 1, ny)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'nstat', DFNT_INT32, 1, nstat)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'ntime', DFNT_INT32, 1, ntime)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'nthresh', DFNT_INT32, 1, nthresh)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'mapproj', DFNT_INT32, 1, mapproj)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'scale', DFNT_FLOAT32, 1, scale)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'trulon', DFNT_FLOAT32, 1, trulon)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'trulat1', DFNT_FLOAT32, 1, trulat1)   
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'trulat2', DFNT_FLOAT32, 1, trulat2)   
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'scswlat', DFNT_FLOAT32, 1, scswlat)   
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'scswlon', DFNT_FLOAT32, 1, scswlon)   
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'mapswlat', DFNT_FLOAT32, 1, mapswlat)   
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'mapswlon', DFNT_FLOAT32, 1, mapswlon)   
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'mapnwlat', DFNT_FLOAT32, 1, mapnwlat)   
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'mapnwlon', DFNT_FLOAT32, 1, mapnwlon)   
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'mapselat', DFNT_FLOAT32, 1, mapselat)   
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'mapselon', DFNT_FLOAT32, 1, mapselon)   
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'mapnelat', DFNT_FLOAT32, 1, mapnelat)   
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'mapnelon', DFNT_FLOAT32, 1, mapnelon)   
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
    
!-----------------------------------------------------------------------
!
! Write the time datasets.
!
!-----------------------------------------------------------------------

  PRINT *, 'wrtqpfstats: Creating time', ntime
  CALL write_sds (sds_id, sd_id, DFNT_INT32, 'time', timesec, 1, ntime, &
                  missing)
  
  string = 'time'
  status = sfscatt(sds_id,'long_name',DFNT_CHAR8,LEN_TRIM(string), &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 's'
  status = sfscatt(sds_id,'units', DFNT_CHAR8, LEN_TRIM(string), &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 'seconds_from_init_time'
  status=sfscatt(sds_id,'description',DFNT_CHAR8,LEN_TRIM(string), &
                 TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')
  
  status = sfendacc(sds_id)
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfendacc')

!-----------------------------------------------------------------------
!
! Write the thresholds
!
!-----------------------------------------------------------------------
  
  PRINT *, 'wrtqpfstats: Creating thresholds', nthresh
  CALL write_sds (sds_id, sd_id, DFNT_INT32, 'thresh', thresh, 1, &
                  nthresh, missing)
  
  string = 'QPF_Threshold'
  status = sfscatt(sds_id,'long_name',DFNT_CHAR8,LEN_TRIM(string), &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')
  
  string = threshunit
  status = sfscatt(sds_id,'units', DFNT_CHAR8, LEN_TRIM(string), &
                   TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 'thresholds_for_QPF_statistics'
  status=sfscatt(sds_id,'description',DFNT_CHAR8,LEN_TRIM(string), &
                 TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')
  
  status = sfendacc(sds_id)
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfendacc')
  
!-----------------------------------------------------------------------
!
! Dump the contingency table scores.
! 
!-----------------------------------------------------------------------

  DO i=1,nstat
    print *, i, ctstatsid(i), ctstatsname(i)
  END DO
  
  dims(2) = nstat
  dims(1) = LEN(ctstatsid(1)) 
  CALL write_sds_char (sds_id, sd_id, 'ctstatsid', ctstatsid, 2, dims)
  dims(1) = LEN(ctstatsname(1))
  CALL write_sds_char (sds_id, sd_id, 'ctstatsnames',ctstatsname, 2, dims)

  dims(1) = ntime
  dims(2) = nthresh
  dims(3) = nstat

  PRINT *, 'wrtqpfstat: Creating ctstats', dims(1:3)
  CALL write_sds (sds_id, sd_id, DFNT_FLOAT32,'ctstats',ctstats,3,dims,&
                  missing)

  status = sfendacc(sds_id)
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfendacc')
  
!-----------------------------------------------------------------------
!
! Now dump the contingency table forecast/observation pairs.
!
!-----------------------------------------------------------------------  

  dims(1) = ntime
  dims(2) = nthresh
  
  PRINT *, 'wrtqpfstat: Creating CTA', dims(1:2)
  CALL write_sds (sds_id, sd_id, DFNT_FLOAT32,'CTA',CTA,2,dims,&
                  missing)
  PRINT *, 'wrtqpfstat: Creating CTB', dims(1:2)
  CALL write_sds (sds_id, sd_id, DFNT_FLOAT32,'CTB',CTB,2,dims,&
                  missing) 
  PRINT *, 'wrtqpfstat: Creating CTC', dims(1:2)
  CALL write_sds (sds_id, sd_id, DFNT_FLOAT32,'CTC',CTC,2,dims,&
                  missing)
  PRINT *, 'wrtqpfstat: Creating CTD', dims(1:2)
  CALL write_sds (sds_id, sd_id, DFNT_FLOAT32,'CTD',CTD,2,dims,&
                  missing)

  status = sfendacc(sds_id)
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfendacc')
  
!-----------------------------------------------------------------------
!
! Close file and exit.
!
!-----------------------------------------------------------------------

  status = sfend(sd_id)    
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfend')
  
  WRITE(6,*)'wrtqpfstats: Done'
  WRITE(6,*)'-----------------------------------------------------------'
  
END SUBROUTINE wrtqpfstats


!########################################################################
!########################################################################
!######                                                            ######
!######                  SUBROUTINE WRTVERIF                       ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################

!TODO: this routine should use wrt_verif_hdr
SUBROUTINE wrtverif ( nx,ny,nlevel,ntime,nvar_p,nvar_sfc, missing,	&
                      filename,model,grid,date,timesec,pressure,        &
                      mapproj, scale, trulat1,trulat2,trulon,           &
                      dx,dy,scswlat,scswlon, &  
                      mapswlat,mapswlon,mapnwlat,mapnwlon,        &
                      mapselat,mapselon,mapnelat,mapnelon,              &
	              flag_p,varid_p, varname_p, varunit_p, var_p,      &
	              flag_sfc,varid_sfc, varname_sfc, varunit_sfc,     &
                      var_sfc)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Writes gridded verification file in HDF format.
!
! AUTHOR:  Eric Kemp, November 1999
!
! Richard Carpenter, January 2000: Extensive revisions
!
! Eric Kemp, February 2000:  Modified code to allow for missing 
! pressure level or surface data.  Also, model and grid are
! now passed in subroutine call.
! 
! Eric Kemp, 31 March 2000:  Added lat/lon coordinates of corner points
! for NCL mapping.
!
!-----------------------------------------------------------------------
!
! Variable declarations
!
!-----------------------------------------------------------------------

  USE hdf_constants

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,nlevel,ntime,nvar_p,nvar_sfc
  INTEGER, INTENT(IN) :: &
      date(6), timesec(ntime), mapproj, pressure(nlevel)
  CHARACTER(LEN=*), INTENT(IN) :: model,grid
  REAL, INTENT(IN) :: scale, trulat1,trulat2,trulon, dx,dy, missing
  REAL, INTENT(IN) :: mapswlat,mapswlon,mapnwlat,mapnwlon, &
                      mapselat,mapselon,mapnelat,mapnelon
  REAL, INTENT(IN) :: scswlat,scswlon
  REAL, INTENT(IN) :: var_p(nx,ny,nlevel,ntime,nvar_p)
  REAL, INTENT(IN) :: var_sfc(nx,ny,ntime,nvar_sfc)
  CHARACTER(LEN=*), INTENT(IN) :: filename
  CHARACTER(LEN=*), DIMENSION(nvar_sfc), INTENT(IN) :: &
      varid_sfc, varname_sfc, varunit_sfc
  CHARACTER(LEN=*), DIMENSION(nvar_p), INTENT(IN) :: &
      varid_p, varname_p, varunit_p
  INTEGER, INTENT(IN) :: flag_p,flag_sfc

!-----------------------------------------------------------------------
!
! HDF parameters and variables
!
!-----------------------------------------------------------------------

  INTEGER,PARAMETER :: maxrank = 5
  INTEGER :: dims(maxrank),start(maxrank),edges(maxrank),stride(maxrank)

  INTEGER :: sd_id,sds_id

!-----------------------------------------------------------------------
!
! Miscellaneous variables
!
!-----------------------------------------------------------------------

  CHARACTER*256 :: string
  INTEGER :: i, j, status

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
! Open HDF file.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'-----------------------------------------------------------' 
  PRINT *, 'wrtverif: Creating HDF file ', TRIM(filename)
  sd_id = sfstart(filename, DFACC_CREATE)
  IF (sd_id == FAIL) CALL hdf_fail (1, sd_id, sd_id, 'sfstart')

!-----------------------------------------------------------------------
!
! Add scalar attributes describing the model data.
!
!-----------------------------------------------------------------------

  PRINT *, 'wrtverif: Adding scalar attributes'
  status = sfsnatt(sd_id, 'start_date', DFNT_INT32, 6, date)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfscatt(sd_id, 'model', DFNT_CHAR8,LEN_TRIM(model), model)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfscatt')
  status = sfscatt(sd_id, 'grid', DFNT_CHAR8,LEN_TRIM(grid), grid)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfscatt')
  status = sfsnatt(sd_id, 'dx', DFNT_FLOAT32, 1, dx)  
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'dy', DFNT_FLOAT32, 1, dy)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'nx', DFNT_INT32, 1, nx)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'ny', DFNT_INT32, 1, ny)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'nlevel', DFNT_INT32, 1, nlevel)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'ntime', DFNT_INT32, 1, ntime)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'nvar_p', DFNT_INT32, 1, nvar_p)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'nvar_sfc', DFNT_INT32, 1, nvar_sfc)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'mapproj', DFNT_INT32, 1, mapproj)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'scale', DFNT_FLOAT32, 1, scale)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'trulon', DFNT_FLOAT32, 1, trulon)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'trulat1', DFNT_FLOAT32, 1, trulat1)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'trulat2', DFNT_FLOAT32, 1, trulat2)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'scswlat', DFNT_FLOAT32, 1, scswlat)   
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'scswlon', DFNT_FLOAT32, 1, scswlon)   
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'mapswlat', DFNT_FLOAT32, 1, mapswlat)   
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'mapswlon', DFNT_FLOAT32, 1, mapswlon)   
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'mapnwlat', DFNT_FLOAT32, 1, mapnwlat)   
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'mapnwlon', DFNT_FLOAT32, 1, mapnwlon)   
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'mapselat', DFNT_FLOAT32, 1, mapselat)   
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'mapselon', DFNT_FLOAT32, 1, mapselon)   
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'mapnelat', DFNT_FLOAT32, 1, mapnelat)   
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'mapnelon', DFNT_FLOAT32, 1, mapnelon)   
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')

!-----------------------------------------------------------------------
!
! Write the time datasets.
!
!-----------------------------------------------------------------------

  PRINT *, 'wrtverif: Creating time', ntime
  CALL write_sds (sds_id, sd_id, DFNT_INT32, 'time', timesec, 1, ntime, missing)

  string = 'time'
  status = sfscatt(sds_id,'long_name', DFNT_CHAR8,LEN_TRIM(string),TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 's'
  status = sfscatt(sds_id,'units', DFNT_CHAR8, LEN_TRIM(string), TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  string = 'seconds_from_init_time'
  status =sfscatt(sds_id,'description',DFNT_CHAR8,LEN_TRIM(string),TRIM(string))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

  status = sfendacc(sds_id)
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfendacc')

!-----------------------------------------------------------------------
!
! Write the pressure datasets.
!
!-----------------------------------------------------------------------

  IF (flag_p == 1) THEN

    PRINT *, 'wrtverif: Creating pressure', nlevel

    CALL write_sds (sds_id, sd_id, DFNT_FLOAT32, 'pressure', &
                  pressure, 1, nlevel, missing)

    string = 'pressure'
    status = sfscatt(sds_id,'long_name',DFNT_CHAR8,LEN_TRIM(string),     &
                     TRIM(string))
    IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

    string = 'Pa'
    status = sfscatt(sds_id,'units', DFNT_CHAR8, LEN_TRIM(string),       &
                     TRIM(string))
    IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfscatt')

    status = sfendacc(sds_id)
    IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfendacc')

!-----------------------------------------------------------------------
!
!   Names, units
!
!-----------------------------------------------------------------------

    DO i=1,nvar_p
      print *, i, varid_p(i), varunit_p(i), varname_p(i)
    END DO
    dims(2) = nvar_p
    dims(1) = LEN(varid_p(1))
    CALL write_sds_char (sds_id, sd_id, 'id_p', varid_p, 2, dims)
    status = sfendacc(sds_id)
    IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfendacc')

    dims(2) = nvar_p
    dims(1) = LEN(varunit_p(1))
    CALL write_sds_char (sds_id, sd_id, 'units_p', varunit_p, 2, dims)
    status = sfendacc(sds_id)
    IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfendacc')

    dims(2) = nvar_p
    dims(1) = LEN(varname_p(1))
    CALL write_sds_char (sds_id, sd_id, 'names_p', varname_p, 2, dims)
    status = sfendacc(sds_id)
    IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfendacc')

    dims(1) = nx
    dims(2) = ny
    dims(3) = nlevel
    dims(4) = ntime
    dims(5) = nvar_p

    PRINT *, 'wrtverif: Creating var_p', dims(1:5)
    CALL write_sds (sds_id, sd_id, DFNT_FLOAT32, 'var_p', var_p, 5, dims,&
                    missing)

    status = sfendacc(sds_id)
    IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfendacc')
  ENDIF

!-----------------------------------------------------------------------
!
! Write the surface datasets.
!
!-----------------------------------------------------------------------

  IF (flag_sfc == 1) THEN 

    DO i=1,nvar_sfc
      print *, i, varid_sfc(i), varunit_sfc(i), varname_sfc(i)
    END DO

    dims(2) = nvar_sfc
    dims(1) = LEN(varid_sfc(1))
    CALL write_sds_char (sds_id, sd_id, 'id_sfc', varid_sfc, 2, dims)
    status = sfendacc(sds_id)
    IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfendacc')

    dims(2) = nvar_sfc
    dims(1) = LEN(varunit_sfc(1))
    CALL write_sds_char (sds_id, sd_id, 'units_sfc', varunit_sfc, 2, dims)
    status = sfendacc(sds_id)
    IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfendacc')

    dims(2) = nvar_sfc
    dims(1) = LEN(varname_sfc(1))
    CALL write_sds_char (sds_id, sd_id, 'names_sfc', varname_sfc, 2, dims)
    status = sfendacc(sds_id)
    IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfendacc')

    dims(1) = nx
    dims(2) = ny
    dims(3) = ntime
    dims(4) = nvar_sfc

    PRINT *, 'wrtverif: Creating var_sfc', dims(1:4)
    CALL write_sds (sds_id, sd_id, DFNT_FLOAT32, 'var_sfc', var_sfc, &
      4, dims, missing)

    status = sfendacc(sds_id)
    IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'sfendacc')
  ENDIF

!-----------------------------------------------------------------------
!
! Close file and exit.
!
!-----------------------------------------------------------------------

  status = sfend(sd_id)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfend')

  PRINT *, 'wrtverif: Done'
  WRITE(6,*)'-----------------------------------------------------------' 

END SUBROUTINE wrtverif


!########################################################################
!########################################################################
!######                                                            ######
!######                   SUBROUTINE RDVERIF                       ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################

SUBROUTINE rdverif ( nx,ny,nlevel,ntime,nvar_p,nvar_sfc,		&
                     filename,model,grid,date,timesec,pressure,	        &
                     mapproj, scale, trulat1,trulat2,trulon,            &
                     dx,dy,scswlat,scswlon,                             &
                     mapswlat,mapswlon,mapnwlat,mapnwlon,               &
                     mapselat,mapselon,mapnelat,mapnelon,               &
   		     flag_p,var_p, varid_p, varname_p, varunit_p,       &
		     flag_sfc,var_sfc, varid_sfc, varname_sfc,          &
                     varunit_sfc)

!-----------------------------------------------------------------------
!
! PURPOSE:
!
! Reads gridded verification file in HDF format.
!
! AUTHOR:  Eric Kemp, November 1999
!
! MODIFICATION HISTORY
!
! Richard Carpenter, January 2000: Extensive revisions
! 
! Eric Kemp, February 2000:  Modified code to allow for missing 
! pressure level or surface data.
!
! Eric Kemp, 31 Marcy 2000:  Added lat/lon coordinates of four corners
! for NCL plotting.
!
!-----------------------------------------------------------------------

  USE hdf_constants

  IMPLICIT NONE

!-----------------------------------------------------------------------
!
! Variable declarations
!
!-----------------------------------------------------------------------

  INTEGER, INTENT(IN) :: nx,ny,nlevel,ntime,nvar_p,nvar_sfc
  INTEGER, INTENT(OUT) :: &
      date(6), timesec(ntime), mapproj
  REAL, INTENT(OUT) :: scale, trulat1,trulat2,trulon, dx,dy, &
      pressure(nlevel)
  REAL, INTENT(IN) :: mapswlat,mapswlon,mapnwlat,mapnwlon, &
                      mapselat,mapselon,mapnelat,mapnelon
  REAL, INTENT(IN) :: scswlat,scswlon

  REAL, INTENT(OUT) :: var_p(nx,ny,nlevel,ntime,nvar_p)
  REAL, INTENT(OUT) :: var_sfc(nx,ny,ntime,nvar_sfc)
  CHARACTER*(*), INTENT(IN) :: filename
  CHARACTER*(*), DIMENSION(nvar_sfc), INTENT(OUT) :: &
      varid_sfc,varname_sfc,varunit_sfc
  CHARACTER*(*), DIMENSION(nvar_p), INTENT(OUT) :: varid_p,varname_p,varunit_p
  CHARACTER*(*), INTENT(OUT) :: model, grid
  INTEGER, INTENT(OUT) :: flag_p,flag_sfc

!-----------------------------------------------------------------------
!
! HDF parameters and variables
!
!-----------------------------------------------------------------------

  INTEGER,PARAMETER :: maxrank = 5
  INTEGER, DIMENSION(maxrank) :: dims, start, edges, stride
  INTEGER :: sd_id,sds_id

!-----------------------------------------------------------------------
!
! Metadata
!
!-----------------------------------------------------------------------

  INTEGER,PARAMETER :: max_metadata = 23
  CHARACTER*10,PARAMETER :: metadata(max_metadata) =                   &
     (/'start_date', 'model     ', 'grid      ', 'dx        ',         &
       'dy        ', 'nx        ', 'ny        ', 'ntime     ',         &
       'mapproj   ', 'scale     ', 'trulon    ', 'trulat1   ',         &
       'trulat2   ', 'scswlat   ', 'scswlon   ', 'mapswlat  ',         &
       'mapswlon  ', 'mapnwlat  ', 'mapnwlon  ', 'mapselat  ',         &
       'mapselon  ', 'mapnelat  ', 'mapnelon  '/)

!-----------------------------------------------------------------------
!
! Miscellaneous variables
!
!-----------------------------------------------------------------------

  CHARACTER*10 :: attr_name

  INTEGER :: attr_index,status,n_values,sd_index,data_type, i, j

  INTEGER :: rank0, dims0, data_type0, n_attrs0, n_datasets, n_file_attrs, &
      nx0,ny0,ntime0 ! Parameters read in from file.

  CHARACTER*32 :: name0*32

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!-----------------------------------------------------------------------
!
! Open HDF file.
!
!-----------------------------------------------------------------------

  WRITE(6,*)'-----------------------------------------------------------'
  PRINT *, 'rdverif: Opening ', TRIM(filename)
  sd_id = sfstart(filename, DFACC_READ)
  IF (sd_id == FAIL) CALL hdf_fail (1, sd_id, sd_id, 'sfstart')

!-----------------------------------------------------------------------
!
! Determine the number of data sets in the file and the number of
! file attributes.
!
!-----------------------------------------------------------------------

  status = sffinfo(sd_id, n_datasets, n_file_attrs)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sffinfo')
  PRINT *, 'rdverif: There are', n_datasets, 'datasets and', &
           n_file_attrs, 'file attributes'

!-----------------------------------------------------------------------
!
! Read file attributes describing the model data.
!
!-----------------------------------------------------------------------

  DO i = 1,max_metadata

    PRINT *, 'rdverif: Getting attribute ', i, TRIM(metadata(i))
    attr_index = sffattr(sd_id,metadata(i))
    status = sfgainfo(sd_id, attr_index, attr_name, data_type, n_values)
    IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfgainfo')

    SELECT CASE (TRIM(attr_name))
      CASE ("start_date")
        status = sfrattr(sd_id, attr_index, date)
      CASE ("model")  
        status = sfrattr(sd_id, attr_index, model)
      CASE ("grid")  
        status = sfrattr(sd_id, attr_index, grid)
      CASE ("dx")  
        status = sfrattr(sd_id, attr_index, dx)
      CASE ("dy")  
        status = sfrattr(sd_id, attr_index, dy)
      CASE ("nx")  
        status = sfrattr(sd_id, attr_index, nx0)
      CASE ("ny")  
        status = sfrattr(sd_id, attr_index, ny0)
      CASE ("ntime")  
        status = sfrattr(sd_id, attr_index, ntime0)
      CASE ("mapproj")  
        status = sfrattr(sd_id, attr_index, mapproj)
      CASE ("scale")  
        status = sfrattr(sd_id, attr_index, scale)
      CASE ("trulon")  
        status = sfrattr(sd_id, attr_index, trulon)
      CASE ("trulat1")  
        status = sfrattr(sd_id, attr_index, trulat1)
      CASE ("trulat2")  
        status = sfrattr(sd_id, attr_index, trulat2)
      CASE ("scswlat")  
        status = sfrattr(sd_id, attr_index, scswlat)
      CASE ("scswlon")  
        status = sfrattr(sd_id, attr_index, scswlon)
      CASE ("mapswlat")  
        status = sfrattr(sd_id, attr_index, mapswlat)
      CASE ("mapswlon")  
        status = sfrattr(sd_id, attr_index, mapswlon)
      CASE ("mapnwlat")  
        status = sfrattr(sd_id, attr_index, mapnwlat)
      CASE ("mapnwlon")  
        status = sfrattr(sd_id, attr_index, mapnwlon)
      CASE ("mapselat")  
        status = sfrattr(sd_id, attr_index, mapselat)
      CASE ("mapselon")  
        status = sfrattr(sd_id, attr_index, mapselon)
      CASE ("mapnelat")  
        status = sfrattr(sd_id, attr_index, mapnelat)
      CASE ("mapnelon")  
        status = sfrattr(sd_id, attr_index, mapnelon)

      CASE DEFAULT
        WRITE(6,*)'rdverif: FATAL: Unknown metadata.  Aborting...'
	PRINT *, 'metadata=', TRIM(attr_name)
        STOP
    END SELECT
  END DO

  IF (ntime0 /= ntime .OR. nx0 /= nx .OR.  ny0 /= ny ) THEN
    PRINT *, 'rdverif: FATAL:  Dimension mismatch.'
    PRINT *, 'File ntime = ',ntime0,' Program ntime = ',ntime
    PRINT *, 'File nx = ',nx0,' Program nx = ',nx
    PRINT *, 'File ny = ',ny0,' Program ny = ',ny
    STOP
  ENDIF 

!-----------------------------------------------------------------------
!
! Read in time data.
!
!-----------------------------------------------------------------------

  sd_index = 0
  PRINT *, 'rdverif: Reading ', sd_index, ' time'
  CALL read_sds (sd_id, sd_index, timesec, 1, ntime)
  PRINT *, 'timesec: ', timesec

  sd_index = sd_index + 1

!-----------------------------------------------------------------------
!
! Read in pressure level data.
!
!-----------------------------------------------------------------------

  sds_id = sfselect(sd_id,sd_index) ! EMK
  status = sfginfo(sds_id, name0, rank0, dims0, data_type0, n_attrs0)
  IF (name0 == 'pressure') THEN
    flag_p = 1
    PRINT *, 'rdverif: Reading ', sd_index, ' pressure'
    CALL read_sds (sd_id, sd_index, pressure, 1, nlevel)
    PRINT *, 'pressure: ', pressure
    sd_index = sd_index + 1

    PRINT *, 'rdverif: Reading ', sd_index, ' varid_p'

    dims(1) = LEN(varid_p(1))
    dims(2) = nvar_p
    CALL read_sds_char (sd_id, sd_index, varid_p, 2, dims)
    DO i=1,nvar_p
      PRINT *, i, TRIM(varid_p(i))
    END DO

    sd_index = sd_index + 1
    dims(1) = LEN(varunit_p(1))
    dims(2) = nvar_p
    CALL read_sds_char (sd_id, sd_index, varunit_p, 2, dims)
    DO i=1,nvar_p
      PRINT *, i, TRIM(varunit_p(i))
    END DO

    sd_index = sd_index + 1
    dims(1) = LEN(varname_p(1))
    dims(2) = nvar_p
    CALL read_sds_char (sd_id, sd_index, varname_p, 2, dims)
    DO i=1,nvar_p
      PRINT *, i, TRIM(varname_p(i))
    END DO

    sd_index = sd_index + 1
    PRINT *, 'rdverif: Reading ', sd_index, ' var_p'
    dims(1) = nx
    dims(2) = ny
    dims(3) = nlevel
    dims(4) = ntime
    dims(5) = nvar_p
    CALL read_sds (sd_id, sd_index, var_p, 5, dims)

    sd_index = sd_index + 1
  ELSE
    flag_p = 0
  ENDIF

!-----------------------------------------------------------------------
!
! Get surface data.
!
!-----------------------------------------------------------------------

  sds_id = sfselect(sd_id,sd_index) ! EMK
  status = sfginfo(sds_id, name0, rank0, dims0, data_type0, n_attrs0)
  IF (name0 == 'id_sfc') THEN
    flag_sfc = 1

    PRINT *, 'rdverif: Reading ', sd_index, ' varid_sfc'
    dims(1) = LEN(varid_sfc(1))
    dims(2) = nvar_sfc
    CALL read_sds_char (sd_id, sd_index, varid_sfc, 2, dims)
    DO i=1,nvar_sfc
      PRINT *, i, TRIM(varid_sfc(i))
    END DO
    sd_index = sd_index + 1

    dims(1) = LEN(varunit_sfc(1))
    dims(2) = nvar_sfc
    CALL read_sds_char (sd_id, sd_index, varunit_sfc, 2, dims)
    DO i=1,nvar_sfc
      PRINT *, i, TRIM(varunit_sfc(i))
    END DO
    
    sd_index = sd_index + 1
    dims(1) = LEN(varname_sfc(1))
    dims(2) = nvar_sfc
    CALL read_sds_char (sd_id, sd_index, varname_sfc, 2, dims)
    DO i=1,nvar_sfc
      PRINT *, i, TRIM(varname_sfc(i))
    END DO
    
    sd_index = sd_index + 1

    dims(1) = nx
    dims(2) = ny
    dims(3) = ntime
    dims(4) = nvar_sfc
    CALL read_sds (sd_id, sd_index, var_sfc, 4, dims)
  ELSE
    flag_sfc = 0
  ENDIF

!-----------------------------------------------------------------------
!
! Close file and exit.
!
!-----------------------------------------------------------------------

  status = sfend(sd_id)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfend')

  n_datasets = n_datasets - 1

  WRITE(6,*)'rdverif: Done'
  WRITE(6,*)'-----------------------------------------------------------' 

  RETURN
END SUBROUTINE rdverif


!########################################################################
!########################################################################
!######                                                            ######
!######                SUBROUTINE WRT_VERIF_STATS_DIFF             ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################

SUBROUTINE wrt_verif_stats_diff (filename, if_diff, diffile, date, &
                     amodel,agrid,adx,anx,any, fmodel,fgrid,fdx, &
                     vdx, vdy, vnx, vny, vmapproj, &
		     vtrulat1, vtrulat2, vtrulon, vsclfct, vctrlat, vctrlon, &
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
! Manages the writing of verification statistics and difference fields to
!  HDF files.
!
! AUTHOR:  Richard Carpenter, January 2000
!
!-----------------------------------------------------------------------  
!                       
! Variable declarations
!                       
!-----------------------------------------------------------------------  

  USE hdf_constants

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: date(6), vnx,vny, anx,any, vmapproj, &
      mastercounter, if_diff
  INTEGER, INTENT(IN) :: nlevel, ntime, nvar_p, nvar_sfc
  INTEGER, INTENT(IN) :: timesec(ntime), pressure(nlevel)
  REAL, INTENT(IN) :: vdx,vdy, adx, fdx, corner_lat(2,2), corner_lon(2,2), &
		     vtrulat1, vtrulat2, vtrulon, vsclfct, vctrlat, vctrlon, &
		     missing
  REAL, INTENT(IN), DIMENSION(nlevel,ntime,nvar_p) :: &
      counter_p, bias_p, rms_p
  REAL, INTENT(IN) :: diff_p(anx,any,nlevel,ntime,nvar_p), &
                    diff_sfc(anx,any,ntime,nvar_p) 
  REAL, INTENT(IN), DIMENSION(ntime,nvar_sfc) :: &
      counter_sfc, bias_sfc, rms_sfc
  CHARACTER(LEN=*), INTENT(IN) :: filename, diffile, &
      amodel, agrid, fmodel, fgrid
  CHARACTER(LEN=*), INTENT(IN), DIMENSION(nvar_p) :: &
      varid_p, varname_p, varunit_p
  CHARACTER(LEN=*), INTENT(IN), DIMENSION(nvar_sfc) :: &
      varid_sfc, varname_sfc, varunit_sfc
  
!-----------------------------------------------------------------------  
!                       
! HDF variables and parameters.
!                       
!-----------------------------------------------------------------------  

  INTEGER,PARAMETER :: maxrank = 5
  INTEGER, DIMENSION(maxrank) :: dims, start, edges, stride

  INTEGER :: sd_id,sds_id
  
!-----------------------------------------------------------------------
!
! Miscellaneous variables
!
!-----------------------------------------------------------------------

  CHARACTER :: string*256

  INTEGER :: status, i

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!                             
! Write the stats
!
!-----------------------------------------------------------------------

  WRITE(6,*)'-----------------------------------------------------------' 
  PRINT *, 'wrt_verif_stats_diff: Stats to file ', TRIM(filename)
print *, 'stats: nvar p, sfc: ', nvar_p, nvar_sfc
  CALL wrt_verif_header ( sd_id, missing, filename, date, &
                     amodel,agrid,adx, fmodel,fgrid,fdx, &
                     vdx, vdy, vnx, vny, vmapproj, &
		     vtrulat1, vtrulat2, vtrulon, vsclfct, vctrlat, vctrlon, &
                     corner_lat, corner_lon,  &
                     nlevel, ntime, nvar_p, nvar_sfc, &
		     timesec, pressure, mastercounter, &
		     counter_p, counter_sfc, &
	             varid_p, varname_p, varunit_p, &
	             varid_sfc, varname_sfc, varunit_sfc)

  PRINT *, 'wrt_verif_stats_diff: Done call wrt_verif_header'
  dims(1) = nlevel
  dims(2) = ntime
  dims(3) = nvar_p
  CALL write_sds (sds_id, sd_id, DFNT_INT32, 'counter_p', counter_p, &
      3, dims, missing)
  CALL write_sds (sds_id, sd_id, DFNT_FLOAT32, 'bias_p', bias_p, &
      3, dims, missing)
  CALL write_sds (sds_id, sd_id, DFNT_FLOAT32, 'rms_p', rms_p, &
      3, dims, missing)

  dims(1) = ntime
  dims(2) = nvar_sfc
  CALL write_sds (sds_id, sd_id, DFNT_INT32, 'counter_sfc', counter_sfc, &
      2, dims, missing)
  CALL write_sds (sds_id, sd_id, DFNT_FLOAT32, 'bias_sfc', bias_sfc, &
      2, dims, missing)
  CALL write_sds (sds_id, sd_id, DFNT_FLOAT32, 'rms_sfc', rms_sfc, &
      2, dims, missing)

  status = sfend(sd_id)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfend')

!-----------------------------------------------------------------------
!                             
! Write the difference fields
!
!-----------------------------------------------------------------------

  IF (if_diff == 0) THEN
    PRINT *, 'wrt_verif_stats_diff: Successful completion'
    RETURN
  END IF

  ! Note that diff_p and diff_sfc have dimensions anx,any

  PRINT *, 'wrt_verif_stats_diff: Differences to file ', TRIM(diffile)
  CALL wrt_verif_header ( sd_id, missing, diffile, date, &
                     amodel,agrid,adx, fmodel,fgrid,fdx, &
                     adx, adx, anx, any, vmapproj, &
		     vtrulat1, vtrulat2, vtrulon, vsclfct, vctrlat, vctrlon, &
                     corner_lat, corner_lon,  &
                     nlevel, ntime, nvar_p, nvar_sfc, &
		     timesec, pressure, mastercounter, &
		     counter_p, counter_sfc, &
	             varid_p, varname_p, varunit_p, &
	             varid_sfc, varname_sfc, varunit_sfc)
  PRINT *, 'wrt_verif_stats_diff: Done call wrt_verif_header'

  dims(1) = anx
  dims(2) = any
  dims(3) = nlevel
  dims(4) = ntime
  dims(5) = nvar_p
  CALL write_sds (sds_id, sd_id, DFNT_FLOAT32, 'diff_p', diff_p, &
      5, dims, missing)

  dims(1) = anx
  dims(2) = any
  dims(3) = ntime
  dims(4) = nvar_sfc
  CALL write_sds (sds_id, sd_id, DFNT_FLOAT32, 'diff_sfc', diff_sfc, &
      4, dims, missing)
 
  status = sfend(sd_id)
    IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfend')

  PRINT *, 'wrt_verif_stats_diff: Successful completion'
  WRITE(6,*)'-----------------------------------------------------------' 

END SUBROUTINE wrt_verif_stats_diff


!########################################################################
!########################################################################
!######                                                            ######
!######                SUBROUTINE WRT_VERIF_HEADER                 ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################

SUBROUTINE wrt_verif_header (sd_id, missing, filename, date, &
                     amodel,agrid,adx, fmodel,fgrid,fdx, &
                     vdx, vdy, vnx, vny, vmapproj, &
		     vtrulat1, vtrulat2, vtrulon, vsclfct, vctrlat, vctrlon, &
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
! Writes header for verification statistics file in HDF format.
!
! AUTHOR:  Richard Carpenter, January 2000
!
!-----------------------------------------------------------------------  
!                       
! Variable declarations
!                       
!-----------------------------------------------------------------------  

  USE hdf_constants

  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: sd_id
  INTEGER, INTENT(IN) :: date(6), vnx,vny, vmapproj, mastercounter
  INTEGER, INTENT(IN) :: nlevel, ntime, nvar_p, nvar_sfc
  INTEGER, INTENT(IN) :: timesec(ntime), pressure(nlevel)
  REAL, INTENT(IN) :: vdx,vdy, adx, fdx, corner_lat(2,2), corner_lon(2,2), &
		     vtrulat1, vtrulat2, vtrulon, vsclfct, vctrlat, vctrlon, &
		     missing
  REAL, INTENT(IN), DIMENSION(nlevel,ntime,nvar_p) :: counter_p
  REAL, INTENT(IN), DIMENSION(ntime,nvar_sfc) :: counter_sfc
  CHARACTER(LEN=*), INTENT(IN) :: filename, amodel, agrid, fmodel, fgrid
  CHARACTER(LEN=*), INTENT(IN), DIMENSION(nvar_p) :: &
      varid_p, varname_p, varunit_p
  CHARACTER(LEN=*), INTENT(IN), DIMENSION(nvar_sfc) :: &
      varid_sfc, varname_sfc, varunit_sfc
  
!-----------------------------------------------------------------------
!
! Miscellaneous variables
!
!-----------------------------------------------------------------------

  INTEGER :: sds_id, status, i, dims(5)
  CHARACTER :: string*256

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
  sd_id = sfstart(filename, DFACC_CREATE)
  IF (status == FAIL) CALL hdf_fail (1, sd_id, sd_id, 'sfstart')

!-----------------------------------------------------------------------
!                             
! Add file attributes describing the verification data.
!
!-----------------------------------------------------------------------
 
  status = sfscatt(sd_id, 'forecast_source', DFNT_CHAR8,LEN_TRIM(fmodel),fmodel)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfscatt')
  status = sfscatt(sd_id, 'forecast_grid', DFNT_CHAR8,LEN_TRIM(fgrid),fgrid)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfscatt')
  status = sfscatt(sd_id, 'analysis_source', DFNT_CHAR8,LEN_TRIM(amodel),amodel)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfscatt')
  status = sfscatt(sd_id, 'analysis_grid', DFNT_CHAR8,LEN_TRIM(agrid),agrid)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfscatt')
  status = sfsnatt(sd_id, 'start_date', DFNT_INT32, 6, date)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'vdx', DFNT_FLOAT32, 1, vdx)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'vdy', DFNT_FLOAT32, 1, vdy)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'vnx', DFNT_INT32, 1, vnx)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'vny', DFNT_INT32, 1, vny)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'vmapproj', DFNT_INT32, 1, vmapproj)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'vtrulat1', DFNT_FLOAT32, 1, vtrulat1)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'vtrulat2', DFNT_FLOAT32, 1, vtrulat2)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'vtrulon', DFNT_FLOAT32, 1, vtrulon)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'vctrlat', DFNT_FLOAT32, 1, vctrlat)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'vctrlon', DFNT_FLOAT32, 1, vctrlon)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'adx', DFNT_FLOAT32, 1, adx)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  status = sfsnatt(sd_id, 'fdx', DFNT_FLOAT32, 1, fdx)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfsnatt')
  
  status = sfsnatt(sd_id, 'mastercounter', DFNT_INT32, 1, mastercounter)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfnatt')
  
!-----------------------------------------------------------------------
!
! Corner lat/lon, time
!
!-----------------------------------------------------------------------
 
  dims(:) = 2
  CALL write_sds (sds_id, sd_id, DFNT_FLOAT32, 'corner_lat', corner_lat, &
      2, dims, missing)
  CALL write_sds (sds_id, sd_id, DFNT_FLOAT32, 'corner_lon', corner_lon, &
      2, dims, missing)
 
  CALL write_sds (sds_id, sd_id, DFNT_FLOAT32, 'pressure', pressure, &
      1, nlevel, missing)
 
  CALL write_sds (sds_id, sd_id, DFNT_INT32, 'time', timesec, &
      1, ntime, missing)
 
  string = 'time'
  status = sfscatt(sds_id, 'code', DFNT_CHAR8, LEN_TRIM(string), string)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfscatt')

  string = 's'
  status = sfscatt(sds_id,'units', DFNT_CHAR8, LEN_TRIM(string), string)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfscatt')

  string = 'seconds_from_init_time'
  status = sfscatt(sds_id,'description', DFNT_CHAR8, LEN_TRIM(string), string)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfscatt')

  status = sfendacc(sds_id)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'sfendacc')

!-----------------------------------------------------------------------
!
! Describe the data
!
!-----------------------------------------------------------------------

print *, 'hdr: nvar p, sfc: ', nvar_p, nvar_sfc
 
  DO i=1,nvar_p
    print *, i, varid_p(i), " ", varunit_p(i), " ", varname_p(i)
  END DO
  dims(2) = nvar_p
  dims(1) = LEN(varid_p(1))
  CALL write_sds_char (sds_id, sd_id, 'id_p', varid_p, 2, dims)
  dims(1) = LEN(varunit_p(1))
  CALL write_sds_char (sds_id, sd_id, 'units_p', varunit_p, 2, dims)
  dims(1) = LEN(varname_p(1))
  CALL write_sds_char (sds_id, sd_id, 'names_p', varname_p, 2, dims)

print *, "varid_sfc: ", varid_sfc
print *, "varunit_sfc: ", varunit_sfc
print *, "varname_sfc: ", varname_sfc

  DO i=1,nvar_sfc
    print *, i, varid_sfc(i), " ", varunit_sfc(i), " ", varname_sfc(i)
  END DO
  dims(2) = nvar_sfc
  dims(1) = LEN(varid_sfc(1))
  CALL write_sds_char (sds_id, sd_id, 'id_sfc', varid_sfc, 2, dims)
  dims(1) = LEN(varunit_sfc(1))
  CALL write_sds_char (sds_id, sd_id, 'units_sfc', varunit_sfc, 2, dims)
  dims(1) = LEN(varname_sfc(1))
  CALL write_sds_char (sds_id, sd_id, 'names_sfc', varname_sfc, 2, dims)

  PRINT *, 'wrt_verif_header: Done'
  WRITE(6,*)'-----------------------------------------------------------' 

END SUBROUTINE wrt_verif_header


!########################################################################
!########################################################################
!######                                                            ######
!######                  SUBROUTINE ALLOC_FAIL                     ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################

SUBROUTINE alloc_fail (status, string)

!-----------------------------------------------------------------------
! 
! PURPOSE:
! 
! Gracefully exits program if there is an array allocation error.
!
! AUTHOR:  Richard Carpenter, January 2000
!
! MODIFICATION HISTORY:
!
! Eric Kemp, March 2000
! Added Documentation.
!
!----------------------------------------------------------------------- 
!
! Variable declarations
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: status
  CHARACTER(LEN=*), INTENT(IN) :: string

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  PRINT *, 'alloc_fail: FATAL: Allocation failure. status=', status
  PRINT *, 'Message: ', TRIM(string)
  STOP
END SUBROUTINE alloc_fail


!########################################################################
!########################################################################
!######                                                            ######
!######                   SUBROUTINE HDF_FAIL                      ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################

SUBROUTINE hdf_fail (disposition, status, id, string)

!-----------------------------------------------------------------------
! 
! PURPOSE:
! 
! Gracefully exits program if there is an HDF library call error.
!
! AUTHOR:  Richard Carpenter, January 2000
!
! MODIFICATION HISTORY:
!
! Eric Kemp, March 2000
! Added Documentation.
!
!----------------------------------------------------------------------- 
!
! Variable declarations
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: disposition, status, id
  CHARACTER(LEN=*), INTENT(IN) :: string
  CHARACTER(LEN=5) :: disposition_string(0:1) = (/ 'ERROR', 'FATAL' /)

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  PRINT *, 'hdf_fail: ', TRIM(disposition_string(disposition)), &
           ': HDF error. status=', status, ', id=', id
  PRINT *, 'hdf_fail: Message: ', TRIM(string)
  IF (disposition > 0) STOP

END SUBROUTINE hdf_fail


!########################################################################
!########################################################################
!######                                                            ######
!######                   SUBROUTINE READ_SDS                      ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################

SUBROUTINE read_sds (sd_id, sd_index, var, rank, dims)

!-----------------------------------------------------------------------
! 
! PURPOSE:
! 
! Reads numerical SDS data from a HDF file. 
!
! AUTHOR:  Richard Carpenter, January 2000
!
! MODIFICATION HISTORY:
!
! Eric Kemp, March 2000
! Added Documentation.
!
!----------------------------------------------------------------------- 
!
! Variable declarations
!
!-----------------------------------------------------------------------

  USE hdf_constants

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: sd_id, sd_index, rank
  INTEGER, INTENT(IN), DIMENSION(rank) :: dims
  REAL, INTENT(IN) :: var

!----------------------------------------------------------------------- 
!
! Local variables
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: maxrank=8

  INTEGER :: status, sds_id, rank0, data_type0, n_attrs0, attr_index, i
  INTEGER, DIMENSION(maxrank) :: start, edges, stride, dims0
  CHARACTER :: name0*32

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!----------------------------------------------------------------------- 
!
! Find SDS in HDF file.
!
!-----------------------------------------------------------------------

  PRINT *, 'read_sds: Reading sd_index', sd_index
  sds_id = sfselect (sd_id, sd_index) 
  IF (sds_id == FAIL) CALL hdf_fail (1, sds_id, sds_id, 'rd1: sfselect')

!----------------------------------------------------------------------- 
!
! Get information on the SDS.
!
!-----------------------------------------------------------------------

  status = sfginfo(sds_id, name0, rank0, dims0, data_type0, n_attrs0)
  IF (attr_index == FAIL) CALL hdf_fail (1, attr_index,sds_id, 'rd1: sfginfo')

  PRINT *, "read_sds: name = ", TRIM(name0)
  PRINT *, "     rank =", rank0, " dims =", dims0(1:rank0)
  PRINT *, "     data_type =", data_type0, " num attrs =", n_attrs0

!----------------------------------------------------------------------- 
!
! Check rank and dimensions
!
!-----------------------------------------------------------------------

  IF (rank /= rank0) THEN
    PRINT *, 'rd1: FATAL: Arrays differ in rank: ', rank, rank0
    STOP
  END IF

  DO i=1,rank
    IF (dims(i) /= dims0(i)) THEN
      PRINT *, 'rd1: FATAL: Dimension ', i, ' differs: ', dims(i), dims0(i)
      STOP
    END IF
  END DO

!----------------------------------------------------------------------- 
!
! Read in SDS data.
!
!-----------------------------------------------------------------------

  start(:) = 0
  edges(:) = dims0(:)
  stride(:) = 1

  status = sfrdata (sds_id, start, stride, edges, var)
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'rd1: sfrdata')

!----------------------------------------------------------------------- 
!
! Close SDS.
!
!-----------------------------------------------------------------------

  status = sfendacc(sds_id)
  IF (status == FAIL) CALL hdf_fail (1, status, sd_id, 'rd1: sfendacc')

  PRINT *, 'read_sds: Read in ', TRIM(name0), dims(1:rank)

END SUBROUTINE read_sds 


!########################################################################
!########################################################################
!######                                                            ######
!######                 SUBROUTINE READ_SDS_CHAR                   ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################

SUBROUTINE read_sds_char (sd_id, sd_index, var, rank, dims)

!-----------------------------------------------------------------------
! 
! PURPOSE:
! 
! Reads character SDS data from a HDF file. 
!
! AUTHOR:  Richard Carpenter, January 2000
!
! MODIFICATION HISTORY:
!
! Eric Kemp, March 2000
! Added Documentation.
!
!----------------------------------------------------------------------- 
!
! Variable declarations
!
!-----------------------------------------------------------------------

  USE hdf_constants

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: sd_id, sd_index, rank
  INTEGER, INTENT(IN), DIMENSION(rank) :: dims
  CHARACTER(LEN=*), INTENT(IN) :: var

!----------------------------------------------------------------------- 
!
! Local variables
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: maxrank=8

  INTEGER :: status, sds_id, rank0, data_type0, n_attrs0, attr_index, i
  INTEGER, DIMENSION(maxrank) :: start, edges, stride, dims0
  CHARACTER :: name0*32

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!----------------------------------------------------------------------- 
!
! Find SDS in HDF file.
!
!-----------------------------------------------------------------------

  PRINT *, 'read_sds_char: Reading sd_index', sd_index
  sds_id = sfselect (sd_id, sd_index) 
  IF (sds_id == FAIL) CALL hdf_fail (0, sds_id, sds_id, 'rd1c: sfselect')

!----------------------------------------------------------------------- 
!
! Get information on the SDS.
!
!-----------------------------------------------------------------------

  status = sfginfo(sds_id, name0, rank0, dims0, data_type0, n_attrs0)
  IF (attr_index == FAIL) CALL hdf_fail (0, attr_index,sds_id, 'rd1c:sfginfo')

  PRINT *, "read_sds_char: name = ", TRIM(name0)
  PRINT *, "     rank =", rank0, " dims =", dims0(1:rank0)
  PRINT *, "     data_type =", data_type0, " num attrs =", n_attrs0

!----------------------------------------------------------------------- 
!
! Check rank and dimensions
!
!-----------------------------------------------------------------------

  IF (rank /= rank0) THEN
    PRINT *, 'rd1c: ERROR: Arrays differ in rank: ', rank, rank0
    RETURN
  END IF

  DO i=1,rank
    IF (dims(i) /= dims0(i)) THEN
      PRINT *, 'rd1c: ERROR: Dimension ', i, ' differs: ', dims(i), dims0(i)
      RETURN
    END IF
  END DO

!----------------------------------------------------------------------- 
!
! Read in SDS data.
!
!-----------------------------------------------------------------------

  start(:) = 0
  edges(:) = dims0(:)
  stride(:) = 1

  status = sfrcdata (sds_id, start, stride, edges, var)
  IF (status == FAIL) CALL hdf_fail (0, status, sds_id, 'rd1c: sfrdata')

!----------------------------------------------------------------------- 
!
! Close SDS.
!
!-----------------------------------------------------------------------

  status = sfendacc(sds_id)
  IF (status == FAIL) CALL hdf_fail (0, status, sd_id, 'rd1c: sfendacc')

  PRINT *, 'read_sds_char: Read in ', TRIM(name0), dims(1:rank)

END SUBROUTINE read_sds_char


!########################################################################
!########################################################################
!######                                                            ######
!######                  SUBROUTINE WRITE_SDS                      ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################

SUBROUTINE write_sds (sds_id, sd_id, type, name, var, rank, dims, missing)

!-----------------------------------------------------------------------
! 
! PURPOSE:
! 
! Writes numerical SDS data to a HDF file. 
!
! AUTHOR:  Richard Carpenter, January 2000
!
! MODIFICATION HISTORY:
!
! Eric Kemp, March 2000
! Added Documentation.
!
! Eric Kemp, 31 March 2000
! Corrected checks for file access errors.
!
!----------------------------------------------------------------------- 
!
! Variable declarations
!
!-----------------------------------------------------------------------

  USE hdf_constants

  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: sds_id
  INTEGER, INTENT(IN) :: sd_id, type, rank
  INTEGER, INTENT(IN), DIMENSION(rank) :: dims
  REAL, INTENT(IN) :: var, missing
  CHARACTER(LEN=*) :: name

!----------------------------------------------------------------------- 
!
! Local variables
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: maxrank=8
  INTEGER :: status
  INTEGER, DIMENSION(maxrank) :: start, edges, stride

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  PRINT *, 'write_sds: ', TRIM(name), dims(1:rank)

!----------------------------------------------------------------------- 
!
! Create SDS
!
!-----------------------------------------------------------------------

  start(:) = 0
  edges(1:rank) = dims(1:rank)
  stride(:) = 1

  sds_id = sfcreate(sd_id, TRIM(name), type, rank, dims)
  IF (sds_id == FAIL) CALL hdf_fail (1, sds_id, sd_id, 'wrt1sds: sfcreate')

!----------------------------------------------------------------------- 
!
! Set compression
!
!-----------------------------------------------------------------------

  status = sfscompress(sds_id, comp_type, comp_prm(1))
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'wrt1sds: sfcompr')
  
!----------------------------------------------------------------------- 
!
! Write the data to the SDS.
!
!-----------------------------------------------------------------------

  status = sfwdata(sds_id, start, stride, edges, var)
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'wrt1sds: sfwdata')

!----------------------------------------------------------------------- 
!
! Set the missing value and exit.
!
!-----------------------------------------------------------------------

  status = sfsnatt(sds_id, '_FillValue', DFNT_FLOAT32, 1, missing)
  IF (status == FAIL) CALL hdf_fail (0, status, sds_id, 'sfsnatt')

END SUBROUTINE write_sds


!########################################################################
!########################################################################
!######                                                            ######
!######                SUBROUTINE WRITE_SDS_CHAR                   ######
!######                                                            ######
!######                      Developed by                          ######
!######         Center for Analysis and Prediction of Storms       ######
!######                   University of Oklahoma                   ######
!######                                                            ######
!########################################################################
!########################################################################

SUBROUTINE write_sds_char (sds_id, sd_id, name, var, rank, dims)

!-----------------------------------------------------------------------
! 
! PURPOSE:
! 
! Writes character SDS data to a HDF file. 
!
! AUTHOR:  Richard Carpenter, January 2000
!
! MODIFICATION HISTORY:
!
! Eric Kemp, March 2000
! Added Documentation.
!
! Eric Kemp, 31 March 2000
! Corrected checks for file access errors.
!
!-----------------------------------------------------------------------
!
! Variable declarations.  Note that the first dimension is the length of
! the strings, and the rank is increased by 1.
!
!-----------------------------------------------------------------------

  USE hdf_constants

  IMPLICIT NONE

  INTEGER, INTENT(OUT) :: sds_id
  INTEGER, INTENT(IN) :: sd_id, rank
  INTEGER, INTENT(IN), DIMENSION(rank) :: dims
  CHARACTER(LEN=*) :: var, name

!-----------------------------------------------------------------------
!
! Local variables
!
!-----------------------------------------------------------------------

  INTEGER, PARAMETER :: maxrank=8
  INTEGER :: status
  INTEGER, DIMENSION(maxrank) :: start, edges, stride

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  PRINT *, 'write_sds_char: ', TRIM(name), dims(1:rank)

!----------------------------------------------------------------------- 
!
! Create SDS
!
!-----------------------------------------------------------------------

  start(:) = 0
  edges(1:rank) = dims(1:rank)
  stride(:) = 1

  sds_id = sfcreate(sd_id, TRIM(name), DFNT_CHAR8, rank, dims)
  IF (sds_id == FAIL) CALL hdf_fail (1, sds_id, sd_id, 'wrt1c: sfcreate')

!----------------------------------------------------------------------- 
!
! Set compression.
!
!-----------------------------------------------------------------------

  !status = sfscompress(sds_id, comp_type, comp_prm(1))
  !  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'wrt1c: sfcompr')
  
!----------------------------------------------------------------------- 
!
! Write the data to the SDS.
!
!-----------------------------------------------------------------------

  status = sfwcdata(sds_id, start, stride, edges, var)
  IF (status == FAIL) CALL hdf_fail (1, status, sds_id, 'wrt1c: sfwcdata')

END SUBROUTINE write_sds_char
