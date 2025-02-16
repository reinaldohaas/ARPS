!
!##################################################################
!##################################################################
!######                                                      ######
!######             SUBROUTINE GET_DIMS_FROM_RAD             ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
SUBROUTINE get_grid_from_rad(radfname,dx,dy,dz,dzmin,strhopt,          &
                  mapproj,ctrlat,ctrlon,trulat1,trulat2,trulon,sclfct, &
                  istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!
!-----------------------------------------------------------------------
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    radfname   file name for radar datasets
!
!  OUTPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!
!    istatus  status indicator
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  CHARACTER (LEN=*) :: radfname
!
  REAL    :: dx,dy,dz,dzmin
  INTEGER :: strhopt,mapproj
  REAL    :: ctrlat,ctrlon
  REAL    :: trulat1,trulat2,trulon
  REAL    :: sclfct
!
  INTEGER :: istatus
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=4) :: stn
  CHARACTER (LEN=80) :: runname
  INTEGER :: ierr,lens,dmpfmt,nchanl
  INTEGER :: ireftim,itime,vcpnum,idummy
  INTEGER :: iradfmt,irngmin,irngmax
  INTEGER :: iyr, imon, idy, ihr, imin, isec
!
!-----------------------------------------------------------------------
!
! hdf variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: sd_id,isource
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  istatus=0
!
!-----------------------------------------------------------------------
!
! Open file
!
!-----------------------------------------------------------------------
!
  CALL asnctl ('NEWLOCAL', 1, ierr)
  CALL asnfile(radfname, '-F f77 -N ieee', ierr)

  lens=LEN_trim(radfname)
  !IF(radfname(lens-3:lens)=='hdf4')THEN
  IF(INDEX(radfname(1:lens),'hdf',.TRUE.) > 0)THEN
    dmpfmt=2
  ELSE
    dmpfmt=1
  ENDIF
  print*,' radfname=',TRIM(radfname),', dmpfmt=',dmpfmt

  IF(dmpfmt==1)THEN
    CALL getunit( nchanl )
    OPEN(UNIT=nchanl,FILE=trim(radfname),ERR=399,                       &
         FORM='unformatted',STATUS='old')
!
!-----------------------------------------------------------------------
!
!  Read radar description variables
!
!-----------------------------------------------------------------------
!
    istatus=1
!
    READ(nchanl) stn
    READ(nchanl) ireftim,itime,vcpnum,idummy,idummy,                    &
                 idummy,idummy,idummy,idummy,idummy
!
    CALL abss2ctim(itime, iyr, imon, idy, ihr, imin, isec )
    WRITE(6,'(/a,i2.2,a,i2.2,a,i4.4,1X,i2.2,a,i2.2,a)')                 &
       'Getting grid info from data recorded at: ',                     &
       imon,'/',idy,'/',iyr,ihr,':',imin,' UTC'
!
    READ(nchanl) runname
    READ(nchanl) iradfmt,strhopt,mapproj,irngmin,irngmax,               &
               idummy,idummy,idummy,idummy,idummy

    READ(nchanl) dx,dy,dz,dzmin,ctrlat,                                 &
           ctrlon,trulat1,trulat2,trulon,sclfct
!
    CLOSE(nchanl)
    CALL retunit( nchanl )

    GO TO 400

    399     CONTINUE
    PRINT*,'Error reading the radar file:',radfname

    400     CONTINUE
  ELSE !HDF4 file 

    CALL hdfopen(trim(radfname), 1, sd_id)
    IF (sd_id < 0) THEN
      WRITE (6,*) "get_grid_from_rad: ERROR opening ",                  &
                trim(radfname)," for reading."
      istatus = 1
      STOP
    END IF
    CALL hdfrdc(sd_id, 4, 'radid', stn, istatus)
    CALL hdfrdi(sd_id, 'ireftim', ireftim, istatus)
    CALL hdfrdi(sd_id, 'itime', itime, istatus)
    CALL hdfrdc(sd_id, 40, 'runname', runname, istatus)
    CALL hdfrdi(sd_id, 'iradfmt', iradfmt, istatus)
    CALL hdfrdi(sd_id, 'strhopt', strhopt, istatus)
    CALL hdfrdi(sd_id, 'mapproj', mapproj, istatus)
    CALL hdfrdr(sd_id, 'dx', dx, istatus)
    CALL hdfrdr(sd_id, 'dy', dy, istatus)
    CALL hdfrdr(sd_id, 'dz', dz, istatus)
    CALL hdfrdr(sd_id, 'dzmin', dzmin, istatus)
    CALL hdfrdr(sd_id, 'ctrlat', ctrlat, istatus)
    CALL hdfrdr(sd_id, 'ctrlon', ctrlon, istatus)
    CALL hdfrdr(sd_id, 'trulat1', trulat1, istatus)
    CALL hdfrdr(sd_id, 'trulat2', trulat2, istatus)
    CALL hdfrdr(sd_id, 'trulon', trulon, istatus)
    CALL hdfrdr(sd_id, 'sclfct', sclfct, istatus)
    CALL abss2ctim(itime, iyr, imon, idy, ihr, imin, isec )
    iyr=MOD(iyr,100)
    WRITE(6,'(/a,i2.2,a,i2.2,a,i2.2,1X,i2.2,a,i2.2,a)')               &
       'Read grid info from remapped raw radar recorded at: ',        &
       imon,'/',idy,'/',iyr,ihr,':',imin,' UTC'

    CALL hdfclose(sd_id,istatus)
    istatus=1
  END IF

  RETURN
!
!  Destination for hdf read error
!

  115   CONTINUE
  WRITE(6,'(/a/)') ' Error reading data in get_grid_from_rad.'

  istatus=-11

  RETURN

END SUBROUTINE get_grid_from_rad
