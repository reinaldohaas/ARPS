!##################################################################
!##################################################################
!######                                                      ######
!######                PROGRAM JOINBIN2HDF                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################

PROGRAM join_bin2hdf
!
!-----------------------------------------------------------------------
!
!  To join together a set of ARPS history or data files produced by the
!  processors of MPP machines with message passing.
!
!  Input data file is in binary format and the output is in HDF4 format
!
!-----------------------------------------------------------------------
!
!  AUTHOR:
!  Yunheng Wang (05/16/2002)
!  Based on joinhdf.f90
!
!  MODIFICATION HISTORY.
!
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INCLUDE 'mp.inc'
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER nxsm,nysm,nz, nzsoil,nstyps
  INTEGER ireturn

  CHARACTER (LEN=256) :: filename
  CHARACTER (LEN=256) :: filename1

  INTEGER :: hdfcmpropt

  INTEGER :: grdbas
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  WRITE(*,FMT='(/a)',ADVANCE='NO') " Enter the filename (base name):  "
  READ(*,'(a)') filename
  WRITE(*,FMT='(a)',ADVANCE='NO')  " Enter nproc_x, nproc_y        :  "
  READ(*,*) nproc_x, nproc_y
  WRITE(*,FMT='(a)',ADVANCE='NO')  " Enter HDF4 compression option :  "
!  WRITE(*,*) "     = 0  no compression;"
!  WRITE(*,*) "     = 1, fast gzip compression;"
!  WRITE(*,*) "     = 2, high gzip compression;"
!  WRITE(*,*) "     = 3, adaptive or skipping Huffman compression;"
!  WRITE(*,*) "     = 4-7, as above plus mapping reals to 16 bit integers."
!  WRITE(*,*) "       Note that only options 0-2 work on Cray platforms."
!  READ(*,*) hdfcmpropt

  hdfcmpropt = 5                  ! hard-coded as Kevin's sugguestion.
  WRITE(*,'(a//)') "(Using hdfcmpropt = 5 as default value)"

  !filename1 = TRIM(filename)//"_0101"
  CALL gtsplitfn(filename,1,1,1,1,1,1,0,0,1,2,filename1,ireturn)

  CALL get_dims_only(1,filename1,nxsm,nysm,nz,nzsoil,nstyps, grdbas, ireturn)

  IF (ireturn /= 0) THEN
    WRITE (6,*) 'JOINFILE: WARNING, error returned from get_dims_from_data', &
       ireturn
  ENDIF

  WRITE (6, '(/a/)') 'Joining files ...'

  CALL joinbin2hdf (filename,nxsm,nysm,nz,nzsoil,nstyps,grdbas,hdfcmpropt)

  WRITE (6, '(/a/)') 'Done joining files ...'

END PROGRAM join_bin2hdf

!
!##################################################################
!##################################################################
!######                                                      ######
!######            SUBROUTINE GET_DIMS_FROM_DATA             ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE get_dims_only(hinfmt,hisfile, nx,ny,nz,nzsoil,nstyps, &
                         grdbas, ireturn)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read in grid dimensions from base state/grid history data.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  7/17/2000
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    hinfmt   The format of the history data dump
!    grdbasfn Name of the grid/base state array file
!
!  OUTPUT:
!
!    nx,ny,nz The dimension of data arrays
!    nzsoil   Number of soil levels
!    nstyps   The number of soil types
!    ireturn  Return status indicator
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: hinfmt            ! The format of the history data dump
  CHARACTER (LEN=*  ) :: hisfile  ! Name of the grid/base state array file
  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of soil levels
  INTEGER :: nstyps            ! Number of soil types
  INTEGER :: ireturn,nchanl,ierr,istat,packed
  INTEGER, INTENT(OUT) :: grdbas
  LOGICAL :: fexist

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  ireturn = 0
!
!-----------------------------------------------------------------------
!
!  Open and read grid and base state data file depending on the
!  values of parameters grdin and basin, which are read in from the
!  time dependent data set. If grdin or basin is zero, the grid and
!  base state arrays have to be read in from a separate file.
!
!-----------------------------------------------------------------------
!
  IF ( hinfmt == 9 ) THEN
    WRITE(6,'(/2(1x,a/))')                                         &
    'Read in from a GrADS file is currently not supported. ',      &
    'Program stopped in get_dims_from_data.'
    CALL arpsstop('arpsstop called from get_dims_from_data',1)
  END IF

  INQUIRE(FILE=trim(hisfile), EXIST = fexist )
  IF( fexist ) GO TO 200

  INQUIRE(FILE=trim(hisfile)//'.Z', EXIST = fexist )
  IF( fexist ) THEN
    CALL uncmprs( trim(hisfile)//'.Z' )
    GO TO 200
  END IF

  INQUIRE(FILE=trim(hisfile)//'.gz', EXIST = fexist )
  IF( fexist ) THEN
    CALL uncmprs( trim(hisfile)//'.gz' )
    GO TO 200
  END IF

  WRITE(6,'(/2(1x,a/))')                                              &
    'File '//trim(hisfile)//                                          &
    ' or its compressed version not found.',                          &
    'Program stopped in get_dims_from_data.'
  CALL arpsstop('arpsstop called from get_dims_from_data',1)

  200     CONTINUE

!
!-----------------------------------------------------------------------
!
!  Read grid and base state fields.
!
!-----------------------------------------------------------------------
!
  IF( hinfmt == 1 .OR. hinfmt == 6 ) THEN

    CALL getunit( nchanl )
!
!-----------------------------------------------------------------------
!
!  Cray routines to force binary data file to be in the IEEE format
!
!-----------------------------------------------------------------------
!
    CALL asnctl ('NEWLOCAL', 1, ierr)
    CALL asnfile(trim(hisfile), '-F f77 -N ieee', ierr)

    OPEN(UNIT=nchanl,FILE=trim(hisfile),                           &
         STATUS='old',FORM='unformatted',IOSTAT=istat)

    IF( istat /= 0 ) GO TO 999

    CALL bin_getdims_only(nchanl, nx,ny,nz,nzsoil,nstyps, grdbas,ireturn)

    CLOSE(UNIT=nchanl)
    CALL retunit( nchanl )

!  ELSE IF( hinfmt == 2 ) THEN
!
!    CALL getunit( nchanl )
!    OPEN(UNIT=nchanl,FILE=trim(grdbasfn),                           &
!         STATUS='old',FORM='formatted',IOSTAT=istat)
!
!    IF( istat /= 0 ) GO TO 999
!
!    CALL asc_getdims(nchanl, nx,ny,nz,nzsoil,nstyps, ireturn)
!
!    CLOSE(UNIT=nchanl)
!    CALL retunit( nchanl )

  ELSE IF( hinfmt == 3 ) THEN

    CALL hdf_getdims(trim(hisfile),nx,ny,nz,nzsoil,nstyps, ireturn)

!  ELSE IF( hinfmt == 4 ) THEN
!
!    CALL getunit( nchanl )
!    OPEN(UNIT=nchanl,FILE=trim(grdbasfn),                           &
!         STATUS='old',FORM='unformatted',IOSTAT=istat)
!
!    IF( istat /= 0 ) GO TO 999
!
!    CALL pak_getdims(nchanl, nx,ny,nz,nzsoil,nstyps, ireturn)
!
!    CLOSE(UNIT=nchanl)
!    CALL retunit( nchanl )

  ELSE IF (hinfmt == 7 .OR. hinfmt == 8) THEN ! NetCDF format

    packed = 0
    CALL netopen(TRIM(hisfile),'R',nchanl)

    CALL net_getdims(nchanl,nx,ny,nz,nzsoil,nstyps,ireturn)

    CALL netclose(nchanl)

!
!-----------------------------------------------------------------------
!
!  Cray routines to force binary data file to be in the IEEE format
!
!-----------------------------------------------------------------------
!
!  ELSE IF( hinfmt == 10 ) THEN
!    CALL asnctl ('NEWLOCAL', 1, ierr)
!    CALL asnfile(trim(grdbasfn), '-F f77 -N ieee', ierr)
!
!    CALL getunit( nchanl )
!
!    OPEN(UNIT=nchanl,FILE=trim(grdbasfn),STATUS='old',              &
!         FORM='unformatted',IOSTAT= istat )
!
!    CALL grib_getdims(nchanl, nx,ny,nz,nzsoil,nstyps,ireturn)
!
!    CLOSE(UNIT=nchanl)
!    CALL retunit( nchanl )

  ELSE

    WRITE(6,'(a,i3,a)')                                                 &
        ' Data format flag had an invalid value ',                      &
          hinfmt ,' program stopped.'
    CALL arpsstop('arpsstop called from get_dims_from_data wrong flag',1)

  END IF

  RETURN

  999   CONTINUE
  WRITE(6,'(1x,a,a,/1x,i3,a)')                                          &
      'Error occured when opening file ',trim(hisfile),                 &
      'using FORTRAN unit ',nchanl,' Program stopped in get_dims_from_data.'

  CALL arpsstop('arpsstop called from get_dims_from_data opening file',1)

END SUBROUTINE get_dims_only

SUBROUTINE bin_getdims_only(inch, nx,ny,nz,nzsoil,nstyps, grdbas,ireturn)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!  Read in grid dimensions from history data.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  7/17/2000.
!
!  MODIFICATION HISTORY:
!
!  09/12/2006 (Y. Wang)
!  Changed to read ARPS history file instead of grid/base file. Beside
!  ARPS domains, it also set the microphysics scalar indices.
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    inch     Channel number for binary reading.
!
!  OUTPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of soil levels in the vertical
!
!    nstyps   Number of soil types
!
!    ireturn  Return status indicator
!             =0, successful read of all data
!             =1, error reading data
!             =2, end-of-file reached during read attempt
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of soil levels
  INTEGER :: nstyps            ! Number of soil types
  INTEGER :: inch              ! Channel number for binary reading
  INTEGER :: ireturn           ! Return status indicator
  INTEGER, INTENT(OUT) :: grdbas
!
! 06/28/2002 Zuwen He
!
! fmtver??: to label each data a version.
! intver??: an integer to allow faster comparison than fmtver??,
!           which are strings.
!
! Verion 5.00: significant change in soil variables since version 4.10.
!
  CHARACTER (LEN=40) :: fmtver320,fmtver400,fmtver410,fmtver500,fmtver530
  INTEGER  :: intver,intver320,intver400,intver410,intver500,intver530

  PARAMETER (fmtver320='003.20 Binary Data',intver320=320)
  PARAMETER (fmtver400='004.00 Binary Data',intver400=400)
  PARAMETER (fmtver410='004.10 Binary Data',intver410=410)
  PARAMETER (fmtver500='005.00 Binary Data',intver500=500)
  PARAMETER (fmtver530='005.30 Binary Data',intver530=530)

  CHARACTER (LEN=40) :: fmtverin

  CHARACTER (LEN=10) :: tmunit
  INTEGER :: i
  REAL    :: time

  INTEGER :: idummy,totin,mstin,icein
  INTEGER :: grdin, basin
  REAL(4) :: rdummy

  INCLUDE   'globcst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  READ(inch,ERR=110,END=120) fmtverin

  IF (fmtverin == fmtver320) THEN
    intver=intver320
  ELSE IF (fmtverin == fmtver400) THEN
    intver=intver400
  ELSE IF (fmtverin == fmtver410) THEN
    intver=intver410
  ELSE IF (fmtverin == fmtver500) THEN
    intver=intver500
  ELSE IF (fmtverin == fmtver530) THEN
    intver=intver530
  ELSE
    WRITE(6,'(/2x,a,a,a/)')                        &
        'Incoming data format, fmtverin=',fmtverin,             &
        ', not found. The Job stopped.'
    CALL arpsstop('arpstop called from bin_getdims. ',1)
  END IF

  WRITE(6,'(/2x,a,a/)') 'Incoming data format, fmtverin=',fmtverin

  READ(inch,ERR=110,END=120) runname
  WRITE(6,'(//''  THE NAME OF THIS RUN IS:  '',A//)') runname

  READ(inch,ERR=110,END=120) nocmnt

  IF( nocmnt > 0 ) THEN
    DO i=1,nocmnt
      READ(inch,ERR=110,END=120)
    END DO
  END IF

  READ(inch,ERR=110,END=120) rdummy,tmunit
  time = rdummy

  IF (intver <= intver410) THEN

    READ(inch,ERR=110,END=120) nx, ny, nz
    nzsoil = 2  ! for version prior to 410, it is a two-layer soil model

  ELSE IF (intver >= intver500) THEN

    READ(inch,ERR=110,END=120) nx, ny, nz,nzsoil

  END IF

  WRITE(6,'(1x,a,4i5,a)') 'nx,ny,nz,nzsoil read in from data are ',     &
                           nx,ny,nz,nzsoil,' respectively.'

  READ(inch,ERR=110,END=120) grdin,basin,idummy, mstin,  icein,         &
                             idummy,idummy,idummy,idummy,  totin,       &
                             idummy,nscalar,idummy,idummy,idummy,       &
                             idummy,idummy,idummy,idummy,idummy

  IF (grdin == 1 .AND. basin == 1) THEN  ! suppose it is time-independent
                                         ! grid and base file
    grdbas = 1
  ELSE
    grdbas = 0
  END IF

  IF (totin /= 0) THEN
    READ(inch,ERR=110,END=120) ! block of 20 REALs ...
    READ(inch,ERR=110,END=120) nstyps,idummy,idummy,idummy,idummy,  &
                               idummy,idummy,idummy,idummy,idummy,  &
                               idummy,idummy,idummy,idummy,idummy,  &
                               idummy,idummy,idummy,idummy,idummy
    WRITE (6,*) " nstyps read in as: ",nstyps
  ELSE
    WRITE (6,*) " nstyps not defined in data, set to 4"
    nstyps = 4
  ENDIF

  ireturn = 0

  RETURN

  110   CONTINUE
  WRITE(6,'(/a/)') ' Error reading data in BIN_GETDIMS'
  ireturn=1
  RETURN

  120   CONTINUE
  WRITE(6,'(/a/)') ' End of file reached in BIN_GETDIMS'
  ireturn=2
  RETURN

END SUBROUTINE bin_getdims_only
