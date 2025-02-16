!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE RDSATFLD                     ######
!######                                                      ######
!######                  Developed by                        ######
!######    Center for Analysis and Prediction of Storms      ######
!######             University of Oklahoma.                  ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE rdsatfld(nx,ny,nfield,                                       &
           sfname,satname,latsat,lonsat,                                &
           itime,isource,fldname,satfld,istatus)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Reads remapped satellite data to a file as one or
!  more 2-d fields.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  09/20/97
!
!  MODIFICATION HISTORY:
!
!  2001/08/03 (Gene Bassett)
!  Read satellite grid parameters into temporary variables and compare
!  to those defined in the common block.
!
!  04/17/01 (Leilei Wang)
!  Added processing for hdf files.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!    nx,ny     horizontal dimensions
!    nfield    number of satellite fields to write
!    sfname    satellite file name (character string)
!    satnam    satellite name (character*6)
!    latsat    sub-satellite latitude (degrees N)
!    lonsat    sub-satellite longitude (degrees E)
!    itime     time, seconds since 1960
!    isource   source number
!                1= GVAR raw 2-byte data file
!                2= IDD 1-byte datafeed
!    fldname   name of variable(s) (character*6 array)
!    satfld    satellite data
!
!  OUTPUT:
!    data are written to file
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: nx,ny
  INTEGER :: nfield
  CHARACTER (LEN=256) :: sfname
  CHARACTER (LEN=6) :: satname
  REAL :: latsat
  REAL :: lonsat
  INTEGER :: itime
  INTEGER :: isource
  CHARACTER (LEN=6) :: fldname(nfield)
  REAL :: satfld(nx,ny,nfield)
  INTEGER :: istatus,istat
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=80) :: dummy_name
  INTEGER :: iunit,iopen
  INTEGER :: idummy
  INTEGER :: nxin,nyin,nfieldin
  REAL :: dxin,dyin,ctrlatin,ctrlonin,trlat1in,trlat2in,trlonin,sclfctin
  INTEGER :: mprojin
  REAL :: rdummy
  INTEGER :: ireturn
  INTEGER :: dmpfmt,sd_id,lens,i
!
!-----------------------------------------------------------------------
!
! hdf temporary arrays
!
!-----------------------------------------------------------------------
!
  INTEGER (KIND=selected_int_kind(4)), ALLOCATABLE :: itmp(:,:,:)
  REAL, ALLOCATABLE :: hmax(:)
  REAL, ALLOCATABLE :: hmin(:)
!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'        ! Grid parameters
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  PRINT *, ' sfname= ',sfname
  PRINT *, ' nx,ny,nfield= ',nx,ny,nfield
!
  CALL getunit(iunit)
!
!-----------------------------------------------------------------------
!
!  Open file for reading
!
!-----------------------------------------------------------------------
!
  lens = LEN(trim(sfname))
  IF(sfname(lens-3:lens)=='hdf4')THEN
    dmpfmt=2
  ELSE
    dmpfmt=1
  END IF

  WRITE(6,'(a,i4,a,a)') ' rdsatfld: dmpfmt= ',dmpfmt,'  sfname=',sfname

  IF (dmpfmt == 1) THEN
    OPEN(iunit,IOSTAT=iopen,FILE=trim(sfname),STATUS='old',             &
       FORM='unformatted')

    IF(iopen == 0) THEN
!
!-----------------------------------------------------------------------
!
!  Read satellite description variables
!
!-----------------------------------------------------------------------
!
      READ (iunit,ERR=200) satname
      READ (iunit,ERR=200) nxin,nyin,nfieldin,itime,idummy,             &
                         idummy,idummy,idummy,idummy,idummy

!-----------------------------------------------------------------------
!
!  Read grid description variables
!  This should provide enough info to uniquely identify the 2-d grid.
!
!-----------------------------------------------------------------------
!

      READ (iunit,ERR=200) dummy_name
      READ (iunit,ERR=200) idummy,idummy,mprojin,idummy,idummy,         &
                 idummy,idummy,idummy,idummy,idummy
      READ (iunit,ERR=200) dxin,dyin,rdummy,rdummy,ctrlatin,            &
                 ctrlonin,trlat1in,trlat2in,trlonin,sclfctin,           &
                 latsat,lonsat,rdummy,rdummy,rdummy
    ELSE
      istatus=iopen
    END IF

  ELSE  ! HDF4 format
    CALL hdfopen(trim(sfname), 1, sd_id)
    IF (sd_id < 0) THEN
      WRITE (6,*) "RDSATFLD: ERROR opening ",                          &
                  trim(sfname)," for reading."
      GO TO 200
    END IF
    IOPEN = 0
    CALL hdfrdc(sd_id, 6, 'satname', satname, istat)
    CALL hdfrdi(sd_id, 'nx', nxin, istat)
    CALL hdfrdi(sd_id, 'ny', nyin, istat)
    CALL hdfrdi(sd_id, 'nfield', nfieldin, istat)
    CALL hdfrdi(sd_id, 'itime', itime, istat)

    CALL hdfrdc(sd_id, 4, 'runname', dummy_name, istat)
    CALL hdfrdi(sd_id, 'hdmpfmt', idummy, istat)
    CALL hdfrdi(sd_id, 'strhopt', idummy, istat)
    CALL hdfrdi(sd_id, 'mapproj', mprojin, istat)
    CALL hdfrdr(sd_id, 'dx', dxin, istat)
    CALL hdfrdr(sd_id, 'dy', dyin, istat)
    CALL hdfrdr(sd_id, 'dz', rdummy, istat)
    CALL hdfrdr(sd_id, 'dzmin', rdummy, istat)
    CALL hdfrdr(sd_id, 'ctrlat', ctrlatin, istat)
    CALL hdfrdr(sd_id, 'ctrlon', ctrlonin, istat)
    CALL hdfrdr(sd_id, 'trulat1', trlat1in, istat)
    CALL hdfrdr(sd_id, 'trulat2', trlat2in, istat)
    CALL hdfrdr(sd_id, 'trulon', trlonin, istat)
    CALL hdfrdr(sd_id, 'sclfct', sclfctin, istat)
    CALL hdfrdr(sd_id, 'latsat', latsat, istat)
    CALL hdfrdr(sd_id, 'lonsat', lonsat, istat)
  ENDIF
!
!-----------------------------------------------------------------------
!
!  Check dimensions and grid parameters of incoming data.
!
!-----------------------------------------------------------------------
!
  IF(iopen==0)THEN
    CALL checkgrid2d(nx,ny,nxin,nyin,                                   &
          dx,dy,ctrlat,ctrlon,                                          &
          mapproj,trulat1,trulat2,trulon,sclfct,                        &
          dxin,dyin,ctrlatin,ctrlonin,                                  &
          mprojin,trlat1in,trlat2in,trlonin,sclfctin,ireturn)

    IF (ireturn /= 0) THEN
      WRITE (6,*) "RDSATFLD: ERROR, grid parameter mismatch in ",       &
                  "satellite date file ",trim(sfname)
      GOTO 200
    END IF

    IF(nfieldin > nfield) THEN
      WRITE(6,'(a,/a,i4,a,i4)')                                        &
        ' Mis-match in number of fields in file.',                     &
        ' nfield in file:',nfieldin,' nfield expected:',nfield
    ELSE IF (nfieldin < nfield) THEN
      WRITE(6,'(a,/a,i4,a,i4,/a)')                                     &
        ' Mis-match in number of fields in file.',                     &
        ' nfield in file:',nfieldin,' nfield expected:',nfield,        &
        ' Rerun mci2arps using version compatible with this version'
      GOTO 200
    END IF
  ELSE
    istatus=iopen
  END IF

!
!-----------------------------------------------------------------------
!
!  Read 2-d fields.
!
!-----------------------------------------------------------------------
!
  IF(dmpfmt==1)THEN
    IF(iopen==0)THEN

      READ(iunit,ERR=200) fldname
      READ(iunit,ERR=200) satfld
!
      CLOSE(iunit)
      CALL retunit(iunit)
    ELSE
      istatus=iopen
    END IF

  ELSE
    DO i = 1,nfield
      CALL hdfrdc(sd_id, 6, 'fldname', fldname(i), istat)
    ENDDO
    ALLOCATE (itmp(nx,ny,nfield),stat=istat)
    IF (istat /= 0) THEN
      WRITE (6,*) "HDFREAD: ERROR allocating itmp, returning"
      CALL arpsstop("Allocate1 failure in RDSATFLD",1)
    END IF
    ALLOCATE (hmax(nfield),stat=istat)
    IF (istat /= 0) THEN
      WRITE (6,*) "HDFREAD: ERROR allocating hmax, returning"
      CALL arpsstop("Allocate2 failure in RDSATFLD",1)
    END IF
    ALLOCATE (hmin(nfield),stat=istat)
    IF (istat /= 0) THEN
      WRITE (6,*) "HDFREAD: ERROR allocating hmin, returning"
      CALL arpsstop("Allocate3 failure in RDSATFLD",1)
    END IF

    CALL hdfrd3d(sd_id,"satfld",nx,ny,nfield,satfld,istat,itmp,          &
                 hmax,hmin)
    IF (istat /= 0) GO TO 200

    CALL hdfclose(sd_id,istat)
    IF (istat /= 0) THEN
      WRITE (6,*) "RDSATFLD: ERROR on closing file ",trim(sfname),       &
                " (status",istat,")"
    END IF

    DEALLOCATE(itmp)
    DEALLOCATE(hmax)
    DEALLOCATE(hmin)

  END IF
!
!-----------------------------------------------------------------------
!
!  Report on what data were read
!
!-----------------------------------------------------------------------
!
  IF(dmpfmt==1)THEN
    IF(iopen==0)THEN
      WRITE(6,'(//a,a,a,a)') ' Read ',fldname(1),' from ',satname
!
      PRINT *, ' satname= ',satname
      PRINT *, ' lat,lon= ',latsat,lonsat
      PRINT *, ' itime= ',itime
      PRINT *, ' isource = ',isource
      PRINT *, ' fldname= ',fldname(1)
      PRINT *, ' satfld(1,1,1)= ',satfld(1,1,1)
      PRINT *, ' satfld(nx,ny,1) = ',satfld(nx,ny,1)

      istatus=0
    ELSE
      istatus=iopen
    END IF
  ELSE
    istatus=0
  END IF

  RETURN

  200 CONTINUE

  istatus=-1

  RETURN
END SUBROUTINE rdsatfld
