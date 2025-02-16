!
!##################################################################
!##################################################################
!######                                                      ######
!######              SUBROUTINE WTSATFLD                     ######
!######                                                      ######
!######                  Developed by                        ######
!######    Center for Analysis and Prediction of Storms      ######
!######             University of Oklahoma.                  ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE wtsatfld(nx,ny,nfield,                                       &
                    sfname,satname,latsat,lonsat,                       &
                    iyr,imon,iday,ihr,imin,isec,isource,                &
                    dmpfmt,hdf4cmpr,                                    &
                    fldname,satfld)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Writes remapped satellite data to a file as one or
!  more 2-d fields.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  09/20/97
!
!  MODIFICATION HISTORY:
!
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
!    iyr       year
!    imon      month
!    iday      day
!    ihr       hour
!    imin      min
!    isec      sec
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
  CHARACTER (LEN=6)   :: satname
  REAL :: latsat
  REAL :: lonsat
  INTEGER :: iyr,imon,iday,ihr,imin,isec
  INTEGER :: isource
  CHARACTER (LEN=6) :: fldname(nfield)
  REAL :: satfld(nx,ny,nfield)
  integer::dmpfmt
  integer::hdf4cmpr
!
!-----------------------------------------------------------------------
!
!  Misc local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: iunit,myr,kyr,itime
  INTEGER :: idummy
  REAL :: rdummy
  INTEGER(2), allocatable :: itmp(:,:,:) ! Temporary array
  REAL, allocatable :: hmax(:), hmin(:) ! Temporary array
  integer::sd_id,i,istat

!
!-----------------------------------------------------------------------
!
!  Include files
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
 IF (hdf4cmpr > 3) THEN
   ALLOCATE (itmp(nx,ny,nfield),stat=istat)
   IF (istat /= 0) THEN
     WRITE (6,*) "HDFDUMP: ERROR allocating itmp, returning"
     RETURN
   END IF
   ALLOCATE (hmax(nfield),stat=istat)
    IF (istat /= 0) THEN
      WRITE (6,*) "HDFDUMP: ERROR allocating hmax, returning"
      RETURN
    END IF
    ALLOCATE (hmin(nfield),stat=istat)
    IF (istat /= 0) THEN
      WRITE (6,*) "HDFDUMP: ERROR allocating hmin, returning"
      RETURN
    END IF
 ENDIF

  myr=1900+iyr
  IF(myr < 1960) myr=myr+100
  CALL ctim2abss(myr,imon,iday,ihr,imin,isec,itime)

  IF(dmpfmt==1)THEN
    CALL getunit(iunit)
!
!-----------------------------------------------------------------------
!
!  Open file for output
!
!-----------------------------------------------------------------------
!
    OPEN(iunit,FILE=TRIM(sfname),STATUS='unknown',                      &
               FORM='unformatted')
!
!-----------------------------------------------------------------------
!
!  Write satellite description variables
!
!-----------------------------------------------------------------------
!
  idummy=0
  rdummy=0.
  WRITE(iunit) satname
  WRITE(iunit) nx,ny,nfield,itime,idummy,                               &
               idummy,idummy,idummy,idummy,idummy
!
!-----------------------------------------------------------------------
!
!  Write grid description variables
!  This should provide enough info to uniquely identify the 2-d grid.
!
!-----------------------------------------------------------------------
!
  WRITE(iunit) runname
  WRITE(iunit) hdmpfmt,strhopt,mapproj,idummy,idummy,                   &
               idummy,idummy,idummy,idummy,idummy
  WRITE(iunit) dx,dy,dz,dzmin,ctrlat,                                   &
               ctrlon,trulat1,trulat2,trulon,sclfct,                    &
               latsat,lonsat,rdummy,rdummy,rdummy
!
!-----------------------------------------------------------------------
!
!  Write 2-d fields.
!
!-----------------------------------------------------------------------
!
  WRITE(iunit) fldname
  WRITE(iunit) satfld
  ELSE     !HDF4 format
    CALL hdfopen(trim(sfname), 2, sd_id)
    IF (sd_id < 0) THEN
      WRITE (6,'(a,a,a)') "WTSATFLD: ERROR opening ",                   &
                  trim(sfname)," for writing."
      istat = 1
      STOP
    END IF
    CALL hdfwrtc(sd_id, 6, 'satname', satname, istat)
    CALL hdfwrti(sd_id, 'nx', nx, istat)
    CALL hdfwrti(sd_id, 'ny', ny, istat)
    CALL hdfwrti(sd_id, 'nfield', nfield, istat)
    CALL hdfwrti(sd_id, 'itime', itime, istat)

    CALL hdfwrtc(sd_id, 4, 'runname', runname, istat)
    CALL hdfwrti(sd_id, 'hdmpfmt', hdmpfmt, istat)
    CALL hdfwrti(sd_id, 'strhopt', strhopt, istat)
    CALL hdfwrti(sd_id, 'mapproj', mapproj, istat)
    CALL hdfwrtr(sd_id, 'dx', dx, istat)
    CALL hdfwrtr(sd_id, 'dy', dy, istat)
    CALL hdfwrtr(sd_id, 'dz', dz, istat)
    CALL hdfwrtr(sd_id, 'dzmin', dzmin, istat)
    CALL hdfwrtr(sd_id, 'ctrlat', ctrlat, istat)
    CALL hdfwrtr(sd_id, 'ctrlon', ctrlon, istat)
    CALL hdfwrtr(sd_id, 'trulat1', trulat1, istat)
    CALL hdfwrtr(sd_id, 'trulat2', trulat2, istat)
    CALL hdfwrtr(sd_id, 'trulon', trulon, istat)
    CALL hdfwrtr(sd_id, 'sclfct', sclfct, istat)
    CALL hdfwrtr(sd_id, 'latsat', latsat, istat)
    CALL hdfwrtr(sd_id, 'lonsat', lonsat, istat)
    DO i = 1,nfield
      WRITE(6,'(a,i4,a,a)') 'wtsatfld fldname(',i,'): ',fldname(i)
      CALL hdfwrtc(sd_id, 6, 'fldname', fldname(i), istat)
    ENDDO
     
    CALL hdfwrt3d(satfld,nx,ny,nfield,sd_id,0,hdf4cmpr,                &
                  'satfld','satfld','',                &
                   itmp,hmax,hmin)
  ENDIF
!
  IF(dmpfmt==1)THEN
    CLOSE(iunit)
    CALL retunit(iunit)
  ELSE
    CALL hdfclose(sd_id,istat)
    IF (istat /= 0) THEN
    WRITE (6,*) "HDFDUMP: ERROR on closing file ",trim(sfname),       &
                " (status",istat,")"
    DEALLOCATE (itmp,stat=istat)
    DEALLOCATE (hmax,stat=istat)
    DEALLOCATE (hmin,stat=istat)
    ENDIF
  ENDIF
  
  kyr = iyr
  if ( kyr > 99 ) kyr = kyr - 100               ! Y2K

!
!-----------------------------------------------------------------------
!
!  Report on what data were written
!
!-----------------------------------------------------------------------
!
  WRITE(6,'(//a,i2.2,i2.2,i2.2,a1,i2.2,a1,i2.2)')                       &
                    ' Wrote satellite fields for time ',                &
                      kyr,imon,iday,' ',ihr,':',imin
!
  RETURN
END SUBROUTINE wtsatfld
