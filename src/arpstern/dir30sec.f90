PROGRAM dir30sec

!##################################################################
!##################################################################
!######                                                      ######
!######                 PROGRAM DIR30SEC                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!
!  PURPOSE:
!
!  A stand alone program which rewrites 30 second terrain elevation data
!  from dma_elev.dat recieved from NCAR DATA SERVICES into a direct access
!  file using either 2 of 8 byte integer format depending upon the type
!  of machine used.  The purpose of writing a direct access data format
!  is to expedite the run time of the ARPSTERN11.F terrain processing
!  program for initializing the ARPS terrain field.
!  In addition, a dir1km.hdr file will be created and needs to be used
!
!  INPUT:
!            dma_elev.dat     (NCAR DATA FILE 30 SECOND RESOLUTION)
!
!            arpstern.input    ARPSTERN11.F input file
!
!
!
!
!  OUTPUT:
!            dir30sec.dat      direct access 30 second unformatted
!                              data file with 1099 seperate records.
!
!            dir30sec.hdr      header file which has record numbers
!                              corresponding to specific blocks or
!                              1 degree by 1 degree lat/lon areas.
!
!  This program reads the ascii version of the 30 second resolution
!  terrain data set for the continental US (NCAR FILE ds756.1).
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  Dan Weber
!           1/12/94
!
!  Modification history:
!
!    Dan Weber 6/22/94.
!       Added documentation.
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!  Variable Declarations
!
!-----------------------------------------------------------------------

  IMPLICIT NONE
!
!  INPUT:
!
  INTEGER :: ndx              !  Number of data points in the longitudinal
                              !  direction
  INTEGER :: ndy              !  Number of data points in the latitudinal
                              !  direction

  PARAMETER (ndx=120,ndy=120)

  CHARACTER (LEN=11      ) :: desc  !  Character string for reading the header
                                    !  portion of the dma_elev.dat file
  INTEGER :: id(5)            !  Array used to read the header portion of
                              !  dma_elev.dat NCAR file.
  INTEGER*2 nlval(121,121) !  array used to read initial NCAR data file
  INTEGER*2 nval(120,120)  !  Array used to write unformatted terrain
                           !  data block to dir30sec.dat file


!
!  COMMON BLOCK VARIABLES:
!


!
!  LOCAL VARIABLES:
!
  INTEGER :: comtype          !  Computer type for program to run on...
                              !  for IBM, comtype = 1
                              !     CRAY,         = 4
  INTEGER :: i,j,n
  INTEGER :: ltdir            !  Dummy variable used for namelist data
  CHARACTER (LEN=80  ) :: tdatadir  !  Directory in which the data file
                                    !  dma_elev.dat is stored.
  CHARACTER (LEN=80   ) :: terndir  !  Directory in which the original ASCII terrain
                                    !  data files are stored.
  INTEGER :: lterndir         !  Length of the non-blank part of string terndir.

  INTEGER :: temi             !  Dummy variable used for namelist data
  INTEGER :: analtype         !  Dummy variable used for namelist data
  INTEGER :: mapproj          !  Dummy variable used for namelist data
  INTEGER :: itertype         !  Dummy variable used for namelist data
  INTEGER :: rmsopt           !  Dummy variable used for namelist data

!
!  OUTPUT:
!

  REAL :: aval(120,120)       !  Array used to write unformatted terrain
                              !  data block to dir30sec.dat file

!-----------------------------------------------------------------------
!
!  Defining namelist:
!
!-----------------------------------------------------------------------

  NAMELIST /terraind/analtype,mapproj,itertype,rmsopt,                  &
      comtype,tdatadir,terndir

!-----------------------------------------------------------------------
!
!  Begin executable portion of the code.
!
!-----------------------------------------------------------------------

  READ(5,terraind)
!
!  open dma_elev.dat, dir30sec.dat and dir30sec.hdr files...
!

  lterndir = LEN(terndir)
  CALL strlnth( terndir, lterndir )

  IF( lterndir == 0 ) THEN
    terndir = '.'
    lterndir=1
  END IF

!
!  open elev.dat, dir1deg.dat and dir1deg.hdr files...
!
  OPEN(10,FILE=terndir(1:lterndir)//'/dma_elev.dat',STATUS='old')

  OPEN(11,FILE='dir30sec.dat',ACCESS='direct',FORM='unformatted',       &
      STATUS='unknown',RECL=28800*comtype)
  OPEN(12,FILE='dir30sec.hdr',FORM='formatted',STATUS='unknown')


!
!  Read in the initial record from the dma_elev.dat file.
!  Setting the record counter n=1.
!

  n=1
  5    READ(10,'(5i6,1x,a11)',END=90)(id(i),i=1,5),desc
  READ(10,'(30i4)')((nlval(i,j),i=1,121),j=1,121)

!
!  Convert from feet to meters...
!  Writing (120,120) from the (121,121) array.
!  The extra data was found to be a source of bogus
!  data after inspection of the data squares which were flagged
!  in a seperate program.
!

  DO i=1,ndx
    DO j=1,ndy
      nval(i,j)=nlval(i,j)*20-4000    ! unpacking the data
      aval(i,j)=nval(i,j)
      aval(i,j)=aval(i,j)*0.3048+0.5  ! converting from feet to meters..
      nval(i,j)=INT(aval(i,j))
      nval(i,j)=(nval(i,j)+4000)/20   ! putting in packed form....
    END DO
  END DO

!
!  Setting the remaining area over Lake Huron to 580 feet
!

  IF(id(1) == 134.AND.id(2) == 277)THEN ! id(1) = lat(N)+90., id(2)=lon (E)
    DO i=1,121
      DO j=1,121
        IF(nlval(i,j) == 0)nlval(i,j)=209   ! 209 meters = 580 feet ....
      END DO
    END DO
  ELSE
  END IF

!
!  write the 1 degree by 1 degree box header to dir30sec.hdr and
!  the 120,120 data array to dir30sec.dat
!

  WRITE(11,REC=n) nval
  WRITE(12,'(i6,i6,i6)') id(1)-90,id(2),n

  PRINT *,'DATA WRITTEN FOR RECORD NUMBER = ',n
  n=n+1
  GO TO 5
  90   CLOSE (10)
  CLOSE (11)
  CLOSE (12)

  STOP
END PROGRAM dir30sec
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE STRLNTH                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE strlnth( string, length )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Return the length of the non-blank part of a character string.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Ming Xue
!  11/20/91
!
!  MODIFICATION HISTORY:
!
!  5/05/92 (M. Xue)
!  Added full documentation.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    string   A character string
!    length   The declared length of the character string 'string'.
!
!  OUTPUT:
!
!    length   The length of the non-blank part of the string.
!
!-----------------------------------------------------------------------
!

!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  CHARACTER (LEN=*     ) :: string ! A character string for the name of
                                   ! this run.
  INTEGER :: length            ! The length of the non-blank part
                               ! of a string.

  INTEGER :: i
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  DO i = length,1,-1

    IF(string(i:i) /= ' ') EXIT

  END DO

!  200   CONTINUE

  length = MAX(1,i)

  RETURN
END SUBROUTINE strlnth
