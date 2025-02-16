PROGRAM dir1deg

!##################################################################
!##################################################################
!######                                                      ######
!######                 PROGRAM DIR1DEG                      ######
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
!  A stand alone program which rewrites 1 degree terrain elevation data
!  from elev.dat recieved from NCAR DATA SERVICES into a direct access
!  file using either 2 of 8 byte integer format depending upon the type
!  of machine used.  The purpose of writing a direct access data format
!  is to expedite the run time of the ARPSTERN11.F terrain processing
!  program for initializing the ARPS terrain field.
!  In addition, a dir1deg.hdr file (header file) will be created and
!  used by arpstern11.f.
!
!  INPUT:
!            elev.dat          (NCAR DATA FILE 5 MINUTE AND 1 DEGREE
!                               RESOLUTION)
!
!            arpstern.input    ARPSTERN11.F input file
!
!
!
!
!  OUTPUT:
!            dir1deg.dat       direct access 1 degree unformatted
!                              data file with xxxx seperate records.
!
!            dir1deg.hdr       header file which has record numbers
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR:  Dan Weber
!           1/12/94
!
!  Modification history:
!
!    Dan Weber 7/06/94.
!       Documentation added.
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
                              !  and latitudinal direction
  INTEGER :: jlim             !  Limit of the number of reads in the
                              !  elev.dat data file
  INTEGER :: klim             !  Number of reads at a given
                              !  longitude (determined by elev.dat = 6)
  INTEGER :: ilim             !  Number of 80 character lines per record
                              !  in elev.dat (=20)
  INTEGER :: ilon             !  Maximum longitude

  PARAMETER (ilim=20,ilon=360,jlim=100000,klim=6,ndx=2)

  INTEGER :: rlat(ilim)       !  Latitude read in from the ascii data file
  INTEGER :: rlon(ilim)       !  Longitude read in from the ascii data file

  INTEGER :: flat(jlim)       !  Latitude corresponding a record number
  INTEGER :: flon(jlim)       !  Longitude corresponding a record number


  INTEGER :: qsize(ilim)      !  Quad size of the block (1= 1 deg by 1 deg and
                              !  30' by 30', 5 = 5' by 5')
  INTEGER :: nlat(ilon)       !  Search Latitude
  INTEGER :: nlon(ilon)       !  Search Longitude

!
!  LOCAL VARIABLES:
!

  INTEGER :: irec,i,ii,j,k,kk,l,num,iloc
  CHARACTER (LEN=80  ) :: tdatadir  !  Directory in which the data file
                                    !  dma_elev.dat is stored.
  CHARACTER (LEN=80   ) :: terndir  !  Directory in which the original ASCII terrain
                                    !  data files are stored.
  INTEGER :: lterndir         !  Length of the non-blank part of string terndir.

  INTEGER :: analtype         !  Dummy variable used for namelist data
  INTEGER :: mapproj          !  Dummy variable used for namelist data
  INTEGER :: itertype         !  Dummy variable used for namelist data
  INTEGER :: rmsopt           !  Dummy variable used for namelist data
  INTEGER :: comtype          !  Computer type for program to run on...
                              !  for IBM, comtype = 1
                              !     CRAY,         = 4
  CHARACTER (LEN=1) :: dum1         !  Dummy variable used to read ascii header
  CHARACTER (LEN=2) :: dum2         !  Dummy variable used to read ascii header
  CHARACTER (LEN=8) :: dum3         !  Dummy variable used to read ascii header
  CHARACTER (LEN=9) :: dum4         !  Dummy variable used to read ascii header
  INTEGER :: units(ilim)      !  Units of measure for the terrain data
                              !   units = 1  meters
                              !         = 2  feet
                              !         = 3  fathoms


  INTEGER*2 TYPE(klim,ilim) ! Type of area the block defines  See
                            !  appendix fro definition of type.
  INTEGER :: ielev(klim,ilim) !  Intermediate terrain block data put in
                              !  terms of increasing record number order.
  INTEGER*2 itp(4,jlim)    !  Intermediate terrain values in record
                           !  number order
!
!  OUTPUT:
!
  INTEGER*2 iw(ndx,ndx)    !  Terrain data written in unformatted
                           !  form for entire 1 degree by 1 degree
                           !  block to dir1deg.f file.


!
!-----------------------------------------------------------------------
!
!  Defining namelist:
!
!-----------------------------------------------------------------------

  NAMELIST /terraind/analtype,mapproj,itertype,rmsopt,comtype,          &
                     tdatadir,terndir

!-----------------------------------------------------------------------
!
!  Begin executable portion of the code.
!
!-----------------------------------------------------------------------

!
!  Read the arpstern.input file for the terraind namelist information.
!  (note: only need comtype for this program to run properly.)

  READ(5,terraind)

  lterndir = LEN(terndir)
  CALL strlnth( terndir, lterndir )

  IF( lterndir == 0 ) THEN
    terndir = '.'
    lterndir=1
  END IF

!
!  open elev.dat, dir1deg.dat and dir1deg.hdr files...
!

  OPEN(10,FILE=terndir(1:lterndir)//'/elev.dat',STATUS='old')

  OPEN(11,FILE='dir1deg.dat',ACCESS='direct',FORM='unformatted',        &
      STATUS='unknown',RECL=8*comtype)
  OPEN(12,FILE='dir1deg.hdr',FORM='formatted',STATUS='unknown')

!
!  Initializing counters....
!

  irec=0
  num=1
  nlat(1)=89    ! latitude of the southwestern part of the northern
                ! most block of data.
  nlon(1)=0

!
!  Starting the main loop which reads the elev.dat data set
!  and check to see if the data is 5 minute or 1 degree (or
!  30 minute) resolution.
!

  DO j=1,jlim

!
!    Reading 1 record at a time. (note: each record consists of
!    20 lines with 80 charactes...)
!

    READ(10,'(20(A2,I5,A1,I5,I1,6(i1,a1,i6),A9,I1,a8))',END=99)         &
        (dum2,rlat(i),dum1,rlon(i),units(i),TYPE(1,i),dum1,ielev(1,i)   &
        ,TYPE(2,i),dum1,ielev(2,i),TYPE(3,i),dum1,ielev(3,i),TYPE(4,i)  &
        ,dum1,ielev(4,i),TYPE(5,i),dum1,ielev(5,i),TYPE(6,i),dum1       &
        ,ielev(6,i),dum4,qsize(i),dum3,i=1,ilim)
!  print *,'ok after read',j
!
!  Test for 1 degree block type
!

    IF(qsize(1) == 1)THEN ! qsize=1 for 1 degree, =5 for 5 minute data

!
!  write 1 degree header to the dir1deg.hdr file and the 1 degree data
!  to the dir1deg.dat file.  Note: the latitude is for the nortwest
!  corner of the block, this is inconsistant with the dir1km.dat storage
!  format.  Thus, 1 will be subtracted from the latitude to provide
!  consistancy. The first record will be 89 south and 0 east. the file
!  will store 89 south and 1 east and so on until 89 south and 359 east
!  is written.  Then the next latitude to the north will be written
!  and so on, all the way up to the north pole.
!
      DO i=1,ilim

        IF(nlat(num) == (rlat(i)/100-1))THEN !testing if rlat is same as
          nlon(num)=1+nlon(num)              ! previous one
        ELSE
          num=num+1
          nlat(num)=rlat(i)/100-1
          nlon(num)=1
        END IF

        irec=irec+1
        flat(irec)=rlat(i)/100-1
        flon(irec)=rlon(i)/100
        IF(TYPE(1,i) == 9)THEN !quad is all ocean and the elev is set to zero

          DO k=2,5
            ielev(k,i)=200       ! and all components of ielev
          END DO

        ELSE IF(TYPE(1,i) == 0.OR.TYPE(1,i) == 4)THEN !quad contains land/ice

          IF(ielev(2,i) == 0.AND.ielev(3,i) == 0.AND.ielev(4,i) == 0    &
                .AND.ielev(5,i) == 0.AND.ielev(1,i) > 0)THEN ! packing data
            DO k=2,5
              ielev(k,i)=(ielev(1,i)+4000)/20
            END DO
          ELSE                ! setting elevation to sea level (0.0 meters)
            DO k=2,5
              IF(ielev(k,i) < 0)THEN
                ielev(k,i)=200
              ELSE
                ielev(k,i)=(ielev(k,i)+4000)/20   ! putting in packed form
              END IF
            END DO
          END IF

        ELSE IF(TYPE(1,i) == 8)THEN ! quad contains some negative land

          DO k=2,5
            ielev(k,i)=(ielev(k,i)+4000)/20     ! putting in packed form
          END DO

        ELSE IF(TYPE(1,i) == 5.OR.TYPE(1,i) == 3.OR.TYPE(1,i) == 2)THEN
                                                ! land/ocean mix and
! ice/ocean mix.

          IF(ielev(2,i) == 0.AND.ielev(3,i) == 0.AND.ielev(4,i) == 0.AND. &
                 ielev(5,i) == 0.AND.ielev(1,i) > 0)THEN
            DO k=2,5
              ielev(k,i)=(ielev(1,i)+4000)/20 ! packing data...
            END DO
          ELSE
            DO k=2,5
              IF(ielev(k,i) < 0)ielev(k,i)=0
              ielev(k,i)=(ielev(k,i)+4000)/20   ! putting in packed form
            END DO
          END IF

        ELSE                               ! type(k,i) out of accepted values
          PRINT *,'type of elevation unknown'
          STOP

        END IF

        DO k=2,5
          itp(k-1,irec)=ielev(k,i)
        END DO

      END DO

    END IF
  END DO

  99   PRINT *,j,jlim,irec    ! print the number of records
  iloc=0
  kk=0
  DO i=num,1,-1
    iloc=iloc+nlon(i)
    PRINT *,'latitude, number of blocks found of the latitude  circle'  &
        ,nlat(i),nlon(i),iloc
    DO j=1,nlon(i)
      kk=kk+1
      ii=irec-iloc+j

!
!   putting lat in the same format as the 30 second file.
!

      iw(1,1)=itp(3,ii)
      iw(2,1)=itp(4,ii)
      iw(1,2)=itp(1,ii)
      iw(2,2)=itp(2,ii)
      WRITE(12,'(i6,i6,i6)') flat(ii),flon(ii),kk  ! writing data to
                                                   ! the dir1deg.hdr header file
      WRITE(11,REC=kk) iw                !  writing terrain data to the
                                         !  dir1deg.dat data file
      DO l=1,2
        DO k=1,2
          IF(iw(k,l) < 200)PRINT *,flat(ii),flon(ii),kk ! printing negative
                                                        !  value locations....
        END DO
      END DO
    END DO
  END DO



  STOP
END PROGRAM dir1deg
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
