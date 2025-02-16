PROGRAM dir5min

!##################################################################
!##################################################################
!######                                                      ######
!######                 PROGRAM DIR5MIN                      ######
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
!  A stand alone program which rewrites 5 minute terrain elevation data
!  from elev.dat recieved from NCAR DATA SERVICES into a direct access
!  file using either 2 of 8 byte integer format depending upon the type
!  of machine used.  The purpose of writing a direct access data format
!  is to expedite the run time of the ARPSTERN11.F terrain processing
!  program for initializing the ARPS terrain field.
!  In addition, a dir5min.hdr file (header file) will be created and
!  used by arpstern11.f.
!
!  INPUT:
!            elev.dat          (NCAR DATA FILE 5 MINUTE RESOLUTION)
!
!            arpstern.input    ARPSTERN11.F input file
!
!
!
!
!  OUTPUT:
!            dir5min.dat       direct access 5 minute unformatted
!                              data file with 4118 seperate records.
!
!            dir5min.hdr       header file which has record numbers
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
  INTEGER :: klim             !  Number of 5 minute reads at a given
                              !  longitude (determined by elev.dat = 6)
  INTEGER :: ilim             !  Number of 80 character lines per record
                              !  in elev.dat (=20)
  INTEGER :: ilon             !  Maximum longitude

  PARAMETER (ilim=20,ilon=360,jlim=100000,klim=6,ndx=12)

  INTEGER :: rlat(ilim)       !  Latitude read in from the ascii data file
  INTEGER :: rlon(ilim)       !  Longitude read in from the ascii data file

  INTEGER :: flat(jlim)       !  Latitude corresponding to the record number
  INTEGER :: flon(jlim)       !  Longitude corresponding to the record number


  INTEGER :: qsize(ilim)      !  Quad size of the block (1= 1 deg by 1 deg and
                              !  30' by 30', 5 = 5' by 5')
  INTEGER :: nlat(ilon)       !  Latitude
  INTEGER :: nlon(ilon)       !  Longitude

  INTEGER :: TYPE(klim,ilim)  ! Type of area the block defines.  See
                              !  appendix for definition of type.

  INTEGER :: relev(klim,ilim) ! Elevation value read in from the elev.dat
                              !  data file

  INTEGER :: ielev(jlim,ndx,ndx) ! Intermediate terrain block data put
                                 !  in terms of increasing record number order

!
!  LOCAL VARIABLES:
!

  INTEGER :: irec,i,ii,j,k,kk,l,num,iloc
  INTEGER :: item1,item2,item3
  INTEGER :: comtype          !  Computer type for program to run on...
                              !  for IBM, comtype = 1
                              !     CRAY,         = 4
  CHARACTER (LEN=80  ) :: tdatadir  !  Directory in which the data file
                                    !  dma_elev.dat is stored.
  CHARACTER (LEN=80   ) :: terndir  !  Directory in which the original ASCII terrain
                                    !  data files are stored.
  INTEGER :: lterndir         !  Length of the non-blank part of string terndir
  INTEGER :: analtype         !  Dummy variable used for namelist data
  INTEGER :: mapproj          !  Dummy variable used for namelist data
  INTEGER :: itertype         !  Dummy variable used for namelist data
  INTEGER :: rmsopt           !  Dummy variable used for namelist data
  CHARACTER (LEN=1) :: dum1         !  Dummy variable used to read ascii header
  CHARACTER (LEN=2) :: dum2         !  Dummy variable used to read ascii header
  CHARACTER (LEN=8) :: dum3         !  Dummy variable used to read ascii header
  CHARACTER (LEN=9) :: dum4         !  Dummy variable used to read ascii header
  INTEGER :: units(ilim)      !  Units of measure for the terrain data
                              !   units = 1  meters
                              !         = 2  feet
                              !         = 3  fathoms
!
!  OUTPUT:
!
  INTEGER*2 iw(ndx,ndx)    !  Terrain data written in unformatted
                           !  form for entire 1 degree by 1 degree
                           !  block to dir5min.f file.

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
!  open elev.dat, dir5min.dat and dir5min.hdr files...
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
  OPEN(10,FILE=terndir(1:lterndir)//'/elev.dat',STATUS='old')

  OPEN(11,FILE='dir5min.dat',ACCESS='direct',FORM='unformatted',        &
      STATUS='unknown',RECL=288*comtype)
  OPEN(12,FILE='dir5min.hdr',FORM='formatted',STATUS='unknown')

!
!  Initializing counters.....
!

  irec=0            ! number of 5 minute data records read
  num=1
  nlat(num)=74      ! no 5 minute data above this latitude (deg north)
  nlon(num)=0       ! starting longitude for each latitude
  item3=0

!
!    Starting the main loop.  Note that the entire data file is read
!    and the data stored and then the data is written to file outside
!    this loop....
!

  DO j=1,jlim
!
!    Reading 1 record at a time. (note: each record consists of
!    20 lines with 80 charactes...)
!

    READ(10,'(20(A2,I5,A1,I5,I1,6(i1,a1,i6),A9,I1,a8))',END=99)         &
        (dum2,rlat(i),dum1,rlon(i),units(i),TYPE(1,i),dum1,relev(1,i)   &
        ,TYPE(2,i),dum1,relev(2,i),TYPE(3,i),dum1,relev(3,i),TYPE(4,i)  &
        ,dum1,relev(4,i),TYPE(5,i),dum1,relev(5,i),TYPE(6,i),dum1       &
        ,relev(6,i),dum4,qsize(i),dum3,i=1,ilim)

!
!  test to see if the block type is 5 minute
!

    IF(qsize(1) /= 1)THEN !qsize =1 for 1 degree, =5 for 5 minute data

!
!  Test to find the character of the data (ocean, land...)
!  and set ocean values to zero elevation.
!

      DO k=1,klim
        DO i=1,ilim
          IF(TYPE(k,i) == 9)THEN ! quad is ocean and the elev is set to zero.
            relev(k,i)=200       ! 200(packed) = 0.0 meters unpacked

          ELSE IF(TYPE(k,i) == 0.OR.TYPE(k,i) == 4)THEN !quad contains land/ice
            relev(k,i)=(relev(k,i)+4000)/20          ! putting in packed form

          ELSE IF(TYPE(k,i) == 8)THEN       ! quad contains some negative land
            relev(k,i)=(relev(k,i)+4000)/20          ! putting in packed form

          ELSE IF(TYPE(k,i) == 5.OR.TYPE(k,i) == 3.OR.                  &
                    TYPE(k,i) == 2)THEN                  !land/ocean mix
            IF(relev(k,i) < 0)THEN
              relev(k,i)=200                         ! ice/ocean mix
            ELSE
              relev(k,i)=(relev(k,i)+4000)/20        ! putting in packed form
            END IF
          ELSE                               ! type(k,i) does not match
                                             ! defined list

            PRINT *,'type of elevation unknown,irec=',irec
            STOP
          END IF
        END DO
      END DO

!
!  Sort the data into southwest to northeast order from northwest to
!  southeast order.
!
!  Note the data is read from the northwest corner of the 1 degree
!  blocks.  This data will be sorted and stored with the (1,1) element
!  as the southwest corner of the block similar to the 30 second data
!  set format.
!

      DO i=1,ilim
        IF(item3 /= rlat(i) ) THEN

        IF( rlat(i) == 99999 ) CYCLE
 
          item1=MOD(rlat(i),100)
          item2=MOD(rlon(i),100)
          IF(item1 == 0.AND.item2 == 0)THEN    ! starting a new block...
            kk=ndx
            irec=irec+1
            flat(irec)=(rlat(i)-100)/100      ! adjusting the latitude...
            flon(irec)=rlon(i)/100            ! setting the longitude...
            DO k=1,ndx/2
              ielev(irec,k,kk)=relev(k,i)
            END DO

            IF(nlat(num) == flat(irec))THEN   ! counting the number of blocks
              nlon(num)=nlon(num)+1           ! for a latitude circle...
            ELSE
              num=num+1
              nlat(num)=flat(irec)
              nlon(num)=1
            END IF
          ELSE IF(item1 == 0.AND.item2 /= 0)THEN ! setting an existing block..
            DO k=ndx/2+1,ndx
              kk=ndx
              ielev(irec,k,kk)=relev(k-6,i)
            END DO
          ELSE IF(item1 /= 0.AND.item2 == 0)THEN ! setting an existing block..
            kk=kk-1
            DO k=1,ndx/2
              ielev(irec,k,kk)=relev(k,i)
            END DO
          ELSE IF(item1 /= 0.AND.item2 /= 0)THEN ! setting an existing block..
            kk=kk-1
            DO k=ndx/2+1,ndx
              ielev(irec,k,kk)=relev(k-6,i)
            END DO
          ELSE
            PRINT *,'lat,lon does not follow format',rlat(i),rlon(i)
          END IF
          item3=rlat(i)
        ELSE
          PRINT *,'duplicate found at',rlat(i),rlon(i)
        END IF
      END DO ! i loop

    END IF

  END DO  ! j loop 

!
!  End of read/sort loop (main loop)
!

!
!  write 5 minute header to the dir5min.hdr file and the 5 minute data
!  to the dir5min.dat file.  Note: the latitude is for the northwest
!  corner of the block, this is inconsistant with the dir30sec.dat storage
!  format.  The 5 minute data will be written following the 30 second
!  format with the southwestern most data block written first and then
!  the next block to the east and so on.  The next latitiude to the north
!  will then be written.
!

  99   PRINT *,'print main loop index, and irec',jlim,irec
  iloc=0
  kk=0

  DO i=num,1,-1
    iloc=iloc+nlon(i)
    PRINT *,'latitude, number of blocks on the latitude circle, ',      &
            'running sum of block written', nlat(i),nlon(i),iloc
    DO j=1,nlon(i)
      kk=kk+1
      ii=irec-iloc+j
      DO l=1,ndx
        DO k=1,ndx
          iw(k,l)=ielev(ii,k,l)
        END DO
      END DO

!
!   Writing the block data to the unformatted file and the header
!   number corresponding to the block location to the header file.
!

      WRITE(12,'(i6,i6,i6)') flat(ii),flon(ii),kk
      WRITE(11,REC=kk) iw

!
!  Testing for negative values...
!

      DO l=1,ndx
        DO k=1,ndx
          IF(iw(k,l) < 200)PRINT *,'negative elev.',flat(ii),flon(ii),kk
        END DO
      END DO

    END DO
  END DO

!
!  Closing all files and exiting the program....
!

  CLOSE(10)
  CLOSE(11)
  CLOSE(12)

  STOP
END PROGRAM dir5min
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
