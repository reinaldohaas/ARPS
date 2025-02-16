!########################################################################
!########################################################################
!#########                                                      #########
!#########              SUBROUTINE get_radial_data              #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE get_radial_data(iheaderdim,itrailerdim,msglen,ival,          &
               isite, isite_lat,isite_lon,isite_elev,                   &
               iyr,imonth,iday,ihr,iminit,isec,                         &
               iprod,maxval1,ivcp,ifield,maxbins,maxradials,            &
               iazmuth,ielev,nbins,nradials,icatt,icats,                &
               maxval2,iyr1,imonth1,iday1,ihr1,iminit1,                 &
               iyr2,imonth2,iday2,ihr2,iminit2,istm_spd,istm_dir,       &
               iheader,len_header,itrailer,len_trailer,icode)

!------------------------------------------------------------------------
!
! PURPOSE:
!
! Decodes a WSR-88D radial graphic product and return pertinent data.
!
! Note 1:  To convert ifield() to physical units, use the transform:
!
!      VALUE = ICATS(IFIELD())
!
! Note 2:  Non-system subroutines used are:
!         w3fs26
!         i2byte
!         i4byte
!
!------------------------------------------------------------------------
!
! AUTHOR: 
!
! NWS Techniques Development Laboratory (April 1998).
!
! MODIFICATION HISTORY:
!
! Yunheng Wang  (July 2001)
! Converted from FORTRAN 77 code.
!
! Eric Kemp, 30 July 2001
! Added USE statement to module n2aconst, minor changes.
! 
! Eric Kemp, 9 August 2001
! Use variable type INTEGER*4 and similar types instead of using
! KIND statements.  Also use dimension parameters from N2ACONST module.
!
! Eric Kemp, 17 August 2001
! Modified code to avoid use of INTEGER*4, INTEGER*1, etc, based on 
! work by Yunheng Wang.  Also, changed variable maxval to maxval1 to
! avoid conflict with Fortran 90 implicit function.  Finally, added
! code to allocate array irad at runtime, with dimension added to
! module N2ACONST.
!
! Keith Brewster, 6 April 2010
! Modifications to accomodate conversion of n2aconst to nidscst.inc
!
!------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------
!
! Declare external functions.
!
!------------------------------------------------------------------------

  INTEGER, EXTERNAL :: i2byte, i4byte

!------------------------------------------------------------------------
!
! Arguments
!
!------------------------------------------------------------------------
!    INPUT:
!                ival = Array containing the graphic product
!
!------------------------------------------------------------------------
!    OUTPUT:
!              msglen = length of input product, bytes
!               isite = wsr-88d site id
! isite_lat,isite_lon = site lat/lon (deg e), deg x 1000
!          isite_elev = site elevation, m msl
!     iyr,imonth,iday = product volume scan date: year, month, day
!     ihr,iminit,isec = product time, utc: hour, minute, second
!               iprod = product code number (output) (eg 19 for z)
!             maxval1 = maximum value in the field,
!                       for velocity products, the greatest inbound
!                       (negative) velocity.
!             maxval2 = for velocity products, the greatest outbound
!                       (positive) velocity
!   istm_spd,istm_dir = for storm-relative velocity products, the velocity
!                       subtracted from base velocity (deg and kt)
!  iyr1,imonth1,iday1 = for variable-period precipitation accumulation
!                       products, starting date of accumulation pd, UTC
!        ihr1,iminit1 = for variable-period precipitation accumulation
!                       products, starting time of accumulation pd, UTC
!  iyr2,imonth2,iday2 = for variable-period precipitation accumulation
!                       products, ending date of accumulation pd, UTC
!        ihr2,iminit2 = for variable-period precipitation accumulation
!                       products, ending time of accumulation pd, UTC
!                ivcp = operating volume coverage pattern
!           ifield(,) = returned values (data levels) for the field
!                       a byte variable, 0-15. dimensioned by maxbins,
!                       maxradials in calling routine.
!
!                       for reflectivity and rainfall:
!                       ifield(,)=0 means 'below minimum threshold',
!                       and the actual value is meaningless. other
!                       values of icats() may be 0 or negative.
!    
!                       for velocity and spectrum width:
!                       category 0 means insufficient backscatter and
!                       that velocity and width cannot be determined.
!                       category 15 indicates range folding (velocity may
!                       be from either primary or second-trip echo) and
!                       again, velocity and width cannot be determined.
!
!           iazmuth() = antenna azimuth directions (deg x 10) for
!                       each radial.  dimensioned by maxradials
!                       in calling routine.
!               ielev = antenna elevation angle, deg x 10
!      nbins,nradials = number of range bins and radials in the graphic
!             icatt() = category-interval values taken directly from
!                       product header.  dimensioned by 16 or 0:15
!                       in calling program.
!             icats() = category-interval values in physical units.
!                       dimensioned by 16 or 0:15 in calling program 
!           iheader() = byte array with copy of header.  should be
!                       dimensioned 160 in calling routine
!          len_header = actual number of bytes of information in iheader()
!
!          itrailer() = character*1 array with text information attached
!                       at end of file (if any). should be dimensioned by
!                       20 000 in calling routine
!         len_trailer = actual number of bytes of data in itrailer().
!                       will be set 0 if no trailer exists
!
!               icode = output return code, conditions as follows:
!                       0 : data returned ok
!                       2 : not a radial graphic product
!                       3 : maxradials or maxbins is too small
!                       4 : product length greater than msglen
!
!------------------------------------------------------------------------

  INTEGER, INTENT(IN) :: iheaderdim,itrailerdim
  INTEGER, INTENT(IN) :: msglen, maxbins, maxradials

  INTEGER, INTENT(IN) :: ival(msglen) 

  INTEGER, INTENT(OUT) :: icode,ivcp,isite,isite_lat,isite_lon,          &
                            isite_elev,iyr,imonth,iday,ihr,iminit,isec,  &
                            iprod,nbins,nradials,maxval2,ihr1,iminit1,   &
                            ihr2,iminit2,istm_spd,istm_dir,len_header,   &
                            len_trailer

  INTEGER, INTENT(OUT) :: iazmuth(maxradials)

  INTEGER, INTENT(OUT) :: icatt(0:15)
  INTEGER, INTENT(OUT) :: icats(0:15)

  INTEGER, INTENT(OUT) :: iyr1,imonth1,iday1,iyr2,imonth2,iday2

  INTEGER, INTENT(OUT) :: maxval1,  ielev
  
  INTEGER, INTENT(OUT) :: iheader(iheaderdim),                           &
                          ifield(maxbins,maxradials)

  CHARACTER(LEN=1), INTENT(OUT) :: itrailer(itrailerdim)

!------------------------------------------------------------------------
!
!  Misc. local variables
!
!------------------------------------------------------------------------
!  LOCAL VARIABLES:
!
!           irad(600) = byte vector hold one radial of data
!         ival0,ival1 = integer*4 words used to hold lower and higher-order
!                       halves of 2-byte integer
!    nex_date,jul_date = nexrad reference date, astronomical julian date, 
!                       used to determine legal date (year-month-day)
!         npos1,npos2 = positions of 1st and last bytes (relative to start
!                       of the file) in the current radial being decoded
!               nruns = number of runs needed to store the present radial
!
!------------------------------------------------------------------------

  INTEGER, ALLOCATABLE :: irad(:)
  INTEGER :: ival0, ival1

  INTEGER :: ipcode

  INTEGER :: nex_date, jul_date, idaywk,idayyr, nex_time
  INTEGER :: npos1, npos2, nruns
  
  INTEGER :: n, k, len1, ibin, ibin1, ibin2, nbin0
  INTEGER :: ilvl, nbyte, ival4, nrun_marker, nr, nc
!
!------------------------------------------------------------------------
!
! Include file
!
!------------------------------------------------------------------------

  INCLUDE 'nidscst.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      
  ALLOCATE(irad(iraddim))

!------------------------------------------------------------------------
!
! Get product type, and check for radial-type packet code
! in bytes 137-138 
!
!------------------------------------------------------------------------

 
  iprod = i2byte(ival(1) ,ival(2))
  ipcode = i2byte(ival(137),ival(138))

  print *, '  iprod = ',iprod
  print *, '  ipcode = ',ipcode

! IF (ipcode /= ipradial1 .AND. ipcode /= ipradial2 .AND. ipcode /= ipradial3 .AND. &
!     ipcode /= ipradial4 .AND. ipcode /= ipradial5 ) THEN

  IF (ipcode /= ipradial1 ) THEN
    icode = 2
    RETURN
  END IF

  isite = 256*ival(13) + ival(14)
  ivcp  = 256*ival(35) + ival(36)

  print *, ' isite =',isite

!------------------------------------------------------------------------
!
! Check internal file-length specification.  Note that if file
! length is shorter than message length we won't proceed.
!
!------------------------------------------------------------------------

  write(6,*) 'ival = ',ival(9),ival(10),ival(11),ival(12)
  len1 = i4byte(ival(9),ival(10),ival(11),ival(12))

  write(0,*) 'msglen = ',msglen, len1
  IF (len1 > msglen) THEN
    icode = 4
    RETURN
  END IF

!------------------------------------------------------------------------
!
! Get header and trailer lengths, initialize arrays.
!
!------------------------------------------------------------------------

  len_header = 150
  DO n=1,150
    iheader(n) = ival(n)
  END DO
  
  len1 = 136 + i4byte(ival(133),ival(134),ival(135),ival(136))
  len_trailer = msglen - len1

  IF (len_trailer > 0) THEN
    k = 0
    DO n=len1+1, msglen
      k = k + 1
      itrailer(k) = ACHAR(ival(n))
    END DO
  END IF
  
  IF (len_trailer < 0) len_trailer = 0

!------------------------------------------------------------------------
!
! Get volume scan date and time
!
!------------------------------------------------------------------------
   
  nex_time = i4byte(ival(43),ival(44),ival(45),ival(46))
  nex_date = i2byte(ival(41),ival(42))
  jul_date = jul_1_1_70 + nex_date
  CALL w3fs26(jul_date,iyr,imonth,iday,idaywk,idayyr)
   
  ihr = nex_time / 3600
  iminit = (nex_time - (3600*ihr)) / 60
  isec = nex_time - (3600*ihr) - (60*iminit)

!------------------------------------------------------------------------
!
! Get site lat/lon, elevation
!
!------------------------------------------------------------------------

  isite_elev = i2byte(ival(29),ival(30))
  isite_lat = i4byte(ival(21),ival(22),ival(23),ival(24))
  isite_lon = i4byte(ival(25),ival(26),ival(27),ival(28))

!------------------------------------------------------------------------
!
! Find number of radials in the product, elevation angle
!
!------------------------------------------------------------------------
 
  nbins = i2byte(ival(141),ival(142))
  nradials = i2byte(ival(149),ival(150))

  IF ((maxbins < nbins) .OR. (maxradials < nradials)) THEN
    print *, ' nbins = ',nbins
    print *, ' nradials = ',nradials
    icode = 3
    RETURN
  END IF

  ielev = i2byte(ival(59),ival(60))
   
!------------------------------------------------------------------------
!
! Assign category threshold values.  Scale up values
! if indicated by comparison of icatt(0) to iscale5.
!
!------------------------------------------------------------------------

  nc = 0
   
  DO n=61,91,2
    icatt(nc) = i2byte(ival(n),ival(n+1))
    ival0 = ival(n+1)
    IF (ival0 < 0) ival0 = ival0 + 256
    ival1 = ival(n)
    IF (MOD(ival1,2) == 0) THEN
      icats(nc) = ival0
    ELSE
      icats(nc) = -1*ival0
    END IF
    nc = nc + 1
  END DO

  IF (icatt(0) >= iscale5) THEN
    DO nc=1,15
      icats(nc) = 5*icats(nc)
    END DO
  END IF

!------------------------------------------------------------------------
!
! Adjust '0' position category back to 0 ('background color') 
!
!------------------------------------------------------------------------
  
  icats(0) = 0

!------------------------------------------------------------------------
!
! Get maximum value(s).  For velocity, maxval1 and maxval2
! are inbound and outbound, respectively.
!
!------------------------------------------------------------------------

  ival0 = ival(94)
  IF (ival0 < 0) ival0 = 256 + ival0
  ival1 = ival(93)
  IF (ival1 < 0) ival1 = 256 + ival1
  maxval1 = i2byte(ival(93),ival(94))
  maxval2 = i2byte(ival(95),ival(96))

!------------------------------------------------------------------------
!
! Get storm velocity for storm-relative products
!
!------------------------------------------------------------------------

  istm_spd = i2byte(ival(101),ival(102))
  istm_dir = i2byte(ival(103),ival(104)) 

!------------------------------------------------------------------------
!
! For user-selected or storm-total precip, derive
! starting and ending date-times for accumulation
!
!------------------------------------------------------------------------

  IF ((iprod == 80) .OR. (iprod == 78) .OR. (iprod == 79)                &
       .OR. (iprod == 81)) THEN
       
    nex_date = i2byte(ival(95),ival(96))
    jul_date = jul_1_1_70 + nex_date
    CALL w3fs26(jul_date,iyr1,imonth1,iday1,idaywk,idayyr)
    nex_time = i2byte(ival(97),ival(98))
    ihr1 = nex_time / 60
    iminit1 = nex_time - (60*ihr1)
    nex_date = i2byte(ival(99),ival(100))
    jul_date = jul_1_1_70 + nex_date
    CALL w3fs26(jul_date,iyr2,imonth2,iday2,idaywk,idayyr)
    nex_time = i2byte(ival(101),ival(102))
    ihr2 = nex_time / 60
    iminit2 = nex_time - (60*ihr2)
  
  END IF

!------------------------------------------------------------------------
!
! Initialize the output array with -127 values.
!
!------------------------------------------------------------------------

  ifield(:,:) = -127

!------------------------------------------------------------------------
!
! Radial by radial, plug run-length encoded values into the
! output array.  For each radial, file has number of 2-byte
! words describing the radial, not bytes as for raster products.
!
!------------------------------------------------------------------------

  nr = 0
  nrun_marker = 151

  DO WHILE (nr < nradials)
    nr = nr + 1
    IF (nr > nradials) EXIT

    ival0 = ival(nrun_marker+1)
    ival1 = ival(nrun_marker)
    IF (ival0 < 0) ival0 = 256 + ival0
    IF (ival1 < 0) ival1 = 256 + ival0
    nruns = 2*i2byte(ival(nrun_marker),ival(nrun_marker+1))
    iazmuth(nr) = i2byte(ival(nrun_marker+2),ival(nrun_marker+3))
    npos1 = nrun_marker + 6
    npos2 = npos1 + nruns - 1
    ibin1 = 1

!------------------------------------------------------------------------
!    
!   Extract run length from 4 high-order bits and run value
!   from 4 low-order bits of this byte, then plug into radial
!   array.
!
!------------------------------------------------------------------------

    DO nbyte=npos1,npos2
      ival4 = ival(nbyte)
      IF (ival4 < 0) ival4 = 256 + ival4 
      nbin0 = ival4/16
      IF (nbin0 <= 0) EXIT
      ilvl = ival4 - (16*nbin0)
      ibin2 = ibin1 + nbin0 - 1
      IF (ibin2 > nbins) ibin2 = nbins
      DO ibin=ibin1,ibin2
        irad(ibin) = ilvl
      END DO
      ibin1 = ibin2 + 1
      IF (ibin2 >= nbins) EXIT
    END DO

!------------------------------------------------------------------------
!
!   Present radial has been decoded.  Now map it into the 
!   full array.  
!
!------------------------------------------------------------------------

    DO ibin=1,nbins
      ifield(ibin,nr) = irad(ibin)
    END DO        

!------------------------------------------------------------------------
!  
!   Finished mapping this radial - go on to the next
!
!------------------------------------------------------------------------

    nrun_marker = npos2 + 1
  END DO
      
  icode = 0
  RETURN
END SUBROUTINE get_radial_data

!########################################################################
!########################################################################
!#########                                                      #########
!#########               SUBROUTINE get_raster_data             #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE get_raster_data(iheaderdim,itrailerdim,msglen,ival,          &
                           isite,isite_lat,isite_lon,isite_elev,        &
                           iyr,imonth,iday,ihr,iminit,isec,             &
                           iprod,maxval1,ivcp,ifield,ncols,nrows,icatt, &
                           icats,ndims_p,iheader,len_header,itrailer,   &
                           len_trailer,icode)

!------------------------------------------------------------------------
!
! PURPOSE:
!
! Decodes WSR-88D raster product and return all pertinent data.
!
! Note 1:  To convert ifield() values to physical units, use the icats()
! array as follows:
!           VALUE = ICATS(IFIELD())
!
! Note 2:  This uses non-system subroutines:
!            W3FS26
!            I2BYTE 
!            I4BYTE
!
!------------------------------------------------------------------------
!
! AUTHOR:
!
! Eric Kemp, 30 July 2001.
! Based on HPUX FORTRAN 77 subroutine written by Kitzmiller of the
! NWS Techniques Development Laboratory (August 1997).
!
! MODIFICATION HISTORY:
!
! Eric Kemp, 9 August 2001
! Use variable type INTEGER*4 and similar types instead of using
! KIND statements.  Also use dimension parameters from N2ACONST module.
!
! Eric Kemp, 17 August 2001
! Modified code to avoid use of INTEGER*4, INTEGER*1, etc., based on
! work by Yunheng Wang.  Also, changed variable maxval to maxval1 to
! avoid conflict with Fortran 90 implicit function.
!    
! Keith Brewster, 6 April 2010
! Modifications to accomodate conversion of n2aconst to nidscst.inc
!
!------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------
!
! External functions
!
!------------------------------------------------------------------------

  INTEGER, EXTERNAL :: i2byte, i4byte

!------------------------------------------------------------------------
!   
! Declare arguments (original documentation):
!
!             IVAL( ) = BYTE OR CHARACTER*1 ARRAY CONTAINING GRAPHIC
!                       MESSAGE, DIMENSIONED MSGLEN (INPUT)
!              MSGLEN = LENGTH OF GRAPHIC MESSAGE IN BYTES (INPUT)
!         NCOLS,NROWS = DIMENSIONS OF IFIELD() (INPUT)
!               ICODE = RETURN CODE, CONDITIONS AS FOLLOWS:
!                       0 : RETURN OK
!                       2 : NOT A RECOGNIZED RASTER PRODUCT
!                       3 : NCOLS/NROWS TOO SMALL FOR PRODUCT
!                       4 : MSGLEN TOO SMALL TO STORE PRODUCT!
!               ISITE = WSR-88D NUMERIC SITE ID (OUTPUT)
! ISITE_LAT,ISITE_LON = SITE LAT/LON (DEG E), DEG X 1000 (OUTPUT)
!          ISITE_ELEV = SITE ELEVATION, M MSL (OUTPUT)
!     IYR,IMONTH,IDAY = PRODUCT DATE (OUTPUT)
!     IHR,IMINIT,ISEC = PRODUCT TIME, UTC (OUTPUT)
!               IPROD = PRODUCT CODE NUMBER (OUTPUT) (EG 57 FOR VIL)
!             MAXVAL1 = MAXIMUM VALUE IN THE RASTER FIELD (OUTPUT)
!                IVCP = OPERATING VOLUME COVERAGE PATTERN (OUTPUT)
! IFIELD(NCOLS,NROWS) = RETURNED VALUES (DATA LEVELS) FOR THE
!                       FIELD (OUTPUT); A BYTE VARIABLE
!             ICATT() = CATEGORY-INTERVAL VALUES DIRECT FROM
!                       PRODUCT HEADER (OUTPUT, INTEGER*2).  SHOULD BE
!                       DIMENSIONED BY 16 OR 0:15 IN CALLING ROUTINE.
!             ICATS() = CATEGORY-INTERVAL VALUES IN PHYSICAL UNITS
!                       (OUTPUT).  SHOULD BE DIMENSIONED BY 16 OR 0:15
!                       IN CALLING ROUTINE.
!             NDIMS_P = ACTUAL NUMBER OF ROWS AND COLUMNS IN OUTPUT
!                       PRODUCT (OUTPUT) (ROWS=COLUMNS IN 88D FILES)
!           IHEADER() = BYTE ARRAY WITH COPY OF HEADER.  SHOULD BE
!                       DIMENSIONED 160 IN CALLING ROUTINE (OUTPUT)
!          LEN_HEADER = ACTUAL NUMBER OF BYTES OF INFORMATION IN IHEADER()
!          ITRAILER() = CHARACTER*1 ARRAY WITH TEXT INFORMATION ATTACHED
!                       AT END OF FILE (IF ANY). SHOULD BE DIMENSIONED BY
!                       20 000 IN CALLING ROUTINE (OUTPUT)
!         LEN_TRAILER = ACTUAL NUMBER OF BYTES OF DATA IN ITRAILER().
!                       WILL BE SET 0 IF NO TRAILER EXISTS (OUTPUT) 
!
!------------------------------------------------------------------------

  INTEGER, INTENT(IN) :: iheaderdim,itrailerdim
  INTEGER, INTENT(IN) :: msglen, nrows, ncols

  INTEGER, INTENT(IN) :: ival(msglen)

  INTEGER, INTENT(OUT) :: icode, ivcp, isite, isite_lat,                &
                                      isite_lon, isite_elev,            &
                                      iyr, imonth, iday, ihr, iminit,   &
                                      isec, iprod, maxval1, len_header, &
                                      len_trailer, ndims_p

  INTEGER, INTENT(OUT) :: icatt(0:15)
  INTEGER, INTENT(OUT) :: icats(0:15)
  
  INTEGER, INTENT(OUT) :: iheader(iheaderdim),ifield(ncols, nrows)

  CHARACTER(LEN=1), INTENT(OUT) :: itrailer(itrailerdim)

!------------------------------------------------------------------------
!
! Declare internal variables
!
!------------------------------------------------------------------------

  INTEGER :: ipcode
  INTEGER :: msglen_int, len1, ival0, ival1

  INTEGER :: n,k,nc

  INTEGER :: nex_time, nex_date, jul_date, idaywk, idayyr
  INTEGER :: nrow, nrun_marker
  INTEGER :: nb
  INTEGER :: npos1,npos2
  INTEGER :: icol,icol1,icol2
  INTEGER :: ival4,ncol0,ilvl
  INTEGER :: nbyte
!
!------------------------------------------------------------------------
!
! Include file
!
!------------------------------------------------------------------------

  INCLUDE 'nidscst.inc'
  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!------------------------------------------------------------------------
!
! Check internal and system-indicated file length for consistency.
!
!------------------------------------------------------------------------

  msglen_int = i4byte(ival(9),ival(10),ival(11),ival(12))

  IF (msglen < msglen_int) THEN
    icode = 4
    RETURN
  END IF

!------------------------------------------------------------------------
!
! Get product ID# and check for raster-type packet code.  Header length
! is constant 158 bytes for raster products.
!
!------------------------------------------------------------------------

  iprod = i2byte(ival(1),ival(2))
  ipcode = i2byte(ival(137),ival(138))

  IF ((ipcode /= ipraster1) .AND. (ipcode /= ipraster2)) THEN
    icode = 2
    RETURN
  END IF

  len_header = 158
  isite = i2byte(ival(13),ival(14))
  ivcp = i2byte(ival(35),ival(36))

!------------------------------------------------------------------------
!
! Set header and trailer arrays.  Text trailer is included in 4-km
! composite reflectivity product (storm cell info).
!
!------------------------------------------------------------------------
  
  DO n = 1,158
    iheader(n) = ival(n)
  END DO

  len1 = 136 + i4byte(ival(133),ival(134),ival(135),ival(136))
  len_trailer = msglen - len1

  IF (len_trailer > 0) THEN
    k = 0
    DO n = len1+1,msglen
      k = k + 1
      itrailer(k) = ACHAR(ival(n))
    END DO
  END IF
  IF (len_trailer < 0) len_trailer = 0

!------------------------------------------------------------------------
!
! Get volume scan date and time.
!
!------------------------------------------------------------------------

  nex_time = i4byte(ival(43),ival(44),ival(45),ival(46))
  nex_date = i2byte(ival(41),ival(42))
  jul_date = jul_1_1_70 + nex_date
  CALL w3fs26(jul_date,iyr,imonth,iday,idaywk,idayyr)
  ihr = nex_time / 3600
  iminit = (nex_time - (3600*ihr))/60
  isec = nex_time - (3600*ihr) - (60*iminit)

!------------------------------------------------------------------------
!
! Get site lat/lon, elevation
!
!------------------------------------------------------------------------

  isite_elev = i2byte(ival(29),ival(30))
  isite_lat = i4byte(ival(21),ival(22),ival(23),ival(24))
  isite_lon = i4byte(ival(25),ival(26),ival(27),ival(28))

!------------------------------------------------------------------------
!
! Confirm number of rows in the product.
!
!------------------------------------------------------------------------

  ndims_p = i2byte(ival(155),ival(156))

  IF ((ncols < ndims_p) .OR. (nrows < ndims_p)) THEN
    icode = 3
    RETURN
  END IF

!------------------------------------------------------------------------
!
! Assign category threshold values.  Scale up the physical values if
! indicated.
!
!------------------------------------------------------------------------

  nc = 0

  DO n = 61, 91, 2
    icatt(nc) = i2byte(ival(n),ival(n+1))
    ival0 = ival(n+1)
    IF (ival0 < 0) ival0 = ival0 + 256
    ival1 = ival(n)
    IF (MOD(ival1,2) == 0) THEN
      icats(nc) = ival0
    ELSE
      icats(nc) = -1*ival0
    END IF
    nc = nc + 1
  END DO

  icats(0) = 0

  IF (icatt(0) >= iscale5) THEN
    DO nc = 1, 15
      icats(nc) = 5*icats(nc)
    END DO
  END IF

  maxval1 = i2byte(ival(93),ival(94))

!------------------------------------------------------------------------
!
! Row by row, plug run-length encoded values into the output array.
! Determine how many bytes store the current row of data, then parse
! out each byte as a run, plug the run of values into the grid.
!
!------------------------------------------------------------------------

  nrow = 0
  nrun_marker= 159

  DO WHILE (nrow < ndims_p)
    nrow = nrow + 1
    nb = i2byte(ival(nrun_marker),ival(nrun_marker+1))
    npos1 = nrun_marker + 2
    npos2 = npos1 + nb - 1
    icol1 = 1
    DO nbyte=npos1,npos2
      ival4 = ival(nbyte)
      IF (ival4 < 0) ival4 = 256 + ival4 
      ncol0 = ival4/16
      ilvl = ival4 - (16*ncol0)
      icol2 = icol1 + ncol0 - 1
      DO icol=icol1,icol2
        ifield(icol,nrow) = ilvl
      END DO
      icol1 = icol2 + 1
    END DO

    nrun_marker = npos2 + 1
  END DO
     
  icode = 0
  RETURN
END SUBROUTINE get_raster_data

!########################################################################
!########################################################################
!#########                                                      #########
!#########                SUBROUTINE get_dpa_data               #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE get_dpa_data(iheaderdim,itrailerdim,msglen,ival,             &
                        isite,isite_lat,isite_lon,isite_elev,           &
                        iyr,imonth,iday,ihr,iminit,isec,                &
                        iprod,maxval1,ivcp,ifield,ncols,nrows,          &
                        ncolsp,nrowsp,iheader,len_header,itrailer,      &
                        len_trailer,icode)

!------------------------------------------------------------------------
!
! PURPOSE:
!
! Decodes a WSR-88D digital precip array (DPA) message and returns all
! pertinent data.  The digital values ifield() are converted to
! rainfall (mm) through the formulae:
!
! dBa = -6.125 + ifield()*0.125
! rainfall (mm) = 10.**(dBa/10.)
!
! ifield() of 0 indicates zero rainfall.
! ifield() of 255 indicates missing data.
!
! Note:  This subroutine uses nonsystem subroutines:
!            W3FS26 
!            I2BYTE
!            I4BYTE
!
!------------------------------------------------------------------------
!
! AUTHOR: Eric Kemp, 2 August 2001.
!         Based on HPUX FORTRAN 77 code developed by Kitzmiller of the
!         Techniques Development Laboratory (August 1997).
!
! MODIFICATION HISTORY:
!
! Eric Kemp, 9 August 2001
! Now uses dimension parameters from N2ACONST module.
!
! Eric Kemp, 17 August 2001
! Modified to avoid use of INTEGER*4, INTEGER*1, etc., based on work
! by Yunheng Wang.  Also, changed variable maxval to maxval1 to avoid
! conflict with Fortran 90 implicit function.
!
! Keith Brewster, 6 April 2010
! Modifications to accomodate conversion of n2aconst to nidscst.inc
!
!------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------
!
! Arguments (original documentation)
!
!              IVAL() = BYTE OR CHARACTER*1 ARRAY CONTAINING DPA
!                       MESSAGE.  MUST BE DIMENSIONED .GE. MSGLEN
!                       (INPUT).
!              MSGLEN = LENGTH OF INPUT MESSAGE IN BYTES (INPUT)
!               ISITE = WSR-88D NUMERIC SITE ID (OUTPUT)
! ISITE_LAT,ISITE_LON = SITE LAT/LON (DEG E), DEG X 1000 (OUTPUT)
!          ISITE_ELEV = SITE ELEVATION, M MSL (OUTPUT)
!     IYR,IMONTH,IDAY = PRODUCT VOLUME DATE: YEAR,MONTH,DAY (OUTPUT)
!     IHR,IMINIT,ISEC = PRODUCT TIME, UTC: HOUR,MINUTE,SECOND (OUTPUT)
!
!               IPROD = PRODUCT CODE NUMBER (OUTPUT)
!             MAXVAL1 = MAXIMUM VALUE IN THE RASTER FIELD (OUTPUT)
!                IVCP = OPERATING VOLUME COVERAGE PATTERN (OUTPUT)
!           IFIELD(,) = RETURNED VALUES (DATA LEVELS) FOR THE
!                       FIELD (OUTPUT); INTEGER*2 VARIABLE.  SHOULD BE
!                       DIMENSIONED BY NCOLS,NROWS IN CALLING ROUTINE.
!         NCOLS,NROWS = DIMENSIONS OF IFIELD(), AT LEAST 131 (INPUT)
!       NCOLSP,NROWSP = ACTUAL NUMBER OF ROWS AND COLUMNS IN OUTPUT
!                       PRODUCT (OUTPUT) (ROWS=COLUMNS IN 88D FILES)
!           IHEADER() = BYTE ARRAY WITH COPY OF HEADER.  SHOULD BE
!                       DIMENSIONED 160 IN CALLING ROUTINE (OUTPUT)
!          LEN_HEADER = ACTUAL NUMBER OF BYTES OF INFORMATION IN 
!                       IHEADER() (OUTPUT)
!          ITRAILER() = CHARACTER*1 ARRAY WITH TEXT INFORMATION ATTACHED
!                       AT END OF FILE (IF ANY). SHOULD BE DIMENSIONED
!                       BY 20 000 IN CALLING ROUTINE (OUTPUT)
!         LEN_TRAILER = ACTUAL NUMBER OF BYTES OF DATA IN ITRAILER().
!                       WILL BE SET 0 IF NO TRAILER EXISTS (OUTPUT) 
!               ICODE = OUTPUT RETURN CODE, CONDITIONS AS FOLLOWS:
!                       0 : DECODED PROPERLY      
!                       2 : NOT A RECOGNIZED DPA PRODUCT
!                       3 : NCOLS/NROWS TOO SMALL FOR PRODUCT
!                       4 : PRODUCT TOO LARGE FOR MSGLEN
!
!------------------------------------------------------------------------

  INTEGER, INTENT(IN) :: iheaderdim,itrailerdim
  INTEGER, INTENT(IN) :: msglen, ncols, nrows
  INTEGER, INTENT(IN) :: ival(msglen)

  INTEGER, INTENT(OUT) :: iheader(iheaderdim)
  INTEGER, INTENT(OUT) :: ifield(ncols,nrows), icode
  INTEGER, INTENT(OUT) :: isite, isite_lat, isite_lon, isite_elev,       &
                            iyr, imonth, iday, ihr, iminit, isec,        &
                            iprod, maxval1, ivcp, ncolsp, nrowsp,        &
                            len_header, len_trailer

  CHARACTER(LEN=1), INTENT(OUT) :: itrailer(itrailerdim)

!------------------------------------------------------------------------
!
! Declare external functions.
!
!------------------------------------------------------------------------

  INTEGER, EXTERNAL :: i2byte, i4byte

!------------------------------------------------------------------------
!
! Internal variables
!
!------------------------------------------------------------------------

  INTEGER :: ipcode
  
  INTEGER :: msglen_int, nrow, nb, len1,nex_time, nex_date, jul_date
  INTEGER :: idaywk, idayyr
  INTEGER :: nrun_marker, npos1, npos2, icol, icol1, icol2, nbyte
  INTEGER :: ncol0, ilvl, n, k
!
!------------------------------------------------------------------------
!
! Include file
!
!------------------------------------------------------------------------

  INCLUDE 'nidscst.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!------------------------------------------------------------------------
!
! Check internal and system-indicated file length for consistency.
!  
!------------------------------------------------------------------------

  msglen_int = i4byte(ival(9),ival(10),ival(11),ival(12))
  
  IF (msglen < msglen_int) THEN
    icode = 4
    RETURN
  END IF

!------------------------------------------------------------------------
!
! Get product type, and check for raster-type packet code.  Header 
! length is constant 146 bytes for DPA product.
!
!------------------------------------------------------------------------

  ipcode = i2byte(ival(137),ival(138))
  iprod = i2byte(ival(1),ival(2))

  IF (iprod /= 81) THEN
    icode = 2
    RETURN
  END IF

  len_header = 146
  isite = i2byte(ival(13),ival(14))
  ivcp = i2byte(ival(35),ival(36))

!------------------------------------------------------------------------
! 
! Set header and trailer arrays.  Text trailer is included in 4-km
! composite reflectivity product (storm cell info).
!
!------------------------------------------------------------------------

  DO n = 1,146
    iheader(n) = ival(n)
  END DO

  len1 = 136 + i4byte(ival(133),ival(134),ival(135),ival(136))
  len_trailer = msglen - len1

  IF (len_trailer > 0) THEN
    k = 0
    DO n = len1+1,msglen
      k = k + 1
      itrailer(k) = ACHAR(ival(n))
    END DO
  END IF
  IF (len_trailer < 0) len_trailer = 0

!------------------------------------------------------------------------
!
! Get volume scan date and time
!
!------------------------------------------------------------------------

  nex_time = i4byte(ival(43),ival(44),ival(45),ival(46))
  nex_date = i2byte(ival(41),ival(42))
  jul_date = jul_1_1_70 + nex_date
  CALL w3fs26(jul_date,iyr,imonth,iday,idaywk,idayyr)
  ihr = nex_time / 3600 
  iminit = (nex_time - (3600*ihr))/60
  isec = nex_time - (3600*ihr) - (60*iminit)

!------------------------------------------------------------------------
!
! Get site lat/lon, elevation.
!
!------------------------------------------------------------------------

  isite_elev = i2byte(ival(29),ival(30))
  isite_lat = i4byte(ival(21),ival(22),ival(23),ival(24))
  isite_lon = i4byte(ival(25),ival(26),ival(27),ival(28))

!------------------------------------------------------------------------
!
! Get maximum value.
!
!------------------------------------------------------------------------

  maxval1 = i2byte(ival(93),ival(94))

!------------------------------------------------------------------------
!
! Confirm number of rows and columns in the product.
!
!------------------------------------------------------------------------

  ncolsp = i2byte(ival(143),ival(144))

  IF (ncols < ncolsp) THEN
    icode = 3
    RETURN
  END IF

  nrowsp = i2byte(ival(145),ival(146))

  IF (nrows < nrowsp) THEN
    icode = 3
    RETURN
  END IF

!------------------------------------------------------------------------
!
! Row by row, plug run-length encoded values into the output array.
! Determine how many bytes store the current row of data, then parse
! out each byte as a run and plug the run of values into the grid.
!
! For DPA product, run length and value are stored in 2-byte words
! rather than one byte as in other gridded products.
!
!------------------------------------------------------------------------

  ifield(:,:) = 255

  nrow = 0
  nrun_marker = 147

  DO WHILE (nrow < nrowsp)
    nrow = nrow + 1
    nb = i2byte(ival(nrun_marker),ival(nrun_marker+1))
    npos1 = nrun_marker + 2
    npos2 = npos1 + nb - 1
    icol1 = 1

    DO nbyte = npos1,npos2,2
      ncol0 = ival(nbyte)
      ilvl = ival(nbyte+1)
      icol2 = icol1 + ncol0 - 1
      DO icol=icol1,icol2
        ifield(icol,nrow) = ilvl
      END DO
      icol1 = icol2 + 1
    END DO ! DO nbyte = npos1,npos2,2

    nrun_marker = npos2 + 1
  END DO ! DO WHILE    

!------------------------------------------------------------------------
!
! The end.
!
!------------------------------------------------------------------------

  icode = 0
  RETURN
END SUBROUTINE get_dpa_data


!########################################################################
!########################################################################
!#########                                                      #########
!#########                SUBROUTINE get_dhr_data               #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE get_dhr_data(ival,msglen,isite,                               &
                        isite_lat,isite_lon,isite_elev,                  &
                        iyr,imonth,iday,                                 &
                        ihr,iminit,iprod,maxval1,ivcp,ifield,maxradials, &
                        maxbins,iazmuth,nbins,nradials,icode)

!------------------------------------------------------------------------
!
! PURPOSE:
!
! Decodes WSR-88D digital hybrid scan (DHR) product.
!
! Note 1:  Raw digital values 0-255 are converted to dBZ by relationship:
!        
!          dBZ = ((ifield-2.)/2.) - 32.
!
! Note 2:  This subroutine uses non-system subroutines:
!            W3FS26 
!            I2BYTE 
!            I4BYTE
!
!------------------------------------------------------------------------
!
! AUTHOR:
!
! Eric Kemp, 1 August 2001
! Based on HPUX FORTRAN 77 subroutine written by Kitzmiller of the
! Techniques Development Laboratory (June 1998).
!
! Eric Kemp, 17 August 2001
! Modified to avoid use of INTEGER*4, INTEGER*1, etc., based on work
! by Yunheng Wang.  Also, changed variable maxval to maxval1 to avoid
! conflict with Fortran 90 implicit function.
!
! Keith Brewster, 6 April 2010
! Modifications to accomodate conversion of n2aconst to nidscst.inc
!
!------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------
!
! Declare arguments (original documentation)
!
!
!             IVAL( ) = BYTE OR CHARACTER*1 VECTOR OF LENGTH MSGLEN,
!                       HOLDING THE DIGITAL HYBRID SCAN MESSAGE (INPUT)
!              MSGLEN = LENGTH OF THE MESSAGE (SHOULD BE 85100) (INPUT)
!               ISITE = WSR-88D SITE ID (OUTPUT)
! ISITE_LAT,ISITE_LON = SITE LAT/LON (DEG E), DEG X 1000 (OUTPUT)
!          ISITE_ELEV = SITE ELEVATION, M MSL (OUTPUT)
!     IYR,IMONTH,IDAY = PRODUCT VOLUME SCAN DATE: YEAR,MONTH, DAY (OUTPUT)
!          IHR,IMINIT = PRODUCT TIME, UTC: HOUR,MINUTE,SECOND (OUTPUT)
!
!               IPROD = PRODUCT CODE NUMBER (OUTPUT)
!             MAXVAL1 = MAXIMUM DIGITAL VALUE IN PRODUCT (OUTPUT)
!                IVCP = OPERATING VOLUME COVERAGE PATTERN (OUTPUT)
!           IFIELD(,) = RETURNED VALUES (DATA LEVELS) FOR THE
!                       FIELD (OUTPUT); AN INTEGER*2 ARRAY.  SHOULD BE
!                       DIMENSIONED BY MAXBINS,MAXRADIALS IN CALLING
!                       ROUTINE.
!
!
!  MAXRADIALS,MAXBINS = DIMENSIONS FOR IFIELD(), IAZMUTH(), INDICATING
!                       NUMBER OF RADIALS AND NUMBER OF RANGE BINS
!                       SHOULD BE .GE. 360 AND 230, RESPECTIVELY (INPUT)
!
!           IAZMUTH() = ANTENNA AZIMUTH DIRECTIONS (DEG X 10) FOR
!                       EACH RADIAL (OUTPUT).  DIMENSIOND BY MAXRADIALS IN
!                       CALLING ROUTINE
!      NBINS,NRADIALS = NUMBER OF RANGE BINS AND RADIALS IN THE 
!                         GRAPHIC (OUTPUT)
!               ICODE = OUTPUT RETURN CODE, CONDITIONS AS FOLLOWS:
!                       0 : MESSAGE DECODED OK  
!                       2 : NOT A DHR (ID CODE 32) PRODUCT
!                       3 : MAXRADIALS OR MAXBINS IS TOO SMALL
!                       4 : PRODUCT LENGTH GREATER THAN MSGLEN
!
!------------------------------------------------------------------------
  
  INTEGER, INTENT(IN) :: msglen,maxradials,maxbins
  INTEGER, INTENT(IN) :: ival(msglen)

  INTEGER, INTENT(OUT) :: ifield(maxbins,maxradials)

  INTEGER, INTENT(OUT) :: icode,ivcp,isite,isite_lat,isite_lon,          &
                            isite_elev,iyr,imonth,iday,ihr,iminit,       &
                            iprod,nbins,nradials

  INTEGER, INTENT(OUT) :: maxval1

  INTEGER, INTENT(OUT) :: iazmuth(maxradials)

!------------------------------------------------------------------------
!
! Declare external functions.
!
!------------------------------------------------------------------------

  INTEGER, EXTERNAL :: i2byte, i4byte

!------------------------------------------------------------------------
!
!     INTERNAL VARIABLES:
!
!         IVAL0,IVAL1 = INTEGER*4 WORDS USED TO HOLD LOWER AND
!                       HIGHER-ORDER HALVES OF 2-BYTE INTEGER
!          JUL_1_1_70 = JULIAN DATE OF 1JAN1970, USED TO CONVERT NEXRAD 
!                       INTERNAL DATE TO YR/MON/DAY (PARAMETER)
!    NEXDATE,JUL_DATE = NEXRAD REFERENCE DATE, ASTRONOMICAL JULIAN DATE, 
!                         USED TO DETERMINE LEGAL DATE (YEAR-MONTH-DAY)
!         NPOS1,NPOS2 = POSITIONS OF 1ST AND LAST BYTES (RELATIVE TO START
!                         OF THE FILE) IN THE CURRENT RADIAL BEING DECODED
!               NRUNS = NUMBER OF RUNS NEEDED TO STORE THE PRESENT RADIAL
!
!------------------------------------------------------------------------

  INTEGER :: ival4
  INTEGER :: nex_time,nex_date,jul_date,idaywk,idayyr,msglen_int

  INTEGER :: npos1,npos2
  INTEGER :: nr, nrun_marker
  INTEGER :: nbyte,ibin
!
!------------------------------------------------------------------------
!
! Include file
!
!------------------------------------------------------------------------

  INCLUDE 'nidscst.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!  
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!------------------------------------------------------------------------
!
! Get product type and check for product code 32 for DHR.
!
!------------------------------------------------------------------------

  iprod = i2byte(ival(1),ival(2))

  IF (iprod /= 32) THEN
    icode = 2
    RETURN
  END IF

  isite = 256*ival(13) + ival(14)
  ivcp = 256*ival(35) + ival(36)

!------------------------------------------------------------------------
!
! Get volume scan date and time.
!
!------------------------------------------------------------------------

  nex_time = i4byte(ival(43),ival(44),ival(45),ival(46))
  nex_date = i2byte(ival(41),ival(42))
  jul_date = jul_1_1_70 + nex_date
  CALL w3fs26(jul_date,iyr,imonth,iday,idaywk,idayyr)
  ihr = nex_time / 3600
  iminit = (nex_time - (3600*ihr)) / 60

!------------------------------------------------------------------------
!
! Get site lat/lon, elevation.
!
!------------------------------------------------------------------------

  isite_elev = i2byte(ival(29),ival(30))
  isite_lat = i4byte(ival(21),ival(22),ival(23),ival(24))
  isite_lon = i4byte(ival(25),ival(26),ival(27),ival(28))

!------------------------------------------------------------------------
!
! Find number of radials in the product.
!
!------------------------------------------------------------------------

  nbins = i2byte(ival(141),ival(142))
  nradials = i2byte(ival(149),ival(150))

  IF ( (maxbins < nbins) .OR. (maxradials < nradials) ) THEN
    icode = 3
    RETURN
  END IF

!------------------------------------------------------------------------
!
! Check internal and user-indicated file length for consistency.
!
!------------------------------------------------------------------------

  msglen_int = i4byte(ival(9),ival(10),ival(11),ival(12))

  IF (msglen < msglen_int) THEN
    icode = 4
    RETURN
  END IF

!------------------------------------------------------------------------
!
! Initialize the reflectivity array with 0 values.
!
!------------------------------------------------------------------------

  ifield(:,:) = 0

!------------------------------------------------------------------------
!
! Radial by radial, read in the data.  There is no run-length encoding
! for DHR product.
!
!------------------------------------------------------------------------

  nr = 0
  nrun_marker = 151
  maxval1 = 0

  DO WHILE (nr < nradials)
    nr = nr + 1
    iazmuth(nr) = i2byte(ival(nrun_marker+2),ival(nrun_marker+3))
    npos1 = nrun_marker + 6
    npos2 = npos1 + 229

    ibin = 0
    DO nbyte = npos1,npos2
      ival4 = ival(nbyte)
      IF (ival4 < 0) ival4 = 256 + ival4
      ibin = ibin + 1
      ifield(ibin,nr) = ival4
      IF (ival4 > maxval1) maxval1 = ival4
    END DO

!------------------------------------------------------------------------
!
!   Finished mapping this radial -- go on to the next.
!
!------------------------------------------------------------------------
    
    nrun_marker = npos2 + 1

  END DO ! DO WHILE (nr < nradials)

!------------------------------------------------------------------------
!
! The end.
!
!------------------------------------------------------------------------
 
  icode = 0
  RETURN
END SUBROUTINE get_dhr_data

!########################################################################
!########################################################################
!#########                                                      #########
!#########                    FUNCTION i2byte                   #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

INTEGER FUNCTION i2byte(ival1,ival0)

!------------------------------------------------------------------------
!
! PURPOSE:
!
! Combines the contents of two 1-byte words into a single 4-byte signed
! integer word.
!
!------------------------------------------------------------------------
!
! AUTHOR: 
!
! MODIFICATIONS:
!
! Yunheng Wang  (July 2001)
! Converted from a FORTRAN 77 code developed by NWS Techniques 
! Development Laboratory in April 1998
!
! Eric Kemp, 30 July 2001
! Added USE statement for module n2aconst, minor changes.
!
! Eric Kemp, 9 August 2001
! Use variable type INTEGER*4 and similar types instead of using
! KIND statements.  
!
! Eric Kemp, 17 August 2001
! Modified to avoid use of INTEGER*4, INTEGER*1, etc., based on work
! by Yunheng Wang.
!
! Keith Brewster, 6 April 2010
! Modifications to accomodate conversion of n2aconst to nidscst.inc
!
!------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------
!
! VARIABLES
!      ival1 = high-order 8 bits of integer word     (byte, INPUT)
!      ival0 = low-order 8 bits of integer word      (byte, INPUT)
!      kval4b = signed sum of the two bytes          (     OUTPUT)
!
!------------------------------------------------------------------------

  INTEGER, INTENT(IN) :: ival1,ival0

!------------------------------------------------------------------------
!
! Misc. local variables
!
!------------------------------------------------------------------------
 
  INTEGER :: ivalh, ivall

!------------------------------------------------------------------------
!
! Include file
!
!------------------------------------------------------------------------

  INCLUDE 'nidscst.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ivalh = ival1;  ivall = ival0

  IF (ival0 < 0) ivall = 256 + ival0
  IF (ival1 < 0) ivalh = 256 + ival1
     
  i2byte = ivall + (256*ivalh)

  RETURN  
END FUNCTION i2byte

!########################################################################
!########################################################################
!#########                                                      #########
!#########                    FUNCTION i4byte     	        #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

INTEGER FUNCTION i4byte(ival3, ival2, ival1,ival0) 

!------------------------------------------------------------------------
!
! PURPOSE:
!
! Combines the contents of four 1-byte words into a single 4-byte signed
! integer word.
!
!------------------------------------------------------------------------
!
! AUTHOR: 
!
! MODIFICATIONS:
!
! Yunheng Wang (July 2001)
! Converted from a FORTRAN 77 code developed by NWS Techniques 
! Development Laboratory in April 1998.
!
! Eric Kemp, 30 July 2001.
! Added USE statement for module n2aconst, minor changes.
!
! Eric Kemp, 9 August 2001
! Use variable type INTEGER*4 and similar types instead of using
! KIND statements.  
!
! Eric Kemp, 17 August 2001
! Modified to avoid use of INTEGER*4, INTEGER*1, etc., based on work
! by Yunheng Wang.
!
! Keith Brewster, 6 April 2010
! Modifications to accomodate conversion of n2aconst to nidscst.inc
!
!------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------
!
! VARIABLES:
!      IVAL3 = HIGH-ORDER 8 BITS OF INTEGER WORD        (byte, INPUT)
!      IVAL2 = 2ND HIGHEST 8 BITS OF INTEGER WORD       (byte, INPUT)
!      IVAL1 = NEXT-TO-LOWEST 8 BITS OF INTEGER WORD    (byte, INPUT)
!      IVAL0 = LOW-ORDER 8 BITS OF INTEGER WORD         (byte, INPUT)
!
!------------------------------------------------------------------------

  INTEGER, INTENT(IN) :: ival0,ival1,ival2,ival3
  
!------------------------------------------------------------------------
!
! Misc. local variables
!
!------------------------------------------------------------------------
 
  INTEGER :: ivalhh, ivalh, ivall, ivalll

!------------------------------------------------------------------------
!
! Include file
!
!------------------------------------------------------------------------

  INCLUDE 'nidscst.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  ivalhh = ival3;   ivalh = ival2
  ivall  = ival1;  ivalll = ival0

  IF (ival3 < 0) ivalhh = 256 + ival3
  IF (ival2 < 0) ivalh  = 256 + ival2
  IF (ival1 < 0) ivall  = 256 + ival1
  IF (ival0 < 0) ivalll = 256 + ival0

  i4byte = ivalll + (256*ivall) + (2**16)*ivalh + (2**24)*ivalhh

  RETURN
END FUNCTION i4byte

!########################################################################
!########################################################################
!#########                                                      #########
!#########                   SUBROUTINE w3fs26                  #########
!#########                                                      #########
!#########                     Developed by                     #########
!#########     Center for Analysis and Prediction of Storms     #########
!#########                University of Oklahoma                #########
!#########                                                      #########
!########################################################################
!########################################################################

SUBROUTINE w3fs26(jldayn,iyear,month,iday,idaywk,idayyr)

!------------------------------------------------------------------------
!
! PURPOSE:
!
! Computes year (4 digits), month, day, day of week, day
! of year from julian day number. This subroutine will work
! from 1583 a.d. to 3300 a.d.
!
!------------------------------------------------------------------------
!
! AUTHOR: 
!
! MODIFICATIONS:
!
! Yunheng Wang  (July 2001)
! Converted from a FORTRAN 77 code developed by Jones, R.E. 
! in March. 29, 1987
!
! Eric Kemp, 30 July 2001
! Added USE statement for module n2aconst, minor changes.
!
! Eric Kemp, 9 August 2001
! Use variable type INTEGER*4 and similar types instead of using
! KIND statements.  
!
! Eric Kemp, 17 August 2001
! Modified to avoid use of INTEGER*4, INTEGER*1, etc., based on work
! by Yunheng Wang.
!
! Keith Brewster, 6 April 2010
! Modifications to accomodate conversion of n2aconst to nidscst.inc
!
!------------------------------------------------------------------------

  IMPLICIT NONE

!------------------------------------------------------------------------
!
! INPUT:
!
!     jldayn:  julian day number
!
!------------------------------------------------------------------------
!
! OUTPUT:
!
!     iyear:    year  (4 digits)
!     month:    month
!     iday:     day
!     idaywk:   day of week (1 is sunday, 7 is sat)
!     idayyr:   day of year (1 to 366)
!
!------------------------------------------------------------------------

  INTEGER, INTENT(IN)  :: jldayn
  INTEGER, INTENT(OUT) :: iyear, month, iday, idaywk, idayyr
  
!------------------------------------------------------------------------
!
! Misc. local variables
!
!------------------------------------------------------------------------

  INTEGER :: l, n, i, j
 
!------------------------------------------------------------------------
!
! Include file
!
!------------------------------------------------------------------------

  INCLUDE 'nidscst.inc'

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
! Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!------------------------------------------------------------------------
!
! REMARKS: 
!
!     A julian day number can be computed by using one of the
!     following statement functions. A day of week can be computed
!     from the julian day number. A day of year can be computed from
!     a julian day number and year.
!
!     JULIAN DAY NUMBER:
!
!        1:   iyear (4 digits)
!
!        jdn(iyear,month,iday) = iday - 32075                          &
!                + 1461 * (iyear + 4800 + (month - 14) / 12) / 4       &
!                + 367 * (month - 2 - (month -14) / 12 * 12) / 12      &
!                - 3 * ((iyear + 4900 + (month - 14) / 12) / 100) / 4
!
!        2:   iyr (4 digits) , idyr(1-366) day of year
!
!        julian(iyr,idyr) = -31739 + 1461 * (iyr + 4799) / 4           &
!                        -3 * ((iyr + 4899) / 100) / 4 + idyr
!
!     DAY OF WEEK FROM JULIAN DAY NUMBER: 1 is sunday, 7 is saturday.
!
!        jdaywk(jldayn) = mod((jldayn + 1),7) + 1
!
!     DAY OF YEAR FROM JULIAN DAY NUMBER AND 4 DIGIT YEAR:
!
!        jdayyr(jldayn,iyear) = jldayn -                               &
!               (-31739+1461*(iyear+4799)/4-3*((iyear+4899)/100)/4)
!
!     The first function was in a letter to the editor communications
!     of the acm  volume 11 / number 10 / october, 1968. The 2nd
!     function was derived from the first. This subroutine was also
!     included in the same letter. Julian day number 1 is
!     jan 1,4713 b.c. A julian day number can be used to replace a
!     day of century, this will take care of the date problem in
!     the year 2000, or reduce program changes to one line change
!     of 1900 to 2000. Julian day numbers can be used for finding
!     record numbers in an archive or day of week, or day of year.
!
!------------------------------------------------------------------------
  
   l      = jldayn + 68569
   n      = 4 * l / 146097
   l      = l - (146097 * n + 3) / 4
   i      = 4000 * (l + 1) / 1461001
   l      = l - 1461 * i / 4 + 31
   j      = 80 * l / 2447
   iday   = l - 2447 * j / 80
   l      = j / 11
   month  = j + 2 - 12 * l
   iyear  = 100 * (n - 49) + i + l
   idaywk = MOD((jldayn + 1),7) + 1
   idayyr = jldayn -                                                  &
       (-31739 +1461 * (iyear+4799) / 4 - 3 * ((iyear+4899)/100)/4)
   RETURN

END SUBROUTINE w3fs26
