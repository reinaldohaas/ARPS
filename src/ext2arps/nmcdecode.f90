!
!##################################################################
!##################################################################
!######                                                      ######
!######                    NMCDECODE                         ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Decode the NMC GRIB files.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: NWS/NMC
!
!  MODIFICATION HISTORY:
!    01/16/1996 (Yuhe Liu)
!    Changed subroutine name of SBYTE(S) and GBYTE(S) to GRBSBYTE(S)
!    and GRBGBYTE(S) to avoid the conflict with the subroutines of
!    same names in NCAR graphics package.
!
!    2003-12-29 (Richard Carpenter)
!    Increased threshold for special (NMC?) handling for large PDS lengths
!    from 50 to 256.
!
!
!-----------------------------------------------------------------------
!
!    Original comments:
!
!    FILE UPDATED ON 11-22-94 BY M. FARLEY.
!   CONTAINS SUBROUTINES NEEDED TO DECODE GRIB MESSAGES:
!     W3FI63 AND W3FI83.
!

SUBROUTINE w3fi63(msga,kpds,kgds,kbms,DATA,kptr,kret)
!$$$  SUBPROGRAM DOCUMENTATION  BLOCK
!             .      .    .                                       .
! SUBPROGRAM:  W3FI63        UNPK GRIB FIELD TO GRIB GRID
!   PRGMMR: FARLEY           ORG: NMC421      DATE:94-11-22
!
! ABSTRACT: UNPACK A GRIB (EDITION 1) FIELD TO THE EXACT GRID
!   SPECIFIED IN THE GRIB MESSAGE, ISOLATE THE BIT MAP, AND MAKE
!   THE VALUES OF THE PRODUCT DESCRIPTON SECTION (PDS) AND THE
!   GRID DESCRIPTION SECTION (GDS) AVAILABLE IN RETURN ARRAYS.
!
!   WHEN DECODING IS COMPLETED, DATA AT EACH GRID POINT HAS BEEN
!       RETURNED IN THE UNITS SPECIFIED IN THE GRIB MANUAL.
!
! PROGRAM HISTORY LOG:
!   91-09-13  CAVANAUGH
!   91-11-12  CAVANAUGH   MODIFIED SIZE OF ECMWF GRIDS 5-8
!   91-12-22  CAVANAUGH   CORRECTED PROCESSING OF MERCATOR PROJECTIONS
!                      IN GRID DEFINITION SECTION (GDS) IN
!                      ROUTINE FI633
!   92-08-05  CAVANAUGH   CORRECTED MAXIMUM GRID SIZE TO ALLOW FOR
!                      ONE DEGREE BY ONE DEGREE GLOBAL GRIDS
!   92-08-27  CAVANAUGH   CORRECTED TYPO ERROR, ADDED CODE TO COMPARE
!                      TOTAL BYTE SIZE FROM SECTION 0 WITH SUM OF
!                      SECTION SIZES.
!   92-10-21  CAVANAUGH   CORRECTIONS WERE MADE (IN FI634) TO REDUCE
!                      PROCESSING TIME FOR INTERNATIONAL GRIDS.
!                      REMOVED A TYPOGRAPHICAL ERROR IN FI635.
!   93-01-07  CAVANAUGH   CORRECTIONS WERE MADE (IN FI635) TO
!                      FACILITATE USE OF THESE ROUTINES ON A PC.
!                      A TYPOGRAPHICAL ERROR WAS ALSO CORRECTED
!   93-01-13  CAVANAUGH   CORRECTIONS WERE MADE (IN FI632) TO
!                      PROPERLY HANDLE CONDITION WHEN
!                      TIME RANGE INDICATOR = 10.
!                      ADDED U.S.GRID 87.
!   93-02-04  CAVANAUGH   ADDED U.S.GRIDS 85 AND 86
!   93-02-26  CAVANAUGH   ADDED GRIDS 2, 3, 37 THRU 44,AND
!                      GRIDS 55, 56, 90, 91, 92, AND 93 TO
!                      LIST OF U.S. GRIDS.
!   93-04-07  CAVANAUGH   ADDED GRIDS 67 THRU 77 TO
!                      LIST OF U.S. GRIDS.
!   93-04-20  CAVANAUGH   INCREASED MAX SIZE TO ACCOMODATE
!                      GAUSSIAN GRIDS.
!   93-05-26  CAVANAUGH   CORRECTED GRID RANGE SELECTION IN FI634
!                      FOR RANGES 67-71 & 75-77
!   93-06-08  CAVANAUGH   CORRECTED FI635 TO ACCEPT GRIB MESSAGES
!                      WITH SECOND ORDER PACKING. ADDED ROUTINE FI636
!                      TO PROCESS MESSAGES WITH SECOND ORDER PACKING.
!   93-09-22  CAVANAUGH   MODIFIED TO EXTRACT SUB-CENTER NUMBER FROM
!                      PDS BYTE 26
!   93-10-13  CAVANAUGH   MODIFIED FI634 TO CORRECT GRID SIZES FOR
!                      GRIDS 204 AND 208
!   93-10-14  CAVANAUGH   INCREASED SIZE OF KGDS TO INCLUDE ENTRIES FOR
!                      NUMBER OF POINTS IN GRID AND NUMBER OF WORDS
!                      IN EACH ROW
!   93-12-08  CAVANAUGH   CORRECTED TEST FOR EDITION NUMBER INSTEAD
!                      OF VERSION NUMBER
!   93-12-15  CAVANAUGH   MODIFIED SECOND ORDER POINTERS TO FIRST ORDER
!                      VALUES AND SECOND ORDER VALUES CORRECTLY
!                      IN ROUTINE FI636
!   94-03-02  CAVANAUGH   ADDED CALL TO W3FI83 WITHIN DECODER.  USER
!                      NO LONGER NEEDS TO MAKE CALL TO THIS ROUTINE
!   94-04-22  CAVANAUGH   MODIFIED FI635, FI636 TO PROCESS ROW BY ROW
!                      SECOND ORDER PACKING, ADDED SCALING CORRECTION
!                      TO FI635, AND CORRECTED TYPOGRAPHICAL ERRORS
!                      IN COMMENT FIELDS IN FI634
!   94-05-17  CAVANAUGH   CORRECTED ERROR IN FI633 TO EXTRACT RESOLUTION
!                      FOR LAMBERT-CONFORMAL GRIDS. ADDED CLARIFYING
!                      INFORMATION TO DOCBLOCK ENTRIES
!   94-05-25  CAVANAUGH   ADDED CODE TO PROCESS COLUMN BY COLUMN AS WELL
!                      AS ROW BY ROW ORDERING OF SECOND ORDER DATA
!   94-06-27  CAVANAUGH   ADDED PROCESSING FOR GRIDS 45, 94 AND 95.
!                      INCLUDES CONSTRUCTION OF SECOND ORDER BIT MAPS
!                      FOR THINNED GRIDS IN FI636.
!   94-07-08  CAVANAUGH   COMMENTED OUT PRINT OUTS USED FOR DEBUGGING
!   94-09-08  CAVANAUGH   ADDED GRIDS 220, 221, 223 FOR FNOC
!   94-11-10  FARLEY      INCREASED MXSIZE FROM 72960 TO 260000
!                      FOR .5 DEGREE SST ANALYSIS FIELDS
!
!
! USAGE:    CALL W3FI63(MSGA,KPDS,KGDS,KBMS,DATA,KPTR,KRET)
!   INPUT ARGUMENT LIST:
!  MSGA     - GRIB FIELD - "GRIB" THRU "7777"   CHAR*1
!                (MESSAGE CAN BE PRECEDED BY JUNK CHARS)
!
!   OUTPUT ARGUMENT LIST:
!  DATA     - ARRAY CONTAINING DATA ELEMENTS
!  KPDS     - ARRAY CONTAINING PDS ELEMENTS.  (EDITION 1)
!       (1)   - ID OF CENTER
!       (2)   - GENERATING PROCESS ID NUMBER
!       (3)   - GRID DEFINITION
!       (4)   - GDS/BMS FLAG (RIGHT ADJ COPY OF OCTET 8)
!       (5)   - INDICATOR OF PARAMETER
!       (6)   - TYPE OF LEVEL
!       (7)   - HEIGHT/PRESSURE , ETC OF LEVEL
!       (8)   - YEAR INCLUDING (CENTURY-1)
!       (9)   - MONTH OF YEAR
!       (10)  - DAY OF MONTH
!       (11)  - HOUR OF DAY
!       (12)  - MINUTE OF HOUR
!       (13)  - INDICATOR OF FORECAST TIME UNIT
!       (14)  - TIME RANGE 1
!       (15)  - TIME RANGE 2
!       (16)  - TIME RANGE FLAG
!       (17)  - NUMBER INCLUDED IN AVERAGE
!       (18)  - VERSION NR OF GRIB SPECIFICATION
!       (19)  - VERSION NR OF PARAMETER TABLE
!       (20)  - NR MISSING FROM AVERAGE/ACCUMULATION
!       (21)  - CENTURY OF REFERENCE TIME OF DATA
!       (22)  - UNITS DECIMAL SCALE FACTOR
!       (23)  - SUBCENTER NUMBER
!    (26-35)  - RESERVED
!    (24-25)  - RESERVED FOR FUTURE USE
!    (36-N)   - CONSECUTIVE BYTES EXTRACTED FROM PROGRAM
!               DEFINITION SECTION (PDS) OF GRIB MESSAGE
!  KGDS     - ARRAY CONTAINING GDS ELEMENTS.
!       (1)   - DATA REPRESENTATION TYPE
!       (19)  - NUMBER OF VERTICAL COORDINATE PARAMETERS
!       (20)  - OCTET NUMBER OF THE LIST OF VERTICAL COORDINATE
!               PARAMETERS
!               OR
!               OCTET NUMBER OF THE LIST OF NUMBERS OF POINTS
!               IN EACH ROW
!               OR
!               255 IF NEITHER ARE PRESENT
!       (21)  - FOR GRIDS WITH PL, NUMBER OF POINTS IN GRID
!       (22)  - NUMBER OF WORDS IN EACH ROW
!    LATITUDE/LONGITUDE GRIDS
!       (2)   - N(I) NR POINTS ON LATITUDE CIRCLE
!       (3)   - N(J) NR POINTS ON LONGITUDE MERIDIAN
!       (4)   - LA(1) LATITUDE OF ORIGIN
!       (5)   - LO(1) LONGITUDE OF ORIGIN
!       (6)   - RESOLUTION FLAG (RIGHT ADJ COPY OF OCTET 17)
!       (7)   - LA(2) LATITUDE OF EXTREME POINT
!       (8)   - LO(2) LONGITUDE OF EXTREME POINT
!       (9)   - DI LATITUDINAL DIRECTION OF INCREMENT
!       (10)  - DJ LONGITUDINAL DIRECTION INCREMENT
!       (11)  - SCANNING MODE FLAG (RIGHT ADJ COPY OF OCTET 28)
!    GAUSSIAN  GRIDS
!       (2)   - N(I) NR POINTS ON LATITUDE CIRCLE
!       (3)   - N(J) NR POINTS ON LONGITUDE MERIDIAN
!       (4)   - LA(1) LATITUDE OF ORIGIN
!       (5)   - LO(1) LONGITUDE OF ORIGIN
!       (6)   - RESOLUTION FLAG  (RIGHT ADJ COPY OF OCTET 17)
!       (7)   - LA(2) LATITUDE OF EXTREME POINT
!       (8)   - LO(2) LONGITUDE OF EXTREME POINT
!       (9)   - DI LATITUDINAL DIRECTION OF INCREMENT
!       (10)  - N - NR OF CIRCLES POLE TO EQUATOR
!       (11)  - SCANNING MODE FLAG (RIGHT ADJ COPY OF OCTET 28)
!       (12)  - NV - NR OF VERT COORD PARAMETERS
!       (13)  - PV - OCTET NR OF LIST OF VERT COORD PARAMETERS
!                          OR
!               PL - LOCATION OF THE LIST OF NUMBERS OF POINTS IN
!                    EACH ROW (IF NO VERT COORD PARAMETERS
!                    ARE PRESENT
!                          OR
!               255 IF NEITHER ARE PRESENT
!    POLAR STEREOGRAPHIC GRIDS
!       (2)   - N(I) NR POINTS ALONG LAT CIRCLE
!       (3)   - N(J) NR POINTS ALONG LON CIRCLE
!       (4)   - LA(1) LATITUDE OF ORIGIN
!       (5)   - LO(1) LONGITUDE OF ORIGIN
!       (6)   - RESOLUTION FLAG  (RIGHT ADJ COPY OF OCTET 17)
!       (7)   - LOV GRID ORIENTATION
!       (8)   - DX - X DIRECTION INCREMENT
!       (9)   - DY - Y DIRECTION INCREMENT
!       (10)  - PROJECTION CENTER FLAG
!       (11)  - SCANNING MODE (RIGHT ADJ COPY OF OCTET 28)
!    SPHERICAL HARMONIC COEFFICIENTS
!       (2)   - J PENTAGONAL RESOLUTION PARAMETER
!       (3)   - K      "          "         "
!       (4)   - M      "          "         "
!       (5)   - REPRESENTATION TYPE
!       (6)   - COEFFICIENT STORAGE MODE
!    MERCATOR GRIDS
!       (2)   - N(I) NR POINTS ON LATITUDE CIRCLE
!       (3)   - N(J) NR POINTS ON LONGITUDE MERIDIAN
!       (4)   - LA(1) LATITUDE OF ORIGIN
!       (5)   - LO(1) LONGITUDE OF ORIGIN
!       (6)   - RESOLUTION FLAG (RIGHT ADJ COPY OF OCTET 17)
!       (7)   - LA(2) LATITUDE OF LAST GRID POINT
!       (8)   - LO(2) LONGITUDE OF LAST GRID POINT
!       (9)   - LATIT - LATITUDE OF PROJECTION INTERSECTION
!       (10)  - RESERVED
!       (11)  - SCANNING MODE FLAG (RIGHT ADJ COPY OF OCTET 28)
!       (12)  - LONGITUDINAL DIR GRID LENGTH
!       (13)  - LATITUDINAL DIR GRID LENGTH
!    LAMBERT CONFORMAL GRIDS
!       (2)   - NX NR POINTS ALONG X-AXIS
!       (3)   - NY NR POINTS ALONG Y-AXIS
!       (4)   - LA1 LAT OF ORIGIN (LOWER LEFT)
!       (5)   - LO1 LON OF ORIGIN (LOWER LEFT)
!       (6)   - RESOLUTION (RIGHT ADJ COPY OF OCTET 17)
!       (7)   - LOV - ORIENTATION OF GRID
!       (8)   - DX - X-DIR INCREMENT
!       (9)   - DY - Y-DIR INCREMENT
!       (10)  - PROJECTION CENTER FLAG
!       (11)  - SCANNING MODE FLAG (RIGHT ADJ COPY OF OCTET 28)
!       (12)  - LATIN 1 - FIRST LAT FROM POLE OF SECANT CONE INTER
!       (13)  - LATIN 2 - SECOND LAT FROM POLE OF SECANT CONE INTER
!  KBMS       - BITMAP DESCRIBING LOCATION OF OUTPUT ELEMENTS.
!                         (ALWAYS CONSTRUCTED)
!  KPTR       - ARRAY CONTAINING STORAGE FOR FOLLOWING PARAMETERS
!       (1)   - TOTAL LENGTH OF GRIB MESSAGE
!       (2)   - LENGTH OF INDICATOR (SECTION  0)
!       (3)   - LENGTH OF PDS       (SECTION  1)
!       (4)   - LENGTH OF GDS       (SECTION  2)
!       (5)   - LENGTH OF BMS       (SECTION  3)
!       (6)   - LENGTH OF BDS       (SECTION  4)
!       (7)   - VALUE OF CURRENT BYTE
!       (8)   - BIT POINTER
!       (9)   - GRIB START BIT NR
!      (10)   - GRIB/GRID ELEMENT COUNT
!      (11)   - NR UNUSED BITS AT END OF SECTION 3
!      (12)   - BIT MAP FLAG (COPY OF BMS OCTETS 5,6)
!      (13)   - NR UNUSED BITS AT END OF SECTION 2
!      (14)   - BDS FLAGS (RIGHT ADJ COPY OF OCTET 4)
!      (15)   - NR UNUSED BITS AT END OF SECTION 4
!  KRET       - FLAG INDICATING QUALITY OF COMPLETION
!
! REMARKS: WHEN DECODING IS COMPLETED, DATA AT EACH GRID POINT HAS BEEN
!       RETURNED IN THE UNITS SPECIFIED IN THE GRIB MANUAL.
!
!       VALUES FOR RETURN FLAG (KRET)
!  KRET = 0 - NORMAL RETURN, NO ERRORS
!       = 1 - 'GRIB' NOT FOUND IN FIRST 100 CHARS
!       = 2 - '7777' NOT IN CORRECT LOCATION
!       = 3 - UNPACKED FIELD IS LARGER THAN 65160
!       = 4 - GDS/ GRID NOT ONE OF CURRENTLY ACCEPTED VALUES
!       = 5 - GRID NOT CURRENTLY AVAIL FOR CENTER INDICATED
!       = 8 - TEMP GDS INDICATED, BUT GDS FLAG IS OFF
!       = 9 - GDS INDICATES SIZE MISMATCH WITH STD GRID
!       =10 - INCORRECT CENTER INDICATOR
!       =11 - BINARY DATA SECTION (BDS) NOT COMPLETELY PROCESSED.
!             PROGRAM IS NOT SET TO PROCESS FLAG COMBINATIONS
!             SHOWN IN OCTETS 4 AND 14.
!       =12 - BINARY DATA SECTION (BDS) NOT COMPLETELY PROCESSED.
!             PROGRAM IS NOT SET TO PROCESS FLAG COMBINATIONS
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 77
!   MACHINE:  NAS9000
!
!$$$
!                                                      4 AUG 1988
!                            W3FI63
!
!
!                    GRIB UNPACKING ROUTINE
!
!
!    THIS ROUTINE WILL UNPACK A 'GRIB' FIELD TO THE EXACT GRID
!  TYPE SPECIFIED IN THE MESSAGE, RETURN A BIT MAP AND MAKE THE
!  VALUES OF THE PRODUCT DEFINITION SEC   (PDS) AND THE GRID
!  DESCRIPTION SEC   (GDS) AVAILABLE IN RETURN ARRAYS.
!  SEE "GRIB - THE WMO FORMAT FOR THE STORAGE OF WEATHER PRODUCT
!  INFORMATION AND THE EXCHANGE OF WEATHER PRODUCT MESSAGES IN
!  GRIDDED BINARY FORM" DATED JULY 1, 1988 BY JOHN D. STACKPOLE
!  DOC, NOAA, NWS, NATIONAL METEOROLOGICAL CENTER.
!
!    THE CALL TO THE GRIB UNPACKING ROUTINE IS AS FOLLOWS:
!
!         CALL W3FI63(MSGA,KPDS,KGDS,LBMS,DATA,KPTR,KRET)
!
!  INPUT:
!
!    MSGA  = CONTAINS THE GRIB MESSAGE TO BE UNPACKED. CHARACTERS
!            "GRIB" MAY BEGIN ANYWHERE WITHIN FIRST 100 BYTES.
!
!  OUTPUT:
!
!    KPDS(100)      INTEGER*4
!            ARRAY TO CONTAIN THE ELEMENTS OF THE PRODUCT
!            DEFINITION SEC  .
!      (VERSION 1)
!         KPDS(1)  - ID OF CENTER
!         KPDS(2)  - MODEL IDENTIFICATION (SEE "GRIB" TABLE 1)
!         KPDS(3)  - GRID IDENTIFICATION (SEE "GRIB" TABLE 2)
!         KPDS(4)  - GDS/BMS FLAG
!                        BIT       DEFINITION
!                         25        0 - GDS OMITTED
!                                   1 - GDS INCLUDED
!                         26        0 - BMS OMITTED
!                                   1 - BMS INCLUDED
!                     NOTE:- LEFTMOST BIT = 1,
!                            RIGHTMOST BIT = 32
!         KPDS(5)  - INDICATOR OF PARAMETER (SEE "GRIB" TABLE 5)
!         KPDS(6)  - TYPE OF LEVEL (SEE "GRIB" TABLES 6 & 7)
!         KPDS(7)  - HEIGHT,PRESSURE,ETC  OF LEVEL
!         KPDS(8)  - YEAR INCLUDING CENTURY
!         KPDS(9)  - MONTH OF YEAR
!         KPDS(10) - DAY OF MONTH
!         KPDS(11) - HOUR OF DAY
!         KPDS(12) - MINUTE OF HOUR
!         KPDS(13) - INDICATOR OF FORECAST TIME UNIT (SEE "GRIB"
!                    TABLE 8)
!         KPDS(14) - TIME 1               (SEE "GRIB" TABLE 8A)
!         KPDS(15) - TIME 2               (SEE "GRIB" TABLE 8A)
!         KPDS(16) - TIME RANGE INDICATOR (SEE "GRIB" TABLE 8A)
!         KPDS(17) - NUMBER INCLUDED IN AVERAGE
!         KPDS(18) - EDITION NR OF GRIB SPECIFICATION
!         KPDS(19) - VERSION NR OF PARAMETER TABLE
!
!    KGDS(13)       INTEGER*4
!          ARRAY CONTAINING GDS ELEMENTS.
!
!         KGDS(1)  - DATA REPRESENTATION TYPE
!
!      LATITUDE/LONGITUDE GRIDS (SEE "GRIB" TABLE 10)
!         KGDS(2)  - N(I) NUMBER OF POINTS ON LATITUDE
!                    CIRCLE
!         KGDS(3)  - N(J) NUMBER OF POINTS ON LONGITUDE
!                    CIRCLE
!         KGDS(4)  - LA(1) LATITUDE OF ORIGIN
!         KGDS(5)  - LO(1) LONGITUDE OF ORIGIN
!         KGDS(6)  - RESOLUTION FLAG
!                        BIT       MEANING
!                         25       0 - DIRECTION INCREMENTS NOT
!                                      GIVEN
!                                  1 - DIRECTION INCREMENTS GIVEN
!         KGDS(7)  - LA(2) LATITUDE OF EXTREME POINT
!         KGDS(8)  - LO(2) LONGITUDE OF EXTREME POINT
!         KGDS(9)  - DI LATITUDINAL DIRECTION INCREMENT
!         KGDS(10) - REGULAR LAT/LON GRID
!                        DJ - LONGITUDINAL DIRECTION
!                             INCREMENT
!                    GAUSSIAN GRID
!                        N  - NUMBER OF LATITUDE CIRCLES
!                             BETWEEN A POLE AND THE EQUATOR
!         KGDS(11) - SCANNING MODE FLAG
!                        BIT       MEANING
!                         25       0 - POINTS ALONG A LATITUDE
!                                      SCAN FROM WEST TO EAST
!                                  1 - POINTS ALONG A LATITUDE
!                                      SCAN FROM EAST TO WEST
!                         26       0 - POINTS ALONG A MERIDIAN
!                                      SCAN FROM NORTH TO SOUTH
!                                  1 - POINTS ALONG A MERIDIAN
!                                      SCAN FROM SOUTH TO NORTH
!                         27       0 - POINTS SCAN FIRST ALONG
!                                      CIRCLES OF LATITUDE, THEN
!                                      ALONG MERIDIANS
!                                      (FORTRAN: (I,J))
!                                  1 - POINTS SCAN FIRST ALONG
!                                      MERIDIANS THEN ALONG
!                                      CIRCLES OF LATITUDE
!                                      (FORTRAN: (J,I))
!
!      POLAR STEREOGRAPHIC GRIDS  (SEE GRIB TABLE 12)
!         KGDS(2)  - N(I) NR POINTS ALONG LAT CIRCLE
!         KGDS(3)  - N(J) NR POINTS ALONG LON CIRCLE
!         KGDS(4)  - LA(1) LATITUDE OF ORIGIN
!         KGDS(5)  - LO(1) LONGITUDE OF ORIGIN
!         KGDS(6)  - RESERVED
!         KGDS(7)  - LOV GRID ORIENTATION
!         KGDS(8)  - DX - X DIRECTION INCREMENT
!         KGDS(9)  - DY - Y DIRECTION INCREMENT
!         KGDS(10) - PROJECTION CENTER FLAG
!         KGDS(11) - SCANNING MODE
!
!      SPHERICAL HARMONIC COEFFICIENTS (SEE "GRIB" TABLE 14)
!         KGDS(2)  - J PENTAGONAL RESOLUTION PARAMETER
!         KGDS(3)  - K PENTAGONAL RESOLUTION PARAMETER
!         KGDS(4)  - M PENTAGONAL RESOLUTION PARAMETER
!         KGDS(5)  - REPRESENTATION TYPE
!         KGDS(6)  - COEFFICIENT STORAGE MODE
!
!    MERCATOR GRIDS
!         KGDS(2)   - N(I) NR POINTS ON LATITUDE CIRCLE
!         KGDS(3)   - N(J) NR POINTS ON LONGITUDE MERIDIAN
!         KGDS(4)   - LA(1) LATITUDE OF ORIGIN
!         KGDS(5)   - LO(1) LONGITUDE OF ORIGIN
!         KGDS(6)   - RESOLUTION FLAG
!         KGDS(7)   - LA(2) LATITUDE OF LAST GRID POINT
!         KGDS(8)   - LO(2) LONGITUDE OF LAST GRID POINT
!         KGDS(9)   - LATIN - LATITUDE OF PROJECTION INTERSECTION
!         KGDS(10)  - RESERVED
!         KGDS(11)  - SCANNING MODE FLAG
!         KGDS(12)  - LONGITUDINAL DIR GRID LENGTH
!         KGDS(13)  - LATITUDINAL DIR GRID LENGTH
!    LAMBERT CONFORMAL GRIDS
!         KGDS(2)   - NX NR POINTS ALONG X-AXIS
!         KGDS(3)   - NY NR POINTS ALONG Y-AXIS
!         KGDS(4)   - LA1 LAT OF ORIGIN (LOWER LEFT)
!         KGDS(5)   - LO1 LON OF ORIGIN (LOWER LEFT)
!         KGDS(6)   - RESOLUTION (RIGHT ADJ COPY OF OCTET 17)
!         KGDS(7)   - LOV - ORIENTATION OF GRID
!         KGDS(8)   - DX - X-DIR INCREMENT
!         KGDS(9)   - DY - Y-DIR INCREMENT
!         KGDS(10)  - PROJECTION CENTER FLAG
!         KGDS(11)  - SCANNING MODE FLAG
!         KGDS(12)  - LATIN 1 - FIRST LAT FROM POLE OF
!                     SECANT CONE INTERSECTION
!         KGDS(13)  - LATIN 2 - SECOND LAT FROM POLE OF
!                     SECANT CONE INTERSECTION
!
!    LBMS(65160)    LOGICAL
!            ARRAY TO CONTAIN THE BIT MAP DESCRIBING THE
!            PLACEMENT OF DATA IN THE OUTPUT ARRAY.  IF A
!            BIT MAP IS NOT INCLUDED IN THE SOURCE MESSAGE,
!            ONE WILL BE GENERATED AUTOMATICALLY BY THE
!            UNPACKING ROUTINE.
!
!
!    DATA(65160)    REAL*4
!            THIS ARRAY WILL CONTAIN THE UNPACKED DATA POINTS.
!
!                   NOTE:- 65160 IS MAXIMUN FIELD SIZE ALLOWABLE
!
!    KPTR(10)       INTEGER*4
!            ARRAY CONTAINING STORAGE FOR THE FOLLOWING
!            PARAMETERS.
!
!              (1)  -    UNUSED
!              (2)  -    UNUSED
!              (3)  -    LENGTH OF PDS (IN BYTES)
!              (4)  -    LENGTH OF GDS (IN BYTES)
!              (5)  -    LENGTH OF BMS (IN BYTES)
!              (6)  -    LENGTH OF BDS (IN BYTES)
!              (7)  -    USED BY UNPACKING ROUTINE
!              (8)  -    NUMBER OF DATA POINTS FOR GRID
!              (9)  -    "GRIB" CHARACTERS START IN BYTE NUMBER
!              (10) -    USED BY UNPACKING ROUTINE
!
!
!    KRET      INTEGER*4
!              THIS VARIABLE WILL CONTAIN THE RETURN INDICATOR.
!
!              0    -    NO ERRORS DETECTED.
!
!              1    -    'GRIB' NOT FOUND IN FIRST 100
!                        CHARACTERS.
!
!              2    -    '7777' NOT FOUND, EITHER MISSING OR
!                        TOTAL OF SEC   COUNTS OF INDIVIDUAL
!                        SECTIONS  IS INCORRECT.
!
!              3    -    UNPACKED FIELD IS LARGER THAN 65160.
!
!              4    -    IN GDS, DATA REPRESENTATION TYPE
!                        NOT ONE OF THE CURRENTLY ACCEPTABLE
!                        VALUES. SEE "GRIB" TABLE 9. VALUE
!                        OF INCORRECT TYPE RETURNED IN KGDS(1).
!
!              5    -    GRID INDICATED IN KPDS(3) IS NOT
!                        AVAILABLE FOR THE CENTER INDICATED IN
!                        KPDS(1) AND NO GDS SENT.
!
!              7    -    EDITION INDICATED IN KPDS(18) HAS NOT
!                        YET BEEN INCLUDED IN THE DECODER.
!
!              8    -    GRID IDENTIFICATION = 255 (NOT STANDARD
!                        GRID) BUT FLAG INDICATING PRESENCE OF
!                        GDS IS TURNED OFF. NO METHOD OF
!                        GENERATING PROPER GRID.
!
!              9    -    PRODUCT OF KGDS(2) AND KGDS(3) DOES NOT
!                        MATCH STANDARD NUMBER OF POINTS FOR THIS
!                        GRID (FOR OTHER THAN SPECTRALS). THIS
!                        WILL OCCUR ONLY IF THE GRID.
!                        IDENTIFICATION, KPDS(3), AND A
!                        TRANSMITTED GDS ARE INCONSISTENT.
!
!             10    -    CENTER INDICATOR WAS NOT ONE INDICATED
!                        IN "GRIB" TABLE 1.  PLEASE CONTACT AD
!                        PRODUCTION MANAGEMENT BRANCH (W/NMC42)
!                                  IF THIS ERROR IS ENCOUNTERED.
!
!             11    -    BINARY DATA SECTION (BDS) NOT COMPLETELY
!                        PROCESSED.  PROGRAM IS NOT SET TO PROCESS
!                        FLAG COMBINATIONS AS SHOWN IN
!                        OCTETS 4 AND 14.
!
!
!  LIST OF TEXT MESSAGES FROM CODE
!
!
!  W3FI63/FI632
!
!         'HAVE ENCOUNTERED A NEW GRID FOR NMC, PLEASE NOTIFY
!         AUTOMATION DIVISION, PRODUCTION MANAGEMENT BRANCH
!         (W/NMC42)'
!
!         'HAVE ENCOUNTERED A NEW GRID FOR ECMWF, PLEASE NOTIFY
!         AUTOMATION DIVISION, PRODUCTION MANAGEMENT BRANCH
!         (W/NMC42)'
!
!         'HAVE ENCOUNTERED A NEW GRID FOR U.K. METEOROLOGICAL
!         OFFICE, BRACKNELL.  PLEASE NOTIFY AUTOMATION DIVISION,
!         PRODUCTION MANAGEMENT BRANCH (W/NMC42)'
!
!         'HAVE ENCOUNTERED A NEW GRID FOR FNOC, PLEASE NOTIFY
!         AUTOMATION DIVISION, PRODUCTION MANAGEMENT BRANCH
!         (W/NMC42)'
!
!
!  W3FI63/FI633
!
!         'POLAR STEREO PROCESSING NOT AVAILABLE'  *
!
!  W3FI63/FI634
!
!         'WARNING - BIT MAP MAY NOT BE ASSOCIATED WITH SPHERICAL
!         COEFFICIENTS'
!
!
!  W3FI63/FI637
!
!         'NO CURRENT LISTING OF FNOC GRIDS'      *
!
!
!  * WILL BE AVAILABLE IN NEXT UPDATE
!  ***************************************************************
!
!                    INCOMING MESSAGE HOLDER
  CHARACTER :: msga(*)
!                    BIT MAP
  LOGICAL*1 :: kbms(*)
!
!                    ELEMENTS OF PRODUCT DESCRIPTION SEC   (PDS)
  INTEGER :: kpds(*)
!                    ELEMENTS OF GRID DESCRIPTION SEC   (PDS)
  INTEGER :: kgds(*)
!
!                    CONTAINER FOR GRIB GRID
  REAL :: DATA(*)
!
!                    ARRAY OF POINTERS AND COUNTERS
  INTEGER :: kptr(*)
!
!  *****************************************************************
  INTEGER :: kkk,jsgn,jexp,ifr,npts
  CHARACTER (LEN=1) :: kk(8)
  REAL :: realkk,fval1,fdiff1
  EQUIVALENCE   (kk(1),kkk)
!  *****************************************************************
!     1.0 LOCATE BEGINNING OF 'GRIB' MESSAGE
!          FIND 'GRIB' CHARACTERS
!     2.0  USE COUNTS IN EACH DESCRIPTION SEC   TO DETERMINE
!          IF '7777' IS IN PROPER PLACE.
!     3.0  PARSE PRODUCT DEFINITION SECTION.
!     4.0  PARSE GRID DESCRIPTION SEC   (IF INCLUDED)
!     5.0  PARSE BIT MAP SEC   (IF INCLUDED)
!     6.0  USING INFORMATION FROM PRODUCT DEFINITION, GRID
!               DESCRIPTION, AND BIT MAP SECTIONS.. EXTRACT
!               DATA AND PLACE INTO PROPER ARRAY.
!  *******************************************************************
!
!                   MAIN DRIVER
!
!  *******************************************************************
  SAVE
  kptr(10) = 0
!               SEE IF PROPER 'GRIB' KEY EXISTS, THEN
!               USING SEC   COUNTS, DETERMINE IF '7777'
!               IS IN THE PROPER LOCATION
!
  CALL fi631(msga,kptr,kpds,kret)
  IF(kret /= 0) THEN
    GO TO 900
  END IF
!  PRINT *,'FI631 KPTR',(KPTR(I),I=1,16)
!
!               PARSE PARAMETERS FROM PRODUCT DESCRIPTION SECTION
!
  CALL fi632(msga,kptr,kpds,kret)
  IF(kret /= 0) THEN
    GO TO 900
  END IF
!  PRINT *,'FI632 KPTR',(KPTR(I),I=1,16)
!
!               IF AVAILABLE, EXTRACT NEW GRID DESCRIPTION
!
  IF (IAND(kpds(4),128) /= 0) THEN
    CALL fi633(msga,kptr,kgds,kret)
    IF(kret /= 0) THEN
      GO TO 900
    END IF
!      PRINT *,'FI633 KPTR',(KPTR(I),I=1,16)
  END IF
!
!               EXTRACT OR GENERATE BIT MAP
!
  CALL fi634(msga,kptr,kpds,kgds,kbms,kret)
  IF(kret /= 0) THEN
    GO TO 900
  END IF
!  PRINT *,'FI634 KPTR',(KPTR(I),I=1,16)
!
!               USING INFORMATION FROM PDS, BMS AND BIT DATA SEC  ,
!               EXTRACT AND SAVE IN GRIB GRID, ALL DATA ENTRIES.
!
  IF (kpds(18) == 1) THEN
    CALL fi635(msga,kptr,kpds,kgds,kbms,DATA,kret)
    IF (kptr(3) >= 256) THEN    ! was 50
      PRINT *, 'SPECIAL HANDLING FOR LARGE PDS LENGTH: ', kptr(3)
!
!                  PDS GREATER THAN 28 BYTES
!                     THEREFORE SOMETHING SPECIAL IS GOING ON
!
!                     IN THIS CASE 2ND DIFFERENCE PACKING
!                             NEEDS TO BE UNDONE.
!
!                EXTRACT FIRST VALUE FROM BYTE 41-44 PDS
!                           KPTR(9) CONTAINS OFFSET TO START OF
!                           GRIB MESSAGE.
!                EXTRACT FIRST FIRST-DIFFERENCE FROM BYTES 45-48 PDS
!
!               AND EXTRACT SCALE FACTOR (E) TO UNDO 2**E
!               THAT WAS APPLIED PRIOR TO 2ND ORDER PACKING
!               AND PLACED IN PDS BYTES 49-51
!               FACTOR IS A SIGNED TWO BYTE INTEGER
!
!               ALSO NEED THE DECIMAL SCALING FROM PDS(27-28)
!               (AVAILABLE IN KPDS(22) FROM UNPACKER)
!               TO UNDO THE DECIMAL SCALING APPLIED TO THE
!               SECOND DIFFERENCES DURING UNPACKING.
!               SECOND DIFFS ALWAYS PACKED WITH 0 DECIMAL SCALE
!               BUT UNPACKER DOESNT KNOW THAT.
!
!          CALL GRBGBYTE  (MSGA,FVAL1,KPTR(9)+384,32)
!
!      NOTE INTEGERS, CHARACTERS AND EQUIVALENCES
!      DEFINED ABOVE TO MAKE THIS KKK EXTRACTION
!      WORK AND LINE UP ON WORD BOUNDARIES
!
      CALL grbgbyte (msga,kkk,kptr(9)+384,32)
!
!    THE NEXT CODE WILL CONVERT THE IBM370 FOATING POINT
!    TO THE FLOATING POINT USED ON YOUR MACHINE.
!
!    1ST TEST TO SEE IN ON 32 OR 64 BIT WORD MACHINE
!    LW = 4 OR 8; IF 8 MAY BE A CRAY
!
      CALL w3fi01(lw)
      IF (lw == 4) THEN
        CALL grbgbyte (kk,jsgn,0,1)
        CALL grbgbyte (kk,jexp,1,7)
        CALL grbgbyte (kk,ifr,8,24)
      ELSE
        CALL grbgbyte (kk,jsgn,32,1)
        CALL grbgbyte (kk,jexp,33,7)
        CALL grbgbyte (kk,ifr,40,24)
      END IF
!
      IF (ifr == 0) THEN
        realkk = 0.0
      ELSE IF (jexp == 0.AND.ifr == 0) THEN
        realkk = 0.0
      ELSE
        realkk = FLOAT(ifr) * 16.0 ** (jexp - 64 - 6)
        IF (jsgn /= 0) realkk = -realkk
      END IF
      fval1 = realkk
!
!          CALL GRBGBYTE  (MSGA,FDIFF1,KPTR(9)+416,32)
!       (REPLACED BY FOLLOWING EXTRACTION)
!
      CALL grbgbyte (msga,kkk,kptr(9)+416,32)
!
!    THE NEXT CODE WILL CONVERT THE IBM370 FOATING POINT
!    TO THE FLOATING POINT USED ON YOUR MACHINE.
!
!    1ST TEST TO SEE IN ON 32 OR 64 BIT WORD MACHINE
!    LW = 4 OR 8; IF 8 MAY BE A CRAY
!
      CALL w3fi01(lw)
      IF (lw == 4) THEN
        CALL grbgbyte (kk,jsgn,0,1)
        CALL grbgbyte (kk,jexp,1,7)
        CALL grbgbyte (kk,ifr,8,24)
      ELSE
        CALL grbgbyte (kk,jsgn,32,1)
        CALL grbgbyte (kk,jexp,33,7)
        CALL grbgbyte (kk,ifr,40,24)
      END IF
!
      IF (ifr == 0) THEN
        realkk = 0.0
      ELSE IF (jexp == 0.AND.ifr == 0) THEN
        realkk = 0.0
      ELSE
        realkk = FLOAT(ifr) * 16.0 ** (jexp - 64 - 6)
        IF (jsgn /= 0) realkk = -realkk
      END IF
      fdiff1 = realkk
!
      CALL grbgbyte  (msga,ISIGN,kptr(9)+448,1)
      CALL grbgbyte  (msga,iscal2,kptr(9)+449,15)
      IF(ISIGN > 0) THEN
        iscal2 = - iscal2
      END IF
!          PRINT *,'DELTA POINT 1-',FVAL1
!          PRINT *,'DELTA POINT 2-',FDIFF1
!          PRINT *,'DELTA POINT 3-',ISCAL2
      npts  = kptr(10)
!          WRITE (6,FMT='(''  2ND DIFF POINTS IN FIELD = '',/,
!    &         10(3X,10F12.2,/))') (DATA(I),I=1,NPTS)
!          PRINT *,'DELTA POINT 4-',KPDS(22)
      CALL w3fi83 (DATA,npts,fval1,fdiff1,                              &
                          iscal2,kpds(22),kpds,kgds)
!          WRITE (6,FMT='(''  2ND DIFF EXPANDED POINTS IN FIELD = '',
!    &            /,10(3X,10F12.2,/))') (DATA(I),I=1,NPTS)
!          WRITE (6,FMT='(''  END OF ARRAY IN FIELD = '',/,
!    &         10(3X,10F12.2,/))') (DATA(I),I=NPTS-5,NPTS)
    END IF
  ELSE
  PRINT *,'FI635 NOT PROGRAMMED FOR EDITION NR',kpds(18)
    kret   = 7
  END IF
!
  900 RETURN
END SUBROUTINE w3fi63

SUBROUTINE fi631(msga,kptr,kpds,kret)
!$$$  SUBPROGRAM DOCUMENTATION  BLOCK
!             .      .    .                                       .
! SUBPROGRAM:    FI631       FIND 'GRIB' CHARS & RESET POINTERS
!   PRGMMR: BILL CAVANAUGH   ORG: W/NMC42    DATE: 91-09-13
!
! ABSTRACT: FIND 'GRIB; CHARACTERS AND SET POINTERS TO THE NEXT
!   BYTE FOLLOWING 'GRIB'. IF THEY EXIST EXTRACT COUNTS FROM GDS AND
!   BMS. EXTRACT COUNT FROM BDS. DETERMINE IF SUM OF COUNTS ACTUALLY
!   PLACES TERMINATOR '7777' AT THE CORRECT LOCATION.
!
! PROGRAM HISTORY LOG:
!   91-09-13  CAVANAUGH
!
! USAGE:    CALL FI631(MSGA,KPTR,KPDS,KRET)
!   INPUT ARGUMENT LIST:
!  MSGA       - GRIB FIELD - "GRIB" THRU "7777"
!  KPTR       - ARRAY CONTAINING STORAGE FOR FOLLOWING PARAMETERS
!       (1)   - TOTAL LENGTH OF GRIB MESSAGE
!       (2)   - LENGTH OF INDICATOR (SECTION  0)
!       (3)   - LENGTH OF PDS       (SECTION  1)
!       (4)   - LENGTH OF GDS       (SECTION  2)
!       (5)   - LENGTH OF BMS       (SECTION  3)
!       (6)   - LENGTH OF BDS       (SECTION  4)
!       (7)   - VALUE OF CURRENT BYTE
!       (8)   - BIT POINTER
!       (9)   - GRIB START BIT NR
!      (10)   - GRIB/GRID ELEMENT COUNT
!      (11)   - NR UNUSED BITS AT END OF SECTION 3
!      (12)   - BIT MAP FLAG
!      (13)   - NR UNUSED BITS AT END OF SECTION 2
!      (14)   - BDS FLAGS
!      (15)   - NR UNUSED BITS AT END OF SECTION 4
!
!   OUTPUT ARGUMENT LIST:      (INCLUDING WORK ARRAYS)
!  KPDS     - ARRAY CONTAINING PDS ELEMENTS.
!       (1)   - ID OF CENTER
!       (2)   - MODEL IDENTIFICATION
!       (3)   - GRID IDENTIFICATION
!       (4)   - GDS/BMS FLAG
!       (5)   - INDICATOR OF PARAMETER
!       (6)   - TYPE OF LEVEL
!       (7)   - HEIGHT/PRESSURE , ETC OF LEVEL
!       (8)   - YEAR OF CENTURY
!       (9)   - MONTH OF YEAR
!       (10)  - DAY OF MONTH
!       (11)  - HOUR OF DAY
!       (12)  - MINUTE OF HOUR
!       (13)  - INDICATOR OF FORECAST TIME UNIT
!       (14)  - TIME RANGE 1
!       (15)  - TIME RANGE 2
!       (16)  - TIME RANGE FLAG
!       (17)  - NUMBER INCLUDED IN AVERAGE
!  KPTR       - SEE INPUT LIST
!  KRET       - ERROR RETURN
!
! REMARKS:
!  ERROR RETURNS
!  KRET  = 1  -  NO 'GRIB'
!          2  -  NO '7777' OR MISLOCATED (BY COUNTS)
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 77
!   MACHINE:  NAS9000
!
!$$$
!
!                    INCOMING MESSAGE HOLDER
  CHARACTER :: msga(*)
!                    ARRAY OF POINTERS AND COUNTERS
  INTEGER :: kptr(*)
!                    PRODUCT DESCRIPTION SECTION DATA.
  INTEGER :: kpds(*)
!
  INTEGER :: kret
!
!  ******************************************************************
  SAVE
  kret = 0
!  -------------------  FIND 'GRIB' KEY
  DO i = 0, 839, 8
    CALL grbgbyte (msga,mgrib,i,32)
    IF (mgrib == 1196575042) THEN
      kptr(9)   = i
      GO TO 60
    END IF
  END DO
  kret  = 1
  RETURN
  60 CONTINUE
!  -------------FOUND 'GRIB'
!                     SKIP GRIB CHARACTERS
!  PRINT *,'FI631 GRIB AT',I
  kptr(8)   = kptr(9) + 32
  CALL grbgbyte (msga,itotal,kptr(8),24)
!                 HAVE LIFTED WHAT MAY BE A MSG TOTAL BYTE COUNT
  ipoint    = kptr(9) + itotal * 8 - 32
  CALL grbgbyte (msga,i7777,ipoint,32)
  IF (i7777 == 926365495) THEN
!              HAVE FOUND END OF MESSAGE '7777' IN PROPER LOCATION
!              MARK AND PROCESS AS GRIB VERSION 1 OR HIGHER
!      PRINT *,'FI631 7777 AT',IPOINT
    kptr(8)   = kptr(8) + 24
    kptr(1)   = itotal
    kptr(2)   = 8
    CALL grbgbyte (msga,kpds(18),kptr(8),8)
    kptr(8)   = kptr(8) + 8
  ELSE
!              CANNOT FIND END OF GRIB EDITION 1 MESSAGE
    kret      = 2
    RETURN
  END IF
!  -------------------  PROCESS SECTION 1
!                EXTRACT COUNT FROM PDS
!  PRINT *,'START OF PDS',KPTR(8)
  CALL grbgbyte (msga,kptr(3),kptr(8),24)
  look      = kptr(8) + 56
!                EXTRACT GDS/BMS FLAG
  CALL grbgbyte (msga,kpds(4),look,8)
  kptr(8)   = kptr(8) + kptr(3) * 8
!  PRINT *,'START OF GDS',KPTR(8)
  IF (IAND(kpds(4),128) /= 0) THEN
!                EXTRACT COUNT FROM GDS
    CALL grbgbyte (msga,kptr(4),kptr(8),24)
    kptr(8)   = kptr(8) + kptr(4) * 8
  ELSE
    kptr(4)   = 0
  END IF
!  PRINT *,'START OF BMS',KPTR(8)
  IF (IAND(kpds(4),64) /= 0) THEN
!                EXTRACT COUNT FROM BMS
    CALL grbgbyte (msga,kptr(5),kptr(8),24)
  ELSE
    kptr(5)   = 0
  END IF
  kptr(8)   = kptr(8) + kptr(5) * 8
!  PRINT *,'START OF BDS',KPTR(8)
!                EXTRACT COUNT FROM BDS
  CALL grbgbyte (msga,kptr(6),kptr(8),24)
!  ---------------  TEST FOR '7777'
!  PRINT *,(KPTR(KJ),KJ=1,10)
  kptr(8)   = kptr(8) + kptr(6) * 8
!                EXTRACT FOUR BYTES FROM THIS LOCATION
!  PRINT *,'FI631 LOOKING FOR 7777 AT',KPTR(8)
  CALL grbgbyte (msga,k7777,kptr(8),32)
  match  = kptr(2) + kptr(3) + kptr(4) + kptr(5) + kptr(6) + 4
  IF (k7777 /= 926365495.OR.match /= kptr(1)) THEN
    kret  = 2
  ELSE
!      PRINT *,'FI631 7777 AT',KPTR(8)
    IF (kpds(18) == 0) THEN
      kptr(1)  = kptr(2) + kptr(3) + kptr(4) + kptr(5) +                &
              kptr(6) + 4
    END IF
  END IF
!  PRINT *,'KPTR',(KPTR(I),I=1,16)
  RETURN
END SUBROUTINE fi631

SUBROUTINE fi632(msga,kptr,kpds,kret)
!$$$  SUBPROGRAM DOCUMENTATION  BLOCK
!             .      .    .                                       .
! SUBPROGRAM:    FI632       GATHER INFO FROM PRODUCT DEFINITION SEC
!   PRGMMR: BILL CAVANAUGH   ORG: W/NMC42    DATE: 91-09-13
!
! ABSTRACT: EXTRACT INFORMATION FROM THE PRODUCT DESCRIPTION
!   SEC  , AND GENERATE LABEL INFORMATION TO PERMIT STORAGE
!   IN OFFICE NOTE 84 FORMAT.
!
! PROGRAM HISTORY LOG:
!   91-09-13  CAVANAUGH
!   93-12-08  CAVANAUGH   CORRECTED TEST FOR EDITION NUMBER INSTEAD
!                      OF VERSION NUMBER
!
! USAGE:    CALL FI632(MSGA,KPTR,KPDS,KRET)
!   INPUT ARGUMENT LIST:
!  MSGA      - ARRAY CONTAINING GRIB MESSAGE
!  KPTR       - ARRAY CONTAINING STORAGE FOR FOLLOWING PARAMETERS
!       (1)   - TOTAL LENGTH OF GRIB MESSAGE
!       (2)   - LENGTH OF INDICATOR (SECTION  0)
!       (3)   - LENGTH OF PDS       (SECTION  1)
!       (4)   - LENGTH OF GDS       (SECTION  2)
!       (5)   - LENGTH OF BMS       (SECTION  3)
!       (6)   - LENGTH OF BDS       (SECTION  4)
!       (7)   - VALUE OF CURRENT BYTE
!       (8)   - BIT POINTER
!       (9)   - GRIB START BIT NR
!      (10)   - GRIB/GRID ELEMENT COUNT
!      (11)   - NR UNUSED BITS AT END OF SECTION 3
!      (12)   - BIT MAP FLAG
!      (13)   - NR UNUSED BITS AT END OF SECTION 2
!      (14)   - BDS FLAGS
!      (15)   - NR UNUSED BITS AT END OF SECTION 4
!
!   OUTPUT ARGUMENT LIST:      (INCLUDING WORK ARRAYS)
!  KPDS     - ARRAY CONTAINING PDS ELEMENTS.
!       (1)   - ID OF CENTER
!       (2)   - MODEL IDENTIFICATION
!       (3)   - GRID IDENTIFICATION
!       (4)   - GDS/BMS FLAG
!       (5)   - INDICATOR OF PARAMETER
!       (6)   - TYPE OF LEVEL
!       (7)   - HEIGHT/PRESSURE , ETC OF LEVEL
!       (8)   - YEAR OF CENTURY
!       (9)   - MONTH OF YEAR
!       (10)  - DAY OF MONTH
!       (11)  - HOUR OF DAY
!       (12)  - MINUTE OF HOUR
!       (13)  - INDICATOR OF FORECAST TIME UNIT
!       (14)  - TIME RANGE 1
!       (15)  - TIME RANGE 2
!       (16)  - TIME RANGE FLAG
!       (17)  - NUMBER INCLUDED IN AVERAGE
!       (18)  -
!       (19)  -
!       (20)  - NUMBER MISSING FROM AVGS/ACCUMULATIONS
!       (21)  - CENTURY
!       (22)  - UNITS DECIMAL SCALE FACTOR
!       (23)  - SUBCENTER
!  KPTR       - ARRAY CONTAINING STORAGE FOR FOLLOWING PARAMETERS
!               SEE INPUT LIST
!  KRET   - ERROR RETURN
!
! REMARKS:
!     ERROR RETURN = 0 - NO ERRORS
!                  = 8 - TEMP GDS INDICATED, BUT NO GDS
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 77
!   MACHINE:  NAS9000
!
!$$$
!
!                    INCOMING MESSAGE HOLDER
  CHARACTER :: msga(*)
!
!                    ARRAY OF POINTERS AND COUNTERS
  INTEGER :: kptr(*)
!                    PRODUCT DESCRIPTION SECTION ENTRIES
  INTEGER :: kpds(*)
!
  INTEGER :: kret
!  -------------------  PROCESS SECTION 1
  SAVE
  kptr(8)  = kptr(9) + kptr(2) * 8 + 24
!  BYTE 4
!                PARAMETER TABLE VERSION NR
  CALL grbgbyte (msga,kpds(19),kptr(8),8)
  kptr(8)   = kptr(8) + 8
!  BYTE 5           IDENTIFICATION OF CENTER
  CALL grbgbyte (msga,kpds(1),kptr(8),8)
  kptr(8)   = kptr(8) + 8
!  BYTE 6
!                    GET GENERATING PROCESS ID NR
  CALL grbgbyte (msga,kpds(2),kptr(8),8)
  kptr(8)   = kptr(8) + 8
!  BYTE 7
!                   GRID DEFINITION
  CALL grbgbyte (msga,kpds(3),kptr(8),8)
  kptr(8)   = kptr(8) + 8
!  BYTE 8
!                   GDS/BMS FLAGS
!  CALL GRBGBYTE (MSGA,KPDS(4),KPTR(8),8)
  kptr(8)   = kptr(8) + 8
!  BYTE 9
!                   INDICATOR OF PARAMETER
  CALL grbgbyte (msga,kpds(5),kptr(8),8)
  kptr(8)   = kptr(8) + 8
!  BYTE 10
!                   TYPE OF LEVEL
  CALL grbgbyte (msga,kpds(6),kptr(8),8)
  kptr(8)   = kptr(8) + 8
!  BYTE 11,12
!                   HEIGHT/PRESSURE
  CALL grbgbyte (msga,kpds(7),kptr(8),16)
  kptr(8)   = kptr(8) + 16
!  BYTE 13
!                   YEAR OF CENTURY
  CALL grbgbyte (msga,kpds(8),kptr(8),8)
  kptr(8)   = kptr(8) + 8
!  BYTE 14
!                   MONTH OF YEAR
  CALL grbgbyte (msga,kpds(9),kptr(8),8)
  kptr(8)   = kptr(8) + 8
!  BYTE 15
!                   DAY OF MONTH
  CALL grbgbyte (msga,kpds(10),kptr(8),8)
  kptr(8)   = kptr(8) + 8
!  BYTE 16
!                   HOUR OF DAY
  CALL grbgbyte (msga,kpds(11),kptr(8),8)
  kptr(8)   = kptr(8) + 8
!  BYTE 17
!                   MINUTE
  CALL grbgbyte (msga,kpds(12),kptr(8),8)
  kptr(8)   = kptr(8) + 8
!  BYTE 18
!                   INDICATOR TIME UNIT RANGE
  CALL grbgbyte (msga,kpds(13),kptr(8),8)
  kptr(8)   = kptr(8) + 8
!  BYTE 19
!                   P1 - PERIOD OF TIME
  CALL grbgbyte (msga,kpds(14),kptr(8),8)
  kptr(8)   = kptr(8) + 8
!  BYTE 20
!                   P2 - PERIOD OF TIME
  CALL grbgbyte (msga,kpds(15),kptr(8),8)
  kptr(8)   = kptr(8) + 8
!  BYTE 21
!                   TIME RANGE INDICATOR
  CALL grbgbyte (msga,kpds(16),kptr(8),8)
  kptr(8)   = kptr(8) + 8
!
!  IF TIME RANGE INDICATOR IS 10, P1 IS PACKED IN
!  PDS BYTES 19-20
!
  IF (kpds(16) == 10) THEN
    kpds(14)  = kpds(14) * 256 + kpds(15)
    kpds(15)  = 0
  END IF
!  BYTE 22,23
!                   NUMBER INCLUDED IN AVERAGE
  CALL grbgbyte (msga,kpds(17),kptr(8),16)
  kptr(8)   = kptr(8) + 16
!  BYTE 24
!                   NUMBER MISSING FROM AVERAGES/ACCUMULATIONS
  CALL grbgbyte (msga,kpds(20),kptr(8),8)
  kptr(8)   = kptr(8) + 8
!  BYTE 25
!                   IDENTIFICATION OF CENTURY
  CALL grbgbyte (msga,kpds(21),kptr(8),8)
  kptr(8)   = kptr(8) + 8
  IF (kptr(3) > 25) THEN
!  BYTE 26              SUB CENTER NUMBER
    CALL grbgbyte (msga,kpds(23),kptr(8),8)
    kptr(8)   = kptr(8) + 8
    IF (kptr(3) >= 28) THEN
!  BYTE 27-28
!                       UNITS DECIMAL SCALE FACTOR
      CALL grbgbyte (msga,ISIGN,kptr(8),1)
      kptr(8)  = kptr(8) + 1
      CALL grbgbyte (msga,idec,kptr(8),15)
      kptr(8)  = kptr(8) + 15
      IF (ISIGN > 0) THEN
        kpds(22)  = - idec
      ELSE
        kpds(22)  = idec
      END IF
      isiz  = kptr(3) - 28
      IF (isiz <= 12) THEN
!  BYTE  29
        CALL grbgbyte (MSGA,KPDS(24),KPTR(8)+8,8)
!  BYTE  30
        CALL grbgbyte (MSGA,KPDS(25),KPTR(8)+16,8)
!  BYTES 31-40                  CURRENTLY RESERVED FOR FUTURE USE
        KPTR(8)  = KPTR(8) + ISIZ * 8
      ELSE
!  BYTE  29
        CALL grbgbyte (MSGA,KPDS(24),KPTR(8)+8,8)
!  BYTE  30
        CALL grbgbyte (MSGA,KPDS(25),KPTR(8)+16,8)
!  BYTES 31-40                  CURRENTLY RESERVED FOR FUTURE USE
        KPTR(8)  = KPTR(8) + 12 * 8
!  BYTES 41 - N                 LOCAL USE DATA
        CALL W3FI01(LW)
        !MWDBIT  = LW * 8
        MWDBIT  = bit_size(KPDS)
        ISIZ    = KPTR(3) - 40
        ITER    = ISIZ / LW
        IF (MOD(ISIZ,LW).NE.0) ITER = ITER + 1
        CALL GBYTESC(MSGA,KPDS(36),KPTR(8),MWDBIT,0,ITER)
        KPTR(8)  = KPTR(8) + ISIZ * 8
!
!!  BYTES 29-40                  CURRENTLY RESERVED FOR FUTURE USE
!        kptr(8)  = kptr(8) + isiz * 8
!      ELSE
!!  BYTES 29-40                  CURRENTLY RESERVED FOR FUTURE USE
!        kptr(8)  = kptr(8) + 12 * 8
!!  BYTES 40 - N                 LOCAL USE DATA
!        CALL w3fi01(lw)
!        mwdbit  = lw * 8
!        isiz  = kptr(3) - 40
!        iter  = isiz / lw
!!FIXME RLC 2004-01-08
!! Causes array out of bounds in grbgbytes if PDS section > 28 bytes
!        CALL grbgbytes (msga,kpds(36),kptr(8),mwdbit,0,iter)
!        kptr(8)  = kptr(8) + isiz * 8
      END IF
    END IF
  END IF
!  ----------- TEST FOR NEW GRID
  IF (IAND(kpds(4),128) /= 0) THEN
    IF (IAND(kpds(4),64) /= 0) THEN
      IF (kpds(3) /= 255) THEN
        IF (kpds(3) >= 21.AND.kpds(3) <= 26)THEN
          RETURN
        ELSE IF (kpds(3) >= 37.AND.kpds(3) <= 44)THEN
          RETURN
        ELSE IF (kpds(3) >= 61.AND.kpds(3) <= 64) THEN
          RETURN
        END IF
        IF (kpds(1) == 7) THEN
          IF (kpds(3) >= 2 .AND. kpds(3) <= 6)THEN
          !ELSE IF (kpds(3) >= 5.AND.kpds(3) <= 6)THEN
          ELSE IF (kpds(3) >= 27.AND.kpds(3) <= 34)THEN
          ELSE IF (kpds(3) == 50) THEN
          ELSE IF (kpds(3) >= 70.AND.kpds(3) <= 77) THEN
          ELSE IF (kpds(3) >= 100.AND.kpds(3) <= 105) THEN
          ELSE IF (kpds(3) >= 201.AND.kpds(3) <= 214) THEN
          ELSE IF (kpds(3) == 218) THEN   ! added 218
          ELSE IF (kpds(3) == 221) THEN   ! added 221
          ELSE
            PRINT *,' HAVE ENCOUNTERED A NEW GRID ',kpds(3),' FOR',     &
                ' NMC WITHOUT A GRID DESCRIPTION SECTION'
            PRINT *,' PLEASE NOTIFY AUTOMATION DIVISION'
            PRINT *,' PRODUCTION MANAGEMENT BRANCH'
            PRINT *,' W/NMC42)'
          END IF
        ELSE IF (kpds(1) == 98) THEN
          IF (kpds(3) >= 1.AND.kpds(3) <= 16) THEN
          ELSE
            PRINT *,' HAVE ENCOUNTERED A NEW GRID FOR',                 &
                ' ECMWF WITHOUT A GRID DESCRIPTION SECTION'
            PRINT *,' PLEASE NOTIFY AUTOMATION DIVISION'
            PRINT *,' PRODUCTION MANAGEMENT BRANCH'
            PRINT *,' W/NMC42)'
          END IF
        ELSE IF (kpds(1) == 74) THEN
          IF (kpds(3) >= 1.AND.kpds(3) <= 12) THEN
          ELSE IF (kpds(3) >= 21.AND.kpds(3) <= 26)THEN
          ELSE IF (kpds(3) >= 61.AND.kpds(3) <= 64) THEN
          ELSE IF (kpds(3) >= 70.AND.kpds(3) <= 77) THEN
          ELSE
            PRINT *,' HAVE ENCOUNTERED A NEW GRID FOR',                 &
                    ' U.K. MET OFFICE, BRACKNELL',                      &
                    ' WITHOUT A GRID DESCRIPTION SECTION'
            PRINT *,' PLEASE NOTIFY AUTOMATION DIVISION'
            PRINT *,' PRODUCTION MANAGEMENT BRANCH'
            PRINT *,' W/NMC42)'
          END IF
        ELSE IF (kpds(1) == 58) THEN
          IF (kpds(3) >= 1.AND.kpds(3) <= 12) THEN
          ELSE
            PRINT *,' HAVE ENCOUNTERED A NEW GRID FOR',                 &
                ' FNOC WITHOUT A GRID DESCRIPTION SECTION'
            PRINT *,' PLEASE NOTIFY AUTOMATION DIVISION'
            PRINT *,' PRODUCTION MANAGEMENT BRANCH'
            PRINT *,' W/NMC42)'
          END IF
        END IF
      END IF
    END IF
  END IF
  RETURN
END SUBROUTINE fi632

SUBROUTINE fi633(msga,kptr,kgds,kret)
!$$$  SUBPROGRAM DOCUMENTATION  BLOCK
!             .      .    .                                       .
! SUBPROGRAM:    FI633       EXTRACT INFO FROM GRIB-GDS
!   PRGMMR: BILL CAVANAUGH   ORG: W/NMC42    DATE: 91-09-13
!
! ABSTRACT: EXTRACT INFORMATION ON UNLISTED GRID TO ALLOW
!   CONVERSION TO OFFICE NOTE 84 FORMAT.
!
! PROGRAM HISTORY LOG:
!   91-09-13  CAVANAUGH
!
! USAGE:    CALL FI633(MSGA,KPTR,KGDS,KRET)
!   INPUT ARGUMENT LIST:
!  MSGA      - ARRAY CONTAINING GRIB MESSAGE
!  KPTR       - ARRAY CONTAINING STORAGE FOR FOLLOWING PARAMETERS
!       (1)   - TOTAL LENGTH OF GRIB MESSAGE
!       (2)   - LENGTH OF INDICATOR (SECTION  0)
!       (3)   - LENGTH OF PDS       (SECTION  1)
!       (4)   - LENGTH OF GDS       (SECTION  2)
!       (5)   - LENGTH OF BMS       (SECTION  3)
!       (6)   - LENGTH OF BDS       (SECTION  4)
!       (7)   - VALUE OF CURRENT BYTE
!       (8)   - BIT POINTER
!       (9)   - GRIB START BIT NR
!      (10)   - GRIB/GRID ELEMENT COUNT
!      (11)   - NR UNUSED BITS AT END OF SECTION 3
!      (12)   - BIT MAP FLAG
!      (13)   - NR UNUSED BITS AT END OF SECTION 2
!      (14)   - BDS FLAGS
!      (15)   - NR UNUSED BITS AT END OF SECTION 4
!
!   OUTPUT ARGUMENT LIST:      (INCLUDING WORK ARRAYS)
!  KGDS     - ARRAY CONTAINING GDS ELEMENTS.
!       (1)   - DATA REPRESENTATION TYPE
!       (19)  - NUMBER OF VERTICAL COORDINATE PARAMETERS
!       (20)  - OCTET NUMBER OF THE LIST OF VERTICAL COORDINATE
!               PARAMETERS
!               OR
!               OCTET NUMBER OF THE LIST OF NUMBERS OF POINTS
!               IN EACH ROW
!               OR
!               255 IF NEITHER ARE PRESENT
!       (21)  - FOR GRIDS WITH PL, NUMBER OF POINTS IN GRID
!       (22)  - NUMBER OF WORDS IN EACH ROW
!    LATITUDE/LONGITUDE GRIDS
!       (2)   - N(I) NR POINTS ON LATITUDE CIRCLE
!       (3)   - N(J) NR POINTS ON LONGITUDE MERIDIAN
!       (4)   - LA(1) LATITUDE OF ORIGIN
!       (5)   - LO(1) LONGITUDE OF ORIGIN
!       (6)   - RESOLUTION FLAG
!       (7)   - LA(2) LATITUDE OF EXTREME POINT
!       (8)   - LO(2) LONGITUDE OF EXTREME POINT
!       (9)   - DI LATITUDINAL DIRECTION OF INCREMENT
!       (10)  - DJ LONGITUDINAL DIRECTION INCREMENT
!       (11)  - SCANNING MODE FLAG
!    POLAR STEREOGRAPHIC GRIDS
!       (2)   - N(I) NR POINTS ALONG LAT CIRCLE
!       (3)   - N(J) NR POINTS ALONG LON CIRCLE
!       (4)   - LA(1) LATITUDE OF ORIGIN
!       (5)   - LO(1) LONGITUDE OF ORIGIN
!       (6)   - RESERVED
!       (7)   - LOV GRID ORIENTATION
!       (8)   - DX - X DIRECTION INCREMENT
!       (9)   - DY - Y DIRECTION INCREMENT
!       (10)  - PROJECTION CENTER FLAG
!       (11)  - SCANNING MODE
!    SPHERICAL HARMONIC COEFFICIENTS
!       (2)   - J PENTAGONAL RESOLUTION PARAMETER
!       (3)   - K      "          "         "
!       (4)   - M      "          "         "
!       (5)   - REPRESENTATION TYPE
!       (6)   - COEFFICIENT STORAGE MODE
!    MERCATOR GRIDS
!       (2)   - N(I) NR POINTS ON LATITUDE CIRCLE
!       (3)   - N(J) NR POINTS ON LONGITUDE MERIDIAN
!       (4)   - LA(1) LATITUDE OF ORIGIN
!       (5)   - LO(1) LONGITUDE OF ORIGIN
!       (6)   - RESOLUTION FLAG
!       (7)   - LA(2) LATITUDE OF LAST GRID POINT
!       (8)   - LO(2) LONGITUDE OF LAST GRID POINT
!       (9)   - LATIN - LATITUDE OF PROJECTION INTERSECTION
!       (10)  - RESERVED
!       (11)  - SCANNING MODE FLAG
!       (12)  - LONGITUDINAL DIR GRID LENGTH
!       (13)  - LATITUDINAL DIR GRID LENGTH
!    LAMBERT CONFORMAL GRIDS
!       (2)   - NX NR POINTS ALONG X-AXIS
!       (3)   - NY NR POINTS ALONG Y-AXIS
!       (4)   - LA1 LAT OF ORIGIN (LOWER LEFT)
!       (5)   - LO1 LON OF ORIGIN (LOWER LEFT)
!       (6)   - RESOLUTION (RIGHT ADJ COPY OF OCTET 17)
!       (7)   - LOV - ORIENTATION OF GRID
!       (8)   - DX - X-DIR INCREMENT
!       (9)   - DY - Y-DIR INCREMENT
!       (10)  - PROJECTION CENTER FLAG
!       (11)  - SCANNING MODE FLAG
!       (12)  - LATIN 1 - FIRST LAT FROM POLE OF SECANT CONE INTER
!       (13)  - LATIN 2 - SECOND LAT FROM POLE OF SECANT CONE INTER
!  KPTR       - ARRAY CONTAINING STORAGE FOR FOLLOWING PARAMETERS
!               SEE INPUT LIST
!  KRET       - ERROR RETURN
!
! REMARKS:
!  KRET = 0
!       = 4   - DATA REPRESENTATION TYPE NOT CURRENTLY ACCEPTABLE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 77
!   MACHINE:  NAS9000
!
!$$$
!  ************************************************************
!                    INCOMING MESSAGE HOLDER
  CHARACTER :: msga(*)
!
!                    ARRAY GDS ELEMENTS
  INTEGER :: kgds(*)
!                    ARRAY OF POINTERS AND COUNTERS
  INTEGER :: kptr(*)
!
  INTEGER :: kret
!  ---------------------------------------------------------------
  SAVE
  kret    = 0
!             PROCESS GRID DEFINITION SECTION (IF PRESENT)
!          MAKE SURE BIT POINTER IS PROPERLY SET
  kptr(8)  = kptr(9) + (kptr(2)*8) + (kptr(3)*8) + 24
  nsave    = kptr(8) - 24
!  BYTE 4
!                NV - NR OF VERT COORD PARAMETERS
  CALL grbgbyte (msga,kgds(19),kptr(8),8)
  kptr(8)  = kptr(8) + 8
!  BYTE 5
!                PV - LOCATION - SEE FM92 MANUAL
  CALL grbgbyte (msga,kgds(20),kptr(8),8)
  kptr(8)  = kptr(8) + 8
!  BYTE 6
!                   DATA REPRESENTATION TYPE
  CALL grbgbyte (msga,kgds(1),kptr(8),8)
  kptr(8)   = kptr(8) + 8
!        BYTES 7-32 ARE GRID DEFINITION DEPENDING ON
!        DATA REPRESENTATION TYPE
  IF (kgds(1) == 0) THEN
    GO TO 1000
  ELSE IF (kgds(1) == 1) THEN
    GO TO 4000
  ELSE IF (kgds(1) == 2 .OR. kgds(1) == 5) THEN
    GO TO 2000
  ELSE IF (kgds(1) == 3) THEN
    GO TO 5000
  ELSE IF (kgds(1) == 4) THEN
    GO TO 1000
!  ELSE IF (KGDS(1).EQ.10) THEN
!  ELSE IF (KGDS(1).EQ.14) THEN
!  ELSE IF (KGDS(1).EQ.20) THEN
!  ELSE IF (KGDS(1).EQ.24) THEN
!  ELSE IF (KGDS(1).EQ.30) THEN
!  ELSE IF (KGDS(1).EQ.34) THEN
  ELSE IF (kgds(1) == 50) THEN
    GO TO 3000
!  ELSE IF (KGDS(1).EQ.60) THEN
!  ELSE IF (KGDS(1).EQ.70) THEN
!  ELSE IF (KGDS(1).EQ.80) THEN
  ELSE
!                   MARK AS GDS/ UNKNOWN DATA REPRESENTATION TYPE
    kret     = 4
    RETURN
  END IF
!  BYTE 33-N   VERTICAL COORDINATE PARAMETERS
!  -----------
!  BYTES 33-42 EXTENSIONS OF GRID DEFINITION FOR ROTATION
!              OR STRETCHING OF THE COORDINATE SYSTEM OR
!              LAMBERT CONFORMAL PROJECTION.
!  BYTE 43-N   VERTICAL COORDINATE PARAMETERS
!  -----------
!  BYTES 33-52 EXTENSIONS OF GRID DEFINITION FOR STRETCHED
!              AND ROTATED COORDINATE SYSTEM
!  BYTE 53-N   VERTICAL COORDINATE PARAMETERS
!  -----------
! ************************************************************
!  ------------------- LATITUDE/LONGITUDE GRIDS
!
!  ------------------- BYTE 7-8     NR OF POINTS ALONG LATITUDE CIRCLE
  1000 CONTINUE
  CALL grbgbyte (msga,kgds(2),kptr(8),16)
  kptr(8)  = kptr(8) + 16
!  ------------------- BYTE 9-10    NR OF POINTS ALONG LONG MERIDIAN
  CALL grbgbyte (msga,kgds(3),kptr(8),16)
  kptr(8)  = kptr(8) + 16
!  ------------------- BYTE 11-13   LATITUDE OF ORIGIN
  CALL grbgbyte (msga,kgds(4),kptr(8),24)
  kptr(8)  = kptr(8) + 24
  IF (IAND(kgds(4),8388608) /= 0) THEN
    kgds(4)  =  IAND(kgds(4),8388607) * (-1)
  END IF
!  ------------------- BYTE 14-16   LONGITUDE OF ORIGIN
  CALL grbgbyte (msga,kgds(5),kptr(8),24)
  kptr(8)  = kptr(8) + 24
  IF (IAND(kgds(5),8388608) /= 0) THEN
    kgds(5)  =  - IAND(kgds(5),8388607)
  END IF
!  ------------------- BYTE 17      RESOLUTION FLAG
  CALL grbgbyte (msga,kgds(6),kptr(8),8)
  kptr(8)  = kptr(8) + 8
!  ------------------- BYTE 18-20   LATITUDE OF LAST GRID POINT
  CALL grbgbyte (msga,kgds(7),kptr(8),24)
  kptr(8)  = kptr(8) + 24
  IF (IAND(kgds(7),8388608) /= 0) THEN
    kgds(7)  =  - IAND(kgds(7),8388607)
  END IF
!  ------------------- BYTE 21-23   LONGITUDE OF LAST GRID POINT
  CALL grbgbyte (msga,kgds(8),kptr(8),24)
  kptr(8)  = kptr(8) + 24
  IF (IAND(kgds(8),8388608) /= 0) THEN
    kgds(8)  =  - IAND(kgds(8),8388607)
  END IF
!  ------------------- BYTE 24-25   LATITUDINAL DIR INCREMENT
  CALL grbgbyte (msga,kgds(9),kptr(8),16)
  kptr(8)  = kptr(8) + 16
!  ------------------- BYTE 26-27   IF REGULAR LAT/LON GRID
!                                    HAVE LONGIT DIR INCREMENT
!                                ELSE IF GAUSSIAN GRID
!                                    HAVE NR OF LAT CIRCLES
!                                    BETWEEN POLE AND EQUATOR
  CALL grbgbyte (msga,kgds(10),kptr(8),16)
  kptr(8)  = kptr(8) + 16
!  ------------------- BYTE 28      SCANNING MODE FLAGS
  CALL grbgbyte (msga,kgds(11),kptr(8),8)
  kptr(8)  = kptr(8) + 8
!  ------------------- BYTE 29-32   RESERVED
!                          SKIP TO START OF BYTE 33
  CALL grbgbyte (msga,kgds(12),kptr(8),32)
  kptr(8)  = kptr(8) + 32
!  -------------------
  GO TO 900
!  ******************************************************************
!         ' POLAR STEREO PROCESSING '
!
!  ------------------- BYTE 7-8     NR OF POINTS ALONG X=AXIS
  2000 CONTINUE
  CALL grbgbyte (msga,kgds(2),kptr(8),16)
  kptr(8)  = kptr(8) + 16
!  ------------------- BYTE 9-10    NR OF POINTS ALONG Y-AXIS
  CALL grbgbyte (msga,kgds(3),kptr(8),16)
  kptr(8)  = kptr(8) + 16
!  ------------------- BYTE 11-13   LATITUDE OF ORIGIN
  CALL grbgbyte (msga,kgds(4),kptr(8),24)
  kptr(8)  = kptr(8) + 24
  IF (IAND(kgds(4),8388608) /= 0) THEN
    kgds(4)  =  - IAND(kgds(4),8388607)
  END IF
!  ------------------- BYTE 14-16   LONGITUDE OF ORIGIN
  CALL grbgbyte (msga,kgds(5),kptr(8),24)
  kptr(8)  = kptr(8) + 24
  IF (IAND(kgds(5),8388608) /= 0) THEN
    kgds(5)  =   - IAND(kgds(5),8388607)
  END IF
!  ------------------- BYTE 17      RESERVED
  CALL grbgbyte (msga,kgds(6),kptr(8),8)
  kptr(8)  = kptr(8) + 8
!  ------------------- BYTE 18-20   LOV ORIENTATION OF THE GRID
  CALL grbgbyte (msga,kgds(7),kptr(8),24)
  kptr(8)  = kptr(8) + 24
  IF (IAND(kgds(7),8388608) /= 0) THEN
    kgds(7)  =  - IAND(kgds(7),8388607)
  END IF
!  ------------------- BYTE 21-23   DX - THE X DIRECTION INCREMENT
  CALL grbgbyte (msga,kgds(8),kptr(8),24)
  kptr(8)  = kptr(8) + 24
  IF (IAND(kgds(8),8388608) /= 0) THEN
    kgds(8)  =  - IAND(kgds(8),8388607)
  END IF
!  ------------------- BYTE 24-26   DY - THE Y DIRECTION INCREMENT
  CALL grbgbyte (msga,kgds(9),kptr(8),24)
  kptr(8)  = kptr(8) + 24
  IF (IAND(kgds(9),8388608) /= 0) THEN
    kgds(9)  =  - IAND(kgds(9),8388607)
  END IF
!  ------------------- BYTE 27      PROJECTION CENTER FLAG
  CALL grbgbyte (msga,kgds(10),kptr(8),8)
  kptr(8)  = kptr(8) + 8
!  ------------------- BYTE 28      SCANNING MODE
  CALL grbgbyte (msga,kgds(11),kptr(8),8)
  kptr(8)  = kptr(8) + 8
!  ------------------- BYTE 29-32   RESERVED
!                          SKIP TO START OF BYTE 33
  CALL grbgbyte (msga,kgds(12),kptr(8),32)
  kptr(8)  = kptr(8) + 32
!
!  -------------------
  GO TO 900
!
!  ******************************************************************
!  ------------------- GRID DESCRIPTION FOR SPHERICAL HARMONIC COEFF.
!
!  ------------------- BYTE 7-8     J PENTAGONAL RESOLUTION PARAMETER
  3000 CONTINUE
  CALL grbgbyte (msga,kgds(2),kptr(8),16)
  kptr(8)  = kptr(8) + 16
!  ------------------- BYTE 9-10    K PENTAGONAL RESOLUTION PARAMETER
  CALL grbgbyte (msga,kgds(3),kptr(8),16)
  kptr(8)  = kptr(8) + 16
!  ------------------- BYTE 11-12   M PENTAGONAL RESOLUTION PARAMETER
  CALL grbgbyte (msga,kgds(4),kptr(8),16)
  kptr(8)  = kptr(8) + 16
!  ------------------- BYTE 13 REPRESENTATION TYPE
  CALL grbgbyte (msga,kgds(5),kptr(8),8)
  kptr(8)  = kptr(8) + 8
!  ------------------- BYTE 14 COEFFICIENT STORAGE MODE
  CALL grbgbyte (msga,kgds(6),kptr(8),8)
  kptr(8)  = kptr(8) + 8
!  -------------------        EMPTY FIELDS - BYTES 15 - 32
!              SET TO START OF BYTE 33
  kptr(8)  = kptr(8) + 18 * 8
  GO TO 900
!  ******************************************************************
!                   PROCESS MERCATOR GRIDS
!
!  ------------------- BYTE 7-8     NR OF POINTS ALONG LATITUDE CIRCLE
  4000 CONTINUE
  CALL grbgbyte (msga,kgds(2),kptr(8),16)
  kptr(8)  = kptr(8) + 16
!  ------------------- BYTE 9-10    NR OF POINTS ALONG LONG MERIDIAN
  CALL grbgbyte (msga,kgds(3),kptr(8),16)
  kptr(8)  = kptr(8) + 16
!  ------------------- BYTE 11-13   LATITUE OF ORIGIN
  CALL grbgbyte (msga,kgds(4),kptr(8),24)
  kptr(8)  = kptr(8) + 24
  IF (IAND(kgds(4),8388608) /= 0) THEN
    kgds(4)  =  - IAND(kgds(4),8388607)
  END IF
!  ------------------- BYTE 14-16   LONGITUDE OF ORIGIN
  CALL grbgbyte (msga,kgds(5),kptr(8),24)
  kptr(8)  = kptr(8) + 24
  IF (IAND(kgds(5),8388608) /= 0) THEN
    kgds(5)  =  - IAND(kgds(5),8388607)
  END IF
!  ------------------- BYTE 17      RESOLUTION FLAG
  CALL grbgbyte (msga,kgds(6),kptr(8),8)
  kptr(8)  = kptr(8) + 8
!  ------------------- BYTE 18-20   LATITUDE OF EXTREME POINT
  CALL grbgbyte (msga,kgds(7),kptr(8),24)
  kptr(8)  = kptr(8) + 24
  IF (IAND(kgds(7),8388608) /= 0) THEN
    kgds(7)  =  - IAND(kgds(7),8388607)
  END IF
!  ------------------- BYTE 21-23   LONGITUDE OF EXTREME POINT
  CALL grbgbyte (msga,kgds(8),kptr(8),24)
  kptr(8)  = kptr(8) + 24
  IF (IAND(kgds(8),8388608) /= 0) THEN
    kgds(8)  =  - IAND(kgds(8),8388607)
  END IF
!  ------------------- BYTE 24-26   LATITUDE OF PROJECTION INTERSECTION
  CALL grbgbyte (msga,kgds(9),kptr(8),24)
  kptr(8)  = kptr(8) + 24
  IF (IAND(kgds(9),8388608) /= 0) THEN
    kgds(9)  =  - IAND(kgds(9),8388607)
  END IF
!  ------------------- BYTE 27   RESERVED
  CALL grbgbyte (msga,kgds(10),kptr(8),8)
  kptr(8)  = kptr(8) + 8
!  ------------------- BYTE 28      SCANNING MODE
  CALL grbgbyte (msga,kgds(11),kptr(8),8)
  kptr(8)  = kptr(8) + 8
!  ------------------- BYTE 29-31   LONGITUDINAL DIR INCREMENT
  CALL grbgbyte (msga,kgds(12),kptr(8),24)
  kptr(8)  = kptr(8) + 24
  IF (IAND(kgds(12),8388608) /= 0) THEN
    kgds(12)  =  - IAND(kgds(12),8388607)
  END IF
!  ------------------- BYTE 32-34   LATITUDINAL DIR INCREMENT
  CALL grbgbyte (msga,kgds(13),kptr(8),24)
  kptr(8)  = kptr(8) + 24
  IF (IAND(kgds(13),8388608) /= 0) THEN
    kgds(13)  =  - IAND(kgds(13),8388607)
  END IF
!  ------------------- BYTE 35-42   RESERVED
!                     SKIP TO START OF BYTE 43
  kptr(8)  = kptr(8) + 8 * 8
!  -------------------
  GO TO 900
!  ******************************************************************
!                   PROCESS LAMBERT CONFORMAL
!
!  ------------------- BYTE 7-8     NR OF POINTS ALONG X-AXIS
  5000 CONTINUE
  CALL grbgbyte (msga,kgds(2),kptr(8),16)
  kptr(8)  = kptr(8) + 16
!  ------------------- BYTE 9-10    NR OF POINTS ALONG Y-AXIS
  CALL grbgbyte (msga,kgds(3),kptr(8),16)
  kptr(8)  = kptr(8) + 16
!  ------------------- BYTE 11-13   LATITUDE OF ORIGIN
  CALL grbgbyte (msga,kgds(4),kptr(8),24)
  kptr(8)  = kptr(8) + 24
  IF (IAND(kgds(4),8388608) /= 0) THEN
    kgds(4)  =  - IAND(kgds(4),8388607)
  END IF
!  ------------------- BYTE 14-16   LONGITUDE OF ORIGIN (LOWER LEFT)
  CALL grbgbyte (msga,kgds(5),kptr(8),24)
  kptr(8)  = kptr(8) + 24
  IF (IAND(kgds(5),8388608) /= 0) THEN
    kgds(5)  = - IAND(kgds(5),8388607)
  END IF
!  ------------------- BYTE 17      RESOLUTION
  CALL grbgbyte (msga,kgds(6),kptr(8),8)
  kptr(8)  = kptr(8) + 8
!  ------------------- BYTE 18-20   LOV -ORIENTATION OF GRID
  CALL grbgbyte (msga,kgds(7),kptr(8),24)
  kptr(8)  = kptr(8) + 24
  IF (IAND(kgds(7),8388608) /= 0) THEN
    kgds(7)  = - IAND(kgds(7),8388607)
  END IF
!  ------------------- BYTE 21-23   DX - X-DIR INCREMENT
  CALL grbgbyte (msga,kgds(8),kptr(8),24)
  kptr(8)  = kptr(8) + 24
!  ------------------- BYTE 24-26   DY - Y-DIR INCREMENT
  CALL grbgbyte (msga,kgds(9),kptr(8),24)
  kptr(8)  = kptr(8) + 24
!  ------------------- BYTE 27       PROJECTION CENTER FLAG
  CALL grbgbyte (msga,kgds(10),kptr(8),8)
  kptr(8)  = kptr(8) + 8
!  ------------------- BYTE 28      SCANNING MODE
  CALL grbgbyte (msga,kgds(11),kptr(8),8)
  kptr(8)  = kptr(8) + 8
!  ------------------- BYTE 29-31   LATIN1 - 1ST LAT FROM POLE
  CALL grbgbyte (msga,kgds(12),kptr(8),24)
  kptr(8)  = kptr(8) + 24
  IF (IAND(kgds(12),8388608) /= 0) THEN
    kgds(12)  =  - IAND(kgds(12),8388607)
  END IF
!  ------------------- BYTE 32-34   LATIN2 - 2ND LAT FROM POLE
  CALL grbgbyte (msga,kgds(13),kptr(8),24)
  kptr(8)  = kptr(8) + 24
  IF (IAND(kgds(13),8388608) /= 0) THEN
    kgds(13)  =  - IAND(kgds(13),8388607)
  END IF
!  ------------------- BYTE 35-37   LATITUDE OF SOUTHERN POLE
  CALL grbgbyte (msga,kgds(14),kptr(8),24)
  kptr(8)  = kptr(8) + 24
  IF (IAND(kgds(14),8388608) /= 0) THEN
    kgds(14)  =  - IAND(kgds(14),8388607)
  END IF
!  ------------------- BYTE 38-40   LONGITUDE OF SOUTHERN POLE
  CALL grbgbyte (msga,kgds(15),kptr(8),24)
  kptr(8)  = kptr(8) + 24
  IF (IAND(kgds(15),8388608) /= 0) THEN
    kgds(15)  =  - IAND(kgds(15),8388607)
  END IF
!  ------------------- BYTE 41-42   RESERVED
  CALL grbgbyte (msga,kgds(16),kptr(8),16)
  kptr(8)  = kptr(8) + 16
!  -------------------
  900 CONTINUE
!
!                     MORE CODE FOR GRIDS WITH PL
!
  IF (kgds(19) == 0.AND.kgds(20) /= 255) THEN
    isum  = 0
    kptr(8)  = nsave + (kgds(20) - 1) * 8
    CALL grbgbytes (msga,kgds(22),kptr(8),16,0,kgds(3))
    DO j = 1, kgds(3)
      isum  = isum + kgds(21+j)
    END DO
    kgds(21)  = isum
  END IF
  RETURN
END SUBROUTINE fi633

SUBROUTINE fi634(msga,kptr,kpds,kgds,kbms,kret)
!$$$  SUBPROGRAM DOCUMENTATION  BLOCK
!             .      .    .                                       .
! SUBPROGRAM:    FI634       EXTRACT OR GENERATE BIT MAP FOR OUTPUT
!   PRGMMR: BILL CAVANAUGH   ORG: W/NMC42    DATE: 91-09-13
!
! ABSTRACT: IF BIT MAP SEC   IS AVAILABLE IN GRIB MESSAGE, EXTRACT
!   FOR PROGRAM USE, OTHERWISE GENERATE AN APPROPRIATE BIT MAP.
!
! PROGRAM HISTORY LOG:
!   91-09-13  CAVANAUGH
!   91-11-12  CAVANAUGH  MODIFIED SIZE OF ECMWF GRIDS 5 - 8.
!
! USAGE:    CALL FI634(MSGA,KPTR,KPDS,KGDS,KBMS,KRET)
!   INPUT ARGUMENT LIST:
!  MSGA       - BUFR MESSAGE
!  KPTR       - ARRAY CONTAINING STORAGE FOR FOLLOWING PARAMETERS
!       (1)   - TOTAL LENGTH OF GRIB MESSAGE
!       (2)   - LENGTH OF INDICATOR (SECTION  0)
!       (3)   - LENGTH OF PDS       (SECTION  1)
!       (4)   - LENGTH OF GDS       (SECTION  2)
!       (5)   - LENGTH OF BMS       (SECTION  3)
!       (6)   - LENGTH OF BDS       (SECTION  4)
!       (7)   - VALUE OF CURRENT BYTE
!       (8)   - BIT POINTER
!       (9)   - GRIB START BIT NR
!      (10)   - GRIB/GRID ELEMENT COUNT
!      (11)   - NR UNUSED BITS AT END OF SECTION 3
!      (12)   - BIT MAP FLAG
!      (13)   - NR UNUSED BITS AT END OF SECTION 2
!      (14)   - BDS FLAGS
!      (15)   - NR UNUSED BITS AT END OF SECTION 4
!  KPDS     - ARRAY CONTAINING PDS ELEMENTS.
!       (1)   - ID OF CENTER
!       (2)   - MODEL IDENTIFICATION
!       (3)   - GRID IDENTIFICATION
!       (4)   - GDS/BMS FLAG
!       (5)   - INDICATOR OF PARAMETER
!       (6)   - TYPE OF LEVEL
!       (7)   - HEIGHT/PRESSURE , ETC OF LEVEL
!       (8)   - YEAR OF CENTURY
!       (9)   - MONTH OF YEAR
!       (10)  - DAY OF MONTH
!       (11)  - HOUR OF DAY
!       (12)  - MINUTE OF HOUR
!       (13)  - INDICATOR OF FORECAST TIME UNIT
!       (14)  - TIME RANGE 1
!       (15)  - TIME RANGE 2
!       (16)  - TIME RANGE FLAG
!       (17)  - NUMBER INCLUDED IN AVERAGE
!
!   OUTPUT ARGUMENT LIST:
!  KBMS       - BITMAP DESCRIBING LOCATION OF OUTPUT ELEMENTS.
!  KPTR       - ARRAY CONTAINING STORAGE FOR FOLLOWING PARAMETERS
!               SEE INPUT LIST
!  KRET       - ERROR RETURN
!
! REMARKS:
!  KRET   = 0 - NO ERROR
!         = 5 - GRID NOT AVAIL FOR CENTER INDICATED
!         =10 - INCORRECT CENTER INDICATOR
!         =12 - BYTES 5-6 ARE NOT ZERO IN BMS, PREDEFINED BIT MAP
!               NOT PROVIDED BY THIS CENTER
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 77
!   MACHINE:  NAS9000
!
!$$$
!
!                    INCOMING MESSAGE HOLDER
  CHARACTER :: msga(*)
!
!                    BIT MAP
  LOGICAL*1 ::  kbms(*)
!
!                    ARRAY OF POINTERS AND COUNTERS
  INTEGER :: kptr(*)
!                    ARRAY OF POINTERS AND COUNTERS
  INTEGER :: kpds(*)
  INTEGER :: kgds(*)
!
  INTEGER :: kret
  INTEGER :: mask(8)
!  ----------------------GRID 21 AND GRID 22 ARE THE SAME
  LOGICAL*1     grd21( 1369)
!  ----------------------GRID 23 AND GRID 24 ARE THE SAME
  LOGICAL*1     grd23( 1369)
  LOGICAL*1     grd25( 1368)
  LOGICAL*1     grd26( 1368)
!  ----------------------GRID 27 AND GRID 28 ARE THE SAME
!  ----------------------GRID 29 AND GRID 30 ARE THE SAME
!  ----------------------GRID 33 AND GRID 34 ARE THE SAME
  LOGICAL*1     grd50( 1188)
!  -----------------------GRID 61 AND GRID 62 ARE THE SAME
  LOGICAL*1     grd61( 4186)
!  -----------------------GRID 63 AND GRID 64 ARE THE SAME
  LOGICAL*1     grd63( 4186)
!  LOGICAL*1     GRD70(16380)/16380*.TRUE./
!  -------------------------------------------------------------
  SAVE
  DATA  grd21 /1333*.true.,36*.false./
  DATA  grd23 /.true.,36*.false.,1332*.true./
  DATA  grd25 /1297*.true.,71*.false./
  DATA  grd26 /.true.,71*.false.,1296*.true./
  DATA  grd50/                                                          &
! LINE 1-4
  7*.false.,22*.true.,14*.false.,22*.true.,                             &
      14*.false.,22*.true.,14*.false.,22*.true.,7*.false.,              &
! LINE 5-8
  6*.false.,24*.true.,12*.false.,24*.true.,                             &
      12*.false.,24*.true.,12*.false.,24*.true.,6*.false.,              &
! LINE 9-12
  5*.false.,26*.true.,10*.false.,26*.true.,                             &
      10*.false.,26*.true.,10*.false.,26*.true.,5*.false.,              &
! LINE 13-16
  4*.false.,28*.true., 8*.false.,28*.true.,                             &
      8*.false.,28*.true., 8*.false.,28*.true.,4*.false.,               &
! LINE 17-20
  3*.false.,30*.true., 6*.false.,30*.true.,                             &
      6*.false.,30*.true., 6*.false.,30*.true.,3*.false.,               &
! LINE 21-24
  2*.false.,32*.true., 4*.false.,32*.true.,                             &
      4*.false.,32*.true., 4*.false.,32*.true.,2*.false.,               &
! LINE 25-28
  .false.,34*.true., 2*.false.,34*.true.,                               &
      2*.false.,34*.true., 2*.false.,34*.true.,  .false.,               &
! LINE 29-33
  180*.true./
  DATA  grd61 /4096*.true.,90*.false./
  DATA  grd63 /.true.,90*.false.,4095*.true./
  DATA  mask  /128,64,32,16,8,4,2,1/

  LOGICAL, PARAMETER :: verbose = .FALSE.

!  PRINT *,'FI634'
  IF (IAND(kpds(4),64) == 64) THEN
!
!                SET UP BIT POINTER
!                       SECTION 0    SECTION 1     SECTION 2
    kptr(8) = kptr(9) + (kptr(2)*8) + (kptr(3)*8) + (kptr(4)*8) + 24
!
!  BYTE 4           NUMBER OF UNUSED BITS AT END OF SECTION 3
!
    CALL grbgbyte (msga,kptr(11),kptr(8),8)
    kptr(8)  = kptr(8) + 8
!
!  BYTE 5,6         TABLE REFERENCE IF 0, BIT MAP FOLLOWS
!
    CALL grbgbyte (msga,kptr(12),kptr(8),16)
    kptr(8)  = kptr(8) + 16
!                IF TABLE REFERENCE = 0, EXTRACT BIT MAP
    IF (kptr(12) == 0) THEN
!                CALCULATE NR OF BITS IN BIT MAP
      ibits   = (kptr(5) - 6) * 8 - kptr(11)
      kptr(10)  = ibits
      IF (kpds(3) == 21.OR.kpds(3) == 22.OR.kpds(3) == 25               &
                .OR.kpds(3) == 61.OR.kpds(3) == 62) THEN
!                 NORTHERN HEMISPHERE  21, 22, 25, 61, 62
        DO i = 1, ibits
          CALL grbgbyte (msga,ichk,kptr(8),1)
          kptr(8)   = kptr(8) + 1
          IF (ichk /= 0) THEN
            kbms(i)   = .true.
          ELSE
            kbms(i)   = .false.
          END IF
        END DO
        IF (kpds(3) == 25) THEN
          kadd     = 71
        ELSE IF (kpds(3) == 61.OR.kpds(3) == 62) THEN
          kadd     = 90
        ELSE
          kadd     = 36
        END IF
        DO i = 1, kadd
          kbms(i+ibits)  = .false.
        END DO
        kptr(10)   = kptr(10) + kadd
        RETURN
      ELSE IF (kpds(3) == 23.OR.kpds(3) == 24.OR.kpds(3) == 26          &
                .OR.kpds(3) == 63.OR.kpds(3) == 64) THEN
!                 SOUTHERN HEMISPHERE  23, 24, 26, 63, 64
        IF (kpds(3) == 26) THEN
          kadd     = 72
        ELSE IF (kpds(3) == 63.OR.kpds(3) == 64) THEN
          kadd     = 91
        ELSE
          kadd     = 37
        END IF
        DO i = 1, kadd
          kbms(i+ibits)  = .false.
        END DO
        DO i = 1, ibits
          CALL grbgbyte (msga,ichk,kptr(8),1)
          kptr(8)   = kptr(8) + 1
          IF (ichk /= 0) THEN
            kbms(i)   = .true.
          ELSE
            kbms(i)   = .false.
          END IF
        END DO
        kptr(10)   = kptr(10) + kadd - 1
        RETURN
      ELSE IF (kpds(3) == 50) THEN
        kpad    = 7
        kin     = 22
        kbits   = 0
        DO i = 1, 7
          DO j = 1, 4
            DO k = 1, kpad
              kbits   = kbits + 1
              kbms(kbits)  = .false.
            END DO
            DO k = 1, kin
              CALL grbgbyte (msga,ichk,kptr(8),1)
              kptr(8)   = kptr(8) + 1
              kbits     = kbits + 1
              IF (ichk /= 0) THEN
                kbms(kbits)  = .true.
              ELSE
                kbms(kbits)  = .false.
              END IF
            END DO
            DO k = 1, kpad
              kbits   = kbits + 1
              kbms(kbits)  = .false.
            END DO
          END DO
          kin    = kin + 2
          kpad   = kpad - 1
        END DO
        DO ii = 1, 5
          DO j = 1, kin
            CALL grbgbyte (msga,ichk,kptr(8),1)
            kptr(8)  = kptr(8) + 1
            kbits    = kbits + 1
            IF (ichk /= 0) THEN
              kbms(kbits)  = .true.
            ELSE
              kbms(kbits)  = .false.
            END IF
          END DO
        END DO
      ELSE
!                     EXTRACT BIT MAP FROM BMS FOR OTHER GRIDS
        DO i = 1, ibits
          CALL grbgbyte (msga,ichk,kptr(8),1)
          kptr(8)  = kptr(8) + 1
          IF (ichk /= 0) THEN
            kbms(i) = .true.
          ELSE
            kbms(i) = .false.
          END IF
        END DO
      END IF
      RETURN
    ELSE
      PRINT *,'FI634-NO PREDEFINED BIT MAP PROVIDED BY THIS CENTER'
      kret = 12
      RETURN
    END IF
!
  END IF
  kret = 0
!  -------------------------------------------------------
!                PROCESS NON-STANDARD GRID
!  -------------------------------------------------------
  IF (kpds(3) == 255) THEN
    IF(verbose) PRINT *,'NON STANDARD GRID, CENTER = ',kpds(1)
    j      = kgds(2) * kgds(3)
    kptr(10) = j
    DO i = 1, j
      kbms(i) = .true.
    END DO
    RETURN
  END IF
!  -------------------------------------------------------
!                CHECK INTERNATIONAL SET
!  -------------------------------------------------------
  IF (kpds(3) == 21.OR.kpds(3) == 22) THEN
!                ----- INT'L GRIDS 21, 22 - MAP SIZE 1369
    j   = 1369
    kptr(10)  = j
    CALL fi637(*820,j,kpds,kgds,kret)
    DO i = 1, 1369
      kbms(i) = grd21(i)
    END DO
    RETURN
  ELSE IF (kpds(3) == 23.OR.kpds(3) == 24) THEN
!                ----- INT'L GRIDS 23, 24 - MAP SIZE 1369
    j   = 1369
    kptr(10)  = j
    CALL fi637(*820,j,kpds,kgds,kret)
    DO i = 1, 1369
      kbms(i) = grd23(i)
    END DO
    RETURN
  ELSE IF (kpds(3) == 25) THEN
!                ----- INT'L GRID 25 - MAP SIZE 1368
    j   = 1368
    kptr(10)  = j
    CALL fi637(*820,j,kpds,kgds,kret)
    DO i = 1, 1368
      kbms(i) = grd25(i)
    END DO
    RETURN
  ELSE IF (kpds(3) == 26) THEN
!               ----- INT'L GRID  26 - MAP SIZE 1368
    j   = 1368
    kptr(10)  = j
    CALL fi637(*820,j,kpds,kgds,kret)
    DO i = 1, 1368
      kbms(i) = grd26(i)
    END DO
    RETURN
  ELSE IF (kpds(3) >= 37.AND.kpds(3) <= 44) THEN
!               ----- INT'L GRID  37-44 - MAP SIZE 3447
    j   = 3447
    GO TO 800
  ELSE IF (kpds(1) == 7.AND.kpds(3) == 50) THEN
!                ----- INT'L GRIDS 50 - MAP SIZE 964
    j     = 1188
    kptr(10)  = j
    CALL fi637(*890,j,kpds,kgds,kret)
    DO i = 1, j
      kbms(i) = grd50(i)
    END DO
    RETURN
  ELSE IF (kpds(3) == 61.OR.kpds(3) == 62) THEN
!                ----- INT'L GRIDS 61, 62 - MAP SIZE 4186
    j     = 4186
    kptr(10)  = j
    CALL fi637(*820,j,kpds,kgds,kret)
    DO i = 1, 4186
      kbms(i) = grd61(i)
    END DO
    RETURN
  ELSE IF (kpds(3) == 63.OR.kpds(3) == 64) THEN
!               ----- INT'L GRIDS 63, 64 - MAP SIZE 4186
    j     = 4186
    kptr(10)  = j
    CALL fi637(*820,j,kpds,kgds,kret)
    DO i = 1, 4186
      kbms(i) = grd63(i)
    END DO
    RETURN
  END IF
!  -------------------------------------------------------
!                CHECK UNITED STATES SET
!  -------------------------------------------------------
  IF (kpds(1) == 7) THEN
    IF (kpds(3) < 100) THEN
      IF (kpds(3) == 1) THEN
!                    ----- U.S. GRID 1 - MAP SIZE 1679
        j   = 1679
        GO TO 800
      END IF
      IF (kpds(3) == 2) THEN
!                    ----- U.S. GRID 2 - MAP SIZE 10512
        j   = 10512
        GO TO 800
      ELSE IF (kpds(3) == 3) THEN
!                    ----- U.S. GRID 3 - MAP SIZE 65160
        j   = 65160
        GO TO 800
      ELSE IF (KPDS(3) == 4) THEN
!                    ----- U.S. GRID 4 - MAP SIZE 259920
        J   = 259920
        GO TO 800
      ELSE IF (kpds(3) == 5) THEN
!                    ----- U.S. GRID 5 - MAP SIZE 3021
        j   = 3021
        GO TO 800
      ELSE IF (kpds(3) == 6) THEN
!                    ----- U.S. GRID 6 - MAP SIZE 2385
        j   = 2385
        GO TO 800
      ELSE IF (kpds(3) == 27.OR.kpds(3) == 28) THEN
!                    ----- U.S. GRIDS 27, 28 - MAP SIZE 4225
        j     = 4225
        GO TO 800
      ELSE IF (kpds(3) == 29.OR.kpds(3) == 30) THEN
!                    ----- U.S. GRIDS 29,30 - MAP SIZE 5365
        j     = 5365
        GO TO 800
      ELSE IF (kpds(3) == 33.OR.kpds(3) == 34) THEN
!                    ----- U.S GRID 33, 34 - MAP SIZE 8326
        j     = 8326
        GO TO 800
      ELSE IF (kpds(3) >= 37.AND.kpds(3) <= 44) THEN
!               -----  U.S. GRID  37-44 - MAP SIZE 3447
        j   = 3447
        GO TO 800
      ELSE IF (kpds(3) == 45) THEN
!               ----- U.S.  GRID  45    - MAP SIZE 41760
        j   = 41760
        GO TO 800
      ELSE IF (kpds(3) == 55.OR.kpds(3) == 56) THEN
!                    ----- U.S GRID 55, 56 - MAP SIZE 6177
        j     = 6177
        GO TO 800
      ELSE IF (kpds(3) >= 67.AND.kpds(3) <= 71) THEN
!                    ----- U.S GRID 67-71 - MAP SIZE 13689
        j     = 13689
        GO TO 800
      ELSE IF (kpds(3) == 72) THEN
!                    ----- U.S GRID    72 - MAP SIZE 406
        j     = 406
        GO TO 800
      ELSE IF (kpds(3) == 73) THEN
!                    ----- U.S GRID    73 - MAP SIZE 13056
        j     = 13056
        GO TO 800
      ELSE IF (kpds(3) == 74) THEN
!                    ----- U.S GRID    74 - MAP SIZE 10800
        j     = 10800
        GO TO 800
      ELSE IF (kpds(3) >= 75.AND.kpds(3) <= 77) THEN
!                    ----- U.S GRID 75-77 - MAP SIZE 12321
        j     = 12321
        GO TO 800
      ELSE IF (kpds(3) == 85.OR.kpds(3) == 86) THEN
!                    ----- U.S GRID 85,86 - MAP SIZE 32400
        j     = 32400
        GO TO 800
      ELSE IF (kpds(3) == 87) THEN
!                    ----- U.S GRID 87     - MAP SIZE 5022
        j     = 5022
        GO TO 800
      ELSE IF (kpds(3) == 90) THEN
!                    ----- U.S GRID 90     - MAP SIZE 12902
        j     = 12902
        GO TO 800
      ELSE IF (kpds(3) == 91) THEN
!                    ----- U.S GRID 91     - MAP SIZE 25803
        j     = 25803
        GO TO 800
      ELSE IF (kpds(3) == 92) THEN
!                    ----- U.S GRID 92     - MAP SIZE 24162
        j     = 24162
        GO TO 800
      ELSE IF (kpds(3) == 93) THEN
!                    ----- U.S GRID 93     - MAP SIZE 48323
        j     = 48323
        GO TO 800
      ELSE IF (kpds(3) == 94) THEN
!                    ----- U.S GRID 94     - MAP SIZE XXXXX
        j     = 48916
        GO TO 800
      ELSE IF (kpds(3) == 95) THEN
!                    ----- U.S GRID 95     - MAP SIZE XXXXX
        j     = 97831
        GO TO 800
      END IF
    ELSE IF (kpds(3) >= 100.AND.kpds(3) < 200) THEN
      IF (kpds(3) == 100) THEN
!                    ----- U.S. GRID 100 - MAP SIZE 6889
        j     = 6889
        GO TO 800
      ELSE IF (kpds(3) == 101) THEN
!                 ----- U.S. GRID 101 - MAP SIZE 10283
        j     = 10283
        GO TO 800
      ELSE IF (kpds(3) == 103) THEN
!                  ----- U.S. GRID 103 - MAP SIZE 3640
        j     = 3640
        GO TO 800
      ELSE IF (kpds(3) == 104) THEN
!                  ----- U.S. GRID 104 - MAP SIZE 16170
        j     = 16170
        GO TO 800
      ELSE IF (kpds(3) == 105) THEN
!              ----- U.S. GRID 105 - MAP SIZE 6889
        j     = 6889
        GO TO 800
      ELSE IF (kpds(3) == 106) THEN
!                  ----- U.S. GRID 106 - MAP SIZE 19305
        j     = 19305
        GO TO 800
      ELSE IF (kpds(3) == 107) THEN
!              ----- U.S. GRID 107 - MAP SIZE 11040
        j     = 11040
        GO TO 800
      ELSE IF (IAND(kpds(4),128) == 128) THEN
!                  ----- U.S. NON-STANDARD GRID
        GO TO 895
      END IF
    ELSE IF (kpds(3) >= 200) THEN
      IF (kpds(3) == 201) THEN
        j = 4225
        GO TO 800
      ELSE IF (kpds(3) == 202) THEN
        j = 2795
        GO TO 800
      ELSE IF (kpds(3) == 203.OR.kpds(3) == 205) THEN
        j = 1755
        GO TO 800
      ELSE IF (kpds(3) == 204) THEN
        j = 6324
        GO TO 800
      ELSE IF (kpds(3) == 206) THEN
        j = 2091
        GO TO 800
      ELSE IF (kpds(3) == 207) THEN
        j = 1715
        GO TO 800
      ELSE IF (kpds(3) == 208) THEN
        j = 783
        GO TO 800
      ELSE IF (kpds(3) == 209) THEN
        j = 8181
        GO TO 800
      ELSE IF (kpds(3) == 210) THEN
        j = 625
        GO TO 800
      ELSE IF (kpds(3) == 211) THEN
        j = 6045
        GO TO 800
      ELSE IF (kpds(3) == 212) THEN
        j = 23865
        GO TO 800
      ELSE IF (kpds(3) == 213) THEN
        j = 10965
        GO TO 800
      ELSE IF (kpds(3) == 214) THEN
        j = 6693
        GO TO 800
      ELSE IF (kpds(3) == 218) THEN   ! 12km NAM data
        GO TO 900
      ELSE IF (IAND(kpds(4),128) == 128) THEN
        GO TO 895
      END IF
      kret  = 5
      RETURN
    END IF
  END IF
!  -------------------------------------------------------
!                CHECK JAPAN METEOROLOGICAL AGENCY SET
!  -------------------------------------------------------
  IF (kpds(1) == 34) THEN
    IF (IAND(kpds(4),128) == 128) THEN
      PRINT *,'JMA MAP IS NOT PREDEFINED, THE GDS WILL'
      PRINT *,'BE USED TO UNPACK THE DATA, MAP = ',kpds(3)
      GO TO 900
    END IF
  END IF
!  -------------------------------------------------------
!                CHECK CANADIAN SET
!  -------------------------------------------------------
  IF (kpds(1) == 54) THEN
    IF (IAND(kpds(4),128) == 128) THEN
      PRINT *,'CANADIAN MAP IS NOT PREDEFINED, THE GDS WILL'
      PRINT *,'BE USED TO UNPACK THE DATA, MAP = ',kpds(3)
      GO TO 900
    END IF
  END IF
!  -------------------------------------------------------
!                CHECK FNOC SET
!  -------------------------------------------------------
  IF (kpds(1) == 58) THEN
    IF (kpds(3) == 220.OR.kpds(3) == 221) THEN
!                   FNOC GRID 220, 221 - MAPSIZE 3969 (63 * 63)
      j  = 3969
      kptr(10)  = j
      DO i = 1, j
        kbms(i)  = .true.
      END DO
      RETURN
    END IF
    IF (kpds(3) == 223) THEN
!                   FNOC GRID 223 - MAPSIZE 10512 (73 * 144)
      j  = 10512
      kptr(10)  = j
      DO i = 1, j
        kbms(i)  = .true.
      END DO
      RETURN
    END IF
    IF (IAND(kpds(4),128) == 128) THEN
      PRINT *,'FNOC MAP IS NOT PREDEFINED, THE GDS WILL'
      PRINT *,'BE USED TO UNPACK THE DATA, MAP = ',kpds(3)
      GO TO 900
    END IF
  END IF
!  -------------------------------------------------------
!                CHECK UKMET SET
!  -------------------------------------------------------
  IF (kpds(1) == 74) THEN
    IF (IAND(kpds(4),128) == 128) THEN
      GO TO 820
    END IF
  END IF
!  -------------------------------------------------------
!                CHECK ECMWF SET
!  -------------------------------------------------------
  IF (kpds(1) == 98) THEN
    IF (kpds(3) >= 1.AND.kpds(3) <= 12) THEN
      IF (kpds(3) >= 5.AND.kpds(3) <= 8) THEN
        j     = 1073
      ELSE
        j     = 1369
      END IF
      kptr(10)  = j
      CALL fi637(*810,j,kpds,kgds,kret)
      DO i = 1, j
        kbms(i) = .true.
      END DO
      RETURN
    ELSE IF (kpds(3) >= 13.AND.kpds(3) <= 16) THEN
      j         = 361
      kptr(10)  = j
      CALL fi637(*810,j,kpds,kgds,kret)
      DO i = 1, j
        kbms(i) = .true.
      END DO
      RETURN
    ELSE IF (IAND(kpds(4),128) == 128) THEN
      GO TO 810
    ELSE
      kret  = 5
      RETURN
    END IF
  ELSE
    PRINT *,'CENTER',kpds(1),' IS NOT DEFINED'
    IF (IAND(kpds(4),128) == 128) THEN
      PRINT *,'GDS WILL BE USED TO UNPACK THE DATA',                    &
                      ' MAP = ',kpds(3)
      GO TO 900
    ELSE
      kret  = 10
      RETURN
    END IF
  END IF
! =======================================
!
  800 CONTINUE
  kptr(10)  = j
  CALL fi637 (*801,j,kpds,kgds,kret)
  DO i = 1, j
    kbms(i)  = .true.
  END DO
  RETURN
  801 CONTINUE
!
!  ----- THE MAP HAS A GDS, BYTE 7 OF THE (PDS) THE GRID IDENTIFICATION
!  ----- IS NOT 255, THE SIZE OF THE GRID IS NOT THE SAME AS THE
!  ----- PREDEFINED SIZES OF THE U.S. GRIDS, OR KNOWN GRIDS OF THE
!  ----- OF THE OTHER CENTERS. THE GRID CAN BE UNKNOWN, OR FROM AN
!  ----- UNKNOWN CENTER, WE WILL USE THE INFORMATION IN THE GDS TO MAKE
!  ----- A BIT MAP.
!
  810 CONTINUE
  PRINT *,'ECMWF PREDEFINED MAP SIZE DOES NOT MATCH, I WILL USE'
  GO TO 895

  820 CONTINUE
  PRINT *,'U.K. MET PREDEFINED MAP SIZE DOES NOT MATCH, I WILL USE'
  GO TO 895

  890 CONTINUE
  PRINT *,'PREDEFINED MAP SIZE DOES NOT MATCH, I WILL USE'
  895 CONTINUE
  PRINT *,'THE GDS TO UNPACK THE DATA, MAP TYPE = ',kpds(3)

  900 CONTINUE
  j      = kgds(2) * kgds(3)
!                 AFOS AFOS AFOS        SPECIAL CASE
!                          INVOLVES NEXT SINGLE STATEMENT ONLY
  IF (kpds(3) == 211) kret = 0
  kptr(10) = j
  DO i = 1, j
    kbms(i) = .true.
  END DO
  !PRINT *,'EXIT FI634'
  RETURN
END SUBROUTINE fi634

SUBROUTINE fi635(msga,kptr,kpds,kgds,kbms,DATA,kret)
!$$$  SUBPROGRAM DOCUMENTATION  BLOCK
!             .      .    .                                       .
! SUBPROGRAM:    FI635         EXTRACT GRIB DATA ELEMENTS FROM BDS
!   PRGMMR: BILL CAVANAUGH   ORG: W/NMC42    DATE: 91-09-13
!
! ABSTRACT: EXTRACT GRIB DATA FROM BINARY DATA SECTION AND PLACE
!        INTO OUTPUT ARRAY IN PROPER POSITION.
!
! PROGRAM HISTORY LOG:
!   91-09-13  CAVANAUGH
!   94-04-01  CAVANAUGH  MODIFIED CODE TO INCLUDE DECIMAL SCALING WHEN
!                     CALCULATING THE VALUE OF DATA POINTS SPECIFIED
!                     AS BEING EQUAL TO THE REFERENCE VALUE
!   94-11-10  FARLEY     INCREASED MXSIZE FROM 72960 TO 260000
!                     FOR .5 DEGREE SST ANALYSIS FIELDS
!
! USAGE:    CALL FI635(MSGA,KPTR,KPDS,KGDS,KBMS,DATA,KRET)
!   INPUT ARGUMENT LIST:
!  MSGA       - ARRAY CONTAINING GRIB MESSAGE
!  KPTR       - ARRAY CONTAINING STORAGE FOR FOLLOWING PARAMETERS
!       (1)   - TOTAL LENGTH OF GRIB MESSAGE
!       (2)   - LENGTH OF INDICATOR (SECTION  0)
!       (3)   - LENGTH OF PDS       (SECTION  1)
!       (4)   - LENGTH OF GDS       (SECTION  2)
!       (5)   - LENGTH OF BMS       (SECTION  3)
!       (6)   - LENGTH OF BDS       (SECTION  4)
!       (7)   - VALUE OF CURRENT BYTE
!       (8)   - BIT POINTER
!       (9)   - GRIB START BIT NR
!      (10)   - GRIB/GRID ELEMENT COUNT
!      (11)   - NR UNUSED BITS AT END OF SECTION 3
!      (12)   - BIT MAP FLAG
!      (13)   - NR UNUSED BITS AT END OF SECTION 2
!      (14)   - BDS FLAGS
!      (15)   - NR UNUSED BITS AT END OF SECTION 4
!  KPDS     - ARRAY CONTAINING PDS ELEMENTS.
!               SEE INITIAL ROUTINE
!  KBMS       - BITMAP DESCRIBING LOCATION OF OUTPUT ELEMENTS.
!
!   OUTPUT ARGUMENT LIST:
!  KBDS       - INFORMATION EXTRACTED FROM BINARY DATA SECTION
!  KBDS(1)  - N1
!  KBDS(2)  - N2
!  KBDS(3)  - P1
!  KBDS(4)  - P2
!  KBDS(5)  - BIT POINTER TO 2ND ORDER WIDTHS
!  KBDS(6)  -  "    "     "   "   "    BIT MAPS
!  KBDS(7)  -  "    "     "  FIRST ORDER VALUES
!  KBDS(8)  -  "    "     "  SECOND ORDER VALUES
!  KBDS(9)  -  "    "     START OF BDS
!  KBDS(10) -  "    "     MAIN BIT MAP
!  KBDS(11) - BINARY SCALING
!  KBDS(12) - DECIMAL SCALING
!  KBDS(13) - BIT WIDTH OF FIRST ORDER VALUES
!  KBDS(14) - BIT MAP FLAG
!              0 = NO SECOND ORDER BIT MAP
!              1 = SECOND ORDER BIT MAP PRESENT
!  KBDS(15) - SECOND ORDER BIT WIDTH
!  KBDS(16) - CONSTANT / DIFFERENT WIDTHS
!              0 = CONSTANT WIDTHS
!              1 = DIFFERENT WIDTHS
!  KBDS(17) - SINGLE DATUM / MATRIX
!              0 = SINGLE DATUM AT EACH GRID POINT
!              1 = MATRIX OF VALUES AT EACH GRID POINT
!    (18-20)- UNUSED
!
!  DATA       - REAL*4 ARRAY OF GRIDDED ELEMENTS IN GRIB MESSAGE.
!  KPTR       - ARRAY CONTAINING STORAGE FOR FOLLOWING PARAMETERS
!               SEE INPUT LIST
!  KRET       - ERROR RETURN
!
! REMARKS:
!  ERROR RETURN
!           3 = UNPACKED FIELD IS LARGER THAN 65160
!           6 = DOES NOT MATCH NR OF ENTRIES FOR THIS GRIB/GRID
!           7 = NUMBER OF BITS IN FILL TOO LARGE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 77
!   MACHINE:  NAS9000
!
!$$$
!  *************************************************************
!    ON A PC THIS CAN BE CHANGED TO A SMALLER SIZE TO BETTER FIT
!    THE DOS MEMORY LIMIT OF 640K BYTES.  YOU COULD DO THIS
!    FOR MICROSOFT 5.0.  A PC 32 BIT FORTRAN COMPILER
!    WOULD NOT NEED THIS CHANGE.  IF NONE OF YOUR GRIB RECORDS
!    IS LARGER THAN 20000, SET MXSIZE TO 20000.
!  *************************************************************
  INTEGER, PARAMETER :: mxsize = 263000  ! increase from 260000 to 263000
                                         ! for #218 data
  CHARACTER :: msga(*)
  CHARACTER (LEN=1) :: kk(8)
  CHARACTER (LEN=1) :: ckref(8)
!
  LOGICAL*1     kbms(*)
!
  INTEGER :: kpds(*)
  INTEGER :: kgds(*)
  INTEGER :: kbds(20)
  INTEGER :: kptr(*)
  INTEGER :: nrbits
  INTEGER :: kref
  INTEGER :: kkk
  INTEGER :: ksave(mxsize)
  INTEGER :: kscale
!
  REAL :: DATA(*)
  REAL :: refnce
  REAL :: scale
  REAL :: realkk
!
  EQUIVALENCE   (ckref(1),kref,refnce)
  EQUIVALENCE   (kk(1),kkk,realkk)
!
!
!  CHANGED HEX VALUES TO DECIMAL TO MAKE CODE MORE PORTABLE
!
!  *************************************************************
  SAVE
!  PRINT *,'ENTER FI635'
!           SET UP BIT POINTER
  kptr(8) = kptr(9) + (kptr(2)*8) + (kptr(3)*8) + (kptr(4)*8)           &
                  + (kptr(5)*8) + 24
!  ------------- EXTRACT FLAGS
!         BYTE 4
  CALL grbgbyte(msga,kptr(14),kptr(8),4)
  kptr(8)  = kptr(8) + 4
!  --------- NR OF UNUSED BITS IN SECTION 4
  CALL grbgbyte(msga,kptr(15),kptr(8),4)
  kptr(8)  = kptr(8) + 4
  kend    = kptr(9) + (kptr(2)*8) + (kptr(3)*8) + (kptr(4)*8)           &
                  + (kptr(5)*8) + kptr(6) * 8 - kptr(15)
!  ------------- GET SCALE FACTOR
!         BYTES 5,6
!                               CHECK SIGN
  CALL grbgbyte (msga,ksign,kptr(8),1)
  kptr(8)  = kptr(8) + 1
!                               GET ABSOLUTE SCALE VALUE
  CALL grbgbyte (msga,kscale,kptr(8),15)
  kptr(8)  = kptr(8) + 15
  IF (ksign > 0) THEN
    kscale  = - kscale
  END IF
  scale = 2.0**kscale
!  ------------ GET REFERENCE VALUE
!         BYTES 7,10
  CALL grbgbyte (msga,kref,kptr(8),32)
  kptr(8)  = kptr(8) + 32
!
!  THE NEXT CODE WILL CONVERT THE IBM370 FLOATING POINT
!  TO THE FLOATING POINT USED ON YOUR COMPUTER.
!
!  1ST TEST TO SEE IN ON 32 OR 64 BIT WORD MACHINE
!  LW = 4 OR 8;  IF 8 MAY BE A CRAY
!
  CALL w3fi01(lw)
  IF (lw == 4) THEN
    CALL grbgbyte (ckref,jsgn,0,1)
    CALL grbgbyte (ckref,jexp,1,7)
    CALL grbgbyte (ckref,ifr,8,24)
  ELSE
    CALL grbgbyte (ckref,jsgn,32,1)
    CALL grbgbyte (ckref,jexp,33,7)
    CALL grbgbyte (ckref,ifr,40,24)
  END IF
!  PRINT *,109,JSGN,JEXP,IFR
! 109 FORMAT (' JSGN,JEXP,IFR = ',3(1X,Z8))
  IF (ifr == 0) THEN
    refnce  = 0.0
  ELSE IF (jexp == 0 .AND. ifr == 0) THEN
    refnce  = 0.0
  ELSE IF (jexp == 0 ) THEN  ! Added for ECMWF data with Shenzhen project
    refnce  = 0.0            ! By Y. Wang on 04/02/2012
  ELSE
    refnce  = FLOAT(ifr) * 16.0 ** (jexp - 64 - 6)
    IF (jsgn /= 0) refnce = - refnce
  END IF
!  PRINT *,'SCALE ',SCALE,' REF VAL ',KREF,REFNCE
!  ------------- NUMBER OF BITS SPECIFIED FOR EACH ENTRY
!         BYTE 11
  CALL grbgbyte (msga,kbits,kptr(8),8)
  kptr(8)  = kptr(8) + 8
  kbds(4)  = kbits
!  KBDS(13) = KBITS
  ibyt12   = kptr(8)
!  ------------------ IF THERE ARE NO EXTENDED FLAGS PRESENT
!                  THIS IS WHERE DATA BEGINS AND AND THE PROCESSING
!                  INCLUDED IN THE FOLLOWING IF...END IF
!                  WILL BE SKIPPED
!  PRINT *,'BASIC FLAGS =',KPTR(14) ,IAND(KPTR(14),1)
  IF (IAND(kptr(14),1) == 0) THEN
!      PRINT *,'NO EXTENDED FLAGS'
  ELSE
!         BYTES 12,13
    CALL grbgbyte (msga,koctet,kptr(8),16)
    kptr(8)  = kptr(8) + 16
!  --------------------------- EXTENDED FLAGS
!         BYTE 14
    CALL grbgbyte (msga,kxflag,kptr(8),8)
!      PRINT *,'HAVE EXTENDED FLAGS',KXFLAG
    kptr(8)  = kptr(8) + 8
    IF (IAND(kxflag,16) == 0) THEN
!                       SECOND ORDER VALUES CONSTANT WIDTHS
      kbds(16)  = 0
    ELSE
!                       SECOND ORDER VALUES DIFFERENT WIDTHS
      kbds(16)  = 1
    END IF
    IF (IAND (kxflag,32) == 0) THEN
!                      NO SECONDARY BIT MAP
      kbds(14)  = 0
    ELSE
!                      HAVE SECONDARY BIT MAP
      kbds(14)  = 1
    END IF
    IF (IAND (kxflag,64) == 0) THEN
!                      SINGLE DATUM AT GRID POINT
      kbds(17)  = 0
    ELSE
!                      MATRIX OF VALUES AT GRID POINT
      kbds(17)  = 1
    END IF
!  ---------------------- NR - FIRST DIMENSION (ROWS) OF EACH MATRIX
!         BYTES 15,16
    CALL grbgbyte (msga,nr,kptr(8),16)
    kptr(8)  = kptr(8) + 16
!  ---------------------- NC - SECOND DIMENSION (COLS) OF EACH MATRIX
!         BYTES 17,18
    CALL grbgbyte (msga,nc,kptr(8),16)
    kptr(8)  = kptr(8) + 16
!  ---------------------- NRV - FIRST DIM COORD VALS
!         BYTE 19
    CALL grbgbyte (msga,nrv,kptr(8),8)
    kptr(8)  = kptr(8) + 8
!  ---------------------- NC1 - NR COEFF'S OR VALUES
!         BYTE 20
    CALL grbgbyte (msga,nc1,kptr(8),8)
    kptr(8)  = kptr(8) + 8
!  ---------------------- NCV - SECOND DIM COORD OR VALUE
!         BYTE 21
    CALL grbgbyte (msga,ncv,kptr(8),8)
    kptr(8)  = kptr(8) + 8
!  ---------------------- NC2 - NR COEFF'S OR VALS
!         BYTE 22
    CALL grbgbyte (msga,nc2,kptr(8),8)
    kptr(8)  = kptr(8) + 8
!  ---------------------- KPHYS1 - FIRST DIM PHYSICAL SIGNIF
!         BYTE 23
    CALL grbgbyte (msga,kphys1,kptr(8),8)
    kptr(8)  = kptr(8) + 8
!  ---------------------- KPHYS2 - SECOND DIM PHYSICAL SIGNIF
!         BYTE 24
    CALL grbgbyte (msga,kphys2,kptr(8),8)
    kptr(8)  = kptr(8) + 8
!         BYTES 25-N
  END IF
  IF (kbits == 0) THEN
!                    HAVE NO BDS ENTRIES, ALL ENTRIES = REFNCE
    scal10  = 10.0 ** kpds(22)
    scal10  = 1.0 / scal10
    refn10  = refnce * scal10
    kentry = kptr(10)
    DO i = 1, kentry
      DATA(i) = 0.0
      IF (kbms(i)) THEN
        DATA(i) = refn10
      END IF
    END DO
    GO TO 900
  END IF
!  PRINT *,'KEND ',KEND,' KPTR(8) ',KPTR(8),'KBITS ',KBITS
  knr     = (kend - kptr(8)) / kbits
!  PRINT *,'NUMBER OF ENTRIES IN DATA ARRAY',KNR
!  --------------------
!    CYCLE THRU BDS UNTIL HAVE USED ALL (SPECIFIED NUMBER)
!    ENTRIES.
!  ------------- UNUSED BITS IN DATA AREA
! NUMBER OF BYTES IN DATA AREA
  nrbyte  = kptr(6) - 11
!  ------------- TOTAL NR OF USABLE BITS
  nrbits  = nrbyte * 8  - kptr(15)
!  ------------- TOTAL NR OF ENTRIES
  kentry  = nrbits / kbits
!                          MAX SIZE CHECK
  IF (kentry > mxsize) THEN
    kret   = 3
    RETURN
  END IF
!
!  IF (IAND(KPTR(14),2).EQ.0) THEN
!     PRINT *,'SOURCE VALUES IN FLOATING POINT'
!  ELSE
!     PRINT *,'SOURCE VALUES IN INTEGER'
!  END IF
!
  IF (IAND(kptr(14),8) == 0) THEN
!     PRINT *,'PROCESSING GRID POINT DATA'
    IF (IAND(kptr(14),4) == 0) THEN
!         PRINT *,'    WITH SIMPLE PACKING'
      IF (IAND(kptr(14),1) == 0) THEN
!             PRINT *,'        WITH NO ADDITIONAL FLAGS'
        GO TO 4000
      ELSE IF (IAND(kptr(14),1) /= 0) THEN
        PRINT *,'        WITH ADDITIONAL FLAGS',kxflag
        IF (kbds(17) == 0) THEN
          PRINT *,'            SINGLE DATUM EACH GRID PT'
          IF (kbds(14) == 0) THEN
            PRINT *,'            NO SEC BIT MAP'
            IF (kbds(16) == 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                  ' VALUES CONSTANT WIDTH'
            ELSE IF (kbds(16) /= 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                   ' VALUES DIFFERENT WIDTHS'
            END IF
          ELSE IF (kbds(14) /= 0) THEN
            PRINT *,'            SEC BIT MAP'
            IF (kbds(16) == 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                    ' VALUES CONSTANT WIDTH'
            ELSE IF (kbds(16) /= 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                    ' VALUES DIFFERENT WIDTHS'
            END IF
          END IF
        ELSE IF (kbds(17) /= 0) THEN
          PRINT *,'            MATRIX OF VALS EACH PT'
          IF (kbds(14) == 0) THEN
            PRINT *,'            NO SEC BIT MAP'
            IF (kbds(16) == 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                  ' VALUES CONSTANT WIDTH'
            ELSE IF (kbds(16) /= 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                     ' VALUES DIFFERENT WIDTHS'
            END IF
          ELSE IF (kbds(14) /= 0) THEN
            PRINT *,'            SEC BIT MAP'
            IF (kbds(16) == 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                    ' VALUES CONSTANT WIDTH'
            ELSE IF (kbds(16) /= 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                     ' VALUES DIFFERENT WIDTHS'
            END IF
          END IF
        END IF
      END IF
    ELSE IF (IAND(kptr(14),4) /= 0) THEN
      PRINT *,'    WITH COMPLEX/SECOND ORDER PACKING'
      IF (IAND(kptr(14),1) == 0) THEN
        PRINT *,'        WITH NO ADDITIONAL FLAGS'
      ELSE IF (IAND(kptr(14),1) /= 0) THEN
        PRINT *,'        WITH ADDITIONAL FLAGS'
        IF (kbds(17) == 0) THEN
          PRINT *,'            SINGLE DATUM AT EACH PT'
          IF (kbds(14) == 0) THEN
            PRINT *,'            NO SEC BIT MAP'
            IF (kbds(16) == 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                    ' VALUES CONSTANT WIDTH'
            ELSE IF (kbds(16) /= 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                     ' VALUES DIFFERENT WIDTHS'
            END IF
!                                    ROW BY ROW - COL BY COL
            CALL fi636 (DATA,msga,kbms,                                 &
                                  refnce,kptr,kpds,kgds)
            GO TO 900
          ELSE IF (kbds(14) /= 0) THEN
            PRINT *,'            SEC BIT MAP'
            IF (kbds(16) == 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                  ' VALUES CONSTANT WIDTH'
            ELSE IF (kbds(16) /= 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                  ' VALUES DIFFERENT WIDTHS'
            END IF
            CALL fi636 (DATA,msga,kbms,                                 &
                                  refnce,kptr,kpds,kgds)
            GO TO 900
          END IF
        ELSE IF (kbds(17) /= 0) THEN
          PRINT *,'            MATRIX OF VALS EACH PT'
          IF (kbds(14) == 0) THEN
            PRINT *,'            NO SEC BIT MAP'
            IF (kbds(16) == 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                   ' VALUES CONSTANT WIDTH'
            ELSE IF (kbds(16) /= 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                     ' VALUES DIFFERENT WIDTHS'
            END IF
          ELSE IF (kbds(14) /= 0) THEN
            PRINT *,'            SEC BIT MAP'
            IF (kbds(16) == 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                   ' VALUES CONSTANT WIDTH'
            ELSE IF (kbds(16) /= 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                  ' VALUES DIFFERENT WIDTHS'
            END IF
          END IF
        END IF
      END IF
    END IF
  ELSE IF (IAND(kptr(14),8) /= 0) THEN
    PRINT *,'PROCESSING SPHERICAL HARMONIC COEFFICIENTS'
    IF (IAND(kptr(14),4) == 0) THEN
      PRINT *,'    WITH SIMPLE PACKING'
      IF (IAND(kptr(14),1) == 0) THEN
        PRINT *,'        WITH NO ADDITIONAL FLAGS'
        GO TO 5000
      ELSE IF (IAND(kptr(14),1) /= 0) THEN
        PRINT *,'        WITH ADDITIONAL FLAGS'
        IF (kbds(17) == 0) THEN
          PRINT *,'            SINGLE DATUM EACH GRID PT'
          IF (kbds(14) == 0) THEN
            PRINT *,'            NO SEC BIT MAP'
            IF (kbds(16) == 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                     ' VALUES CONSTANT WIDTH'
            ELSE IF (kbds(16) /= 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                     ' VALUES DIFFERENT WIDTHS'
            END IF
          ELSE IF (kbds(14) /= 0) THEN
            PRINT *,'            SEC BIT MAP'
            IF (kbds(16) == 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                     ' VALUES CONSTANT WIDTH'
            ELSE IF (kbds(16) /= 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                   ' VALUES DIFFERENT WIDTHS'
            END IF
          END IF
        ELSE IF (kbds(17) /= 0) THEN
          PRINT *,'            MATRIX OF VALS EACH PT'
          IF (kbds(14) == 0) THEN
            PRINT *,'            NO SEC BIT MAP'
            IF (kbds(16) == 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                     ' VALUES CONSTANT WIDTH'
            ELSE IF (kbds(16) /= 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                    ' VALUES DIFFERENT WIDTHS'
            END IF
          ELSE IF (kbds(14) /= 0) THEN
            PRINT *,'            SEC BIT MAP'
            IF (kbds(16) == 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                     ' VALUES CONSTANT WIDTH'
            ELSE IF (kbds(16) /= 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                    ' VALUES DIFFERENT WIDTHS'
            END IF
          END IF
        END IF
      END IF
    ELSE IF (IAND(kptr(14),4) /= 0) THEN
!                               COMPLEX/SECOND ORDER PACKING
      PRINT *,'    WITH COMPLEX/SECOND ORDER PACKING'
      IF (IAND(kptr(14),1) == 0) THEN
        PRINT *,'        WITH NO ADDITIONAL FLAGS'
      ELSE IF (IAND(kptr(14),1) /= 0) THEN
        PRINT *,'        WITH ADDITIONAL FLAGS'
        IF (kbds(17) == 0) THEN
          PRINT *,'            SINGLE DATUM EACH GRID PT'
          IF (kbds(14) == 0) THEN
            PRINT *,'            NO SEC BIT MAP'
            IF (kbds(16) == 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                    ' VALUES CONSTANT WIDTH'
            ELSE IF (kbds(16) /= 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                     ' VALUES DIFFERENT WIDTHS'
            END IF
          ELSE IF (kbds(14) /= 0) THEN
            PRINT *,'            SEC BIT MAP'
            IF (kbds(16) == 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                     ' VALUES CONSTANT WIDTH'
            ELSE IF (kbds(16) /= 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                     ' VALUES DIFFERENT WIDTHS'
            END IF
          END IF
        ELSE IF (kbds(17) /= 0) THEN
          PRINT *,'            MATRIX OF VALS EACH PT'
          IF (kbds(14) == 0) THEN
            PRINT *,'            NO SEC BIT MAP'
            IF (kbds(16) == 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                   ' VALUES CONSTANT WIDTH'
            ELSE IF (kbds(16) /= 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                     ' VALUES DIFFERENT WIDTHS'
            END IF
          ELSE IF (kbds(14) /= 0) THEN
            PRINT *,'            SEC BIT MAP'
            IF (kbds(16) == 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                     ' VALUES CONSTANT WIDTH'
            ELSE IF (kbds(16) /= 0) THEN
              PRINT *,'            SECOND ORDER',                       &
                     ' VALUES DIFFERENT WIDTHS'
            END IF
          END IF
        END IF
      END IF
    END IF
  END IF
  PRINT *,' NOT PROCESSED - NOT PROCESSED - NOT PROCESSED'
  kret   = 11
  RETURN
  4000 CONTINUE
!  ****************************************************************
!
! GRID POINT DATA, SIMPLE PACKING, FLOATING POINT, NO ADDN'L FLAGS
!
  scal10  = 10.0 ** kpds(22)
  scal10  = 1.0 / scal10
  IF (kpds(3) == 23.OR.kpds(3) == 24.OR.kpds(3) == 26                   &
               .OR.kpds(3) == 63.OR.kpds(3) == 64) THEN
    IF (kpds(3) == 26) THEN
      kadd    = 72
    ELSE IF (kpds(3) == 63.OR.kpds(3) == 64) THEN
      kadd    = 91
    ELSE
      kadd    = 37
    END IF
    CALL grbgbytes (msga,ksave,kptr(8),kbits,0,knr)
    kptr(8)   = kptr(8) + kbits * knr
    ii        = 1
    kentry    = kptr(10)
    DO i = 1, kentry
      IF (kbms(i)) THEN
        DATA(i)   = (refnce+FLOAT(ksave(ii))*scale)*scal10
        ii        = ii + 1
      ELSE
        DATA(i)   = 0.0
      END IF
    END DO
    DO i = 2, kadd
      DATA(i)   = DATA(1)
    END DO
  ELSE IF (kpds(3) == 21.OR.kpds(3) == 22.OR.kpds(3) == 25              &
               .OR.kpds(3) == 61.OR.kpds(3) == 62) THEN
    CALL grbgbytes (msga,ksave,kptr(8),kbits,0,knr)
    ii    = 1
    kentry = kptr(10)
    DO i = 1, kentry
      IF (kbms(i)) THEN
        DATA(i) = (refnce + FLOAT(ksave(ii)) * scale) * scal10
        ii  = ii + 1
      ELSE
        DATA(i) = 0.0
      END IF
    END DO
    IF (kpds(3) == 25) THEN
      kadd    = 71
    ELSE IF (kpds(3) == 61.OR.kpds(3) == 62) THEN
      kadd    = 90
    ELSE
      kadd    = 36
    END IF
    lastp   = kentry - kadd
    DO i = lastp+1, kentry
      DATA(i) = DATA(lastp)
    END DO
  ELSE
    CALL grbgbytes (msga,ksave,kptr(8),kbits,0,knr)
    ii    = 1
    kentry = kptr(10)
    DO i = 1, kentry
      IF (kbms(i)) THEN
        DATA(i) = (refnce + FLOAT(ksave(ii)) * scale) * scal10
        ii  = ii + 1
      ELSE
        DATA(i) = 0.0
      END IF
    END DO
  END IF
  GO TO 900
!  ------------- PROCESS SPHERICAL HARMONIC COEFFICIENTS,
!            SIMPLE PACKING, FLOATING POINT, NO ADDN'L FLAGS
  5000 CONTINUE
!  PRINT *,'CHECK POINT SPECTRAL COEFF'
  kptr(8)  = ibyt12
  CALL grbgbyte (msga,kkk,kptr(8),32)
  kptr(8)  = kptr(8) + 32
!
!  THE NEXT CODE WILL CONVERT THE IBM370 FOATING POINT
!  TO THE FLOATING POINT USED ON YOUR MACHINE.
!
!  1ST TEST TO SEE IN ON 32 OR 64 BIT WORD MACHINE
!  LW = 4 OR 8;  IF 8 MAY BE A CRAY
!
  CALL w3fi01(lw)
  IF (lw == 4) THEN
    CALL grbgbyte (kk,jsgn,0,1)
    CALL grbgbyte (kk,jexp,1,7)
    CALL grbgbyte (kk,ifr,8,24)
  ELSE
    CALL grbgbyte (kk,jsgn,32,1)
    CALL grbgbyte (kk,jexp,33,7)
    CALL grbgbyte (kk,ifr,40,24)
  END IF
!
  IF (ifr == 0) THEN
    realkk  = 0.0
  ELSE IF (jexp == 0.AND.ifr == 0) THEN
    realkk  = 0.0
  ELSE
    realkk  = FLOAT(ifr) * 16.0 ** (jexp - 64 - 6)
    IF (jsgn /= 0) realkk  = -realkk
  END IF
  DATA(1)  = realkk
  CALL grbgbytes (msga,ksave,kptr(8),kbits,0,knr)
!  --------------
  DO i = 1, kentry
    DATA(i+1)  = refnce + FLOAT(ksave(i)) * scale
  END DO
  900 CONTINUE
!  PRINT *,'EXIT FI635'
  RETURN
END SUBROUTINE fi635

SUBROUTINE fi636 (DATA,msga,kbms,refnce,kptr,kpds,kgds)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!             .      .    .                                       .
! SUBPROGRAM:    FI636       PROCESS SECOND ORDER PACKING
!   PRGMMR: CAVANAUGH        ORG: W/NMC42    DATE: 92-09-22
!
! ABSTRACT: PROCESS SECOND ORDER PACKING FROM THE BINARY DATA SECTION
!   (BDS) FOR SINGLE DATA ITEMS GRID POINT DATA
!
! PROGRAM HISTORY LOG:
!   93-06-08  CAVANAUGH
!   93-12-15  CAVANAUGH   MODIFIED SECOND ORDER POINTERS TO FIRST ORDER
!                      VALUES AND SECOND ORDER VALUES CORRECTLY.
!
! USAGE:    CALL FI636 (DATA,MSGA,KBMS,REFNCE,KPTR,KPDS,KGDS)
!   INPUT ARGUMENT LIST:
!
!  MSGA     - ARRAY CONTAINING GRIB MESSAGE
!  REFNCE   - REFERENCE VALUE
!  KPTR     - WORK ARRAY
!
!   OUTPUT ARGUMENT LIST:      (INCLUDING WORK ARRAYS)
!  DATA     - LOCATION OF OUTPUT ARRAY
!           WORKING ARRAY
!  KBDS(1)  - N1
!  KBDS(2)  - N2
!  KBDS(3)  - P1
!  KBDS(4)  - P2
!  KBDS(5)  - BIT POINTER TO 2ND ORDER WIDTHS
!  KBDS(6)  -  "    "     "   "   "    BIT MAPS
!  KBDS(7)  -  "    "     "  FIRST ORDER VALUES
!  KBDS(8)  -  "    "     "  SECOND ORDER VALUES
!  KBDS(9)  -  "    "     START OF BDS
!  KBDS(10) -  "    "     MAIN BIT MAP
!  KBDS(11) - BINARY SCALING
!  KBDS(12) - DECIMAL SCALING
!  KBDS(13) - BIT WIDTH OF FIRST ORDER VALUES
!  KBDS(14) - BIT MAP FLAG
!              0 = NO SECOND ORDER BIT MAP
!              1 = SECOND ORDER BIT MAP PRESENT
!  KBDS(15) - SECOND ORDER BIT WIDTH
!  KBDS(16) - CONSTANT / DIFFERENT WIDTHS
!              0 = CONSTANT WIDTHS
!              1 = DIFFERENT WIDTHS
!  KBDS(17) - SINGLE DATUM / MATRIX
!              0 = SINGLE DATUM AT EACH GRID POINT
!              1 = MATRIX OF VALUES AT EACH GRID POINT
!    (18-20)- UNUSED
!
! REMARKS:
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 77
!   MACHINE:  NAS, CYBER
!
!$$$
  REAL :: DATA(*)
  REAL :: refn
  REAL :: refnce
!
  INTEGER :: kbds(20)
  INTEGER :: kptr(*)
  INTEGER :: jref,bmap2(12500)
  INTEGER :: i,ibds
  INTEGER :: kbit,ifoval,isoval
  INTEGER :: kpds(*),kgds(*)
!
  LOGICAL*1    kbms(*)
!
  CHARACTER :: msga(*)
!
  EQUIVALENCE  (jref,refn)
!  *******************     SETUP     ******************************
  SAVE
!  PRINT *,'ENTER FI636'
  ij  = 0
!  ========================================================
  IF (kbds(14) == 0) THEN
!                        NO BIT MAP, MUST CONSTRUCT ONE
    IF (kgds(2) == 65535) THEN
      IF (kgds(20) == 255) THEN
        PRINT *,'CANNOT BE USED HERE'
      ELSE
!                            POINT TO PL
        lp  = kptr(9) + kptr(2)*8 + kptr(3)*8 + kgds(20)*8 - 8
        jt  = 0
        DO jz = 1, kgds(3)
!                            GET NUMBER IN CURRENT ROW
          CALL grbgbyte (msga,NUMBER,lp,16)
!                            INCREMENT TO NEXT ROW NUMBER
          lp  = lp + 16
!                  PRINT *,'NUMBER IN ROW',JZ,' = ',NUMBER
          DO jq = 1, NUMBER
            jt  = jt + 1
            IF (jq == 1) THEN
              bmap2(jt) = 1
            ELSE
              bmap2(jt) = 0
            END IF
          END DO
        END DO
      END IF
    ELSE
      IF (IAND(kgds(11),32) == 0) THEN
!                        ROW BY ROW
!              PRINT *,'     ROW BY ROW'
        kout  = kgds(3)
        kin   = kgds(2)
      ELSE
!                        COL BY COL
!              PRINT *,'     COL BY COL'
        kin   = kgds(3)
        kout  = kgds(2)
      END IF
!          PRINT *,'KIN=',KIN,' KOUT= ',KOUT
      DO i = 1, kout
        DO j = 1, kin
          IF (j == 1) THEN
            CALL grbsbyte (bmap2,1,ij,1)
          ELSE
            CALL grbsbyte (bmap2,0,ij,1)
          END IF
          ij  = ij + 1
        END DO
      END DO
    END IF
  END IF
!  ========================================================
!  CALL BINARY (BMAP2,2)
!                             START OF BMS (BIT POINTER)
!  KBDS(10)  = 0
!             BYTE START OF BDS
  ibds  = kptr(2) + kptr(3) + kptr(4) + kptr(5)
!  PRINT *,'KPTR(2-5) ',KPTR(2),KPTR(3),KPTR(4),KPTR(5)
!             BIT START OF BDS
  jptr  = ibds * 8
!  PRINT *,'JPTR ',JPTR
  kbds(9) = jptr
!  PRINT *,'START OF BDS         ',KBDS(9)
!                 BINARY SCALE VALUE  BDS BYTES 5-6
  CALL grbgbyte (msga,ISIGN,jptr+32,1)
  CALL grbgbyte (msga,kbds(11),jptr+33,15)
  IF (ISIGN > 0) THEN
    kbds(11)  = - kbds(11)
  END IF
!  PRINT *,'BINARY SCALE VALUE =',KBDS(11)
!               EXTRACT REFERENCE VALUE
  CALL grbgbyte(msga,jref,jptr+48,32)
!  PRINT *,'DECODED REFERENCE VALUE =',REFN,REFNCE
!             F O BIT WIDTH
  CALL grbgbyte(msga,kbds(13),jptr+80,8)
  jptr  = jptr + 88
!           AT START OF BDS BYTE 12
!             EXTRACT N1
  CALL grbgbyte (msga,kbds(1),jptr,16)
!  PRINT *,'N1  = ',KBDS(1)
  jptr  = jptr + 16
!              EXTENDED FLAGS
  CALL grbgbyte (msga,kflag,jptr,8)
!              ISOLATE BIT MAP FLAG
  kbds(14)  = IAND(kflag,32)
  jptr  = jptr + 8
!             EXTRACT N2
  CALL grbgbyte (msga,kbds(2),jptr,16)
!  PRINT *,'N2  = ',KBDS(2)
  jptr  = jptr + 16
!             EXTRACT P1
  CALL grbgbyte (msga,kbds(3),jptr,16)
!  PRINT *,'P1  = ',KBDS(3)
  jptr  = jptr + 16
!             EXTRACT P2
  CALL grbgbyte (msga,kbds(4),jptr,16)
!  PRINT *,'P2  = ',KBDS(4)
  jptr  = jptr + 16
!              SKIP RESERVED BYTE
  jptr    = jptr + 8
!             START OF SECOND ORDER BIT WIDTHS
  kbds(5) = jptr
!             COMPUTE START OF SECONDARY BIT MAP
  IF (kbds(14) /= 0) THEN
!                        FOR INCLUDED SECONDARY BIT MAP
    jptr    = jptr + (kbds(3) * 8)
    kbds(6) = jptr
  ELSE
!                        FOR CONSTRUCTED SECONDARY BIT MAP
    kbds(6)  = 0
  END IF
!             CREATE POINTER TO START OF FIRST ORDER VALUES
  kbds(7) =  kbds(9) + kbds(1) * 8 - 8
  PRINT *,'BIT POINTER TO START OF FOVALS',kbds(7)
!             CREATE POINTER TO START OF SECOND ORDER VALUES
!  KBDS(8) =  KBDS(9) + KBDS(2) * 8 - 8
!  PRINT *,'BIT POINTER TO START OF SOVALS',KBDS(8)
!  PRINT *,'KBDS( 1) - N1                         ',KBDS( 1)
!  PRINT *,'KBDS( 2) - N2                         ',KBDS( 2)
!  PRINT *,'KBDS( 3) - P1                         ',KBDS( 3)
!  PRINT *,'KBDS( 4) - P2                         ',KBDS( 4)
!  PRINT *,'KBDS( 5) - BIT PTR - 2ND ORDER WIDTHS ',KBDS( 5)
!  PRINT *,'KBDS( 6) -  "   "     "   " BIT MAPS  ',KBDS( 6)
!  PRINT *,'KBDS( 7) -  "   "     F O VALS        ',KBDS( 7)
!  PRINT *,'KBDS( 8) -  "   "     S O VALS        ',KBDS( 8)
!  PRINT *,'KBDS( 9) -  "   "    START OF BDS     ',KBDS( 9)
!  PRINT *,'KBDS(10) -  "   "    MAIN BIT MAP     ',KBDS(10)
!  PRINT *,'KBDS(11) - BINARY SCALING             ',KBDS(11)
!  PRINT *,'KPDS(22) - DECIMAL SCALING            ',KPDS(22)
!  PRINT *,'KBDS(13) - FO BIT WIDTH               ',KBDS(13)
!  PRINT *,'KBDS(14) - 2ND ORDER BIT MAP FLAG     ',KBDS(14)
!  PRINT *,'KBDS(15) - 2ND ORDER BIT WIDTH        ',KBDS(15)
!  PRINT *,'KBDS(16) - CONSTANT/DIFFERENT WIDTHS  ',KBDS(16)
!  PRINT *,'KBDS(17) - SINGLE DATUM/MATRIX        ',KBDS(17)
!  PRINT *,'REFNCE VAL                            ',REFNCE
!  ************************* PROCESS DATA  **********************
!             FOR EACH GRID POINT ENTRY
!
  scale2  = 2.0**kbds(11)
  scal10  = 10.0**kpds(22)
!  PRINT *,'SCALE VALUES - ',SCALE2,SCAL10
  DO i = 1, kptr(10)
!                 GET NEXT MASTER BIT MAP BIT POSITION
!                 IF NEXT MASTER BIT MAP BIT POSITION IS 'ON' (1)
    IF (kbms(i)) THEN
!          WRITE(6,900)I,KBMS(I)
! 900         FORMAT (1X,I4,3X,14HMAIN BIT IS ON,3X,L4)
      IF (kbds(14) /= 0) THEN
        CALL grbgbyte (msga,kbit,kbds(6),1)
      ELSE
        CALL grbgbyte (bmap2,kbit,kbds(6),1)
      END IF
!          PRINT *,'KBDS(6) =',KBDS(6),' KBIT =',KBIT
      kbds(6)  = kbds(6) + 1
      IF (kbit /= 0) THEN
!              PRINT *,'          SOB ON'
!                               GET NEXT FIRST ORDER PACKED VALUE
        CALL grbgbyte (msga,ifoval,kbds(7),kbds(13))
        kbds(7)  = kbds(7) + kbds(13)
!              PRINT *,'FOVAL =',IFOVAL
!                                GET SECOND ORDER BIT WIDTH
        CALL grbgbyte (msga,kbds(15),kbds(5),8)
        kbds(5)  = kbds(5) + 8
!             PRINT *,KBDS(7)-KBDS(13),' FOVAL =',IFOVAL,' KBDS(5)=',
!    *                           ,KBDS(5), 'ISOWID =',KBDS(15)
      ELSE
!              PRINT *,'          SOB NOT ON'
      END IF
      isoval  = 0
      IF (kbds(15) == 0) THEN
!                     IF SECOND ORDER BIT WIDTH = 0
!                          THEN SECOND ORDER VALUE IS 0
!                         SO CALCULATE DATA VALUE FOR THIS POINT
!              DATA(I) = (REFNCE + (FLOAT(IFOVAL) * SCALE2)) / SCAL10
      ELSE
        CALL grbgbyte (msga,isoval,kbds(8),kbds(15))
        kbds(8)  = kbds(8) + kbds(15)
      END IF
      DATA(i) = (refnce + (FLOAT(ifoval + isoval) *                     &
                       scale2)) / scal10
!          PRINT *,I,DATA(I),REFNCE,IFOVAL,ISOVAL,SCALE2,SCAL10
    ELSE
!          WRITE(6,901) I,KBMS(I)
! 901         FORMAT (1X,I4,3X,15HMAIN BIT NOT ON,3X,L4)
      DATA(i)  = 0.0
    END IF
!      PRINT *,I,DATA(I),IFOVAL,ISOVAL,KBDS(5),KBDS(15)
  END DO
!  **************************************************************
!  PRINT *,'EXIT FI636'
  RETURN
END SUBROUTINE fi636

SUBROUTINE fi637(*,j,kpds,kgds,kret)
!$$$  SUBPROGRAM DOCUMENTATION  BLOCK
!             .      .    .                                       .
! SUBPROGRAM:    FI637       GRIB GRID/SIZE TEST
!   PRGMMR: CAVANAUGH        ORG: W/NMC42    DATE: 91-09-13
!
! ABSTRACT: TO TEST WHEN GDS IS AVAILABLE TO SEE IF SIZE MISMATCH
!   ON EXISTING GRIDS (BY CENTER) IS INDICATED
!
! PROGRAM HISTORY LOG:
!   91-09-13  CAVANAUGH
!
! USAGE:    CALL FI637(*,J,KPDS,KGDS,KRET)
!   INPUT ARGUMENT LIST:
!  J        - SIZE FOR INDICATED GRID
!  KPDS     -
!  KGDS     -
!
!   OUTPUT ARGUMENT LIST:      (INCLUDING WORK ARRAYS)
!  KRET     - ERROR RETURN
!
! REMARKS:
!  KRET     -
!       = 9 - GDS INDICATES SIZE MISMATCH WITH STD GRID
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 77
!   MACHINE:  NAS
!
!$$$
  INTEGER :: kpds(*)
  INTEGER :: kgds(*)
  INTEGER :: j
  INTEGER :: i
!  ---------------------------------------
  SAVE
!  ---------------------------------------
!        IF GDS NOT INDICATED, RETURN
!  ----------------------------------------
  IF (IAND(kpds(4),128) == 0) RETURN
!  ---------------------------------------
!         GDS IS INDICATED, PROCEED WITH TESTING
!  ---------------------------------------
  IF (kgds(2) == 65535) THEN
    RETURN
  END IF
  i     = kgds(2) * kgds(3)
!  ---------------------------------------
!         INTERNATIONAL SET
!  ---------------------------------------
  IF (kpds(3) >= 21.AND.kpds(3) <= 26) THEN
    IF (i /= j) THEN
      RETURN 1
    END IF
  ELSE IF (kpds(3) >= 37.AND.kpds(3) <= 44) THEN
    IF (i /= j) THEN
      RETURN 1
    END IF
  ELSE IF (kpds(3) == 50) THEN
    IF (i /= j) THEN
      RETURN 1
    END IF
  ELSE IF (kpds(3) >= 61.AND.kpds(3) <= 64) THEN
    IF (i /= j) THEN
      RETURN 1
    END IF
!  ---------------------------------------
!         TEST ECMWF CONTENT
!  ---------------------------------------
  ELSE IF (kpds(1) == 98) THEN
    kret  = 9
    IF (kpds(3) >= 1.AND.kpds(3) <= 16) THEN
      IF (i /= j) THEN
        RETURN 1
      END IF
    ELSE
      kret  = 5
      RETURN 1
    END IF
!  ---------------------------------------
!        U.K. MET OFFICE, BRACKNELL
!  ---------------------------------------
  ELSE IF (kpds(1) == 74) THEN
    kret  = 9
    IF (kpds(3) >= 25.AND.kpds(3) <= 26) THEN
      IF (i /= j) THEN
        RETURN 1
      END IF
    ELSE
      kret  = 5
      RETURN 1
    END IF
!  ---------------------------------------
!        CANADA
!  ---------------------------------------
  ELSE IF (kpds(1) == 54) THEN
    PRINT *,' NO CURRENT LISTING OF CANADIAN GRIDS'
    RETURN 1
!  ---------------------------------------
!        JAPAN METEOROLOGICAL AGENCY
!  ---------------------------------------
  ELSE IF (kpds(1) == 34) THEN
    PRINT *,' NO CURRENT LISTING OF JMA GRIDS'
    RETURN 1
!  ---------------------------------------
!        NAVY - FNOC
!  ---------------------------------------
  ELSE IF (kpds(1) == 58) THEN
    IF (kpds(3) >= 220.AND.kpds(3) <= 221) THEN
      IF (i /= j) THEN
        RETURN 1
      END IF
    ELSE IF (kpds(3) == 223) THEN
      IF (i /= j) THEN
        RETURN 1
      END IF
    END IF
!  ---------------------------------------
!              U.S. GRIDS
!  ---------------------------------------
  ELSE IF (kpds(1) == 7) THEN
    kret  = 9
    IF (kpds(3) >= 1.AND.kpds(3) <= 6) THEN
      IF (i /= j) THEN
        RETURN 1
      END IF
    !ELSE IF (kpds(3) == 5.OR.kpds(3) == 6) THEN
    !  IF (i /= j) THEN
    !    RETURN 1
    !  END IF
    ELSE IF (kpds(3) >= 27.AND.kpds(3) <= 30) THEN
      IF (i /= j) THEN
        RETURN 1
      END IF
    ELSE IF (kpds(3) >= 33.AND.kpds(3) <= 34) THEN
      IF (i /= j) THEN
        RETURN 1
      END IF
    ELSE IF (kpds(3) >= 37.AND.kpds(3) <= 45) THEN
      IF (i /= j) THEN
        RETURN 1
      END IF
    ELSE IF (kpds(3) >= 55.AND.kpds(3) <= 56) THEN
      IF (i /= j) THEN
        RETURN 1
      END IF
    ELSE IF (kpds(3) >= 67.AND.kpds(3) <= 77) THEN
      IF (i /= j) THEN
        RETURN 1
      END IF
    ELSE IF (kpds(3) >= 85.AND.kpds(3) <= 86) THEN
      IF (i /= j) THEN
        RETURN 1
      END IF
    ELSE IF (kpds(3) == 87) THEN
      IF (i /= j) THEN
        RETURN 1
      END IF
    ELSE IF (kpds(3) >= 90.AND.kpds(3) <= 95) THEN
      IF (i /= j) THEN
        RETURN 1
      END IF
    ELSE IF (kpds(3) == 100.OR.kpds(3) == 101) THEN
      IF (i /= j) THEN
        RETURN 1
      END IF
    ELSE IF (kpds(3) >= 103.AND.kpds(3) <= 107) THEN
      IF (i /= j) THEN
        RETURN 1
      END IF
    ELSE IF (kpds(3) >= 201.AND.kpds(3) <= 214) THEN
      IF (i /= j) THEN
        RETURN 1
      END IF
    ELSE
      kret  = 5
      RETURN 1
    END IF
  ELSE
    kret  = 10
    RETURN 1
  END IF
!  ------------------------------------
!                 NORMAL EXIT
!  ------------------------------------
  kret  = 0
  RETURN
END SUBROUTINE fi637

SUBROUTINE w3fi83 (DATA,npts,fval1,fdiff1,iscal2,                       &
           isc10,kpds,kgds)
!$$$  SUBPROGRAM DOCUMENTATION  BLOCK
!             .      .    .                                       .
! SUBPROGRAM:  W3FI83        RESTORE DELTA PACKED DATA TO ORIGINAL
!   PRGMMR: CAVANAUGH        ORG: NMC421      DATE:94-03-09
!
! ABSTRACT: RESTORE DELTA PACKED DATA TO ORIGINAL VALUES
!        RESTORE FROM BOUSTREPHEDONIC ALIGNMENT
!
! PROGRAM HISTORY LOG:
!   93-07-14  CAVANAUGH
!   93-07-22  STACKPOLE      ADDITIONS TO FIX SCALING
!   94-01-27  CAVANAUGH   ADDED REVERSAL OF EVEN NUMBERED ROWS
!                      (BOUSTROPHEDONIC PROCESSING) TO RESTORE
!                      DATA TO ORIGINAL SEQUENCE.
!   94-03-02  CAVANAUGH   CORRECTED REVERSAL OF EVEN NUMBERED ROWS
!
! USAGE:    CALL W3FI83(DATA,NPTS,FVAL1,FDIFF1,ISCAL2,
!    *                                ISC10,KPDS,KGDS)
!   INPUT ARGUMENT LIST:
!  DATA     - SECOND ORDER DIFFERENCES
!  NPTS     - NUMBER OF POINTS IN ARRAY
!  FVAL1    - ORIGINAL FIRST ENTRY IN ARRAY
!  FDIFF1   - ORIGINAL FIRST FIRST-DIFFERENCE
!  ISCAL2   - POWER-OF-TWO EXPONENT FOR UNSCALING
!  ISC10    - POWER-OF-TEN EXPONENT FOR UNSCALING
!  KPDS     - ARRAY OF INFORMATION FOR PDS
!  KGDS     - ARRAY OF INFORMATION FOR GDS
!
!   OUTPUT ARGUMENT LIST:
!  DATA     - EXPANDED ORIGINAL DATA VALUES
!
! REMARKS:
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 77
!   MACHINE:  NAS
!
!$$$
  REAL :: fval1,fdiff1
  REAL :: DATA(*),boust(200)
  INTEGER :: npts,nrow,ncol,kpds(*),kgds(*),isc10
!  ---------------------------------------
  SAVE
!
!  REMOVE DECIMAL UN-SCALING INTRODUCED DURING UNPACKING
!
  dscal = 10.0 ** isc10
  IF (dscal == 0.0) THEN
    DO i=1,npts
      DATA(i) = 1.0
    END DO
  ELSE IF (dscal == 1.0) THEN
  ELSE
    DO i=1,npts
      DATA(i) = DATA(i) * dscal
    END DO
  END IF
!
  DATA(1)  = fval1
  DATA(2)  = fdiff1
  DO j = 3,2,-1
    DO k = j, npts
      DATA(k)  = DATA(k) + DATA(k-1)
    END DO
  END DO
!
!  NOW REMOVE THE BINARY SCALING FROM THE RECONSTRUCTED FIELD
!  AND THE DECIMAL SCALING TOO
!
  IF (dscal == 0) THEN
    scale  = 0.0
  ELSE
    scale =(2.0**iscal2)/dscal
  END IF
  DO i=1,npts
    DATA(i) = DATA(i) * scale
  END DO
!  ==========================================================
  IF (IAND(kpds(4),128) /= 0) THEN
    nrow  = kgds(3)
    ncol  = kgds(2)
!
!   DATA LAID OUT BOUSTROPHEDONIC STYLE
!
!
!      PRINT*, '  REVERSE BOUSTROPHEDON'
    DO i = 2, nrow, 2
!
!       REVERSE THE EVEN NUMBERED ROWS
!
      DO j = 1, ncol
        npos  = i * ncol - j + 1
        boust(j) = DATA(npos)
      END DO
      DO j = 1, ncol
        npos  = ncol * (i-1) + j
        DATA(npos)  = boust(j)
      END DO
    END DO
!
!
  END IF
!  =================================================================
  RETURN
END SUBROUTINE w3fi83


SUBROUTINE w3fi01(lw)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!             .      .    .                                       .
! SUBPROGRAM:  W3FI01        DETERMINES MACHINE WORD LENGTH IN BYTES
!   PRGMMR: KEYSER           ORG: NMC22       DATE:92-07-30
!
! ABSTRACT: DETERMINES THE NUMBER OF BYTES IN A FULL WORD FOR THE
!   PARTICULAR MACHINE (IBM, CRAY, WORKSTATION, PC).
!
! PROGRAM HISTORY LOG:
!   92-01-10  R. KISTLER (W/NMC23)
!   92-05-22  D. A. KEYSER -- DOCBLOCKED/COMMENTED
!   93-03-29  R.E.JONES    -- ADD SAVE STATEMENT
!
! USAGE:    CALL W3FI01(LW)
!   OUTPUT ARGUMENT LIST:      (INCLUDING WORK ARRAYS)
!  LW       - MACHINE WORD LENGTH IN BYTES
!
!   OUTPUT FILES:
!  FT06F001 - PRINTOUT
!
! REMARKS: NONE.
!
! ATTRIBUTES:
!   LANGUAGE: IBM AIX XL FORTRAN Compiler/6000
!   MACHINE:  IBM RS6000 model 530
!
!$$$
!
  CHARACTER (LEN=8) :: ctest1,ctest2
!
  INTEGER :: itest1,itest2
!
  EQUIVALENCE (ctest1,itest1),(ctest2,itest2)
!
  SAVE
!
  DATA  ctest1/'12345678'/
!
  itest2 = itest1
!CCCC PRINT *,' CTEST1 = ',CTEST1,' CTEST2 = ',CTEST2
  IF (ctest1 == ctest2) THEN
!CCCC   PRINT*,' MACHINE WORD LENGTH IS 8 BYTES'
    lw = 8
  ELSE
!CCCC   PRINT*,' MACHINE WORD LENGTH IS 4 BYTES'
    lw = 4
  END IF
  RETURN
END SUBROUTINE w3fi01

SUBROUTINE w3fi04(iendn,itypec,lw)
!$$$   SUBPROGRAM  DOCUMENTATION  BLOCK
!
! SUBPROGRAM: W3FI04         FIND WORD SIZE, ENDIAN, CHARACTER SET
!   PRGMNR: JONES,R.E.       ORG: W/NMC42       DATE: 94-10-07
!
! ABSTRACT: SUBROUTINE COMPUTES WORD SIZE, THE TYPE OF CHARACTER
!   SET, ASCII OR EBCDIC, AND IF THE COMPUTER IS BIG-ENDIAN, OR
!   LITTLE-ENDIAN.
!
! PROGRAM HISTORY LOG:
!   94-10-07  R.E.JONES
!
! USAGE:  CALL W3FI04 (IENDN, ITYPEC, LW)
!
!   OUTPUT ARGUMENT LIST:
!  IENDN     -  INTEGER FOR BIG-ENDIAN OR LITTLE-ENDIAN
!               = 0   BIG-ENDIAN
!               = 1   LITTLE-ENDIAN
!               = 2   CANNOT COMPUTE
!  ITYPEC    -  INTEGER FOR TYPE OF CHARACTER SET
!               = 0   ASCII  CHARACTER SET
!               = 1   EBCDIC CHARACTER SET
!               = 2   NOT ASCII OR EBCDIC
!  LW        -  INTEGER FOR WORDS SIZE OF COMPUTER IN BYTES
!               = 4   FOR 32 BIT COMPUTERS
!               = 8   FOR 64 BIT COMPUTERS
!
! ATTRIBUTES:
!   LANGUAGE: SiliconGraphics 3.5 FORTRAN 77
!   MACHINE:  SiliconGraphics IRIS-4D/25, 35, INDIGO, INDY
!
!$$$
!
  INTEGER :: itest1
  INTEGER :: itest2
  INTEGER :: itest3
  INTEGER :: iendn
  INTEGER :: itypec
  INTEGER :: lw
!
  CHARACTER (LEN=8) :: ctest1
  CHARACTER (LEN=8) :: ctest2
  CHARACTER (LEN=1) :: ctest3(8)
  CHARACTER (LEN=1) :: BLANK
!
  EQUIVALENCE   (ctest1,itest1),(ctest2,itest2)
!
  EQUIVALENCE   (itest3,ctest3(1))
!
  DATA  ctest1/'12345678'/
  DATA  itest3/z'01020304'/
  DATA  BLANK /' '/
!
  SAVE
!
!  TEST FOR TYPE OF CHARACTER SET
!  BLANK IS 32 (20 HEX) IN ASCII, 64 (40 HEX) IN EBCDEC
!
  IF (ICHAR(BLANK) == 32) THEN
    itypec = 0
  ELSE IF (ICHAR(BLANK) == 64) THEN
!
!  COMPUTER IS PROBABLY AN IBM360, 370, OR 390 WITH
!  A 32 BIT WORD SIZE, AND BIG-ENDIAN.
!
    itypec = 1
  ELSE
    itypec = 2
  END IF
!
!  TEST FOR WORD SIZE, SET LW TO 4 FOR 32 BIT COMPUTER,
!  8 FOR FOR 64 BIT COMPUTERS
!
  itest2 = itest1
  IF (ctest1 == ctest2) THEN
!
!  COMPUTER MAY BE A CRAY, OR COULD BE DEC VAX ALPHA
!  OR SGI WITH R4000, R4400, R8800 AFTER THEY CHANGE
!  FORTRAN COMPILERS FOR 64 BIT INTEGER.
!
    lw = 8
  ELSE
    lw = 4
  END IF
!
!  IF CHARACTER SET IS NOT ASCII OR EBCDIC SET IENDN = 2
!  CAN NOT TEST FOR ENDIAN TYPE
!
  IF (itypec == 2) THEN
    iendn = 2
    RETURN
  END IF
!
!  USING ITEST3 WITH Z'01020304' EQUIVALNCED TO CTEST3
!  ON A 32 BIT BIG-ENDIAN COMPUTER 03 IS IN THE 3RD
!  BYTE OF A 4 BYTE WORD. ON A 32 BIT LITTLE-ENDIAN
!  COMPUTER IT IS IN 2ND BYTE.
!  ON A 64 BIT COMPUTER Z'01020304' IS RIGHT ADJUSTED IN
!  A 64 BIT WORD, 03 IS IN THE 7TH BYTE.  ON A LITTLE-
!  ENDIAN 64 BIT COMPUTER IT IS IN THE 2ND BYTE.
!
  IF (lw == 4) THEN
    IF (ICHAR(ctest3(3)) == 3) THEN
      iendn = 0
    ELSE IF (ICHAR(ctest3(3)) == 2) THEN
      iendn = 1
    ELSE
      iendn = 2
    END IF
  ELSE IF (lw == 8) THEN
    IF (ICHAR(ctest3(7)) == 3) THEN
      iendn = 0
    ELSE IF (ICHAR(ctest3(2)) == 3) THEN
      iendn = 1
    ELSE
      iendn = 2
    END IF
  ELSE
    iendn = 2
  END IF
!
  RETURN
END SUBROUTINE w3fi04

!######################################################################
!
! Copy here from gbytes_char.f in w3lib-1.8 because the origional one from
! grbsbytes from lnxlib.F casue an error when PDS section > 28 bytes.
!
!######################################################################

      SUBROUTINE GBYTESC(IN,IOUT,ISKIP,NBYTE,NSKIP,N)
!C          Get bytes - unpack bits:  Extract arbitrary size values from a
!C          packed bit string, right justifying each value in the unpacked
!C          array.
!C            IN    = character*1 array input
!C            IOUT  = unpacked array output
!C            ISKIP = initial number of bits to skip
!C            NBYTE = number of bits to take
!C            NSKIP = additional number of bits to skip on each iteration
!C            N     = number of iterations
!C v1.1
!C
      character*1 in(*)
      integer iout(*)
      integer ones(8), tbit, bitcnt
      save ones
      data ones/1,3,7,15,31,63,127,255/

!c     nbit is the start position of the field in bits
      nbit = iskip
      do i = 1, n
         bitcnt = nbyte
         index=nbit/8+1
         ibit=mod(nbit,8)
         nbit = nbit + nbyte + nskip

!c        first byte
         tbit = min(bitcnt,8-ibit)
         itmp = iand(mova2i(in(index)),ones(8-ibit))
         if (tbit.ne.8-ibit) itmp = ishft(itmp,tbit-8+ibit)
         index = index + 1
         bitcnt = bitcnt - tbit

!c        now transfer whole bytes
         do while (bitcnt.ge.8)
             itmp = ior(ishft(itmp,8),mova2i(in(index)))
             bitcnt = bitcnt - 8
             index = index + 1
         enddo

!c        get data from last byte
         if (bitcnt.gt.0) then
             itmp = ior(ishft(itmp,bitcnt),iand(ishft(mova2i(in(index)), &
     &          -(8-bitcnt)),ones(bitcnt)))
         endif

         iout(i) = itmp
      enddo

      RETURN
      END
