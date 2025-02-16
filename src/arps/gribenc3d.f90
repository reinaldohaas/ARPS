!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE GRIBENC                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE gribenc(itype,wdlen,grbpkbit,npts,fvar,ivar,                 &
           ipds,igds,ibdshd,pds,gds,                                    &
           nbufsz,bds,ibds,msglen,mgrib,                                &
           tem1,item1)

!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Make a GRIB message for variable fvar or ivar for given PDS and
!  GDS.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  11/05/1995
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    itype    Data type: 0 for floating array and 1 for integer array
!    npts     Number of grid points
!    fvar     Floating array to be written into GRIB file
!    ivar     Integer  array to be written into GRIB file
!    ipds     Integer array of PDS
!    igds     Integer array of GDS
!    pds      PDS (GRIB Section 1).
!    gds      GDS (GRIB Section 3).
!
!    nbusz    Buffer size of a GRIB message array
!
!  OUTPUT:
!
!    bds      BDS except the first 11 octets header
!    mgrib    Buffer carrying the GRIB message
!    msglen   Length of the GRIB message
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: itype              ! Data type: 0 - floating, 1 - integer
  INTEGER :: wdlen              ! Length of machine word in bytes
  INTEGER :: grbpkbit           ! Bit length in packing GRIB data
  INTEGER :: npts               ! Number of grid points

  REAL :: fvar(npts)         ! Floating array
  INTEGER :: ivar(npts)         ! Integer  array

  INTEGER :: ipds(25)
  INTEGER :: igds(25)
  INTEGER :: ibdshd(4)

  CHARACTER (LEN=1) :: pds(28)        ! PDS
  CHARACTER (LEN=1) :: gds(42)        ! GDS
  CHARACTER (LEN=1) :: bdshd(11)      ! BDS header

  INTEGER :: nbufsz             ! Size of GRIB message array
  INTEGER :: ibds(nbufsz/4)     ! Identical to BDS
  CHARACTER (LEN=1) :: bds(nbufsz)    ! BDS

  INTEGER :: msglen             ! Length of GRIB message
  CHARACTER (LEN=1) :: mgrib(nbufsz)  ! GRIB message array

  REAL :: tem1(npts)            ! working array

  INTEGER :: item1(npts)        ! working array
!
!-----------------------------------------------------------------------
!
!  Local GRIB parameters
!
!-----------------------------------------------------------------------
!
  INTEGER :: gribid(4)            ! GRIB indicator 'GRIB'
  DATA gribid/ 71,82,73,66/    ! ASCII code for 'G', 'R', 'I', 'B'
  INTEGER :: gribend(4)
  DATA gribend/ 55,55,55,55/   ! ASCII code for '7', '7', '7', '7'

  INTEGER :: bdslen
  INTEGER :: pdslen
  INTEGER :: gdslen
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i
  INTEGER :: npos,npts1
!
!-----------------------------------------------------------------------
!
!  External function
!
!-----------------------------------------------------------------------

  INTEGER :: char2i

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

!
!-----------------------------------------------------------------------
!
!  Making PDS - GRIB Section 1
!
!-----------------------------------------------------------------------
!
  CALL mkpds(ipds,pds)
!
!-----------------------------------------------------------------------
!
!  Making GDS - GRIB Section 2
!
!-----------------------------------------------------------------------
!
  CALL mkgds(igds,gds,npts1)

  IF ( npts /= npts1 ) THEN
    WRITE (6,'(a/a,i8,a,i8/a)')                                         &
        'The point number was not correct',                             &
        'npts = ',npts1, 'while nx*ny = ',npts,                         &
        'Model stopped in GRIBENC'
    CALL arpsstop(' ',1)
  END IF

  pdslen = char2i(pds(1))*65536 + char2i(pds(2))*256 + char2i(pds(3))
  gdslen = char2i(gds(1))*65536 + char2i(gds(2))*256 + char2i(gds(3))
!
!-----------------------------------------------------------------------
!
!  Making BDS - GRIB Section 4
!
!-----------------------------------------------------------------------
!
  CALL mkbds(itype,wdlen,grbpkbit,npts,fvar,ivar,                       &
             pds,ibdshd,bdshd,                                          &
             nbufsz,ibds,bdslen)

  msglen = 8 + pdslen + gdslen + bdslen + 11 + 4   ! no bit map section

!  write (6,'(a,i3,a,i3,a,i8,a,i8)')
!    : ' pdslen = ', pdslen,    ' gdslen = ', gdslen,
!    : ' bdslen = ', bdslen+11, ' msglen = ', msglen
!
!-----------------------------------------------------------------------
!
!  Output sections to the buffer.
!
!  Move SECTION 0 (8 BYTES) into mgrib
!
!-----------------------------------------------------------------------
!
  npos = 0

  DO i = 1,4
    mgrib(i) = CHAR(gribid(i))
  END DO

  mgrib(5) = CHAR(MOD(msglen/65536,256))      ! Total length of
  mgrib(6) = CHAR(MOD(msglen/256,  256))      ! the GRIB message
  mgrib(7) = CHAR(MOD(msglen,      256))      ! in three octets
  mgrib(8) = CHAR(01)                          ! Edition 1

  npos  = npos + 8
!
!-----------------------------------------------------------------------
!
!  Move SECTION 1 - 'PDS' into mgrib
!
!-----------------------------------------------------------------------
!
  IF ( pdslen > 0) THEN
    DO i=1,pdslen
      mgrib(npos+i) = pds(i)
    END DO
  ELSE
    WRITE (6,'(a)') 'PDS length <= 0, pdslen = ', pdslen
  END IF

  npos  = npos + pdslen
!
!-----------------------------------------------------------------------
!
!  Move SECTION 2 - 'GDS' into mgrib
!
!-----------------------------------------------------------------------
!
  IF ( gdslen > 0) THEN
    DO i=1,gdslen
      mgrib(npos+i) = gds(i)
    END DO
  ELSE
    WRITE (6,'(a)') 'GDS length <= 0, gdslen = ', gdslen
  END IF

  npos  = npos + gdslen
!
!-----------------------------------------------------------------------
!
!  No SECTION 3 (BMS) to move
!
!-----------------------------------------------------------------------
!
!  IF ( bmslen.gt.0) THEN
!    DO 230 i=1,bmslen
!      mgrib(npos+i) = bms(i)
!230     CONTINUE
!  END IF
!
!  npos  = npos + bmslen
!
!-----------------------------------------------------------------------
!
!  Move SECTION 4 - 'BDS' into mgrib
!
!-----------------------------------------------------------------------
!
  DO i=1,11
    mgrib(npos+i) = bdshd(i)
  END DO

  npos = npos + 11

  IF ( bdslen > 0) THEN
    DO i=1,bdslen
      mgrib(npos+i) = bds(i)
    END DO
  END IF

  npos = npos + bdslen
!
!-----------------------------------------------------------------------
!
!  End with SECTION 5 ('7777' where '7's ASCII code is 55)
!
!-----------------------------------------------------------------------
!
  DO i=1,4
    mgrib(npos+i) = CHAR(gribend(i))
  END DO

  npos = npos + 4

  IF ( npos /= msglen ) THEN
    WRITE (6,'(a/a,i8,a,i8/a)')                                         &
        'The GRIB message length was incorrect.',                       &
        'npos = ',npos, ' msglen = ',msglen,                            &
        'Program stopped in GRIBENC'
    CALL arpsstop(' ',1)
  END IF

  RETURN
END SUBROUTINE gribenc
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE MKPDS                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mkpds(ipds,pds)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Make PDS from ipds
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  11/05/1995
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    ipds     Integer array of PDS
!
!  OUTPUT:
!
!    pds      PDS (GRIB Section 1).
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: ipds(25)

  CHARACTER (LEN=1) :: pds(28)     ! PDS

  INTEGER :: k,levflg,level,dscl
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  pds(1)  = CHAR(MOD(ipds(1)/65536,256))
  pds(2)  = CHAR(MOD(ipds(1)/256,256))
  pds(3)  = CHAR(MOD(ipds(1),256))
  pds(4)  = CHAR(ipds(2))
  pds(5)  = CHAR(ipds(3))
  pds(6)  = CHAR(ipds(4))
  pds(7)  = CHAR(ipds(5))
  pds(8)  = CHAR(ior(ishft(ipds(6),7),ishft(ipds(7),6)))
  pds(9)  = CHAR(ipds(8))
  pds(10) = CHAR(ipds(9))
  levflg  = ipds(9)
!
!-----------------------------------------------------------------------
!
!  check if level is in two octets or one
!
!-----------------------------------------------------------------------
!
  IF ( (levflg >= 1.AND.levflg <= 100).OR.levflg == 102.OR.             &
          levflg == 103.OR.levflg == 105.OR.levflg == 107.OR.           &
          levflg == 109.OR.levflg == 111.OR.levflg == 113.OR.           &
          levflg == 115.OR.levflg == 125.OR.levflg == 160.OR.           &
          levflg == 200.OR.levflg == 201 ) THEN
    level   = ipds(11)      ! one level stored in pds(11) & pds(12)
    IF ( level < 0 ) THEN
      level = - level          ! positive integer for the abs(level)
      level = ior(level,32768) ! the first bit of 16 bits for sign
    END IF
    pds(11) = CHAR(MOD(level/256,256))
    pds(12) = CHAR(MOD(level,256))
  ELSE
    pds(11) = CHAR(ipds(10)) ! level 1 stored in pds(11)
    pds(12) = CHAR(ipds(11)) ! level 2 stored in pds(12)
  END IF

  pds(13) = CHAR(ipds(12))
  pds(14) = CHAR(ipds(13))
  pds(15) = CHAR(ipds(14))
  pds(16) = CHAR(ipds(15))
  pds(17) = CHAR(ipds(16))
  pds(18) = CHAR(ipds(17))
!
!-----------------------------------------------------------------------
!
!  Check if time P1 is in two octets or one
!
!-----------------------------------------------------------------------
!
  IF (ipds(20) == 10) THEN
    pds(19) = CHAR(MOD(ipds(18)/256,256))
    pds(20) = CHAR(MOD(ipds(18),256))
  ELSE
    pds(19) = CHAR(ipds(18))
    pds(20) = CHAR(ipds(19))
  END IF

  pds(21) = CHAR(ipds(20))
  pds(22) = CHAR(MOD(ipds(21)/256,256))
  pds(23) = CHAR(MOD(ipds(21),256))
  pds(24) = CHAR(ipds(22))
  pds(25) = CHAR(ipds(23))
  pds(26) = CHAR(ipds(24))
  dscl = ipds(25)
  IF (dscl < 0) THEN
    dscl = -dscl            ! If D-scale less than 0,
    dscl = ior(dscl,32768)  ! set the highest bit of two bytes to 1
  END IF

  pds(27) = CHAR(MOD(dscl/256,256))
  pds(28) = CHAR(MOD(dscl,    256))
!
!-----------------------------------------------------------------------
!
!  Set the rest bytes of PDS to zero
!
!-----------------------------------------------------------------------
!
  IF (ipds(1) > 28) THEN
    k = ipds(1)
    DO i = 29,k
      pds(i) = CHAR(0)
    END DO
  ELSE IF ( ipds(1) < 28 ) THEN
    WRITE (6,'(a,i3/a)')                                                &
        'The PDS length must be no less than 28 octets, ipds(1) = ',    &
        ipds(1), 'Program stopped in MKPDS'
  END IF

  RETURN
END SUBROUTINE mkpds
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE MKGDS                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mkgds(igds,gds,npts)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Make GDS from igds
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  11/05/1995
!
!  MODIFICATION HISTORY:
!  2003-12-02  Allow trulats to be negative for Lamb Conf. (Richard Carpenter)
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    igds     Integer array of GDS
!
!  OUTPUT:
!
!    gds      GDS (GRIB Section 2)
!    npts     Number of grid points
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: igds(25)

  CHARACTER (LEN=1) :: gds(42)     ! GDS
  INTEGER :: npts
!
!-----------------------------------------------------------------------
!
!  Local GRIB parameters
!
!-----------------------------------------------------------------------
!
  INTEGER :: gdslen
  INTEGER :: lat,lon
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!-----------------------------------------------------------------------
!
!  ARPS only support 3 types of projection
!
!-----------------------------------------------------------------------
!
  IF ( igds(3) == 5 ) THEN
    gdslen = 32
  ELSE IF (igds(3) == 1 .OR. igds(3) == 3 ) THEN
    gdslen = 42
  ELSE             ! grid type not valid
    WRITE (6,'(a)')                                                     &
        'ARPS does not support this type of projection grid: ',         &
        igds(3), 'Program stopped in MKGDS'
    CALL arpsstop(' ',1)
  END IF
!
!-----------------------------------------------------------------------
!
!  Number of grid points
!
!-----------------------------------------------------------------------
!
  npts = igds(4)*igds(5)   ! npts = nx * ny
!
!-----------------------------------------------------------------------
!
!  Put GDS length in bytes 1,2,3
!
!-----------------------------------------------------------------------
!
  gds(1)  = CHAR(MOD(gdslen/65536,256))
  gds(2)  = CHAR(MOD(gdslen/  256,256))
  gds(3)  = CHAR(MOD(gdslen      ,256))

  gds(4)  = CHAR(igds(1))    ! nz use one byte
  gds(5)  = CHAR(igds(2))    ! = 255, N/A for PV & PL
  gds(6)  = CHAR(igds(3))    ! data representation
  gds(17) = CHAR(igds(8))    ! resolution and u&v component flag
!
!-----------------------------------------------------------------------
!
!  Mercator projection
!
!-----------------------------------------------------------------------
!
  IF ( igds(3) == 1 ) THEN
    gds(7)  = CHAR(MOD(igds(4)/256,256)) ! nx uses two bytes
    gds(8)  = CHAR(MOD(igds(4),    256))
    gds(9)  = CHAR(MOD(igds(5)/256,256)) ! ny uses two bytes
    gds(10) = CHAR(MOD(igds(5),    256))

    lat = igds(6)                 ! swlat
    IF ( lat < 0 ) THEN
      lat = -lat
      lat = ior(lat,8388608)      ! minus sign put in the highest bit
    END IF
    gds(11) = CHAR(MOD(lat/65536,256))
    gds(12) = CHAR(MOD(lat/256,  256))
    gds(13) = CHAR(MOD(lat,      256))

    lon = igds(7)                 ! swlon
    IF ( lon < 0 ) THEN
      lon = -lon
      lon = ior(lon,8388608)      ! minus sign put in the highest bit
    END IF
    gds(14) = CHAR(MOD(lon/65536,256))
    gds(15) = CHAR(MOD(lon/256,  256))
    gds(16) = CHAR(MOD(lon,      256))

    lat = igds(9)                 ! latitude  increment
    IF ( lat < 0 ) THEN
      lat = -lat
      lat = ior(lat,8388608)      ! minus sign put in the highest bit
    END IF
    gds(18) = CHAR(MOD(lat/65536,256))
    gds(19) = CHAR(MOD(lat/  256,256))
    gds(20) = CHAR(MOD(lat      ,256))

    lon = igds(10)                ! longitude increment
    IF ( lon < 0 ) THEN
      lon = -lon
      lon = ior(lon,8388608)      ! minus sign put in the highest bit
    END IF
    gds(21) = CHAR(MOD(lon/65536,256))
    gds(22) = CHAR(MOD(lon/256,  256))
    gds(23) = CHAR(MOD(lon,      256))

    gds(24) = CHAR(MOD(igds(13)/65536,256))  ! true latitude
    gds(25) = CHAR(MOD(igds(13)/256,  256))
    gds(26) = CHAR(MOD(igds(13),      256))
    gds(27) = CHAR(00)
    gds(28) = CHAR(igds(14))
    gds(29) = CHAR(MOD(igds(12)/65536,256))
    gds(30) = CHAR(MOD(igds(12)/256,  256))
    gds(31) = CHAR(MOD(igds(12),      256))
    gds(32) = CHAR(MOD(igds(11)/65536,256))
    gds(33) = CHAR(MOD(igds(11)/256,  256))
    gds(34) = CHAR(MOD(igds(11),      256))
    gds(35) = CHAR(00)
    gds(36) = CHAR(00)
    gds(37) = CHAR(00)
    gds(38) = CHAR(00)
    gds(39) = CHAR(00)
    gds(40) = CHAR(00)
    gds(41) = CHAR(00)
    gds(42) = CHAR(00)
!
!-----------------------------------------------------------------------
!
!  Lambert Conformal projection
!
!-----------------------------------------------------------------------
!
  ELSE IF ( igds(3) == 3 ) THEN
    gds( 7) = CHAR(MOD(igds(4)/256,256))   ! nx
    gds( 8) = CHAR(MOD(igds(4),    256))
    gds( 9) = CHAR(MOD(igds(5)/256,256))   ! ny
    gds(10) = CHAR(MOD(igds(5),    256))

    lat = igds(6)
    IF ( lat < 0 ) THEN
      lat = -lat
      lat = ior(lat,8388608)
    END IF
    gds(11) = CHAR(MOD(lat/65536,256))
    gds(12) = CHAR(MOD(lat/256,  256))
    gds(13) = CHAR(MOD(lat,      256))

    lon = igds(7)
    IF ( lon < 0 ) THEN
      lon = -lon
      lon = ior(lon,8388608)
    END IF
    gds(14) = CHAR(MOD(lon/65536,256))
    gds(15) = CHAR(MOD(lon/256,  256))
    gds(16) = CHAR(MOD(lon,      256))

    lon = igds(9)
    IF ( lon < 0 ) THEN
      lon = -lon
      lon = ior(lon,8388608)
    END IF
    gds(18) = CHAR(MOD(lon/65536,256))
    gds(19) = CHAR(MOD(lon/256,  256))
    gds(20) = CHAR(MOD(lon,      256))

    gds(21) = CHAR(MOD(igds(10)/65536,256))
    gds(22) = CHAR(MOD(igds(10)/  256,256))
    gds(23) = CHAR(MOD(igds(10),      256))
    gds(24) = CHAR(MOD(igds(11)/65536,256))
    gds(25) = CHAR(MOD(igds(11)/256,  256))
    gds(26) = CHAR(MOD(igds(11),      256))
    gds(27) = CHAR(igds(12))
    gds(28) = CHAR(igds(13))

    lat = igds(15)              ! trulat1 = Latin1
    IF ( lat < 0 ) THEN
      lat = -lat
      lat = ior(lat,8388608)
    END IF
    gds(29) = CHAR(MOD(lat/65536,256))
    gds(30) = CHAR(MOD(lat/256,  256))
    gds(31) = CHAR(MOD(lat,      256))

    lat = igds(16)              ! trulat2 = Latin2
    IF ( lat < 0 ) THEN
      lat = -lat
      lat = ior(lat,8388608)
    END IF
    gds(32) = CHAR(MOD(lat/65536,256))
    gds(33) = CHAR(MOD(lat/256,  256))
    gds(34) = CHAR(MOD(lat,      256))

    lat = igds(17)
    IF ( lat < 0 ) THEN
      lat = -lat
      lat = ior(lat,8388608)
    END IF
    gds(35) = CHAR(MOD(lat/65536,256))
    gds(36) = CHAR(MOD(lat/256,  256))
    gds(37) = CHAR(MOD(lat,      256))

    lon = igds(18)
    IF ( lon < 0 ) THEN
      lon = -lon
      lon = ior(lon,8388608)
    END IF
    gds(38) = CHAR(MOD(lon/65536,256))
    gds(39) = CHAR(MOD(lon/256,  256))
    gds(40) = CHAR(MOD(lon,      256))

    gds(41) = CHAR(00)
    gds(42) = CHAR(00)
!
!-----------------------------------------------------------------------
!
!  Polar stereographic projection
!
!-----------------------------------------------------------------------
!
  ELSE IF ( igds(3) == 5 ) THEN
    gds( 7) = CHAR(MOD(igds(4)/256,256))
    gds( 8) = CHAR(MOD(igds(4),    256))
    gds( 9) = CHAR(MOD(igds(5)/256,256))
    gds(10) = CHAR(MOD(igds(5),    256))

    lat = igds(6)
    IF ( lat < 0 ) THEN
      lat = -lat
      lat = ior(lat,8388608)
    END IF
    gds(11) = CHAR(MOD(lat/65536,256))
    gds(12) = CHAR(MOD(lat/256,  256))
    gds(13) = CHAR(MOD(lat,      256))

    lon = igds(7)
    IF ( lon < 0 ) THEN
      lon = -lon
      lon = ior(lon,8388608)
    END IF
    gds(14) = CHAR(MOD(lon/65536,256))
    gds(15) = CHAR(MOD(lon/256,  256))
    gds(16) = CHAR(MOD(lon,      256))

    lon = igds(9)
    IF ( lon < 0 ) THEN
      lon = -lon
      lon = ior(lon,8388608)
    END IF
    gds(18) = CHAR(MOD(lon/65536,256))
    gds(19) = CHAR(MOD(lon/256,   256))
    gds(20) = CHAR(MOD(lon,       256))

    gds(21) = CHAR(MOD(igds(10)/65536,256))
    gds(22) = CHAR(MOD(igds(10)/256,  256))
    gds(23) = CHAR(MOD(igds(10),      256))
    gds(24) = CHAR(MOD(igds(11)/65536,256))
    gds(25) = CHAR(MOD(igds(11)/256,  256))
    gds(26) = CHAR(MOD(igds(11),      256))
    gds(27) = CHAR(igds(12))
    gds(28) = CHAR(igds(13))
    gds(29) = CHAR(00)
    gds(30) = CHAR(00)
    gds(31) = CHAR(00)
    gds(32) = CHAR(00)
  END IF

  RETURN
END SUBROUTINE mkgds
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE MKBDS                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE mkbds(itype,wdlen,grbpkbit,npts,fvar,ivar,                   &
           pds,ibdshd,bdshd,                                            &
           nbufsz,ibds,bdslen)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Pack the data and make the GRIB BDS Section.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  11/05/1995
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    itype    Data type: 0 for floating array and 1 for integer array
!    wdlen    Length of a computer word in bits
!    grbpkbit DATA packing length in bits
!
!    npts     Number of grid points
!    fvar     Floating array to be written into GRIB file
!    ivar     Integer  array to be written into GRIB file
!
!    ibdshd   Integer array for BDS header
!
!    pds      PDS (GRIB Section 1).
!
!    nbfusz   Buffer size of a GRIB message array
!
!  OUTPUT:
!
!    ibds     Identical to BDS
!    bdshd    BDS first 11 octets header
!    bdslen   Length of bds
!
!  TEMPORARY:
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: itype              ! Data type: 0 - floating, 1 - integer

  INTEGER :: grbpkbit           ! DATA packing length in bits
  INTEGER :: wdlen              ! Length of a computer word in bits

  INTEGER :: npts               ! Array size

  REAL :: fvar(npts)         ! Floating array
  INTEGER :: ivar(npts)         ! Integer  array

  CHARACTER (LEN=1) :: pds(28)        ! PDS

  INTEGER :: ibdshd(4)          ! Integer array for BDS header
  CHARACTER (LEN=1) :: bdshd(11)      ! BDS header

  INTEGER :: nbufsz             ! Size of GRIB message array
  INTEGER :: ibds(nbufsz/wdlen) ! Identical to BDS

  INTEGER :: bdslen             ! Length of BDS
!
!-----------------------------------------------------------------------
!
!  Local GRIB parameters
!
!-----------------------------------------------------------------------
!
  INTEGER :: dscl, bscl

  REAL :: scale, fmax,fmin
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: bwidth, nbits, iskip1, iskip2

  INTEGER :: i, j
  REAL :: tema
!
!-----------------------------------------------------------------------
!
!  External function
!
!-----------------------------------------------------------------------
!
  INTEGER :: char2i

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  iskip1 = 0
  iskip2 = 0

  ibdshd(3) = itype      ! set the floating or integer data flag
!
!-----------------------------------------------------------------------
!
!  Perform decimal scale to either floating or integet data array
!
!-----------------------------------------------------------------------
!
  dscl = char2i(pds(27))*256 + char2i(pds(28))
  IF ( IAND(dscl,32768) /= 0 ) THEN
    dscl = - IAND(dscl,32761)
  END IF

  scale = 10.0**dscl
  IF ( itype == 0 ) THEN
    DO i = 1,npts
      fvar(i) = fvar(i)*scale
    END DO
  ELSE IF ( itype == 1 ) THEN
    DO i = 1,npts
      fvar(i) = FLOAT(ivar(i))*scale
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Find the minimum and maximum value of the array and subtract the
!  minimum
!
!-----------------------------------------------------------------------
!
  CALL a3dmax0(fvar,1,npts,1,npts,1,1,1,1,1,1,1,1, fmax,fmin)

!wdt 2001/05/24 Eliminate constant field hack. RLC
!!tmp arpscntl - try having no contsant fields to fix ncl:
!  IF (fmax.eq.fmin) THEN
!    fmax = fmin + 1.0e-10
!    fvar(1) = fmax
!  ENDIF

  DO i=1,npts
    fvar(i) = fvar(i) - fmin
  END DO

  IF ( (fmax-fmin) == 0.0 ) THEN
    bwidth = 0
    bscl = 0
    scale = 1.0
  ELSE IF ( itype /= 0 ) THEN
    bwidth = 8
    bscl = 0
    scale = 1.0
  ELSE
    bwidth = grbpkbit            ! use the value specified by user
    IF (bwidth == 32) THEN
      tema = (fmax-fmin)/(2.**(bwidth)-1.)
    ELSE
      tema = (fmax-fmin)/(2.**(bwidth+1)-1.)
    END IF
    IF ( tema /= 0.0 ) THEN
      tema = LOG(tema)/LOG(2.0) + 2
    END IF
    bscl = MIN( INT(tema), INT(tema+SIGN(1.0,tema)) )
    scale = 2.0**bscl
  END IF
!
!-----------------------------------------------------------------------
!
!  Do binary scale to the data
!
!-----------------------------------------------------------------------
!
  DO i=1,npts
    ivar(i) = nint( fvar(i)/scale )
  END DO

!
!-----------------------------------------------------------------------
!
!  Calculate the length of BDS which must be an even number of bytes
!
!-----------------------------------------------------------------------
!
  nbits = npts*bwidth + 11*8
  j = MOD( nbits, 16 )
  IF ( j /= 0 ) THEN
    nbits = nbits + 16 -j
  END IF

  bdslen = nbits/8
!
!-----------------------------------------------------------------------
!
!  Make the header of BDS
!
!-----------------------------------------------------------------------
!
  bdshd(1)  = CHAR(MOD(bdslen/65536,256))
  bdshd(2)  = CHAR(MOD(bdslen/256,  256))
  bdshd(3)  = CHAR(MOD(bdslen,      256))

  bdshd(4)  = CHAR( ibdshd(1)*128 + ibdshd(2)*64                        &
                  + ibdshd(3)*32  + ibdshd(4)*16 )

  IF ( bscl < 0 ) THEN
    bscl = 32768 - bscl
  END IF

  bdshd(5)  = CHAR(MOD(bscl/256,256))
  bdshd(6)  = CHAR(MOD(bscl,    256))
!
!-----------------------------------------------------------------------
!
!  Convert the reference floating value to IBM GRIB representation
!
!-----------------------------------------------------------------------
!
  CALL flt2ibm( fmin, bdshd(7), wdlen*8 )

!  write (6,'(a,i6,a,e16.10,a,i3,a,i4,a,i8,a,i8)')
!    :' E = ',bscl, ' R = ',fmin,  ' L = ',bwidth,
!    :' A = ',iexp, ' B = ',imant, ' M = ',nint((fmax-fmin)/scale)

  bdshd(11) = CHAR(MOD(bwidth,256))

  IF ( bwidth /= 0 ) THEN
    CALL grbsbytes( ibds,ivar,iskip1,bwidth,iskip2,npts )
  ELSE
    ibds(1) = 0
  END IF

  bdslen = bdslen - 11

  RETURN
END SUBROUTINE mkbds
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE FLT2IBM                     ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE flt2ibm(rfloat,chribm,kbits)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Convert floating point number from machine
!  representation to GRIB representation.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: John Hennessy, ECMWF
!  06/18/1991
!
!  MODIFICATIONS:
!
!  11/14/1995 (Yuhe Liu)
!  Changed name, converted to ARPS standard format, and added some
!  document
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  rfloat    Floating point number to be converted.
!
!  kbits     Number of bits in computer word.
!
!  OUTPUT:
!
!  chribm    32 bits IBM floating format in 4-byte character string
!
!-----------------------------------------------------------------------
!
!  METHOD:
!
!  Floating point number represented as 8 bit signed
!  exponent and 24 bit mantissa in integer values.
!
!  REFERENCE:
!
!  WMO Manual on Codes re GRIB representation.
!
!  Common block variables.
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  REAL :: rfloat
  CHARACTER (LEN=1) :: chribm(4)
  INTEGER :: kbits

  INTEGER :: iexp
  INTEGER :: ISIGN

  INTEGER :: kexp
  INTEGER :: kmant

  REAL :: zeps
  REAL :: zref
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF ( kbits == 32 ) THEN
    zeps = 1.0E-8
  ELSE
    zeps = 1.0E-12
  END IF
  zref = rfloat
!
!-----------------------------------------------------------------------
!
!  Sign of value.
!
!-----------------------------------------------------------------------
!
  ISIGN = 0
  IF (zref < 0.) THEN
    ISIGN = 128
    zref  = - zref
  END IF

  IF ( zref < 1.0E-20 ) zref = 0.0
!
!-----------------------------------------------------------------------
!
!  Convert value of 0.0.
!
!-----------------------------------------------------------------------
!
  IF (zref == 0.0) THEN
    kexp  = 0
    kmant = 0
    iexp  = 0
  ELSE
!
!-----------------------------------------------------------------------
!
!  Convert other values.
!
!  Exponent first
!
!-----------------------------------------------------------------------
!
    iexp = INT(ALOG(zref)*(1.0/ALOG(16.0))+64.0+1.0+zeps)

    IF ( iexp < 0   ) iexp = 0
    IF ( iexp > 127 ) iexp = 127
!
!-----------------------------------------------------------------------
!
!  Mantissa.
!
!-----------------------------------------------------------------------
!
    kmant = nint (zref/16.0**(iexp-70))
!
!-----------------------------------------------------------------------
!
!  Check that mantissa value does not exceed 24 bits.
!  16777215 = 2**24 - 1
!
!-----------------------------------------------------------------------
!
    IF (kmant > 16777215) THEN
      iexp = iexp + 1
      kmant = nint (zref/16.0**(iexp-70))
!
!-----------------------------------------------------------------------
!
!    Check for mantissa overflow. If so, set mantissa to zero
!
!-----------------------------------------------------------------------
!
      IF (kmant > 16777215) THEN
        WRITE (6,'(a,e20.12)') 'Bad mantissa value: overflow',          &
                               rfloat
        WRITE (6,'(a,i2,a,i10,a,i10)')                                  &
            'isign = ', ISIGN, ' iexp = ', iexp, ' kmant = ', kmant
        kmant = 0
      END IF
    END IF
!
!-----------------------------------------------------------------------
!
!  Add sign bit to exponent.
!
!-----------------------------------------------------------------------
!
  END IF

  kexp = iexp + ISIGN

  chribm(1) = CHAR(MOD(kexp,256))
  chribm(2) = CHAR(MOD(kmant/65536,256))
  chribm(3) = CHAR(MOD(kmant/256,  256))
  chribm(4) = CHAR(MOD(kmant,      256))

  RETURN
END SUBROUTINE flt2ibm
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE WRTHISHD                   ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE wrthishd(nx,ny,nz,nzsoil,nchanl,grdbas, fmtver, x,y,z,        &
           wdlen,nbufsz,mgrib,nbytes)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Write the header of ARPS history dump by using GRIB representation
!  to floating point numbers.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  11/27/1995
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nx       Number of grid points in the x-direction (east/west)
!    ny       Number of grid points in the y-direction (north/south)
!    nz       Number of grid points in the vertical
!    nzsoil   Number of soil levels
!
!    nchanl   I/O unit of history dump file
!
!    grdbas   Flag for grid and base state dump
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!
!    fmtver   Format and version identical
!
!    wdlen    Length in bytes of a machine word
!    nbufsz   Buffer size of GRIB message
!
!  OUTPUT:
!
!    mgrib    String containing the header of ARPS history dump
!    nbytes   Length in bytes of the header
!
!  WORK ARRAY:
!
!-----------------------------------------------------------------------
!
!  The following parameters are passed into this subroutine through
!  a common block in globcst.inc, and they determine which
!  variables are output.
!
!  grdout =0 or 1. If grdout=0, grid variables are not dumped.
!  basout =0 or 1. If basout=0, base state variables are not dumped.
!  varout =0 or 1. If varout=0, model perturbation variables are not dumped.
!  mstout =0 or 1. If mstout=0, water variables are not dumped.
!  iceout =0 or 1. If iceout=0, qi, qs and qh are not dumped.
!  rainout=0 or 1. If rainout=0, rain variables are not dumped.
!  prcout =0 or 1. If prcout=0, precipitation rates are not dumped.
!  tkeout =0 or 1. If tkeout=0, tke is not dumped.
!  trbout =0 or 1. If trbout=0, turbulence parameter km is not dumped.
!  sfcout =0 or 1. If sfcout=0, surface variables are not dumped.
!  landout=0 or 1. If landout=0, surface propertty arrays are not dumped.
!  radout =0 or 1. If radout=0, radiation arrays are not dumped
!  flxout =0 or 1. If flxout=0, surface fluxes are not dumped.
!
!  These following parameters are also passed in through common
!  blocks in globcst.inc.
!
!  runname,curtim,umove,vmove,xgrdorg,ygrdorg
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nx,ny,nz          ! Number of grid points in 3 directions
  INTEGER :: nzsoil            ! Number of soil levels

  INTEGER :: nchanl            ! I/O unit

  INTEGER :: grdbas            ! Flag for grid and base state dump

  REAL :: x     (nx)           ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y     (ny)           ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z     (nz)           ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.

  CHARACTER (LEN=40) :: fmtver

  INTEGER :: wdlen             ! Length in bytes of machine word
  INTEGER :: nbufsz
  INTEGER :: nbytes            ! Length in bytes of the header
  CHARACTER (LEN=1) :: mgrib(nbufsz) ! Header of ARPS GRIB file
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i,j,k
  INTEGER :: idummy
  INTEGER :: r0exp,r0mant

  CHARACTER (LEN=1) :: chrtmp(4)
!
!-----------------------------------------------------------------------
!
!  Include files:
!
!-----------------------------------------------------------------------
!
  INCLUDE 'globcst.inc'
  INCLUDE 'grid.inc'          ! Grid & map parameters.
  INCLUDE 'phycst.inc'
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  idummy = 0
  r0exp  = 0                   ! for rdummy = 0.0
  r0mant = 0                   ! for rdummy = 0.0

  nbytes = 3

  DO i=1,40
    mgrib(nbytes+i) = fmtver(i:i)
  END DO
  nbytes = nbytes + 40

  DO i=1,80
    mgrib(nbytes+i) = runname(i:i)
  END DO
  nbytes = nbytes + 80

  mgrib(nbytes+1) = CHAR(nocmnt)
  nbytes = nbytes + 1

  IF( nocmnt > 0 ) THEN
    DO j=1,nocmnt
      DO  i=1,80
        mgrib(nbytes+i) = cmnt(j)(i:i)
      END DO
      nbytes = nbytes + 80
    END DO
  END IF

  CALL flt2ibm(curtim,chrtmp,wdlen*8)
  mgrib(nbytes+1) = chrtmp(1)
  mgrib(nbytes+2) = chrtmp(2)
  mgrib(nbytes+3) = chrtmp(3)
  mgrib(nbytes+4) = chrtmp(4)
  nbytes = nbytes + 4

  mgrib(nbytes+1) = CHAR(MOD(nx/256,256))
  mgrib(nbytes+2) = CHAR(MOD(nx,    256))
  nbytes = nbytes + 2

  mgrib(nbytes+1) = CHAR(MOD(ny/256,256))
  mgrib(nbytes+2) = CHAR(MOD(ny,    256))
  nbytes = nbytes + 2

  mgrib(nbytes+1) = CHAR(MOD(nz/256,256))
  mgrib(nbytes+2) = CHAR(MOD(nz,    256))
  nbytes = nbytes + 2

  mgrib(nbytes+1) = CHAR(MOD(nzsoil/256,256))
  mgrib(nbytes+2) = CHAR(MOD(nzsoil,    256))
  nbytes = nbytes + 2

  IF( grdbas == 1 ) THEN

    mgrib(nbytes+1) = CHAR(1)
    mgrib(nbytes+2) = CHAR(1)
    mgrib(nbytes+3) = CHAR(0)
    mgrib(nbytes+4) = CHAR(mstout)
    mgrib(nbytes+5) = CHAR(0)
    mgrib(nbytes+6) = CHAR(0)
    mgrib(nbytes+7) = CHAR(0)
    mgrib(nbytes+8) = CHAR(0)
    mgrib(nbytes+9) = CHAR(landout)
    mgrib(nbytes+10) = CHAR(totout)
    mgrib(nbytes+11) = CHAR(0)
    mgrib(nbytes+12) = CHAR(idummy)
    mgrib(nbytes+13) = CHAR(idummy)
    mgrib(nbytes+14) = CHAR(mapproj)
    mgrib(nbytes+15) = CHAR(month)
    mgrib(nbytes+16) = CHAR(day)
    mgrib(nbytes+17) = CHAR(MOD(year/256,256))
    mgrib(nbytes+18) = CHAR(MOD(year,    256))
    mgrib(nbytes+19) = CHAR(hour)
    mgrib(nbytes+20) = CHAR(minute)
    mgrib(nbytes+21) = CHAR(second)

    nbytes = nbytes + 21

  ELSE

    mgrib(nbytes+1) = CHAR(grdout)
    mgrib(nbytes+2) = CHAR(basout)
    mgrib(nbytes+3) = CHAR(varout)
    mgrib(nbytes+4) = CHAR(mstout)
    mgrib(nbytes+5) = CHAR(iceout)
    mgrib(nbytes+6) = CHAR(trbout)
    mgrib(nbytes+7) = CHAR(sfcout)
    mgrib(nbytes+8) = CHAR(rainout)
    mgrib(nbytes+9) = CHAR(landout)
    mgrib(nbytes+10) = CHAR(totout)
    mgrib(nbytes+11) = CHAR(tkeout)
    mgrib(nbytes+12) = CHAR(idummy)
    mgrib(nbytes+13) = CHAR(idummy)
    mgrib(nbytes+14) = CHAR(mapproj)
    mgrib(nbytes+15) = CHAR(month)
    mgrib(nbytes+16) = CHAR(day)
    mgrib(nbytes+17) = CHAR(MOD(year/256,256))
    mgrib(nbytes+18) = CHAR(MOD(year,    256))
    mgrib(nbytes+19) = CHAR(hour)
    mgrib(nbytes+20) = CHAR(minute)
    mgrib(nbytes+21) = CHAR(second)

    nbytes = nbytes + 21

  END IF

  CALL flt2ibm(umove,chrtmp,wdlen*8)
  mgrib(nbytes+1) = chrtmp(1)
  mgrib(nbytes+2) = chrtmp(2)
  mgrib(nbytes+3) = chrtmp(3)
  mgrib(nbytes+4) = chrtmp(4)
  nbytes = nbytes + 4

  CALL flt2ibm(vmove,chrtmp,wdlen*8)
  mgrib(nbytes+1) = chrtmp(1)
  mgrib(nbytes+2) = chrtmp(2)
  mgrib(nbytes+3) = chrtmp(3)
  mgrib(nbytes+4) = chrtmp(4)
  nbytes = nbytes + 4

  CALL flt2ibm(xgrdorg,chrtmp,wdlen*8)
  mgrib(nbytes+1) = chrtmp(1)
  mgrib(nbytes+2) = chrtmp(2)
  mgrib(nbytes+3) = chrtmp(3)
  mgrib(nbytes+4) = chrtmp(4)
  nbytes = nbytes + 4

  CALL flt2ibm(ygrdorg,chrtmp,wdlen*8)
  mgrib(nbytes+1) = chrtmp(1)
  mgrib(nbytes+2) = chrtmp(2)
  mgrib(nbytes+3) = chrtmp(3)
  mgrib(nbytes+4) = chrtmp(4)
  nbytes = nbytes + 4

  CALL flt2ibm(trulat1,chrtmp,wdlen*8)
  mgrib(nbytes+1) = chrtmp(1)
  mgrib(nbytes+2) = chrtmp(2)
  mgrib(nbytes+3) = chrtmp(3)
  mgrib(nbytes+4) = chrtmp(4)
  nbytes = nbytes + 4

  CALL flt2ibm(trulat2,chrtmp,wdlen*8)
  mgrib(nbytes+1) = chrtmp(1)
  mgrib(nbytes+2) = chrtmp(2)
  mgrib(nbytes+3) = chrtmp(3)
  mgrib(nbytes+4) = chrtmp(4)
  nbytes = nbytes + 4

  CALL flt2ibm(trulon,chrtmp,wdlen*8)
  mgrib(nbytes+1) = chrtmp(1)
  mgrib(nbytes+2) = chrtmp(2)
  mgrib(nbytes+3) = chrtmp(3)
  mgrib(nbytes+4) = chrtmp(4)
  nbytes = nbytes + 4

  CALL flt2ibm(sclfct,chrtmp,wdlen*8)
  mgrib(nbytes+1) = chrtmp(1)
  mgrib(nbytes+2) = chrtmp(2)
  mgrib(nbytes+3) = chrtmp(3)
  mgrib(nbytes+4) = chrtmp(4)
  nbytes = nbytes + 4

  CALL flt2ibm(n0rain,chrtmp,wdlen*8)
  mgrib(nbytes+1) = chrtmp(1)
  mgrib(nbytes+2) = chrtmp(2)
  mgrib(nbytes+3) = chrtmp(3)
  mgrib(nbytes+4) = chrtmp(4)
  nbytes = nbytes + 4

  CALL flt2ibm(n0snow,chrtmp,wdlen*8)
  mgrib(nbytes+1) = chrtmp(1)
  mgrib(nbytes+2) = chrtmp(2)
  mgrib(nbytes+3) = chrtmp(3)
  mgrib(nbytes+4) = chrtmp(4)
  nbytes = nbytes + 4

  CALL flt2ibm(n0hail,chrtmp,wdlen*8)
  mgrib(nbytes+1) = chrtmp(1)
  mgrib(nbytes+2) = chrtmp(2)
  mgrib(nbytes+3) = chrtmp(3)
  mgrib(nbytes+4) = chrtmp(4)
  nbytes = nbytes + 4

  CALL flt2ibm(rhosnow,chrtmp,wdlen*8)
  mgrib(nbytes+1) = chrtmp(1)
  mgrib(nbytes+2) = chrtmp(2)
  mgrib(nbytes+3) = chrtmp(3)
  mgrib(nbytes+4) = chrtmp(4)
  nbytes = nbytes + 4

  CALL flt2ibm(rhohail,chrtmp,wdlen*8)
  mgrib(nbytes+1) = chrtmp(1)
  mgrib(nbytes+2) = chrtmp(2)
  mgrib(nbytes+3) = chrtmp(3)
  mgrib(nbytes+4) = chrtmp(4)
  nbytes = nbytes + 4

  mgrib(nbytes+1) = CHAR(r0exp)
  mgrib(nbytes+2) = CHAR(r0mant)
  mgrib(nbytes+3) = CHAR(r0mant)
  mgrib(nbytes+4) = CHAR(r0mant)
  nbytes = nbytes + 4

  mgrib(nbytes+1) = CHAR(r0exp)
  mgrib(nbytes+2) = CHAR(r0mant)
  mgrib(nbytes+3) = CHAR(r0mant)
  mgrib(nbytes+4) = CHAR(r0mant)
  nbytes = nbytes + 4

  CALL flt2ibm(tstop,chrtmp,wdlen*8)
  mgrib(nbytes+1) = chrtmp(1)
  mgrib(nbytes+2) = chrtmp(2)
  mgrib(nbytes+3) = chrtmp(3)
  mgrib(nbytes+4) = chrtmp(4)
  nbytes = nbytes + 4

  CALL flt2ibm(thisdmp,chrtmp,wdlen*8)
  mgrib(nbytes+1) = chrtmp(1)
  mgrib(nbytes+2) = chrtmp(2)
  mgrib(nbytes+3) = chrtmp(3)
  mgrib(nbytes+4) = chrtmp(4)
  nbytes = nbytes + 4

  CALL flt2ibm(latitud,chrtmp,wdlen*8)
  mgrib(nbytes+1) = chrtmp(1)
  mgrib(nbytes+2) = chrtmp(2)
  mgrib(nbytes+3) = chrtmp(3)
  mgrib(nbytes+4) = chrtmp(4)
  nbytes = nbytes + 4

  CALL flt2ibm(ctrlat,chrtmp,wdlen*8)
  mgrib(nbytes+1) = chrtmp(1)
  mgrib(nbytes+2) = chrtmp(2)
  mgrib(nbytes+3) = chrtmp(3)
  mgrib(nbytes+4) = chrtmp(4)
  nbytes = nbytes + 4

  CALL flt2ibm(ctrlon,chrtmp,wdlen*8)
  mgrib(nbytes+1) = chrtmp(1)
  mgrib(nbytes+2) = chrtmp(2)
  mgrib(nbytes+3) = chrtmp(3)
  mgrib(nbytes+4) = chrtmp(4)
  nbytes = nbytes + 4

  IF ( totout /= 0 ) THEN
    mgrib(nbytes+1) = CHAR(nstyp)
    nbytes = nbytes + 1
!
!-----------------------------------------------------------------------
!
!  Reserve 20 integers for future compatibility.
!
!-----------------------------------------------------------------------
!
    mgrib(nbytes+1) = CHAR(prcout)
    nbytes = nbytes + 1

    mgrib(nbytes+1) = CHAR(radout)
    nbytes = nbytes + 1

    mgrib(nbytes+1) = CHAR(flxout)
    nbytes = nbytes + 1

    mgrib(nbytes+1) = CHAR(0)  ! snowcvr not output anymore
    nbytes = nbytes + 1

    mgrib(nbytes+1) = CHAR(snowout)
    nbytes = nbytes + 1

    mgrib(nbytes+1) = CHAR(nscalar)
    nbytes = nbytes + 1

    IF (P_QC > 0) THEN
      mgrib(nbytes+1) = CHAR(P_QC)
    ELSE
      mgrib(nbytes+1) = CHAR(0)
    END IF
    nbytes = nbytes + 1

    IF (P_QR > 0) THEN
      mgrib(nbytes+1) = CHAR(P_QR)
    ELSE
      mgrib(nbytes+1) = CHAR(0)
    END IF
    nbytes = nbytes + 1

    IF (P_QI > 0) THEN
      mgrib(nbytes+1) = CHAR(P_QI)
    ELSE
      mgrib(nbytes+1) = CHAR(0)
    END IF
    nbytes = nbytes + 1

    IF (P_QS > 0) THEN
      mgrib(nbytes+1) = CHAR(P_QS)
    ELSE
      mgrib(nbytes+1) = CHAR(0)
    END IF
    nbytes = nbytes + 1

    IF (P_QH > 0) THEN
      mgrib(nbytes+1) = CHAR(P_QH)
    ELSE
      mgrib(nbytes+1) = CHAR(0)
    END IF
    nbytes = nbytes + 1

    DO i=12,20
      mgrib(nbytes+1) = CHAR(idummy)
      nbytes = nbytes + 1
    END DO

    DO i=0,16,4
      mgrib(nbytes+1) = CHAR(r0exp)
      mgrib(nbytes+2) = CHAR(r0mant)
      mgrib(nbytes+3) = CHAR(r0mant)
      mgrib(nbytes+4) = CHAR(r0mant)
      nbytes = nbytes + 4
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  If grdout=1 or grdbas=1, write out grid variables
!
!-----------------------------------------------------------------------
!
  IF(grdout == 1 .OR. grdbas == 1 ) THEN

    DO i=1,nx
      CALL flt2ibm(x(i),chrtmp,wdlen*8)
      mgrib(nbytes+1) = chrtmp(1)
      mgrib(nbytes+2) = chrtmp(2)
      mgrib(nbytes+3) = chrtmp(3)
      mgrib(nbytes+4) = chrtmp(4)
      nbytes = nbytes + 4
    END DO

    DO j=1,ny
      CALL flt2ibm(y(j),chrtmp,wdlen*8)
      mgrib(nbytes+1) = chrtmp(1)
      mgrib(nbytes+2) = chrtmp(2)
      mgrib(nbytes+3) = chrtmp(3)
      mgrib(nbytes+4) = chrtmp(4)
      nbytes = nbytes + 4
    END DO

    DO k=1,nz
      CALL flt2ibm(z(k),chrtmp,wdlen*8)
      mgrib(nbytes+1) = chrtmp(1)
      mgrib(nbytes+2) = chrtmp(2)
      mgrib(nbytes+3) = chrtmp(3)
      mgrib(nbytes+4) = chrtmp(4)
      nbytes = nbytes + 4
    END DO

  END IF

  mgrib(1) = CHAR(MOD(nbytes/65536,256))
  mgrib(2) = CHAR(MOD(nbytes/256,  256))
  mgrib(3) = CHAR(MOD(nbytes,      256))

  WRITE (nchanl) (mgrib(i),i=1,nbytes)

  RETURN
END SUBROUTINE wrthishd
