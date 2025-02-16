!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE GRIBDEC                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE gribdec(npts,nlev,                                           &
           wdlen,ipds,igds,ibdshd,nbufsz,mgrib,                         &
           fvar,ivar, bds,ibds)
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
!    npts     Number of grid points
!    nlev     Number of levels
!
!    wdlen    Length of machine word
!
!    ipds     Integer array of PDS
!    igds     Integer array of GDS
!    ibdshd   Integer array of BDS header
!
!    nbufsz   Buffer size of GRIB message
!    msglen   Length of the GRIB message
!    mgrib    Buffer carrying the GRIB message
!
!  OUTPUT:
!
!    pds      PDS (GRIB Section 1).
!    gds      GDS (GRIB Section 3).
!    bds      BDS except the first 11 octets header
!    ibds     BDS's integer array
!
!    fvar     Floating array to be written into GRIB file
!    ivar     Integer  array to be written into GRIB file
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

  INTEGER :: npts               ! Number of grid points
  INTEGER :: nlev               ! Number of levels

  INTEGER :: wdlen              ! Length of machine word in bytes

  INTEGER :: ipds(25)
  INTEGER :: igds(25)
  INTEGER :: ibdshd(4)

  INTEGER :: nbufsz             ! Size of GRIB message array
  CHARACTER (LEN=1) :: mgrib(nbufsz)  ! GRIB message array

  REAL :: fvar(npts)         ! Floating array
  INTEGER :: ivar(npts)         ! Integer  array

  INTEGER :: ibds(nbufsz/4)     ! Identical to BDS
  CHARACTER (LEN=1) :: bds(nbufsz)    ! BDS
!
!-----------------------------------------------------------------------
!
!  Local GRIB parameters
!
!-----------------------------------------------------------------------
!
  INTEGER :: bdslen
  INTEGER :: pdslen
  INTEGER :: gdslen
  INTEGER :: msglen             ! Length of GRIB message

  CHARACTER (LEN=1) :: pds(28)        ! PDS
  CHARACTER (LEN=1) :: gds(42)        ! GDS
  CHARACTER (LEN=1) :: bdshd(11)      ! BDS header

  INTEGER :: ipdsin(25)
  INTEGER :: igdsin(25)
  INTEGER :: ibdshdin(4)

  INTEGER :: iver,bwidth,bscl
  INTEGER :: itema

  REAL :: vref, scale
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i
  INTEGER :: npos, iskip1, iskip2
!
!-----------------------------------------------------------------------
!
!  External function
!
!-----------------------------------------------------------------------
!
  INTEGER :: char2i
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  npos = 0

  iskip1 = 0
  iskip2 = 0
!
!-----------------------------------------------------------------------
!
!  Check if the record is a GRIB message
!
!-----------------------------------------------------------------------
!
  iver = char2i(mgrib(8))

  IF ( mgrib(1) /= 'G' .OR. mgrib(2) /= 'R' .OR.                        &
         mgrib(3) /= 'I' .OR. mgrib(4) /= 'B' ) THEN
    WRITE (6,'(a/a)')                                                   &
        'Error: this record is not GRIB message',                       &
        'Program stopped in GRIBDEC'
    CALL arpsstop('arpsstop called from gribdec mgrib error',1)
  ELSE IF ( iver /= 1 ) THEN
    WRITE (6,'(a,i1/a)')                                                &
        'Error: GRIB version ID was not correct, iver = ',iver,         &
        'Program stopped in GRIBDEC'
    CALL arpsstop('arpsstop called from gribdec iver error ',1)
  END IF

  msglen = char2i(mgrib(5))*65536                                        &
         + char2i(mgrib(6))*256                                          &
         + char2i(mgrib(7))

  npos = npos + 8
!
!-----------------------------------------------------------------------
!
!  Get IPDS - GRIB Section 1
!
!-----------------------------------------------------------------------
!
  pdslen = char2i(mgrib(npos+1))*65536                                   &
         + char2i(mgrib(npos+2))*256                                     &
         + char2i(mgrib(npos+3))

  IF ( pdslen /= 28 ) THEN
    WRITE (6,'(a,i5/a)')                                                &
        'PDS length is not equal to 28, pdslen = ',pdslen,              &
        'Program stopped in GRIBDEC.'
    CALL arpsstop('arpsstop called from gribdec pdslen error',1)
  END IF

  DO i=1,pdslen
    pds(i) = mgrib(npos+i)
  END DO

  npos = npos + pdslen

  CALL gtipds( pds, ipdsin )

  IF ( ipdsin(8) /= ipds(8) .OR. ipdsin(9) /= ipds(9) ) THEN
    WRITE (6,'(a)') 'ERROR, variable mismatch (ipds 8 or 9).'
    GO TO 150
  END IF

  IF ( ipds(25) /= 0 ) THEN
    WRITE (6,'(a,i6,a)')                                                &
        'Incorrect decimal scale factor, D-scale = ',ipds(25),          &
        'ARPS GRIB dump should always has zero D-scale.'
    GO TO 150
  END IF

  DO i=1,25
    IF ( ipdsin(i) /= ipds(i) ) THEN
      WRITE (6,'(a,i2)')                                                &
          'WARNING: A difference was found in IPDS(i): i = ',i
    END IF
  END DO

  GO TO 170

  150   CONTINUE
  WRITE (6,'(a)')                                                       &
      ' i        ipds      ipdsin'
  DO i=1,25
    WRITE (6,'(i2,4x,i8,4x,i8)') i, ipds(i),ipdsin(i)
  END DO
  WRITE (6,*) 'ERROR, critical mismatch found in ',                     &
              'grib file, subroutine GRIBDEC, ABORTING.'
  CALL arpsstop('arpsstop called from gribdec file mismatch error',1)

  170   CONTINUE
!
!-----------------------------------------------------------------------
!
!  Get IGDS - GRIB Section 2
!
!-----------------------------------------------------------------------
!
  gdslen = char2i(mgrib(npos+1))*65536                                   &
         + char2i(mgrib(npos+2))*256                                     &
         + char2i(mgrib(npos+3))

  DO i=1,gdslen
    gds(i) = mgrib(npos+i)
  END DO

  npos = npos + gdslen

  CALL gtigds( gds,igdsin )

!  DO 210 i=1,5
!    IF ( igdsin(i).ne.igds(i) ) THEN
!      write (6,'(a,i2)')
!    :    'Error: A difference was found in IGDS(i): i = ',i
!      GOTO 230
!    ENDIF
!210  CONTINUE

  GO TO 250

!  230   CONTINUE

!  DO i=1,25
!    WRITE (6,'(i2,4x,i8,4x,i8)') i, igds(i),igdsin(i)
!  END DO

!  WRITE (6,'(a)') 'Program stopped in GRIBDEC.'
!  CALL arpsstop("arpsstop called from gribdec at 230",1)

  250   CONTINUE

!  DO 220 i=6,25
!    IF ( igdsin(i).ne.igds(i) ) THEN
!      write (6,'(a/a)')
!    :    'Warning: Data read in may contain errors.',
!    :    'A difference was found between IGDS and read-in IGDS:'
!      write (6,'(i2,4x,i8,4x,i8)') i, igds(i),igdsin(i)
!    ENDIF
!220  CONTINUE

!
!-----------------------------------------------------------------------
!
!  Making BDS - GRIB Section 4
!
!-----------------------------------------------------------------------
!
  DO i=1,11
    bdshd(i) = mgrib(npos+i)
  END DO

  npos = npos + 11

  bdslen = char2i(bdshd(1))*65536                                        &
         + char2i(bdshd(2))*256                                          &
         + char2i(bdshd(3))

  itema = char2i(bdshd(4))
  ibdshdin(1) = ishft(IAND(itema,128),-7)
  ibdshdin(2) = ishft(IAND(itema, 64),-6)
  ibdshdin(3) = ishft(IAND(itema, 32),-5)
  ibdshdin(4) = ishft(IAND(itema, 16),-4)

  DO i=1,4
    IF ( ibdshdin(i) /= ibdshd(i) ) THEN
      WRITE (6,'(a,i2.2,a,i4,a,i2.2,a,i4/a)')                           &
          'Error in read-in IBDSHD: ibdshdin(',i,') = ', ibdshdin(i),   &
          ', ibdshd(',i,') = ', ibdshd(i),                              &
          'Program stopped in GRIBDEC'
      CALL arpsstop('arpsstop called from gribdec in reading ibdshd',1)
    END IF
  END DO

  bscl = char2i(bdshd(5))*256 + char2i(bdshd(6))
  IF ( bscl > 32768 ) bscl = 32768 - bscl

  IF ( ibdshd(3) == 1 .AND. bscl /= 0 ) THEN
    WRITE (6,'(a,a,i2,a,i6)')                                           &
        'The binary scale factor should be zero for ',                  &
        'original integer data. itype = ',ibdshd(3),', b-scale = ',bscl, &
        'Program stopped in GRIBDEC'
    CALL arpsstop('arpsstop called from gribdec with bscl ',1)
  END IF

  scale = 2.**bscl

  CALL ibm2flt( bdshd(7), vref )
!  write (6,'(a,e20.12)')
!    :'The reference value is Vref = ',vref

  bwidth = char2i(bdshd(11))

  DO i=1,bdslen-11
    bds(i) = mgrib(npos+i)
  END DO

  npos = npos + bdslen - 11

  IF ( bwidth /= 0 ) THEN
    CALL grbgbytes( ibds,ivar,iskip1,bwidth,iskip2, npts )
  ELSE
    DO i=1,npts
      ivar(i) = 0
    END DO
  END IF

  IF ( ibdshd(3) == 0 ) THEN
    DO i=1,npts
      fvar(i) = FLOAT(ivar(i)) * scale + vref
    END DO
  ELSE
    DO i=1,npts
      ivar(i) = ivar(i) + nint(vref)
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  Check the end of GRIB message: '7777'
!
!-----------------------------------------------------------------------
!
  DO i=1,4
    IF ( mgrib(npos+i) /= '7' ) THEN
      WRITE (6,'(a)')                                                   &
          'Incorrect GRIB end flag: ',mgrib(npos+i),                    &
          'It should be alway character ''7''.',                        &
          'Program stopped in GRIBDEC'
      CALL arpsstop('arpsstop called from gribdec message 7777',1)
    END IF
  END DO

  npos = npos + 4

  IF ( npos /= msglen ) THEN
    WRITE (6,'(a/a,i8,a,i8)')                                           &
        'The GRIB message length was incorrect.',                       &
        'npos = ',npos, ' msglen = ',msglen
    CALL arpsstop('arpsstop called from gribdec with msglen',1)
  END IF

  RETURN
END SUBROUTINE gribdec
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE RDHISHD                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE rdhishd(nx,ny,nz,nzsoil,inch,grdbas, fmtverin, time, x,y,z,     &
           nbufsz,mgrib)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Read the header of ARPS history dump by using GRIB representation
!  to floating point numbers.
!
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  11/29/1995
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
!    inch     I/O unit of the history file
!
!    grdbas   Flag for grid and base state dump
!
!    x        x coordinate of grid points in physical/comp. space (m)
!    y        y coordinate of grid points in physical/comp. space (m)
!    z        z coordinate of grid points in computational space (m)
!
!    fmtver   Format and version identical
!
!    nbufsz   Buffer size of the GRIB message
!
!  OUTPUT:
!
!    mgrib    String containing the header of ARPS history dump
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
!  rainout=0 or 1. If rainout=0, rain variables are not dumped.
!  prcout =0 or 1. If prcout=0, precipitation rates are not dumped.
!  iceout =0 or 1. If iceout=0, qi, qs and qh are not dumped.
!  tkeout =0 or 1. If tkeout=0, tke is not dumped.
!  trbout =0 or 1. If trbout=0, kmh and kmv are not dumped.
!  sfcout =0 or 1. If sfcout=0, surface variables are not dumped.
!  landout=0 or 1. If landout=0, surface propertty arrays are not dumped.
!  ...
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

  INTEGER :: inch              ! I/O unit

  INTEGER :: grdbas            ! Flag for grid and base state dump

  CHARACTER (LEN=40), INTENT(OUT) :: fmtverin

  REAL :: x(nx)                ! The x-coord. of the physical and
                               ! computational grid. Defined at u-point.
  REAL :: y(ny)                ! The y-coord. of the physical and
                               ! computational grid. Defined at v-point.
  REAL :: z(nz)                ! The z-coord. of the computational grid.
                               ! Defined at w-point on the staggered grid.

!  CHARACTER (LEN=40) :: fmtver

  INTEGER :: nbufsz
  CHARACTER (LEN=1) :: mgrib(nbufsz) ! Header of ARPS GRIB file

  REAL :: time
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nxin,nyin,nzin
  INTEGER :: nzsoilin

  INTEGER :: npos,hhdlen

  INTEGER :: i,j,k

!  integer idummy
!  real    rdummy
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
  INCLUDE 'indtflg.inc'
!
!-----------------------------------------------------------------------
!
!  External function
!
!-----------------------------------------------------------------------
!
  INTEGER :: char2i
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  READ (inch) (mgrib(i),i=1,3)

  npos = 0

  hhdlen = char2i(mgrib(npos+1))*65536                                   &
         + char2i(mgrib(npos+2))*256                                     &
         + char2i(mgrib(npos+3))

  BACKSPACE (inch)

  READ (inch) (mgrib(i),i=1,hhdlen)

  npos = npos + 3

  DO i=1,40
    fmtverin(i:i) = mgrib(npos+i)
  END DO

!  IF ( fmtverin /= fmtver ) THEN
!    WRITE (6,'(/1x,a,/1x,2a,/1x,3a)')                                   &
!        'Data format incompatible with the data reader.',               &
!        'Format of data is ',fmtverin,' Format of reader is ',fmtver,   &
!        '. Job stopped.'
!    CALL arpsstop('arpsstop called from rdhishd with fmtver',1)
!  END IF

  npos = npos + 40

  DO i=1,80
    runname(i:i) = mgrib(npos+i)
  END DO

  npos = npos + 80

  nocmnt = char2i(mgrib(npos+1))
  npos = npos + 1

  IF( nocmnt > 0 ) THEN
    DO j=1,nocmnt
      DO  i=1,80
        cmnt(j)(i:i) = mgrib(npos+i)
      END DO
      npos = npos + 80
    END DO
  END IF

  CALL ibm2flt( mgrib(npos+1), time )
  npos = npos + 4

  nxin = char2i(mgrib(npos+1))*256 + char2i(mgrib(npos+2))
  npos = npos + 2
  nyin = char2i(mgrib(npos+1))*256 + char2i(mgrib(npos+2))
  npos = npos + 2
  nzin = char2i(mgrib(npos+1))*256 + char2i(mgrib(npos+2))
  npos = npos + 2
  nzsoilin = char2i(mgrib(npos+1))*256 + char2i(mgrib(npos+2))
  npos = npos + 2

  IF( nxin /= nx .OR. nyin /= ny .OR. nzin /= nz .OR.                   &
      nzsoilin /= nzsoil) THEN
    WRITE (6,'(1x,a/1x,4(a,I4)/1x,4(a,I4)/1x,a)')                       &
        'Dimensions in RDHISHD inconsistent with data.',                &
        'nx   = ',nx,' ny   = ',ny,' nz   = ',nz,' nzsoil   = ', nzsoil,&
        'nxin = ',nxin,' nyin = ',nyin,' nzin = ',nzin,                 &
        ' nzsoilin   = ', nzsoilin,                                     &
        'Program aborted in RDHISHD.'
    CALL arpsstop('arpsstop called from rdhishd with nxin...',1)
  END IF

  grdin = char2i(mgrib(npos+1))
  basin = char2i(mgrib(npos+2))
  varin = char2i(mgrib(npos+3))
  mstin = char2i(mgrib(npos+4))
  icein = char2i(mgrib(npos+5))
  trbin = char2i(mgrib(npos+6))
  sfcin = char2i(mgrib(npos+7))
  rainin = char2i(mgrib(npos+8))
  landin = char2i(mgrib(npos+9))
  totin = char2i(mgrib(npos+10))
  tkein = char2i(mgrib(npos+11))

  mapproj = char2i(mgrib(npos+14))
  month = char2i(mgrib(npos+15))
  day = char2i(mgrib(npos+16))
  year = char2i(mgrib(npos+17))*256 + char2i(mgrib(npos+18))
  hour = char2i(mgrib(npos+19))
  minute = char2i(mgrib(npos+20))
  second = char2i(mgrib(npos+21))

  npos = npos + 21

  CALL ibm2flt( mgrib(npos+1), umove )
  npos = npos + 4

  CALL ibm2flt( mgrib(npos+1), vmove )
  npos = npos + 4

  CALL ibm2flt( mgrib(npos+1), xgrdorg )
  npos = npos + 4

  CALL ibm2flt( mgrib(npos+1), ygrdorg )
  npos = npos + 4

  CALL ibm2flt( mgrib(npos+1), trulat1 )
  npos = npos + 4

  CALL ibm2flt( mgrib(npos+1), trulat2 )
  npos = npos + 4

  CALL ibm2flt( mgrib(npos+1), trulon )
  npos = npos + 4

  CALL ibm2flt( mgrib(npos+1), sclfct )
  npos = npos + 4

  CALL ibm2flt( mgrib(npos+1), ntcloud )
  npos = npos + 4

  CALL ibm2flt( mgrib(npos+1), n0rain )
  npos = npos + 4

  CALL ibm2flt( mgrib(npos+1), n0snow )
  npos = npos + 4

  CALL ibm2flt( mgrib(npos+1), n0grpl )
  npos = npos + 4

  CALL ibm2flt( mgrib(npos+1), n0hail )
  npos = npos + 4

  CALL ibm2flt( mgrib(npos+1), rhoice )
  npos = npos + 4

  CALL ibm2flt( mgrib(npos+1), rhosnow )
  npos = npos + 4

  CALL ibm2flt( mgrib(npos+1), rhogrpl )
  npos = npos + 4

  CALL ibm2flt( mgrib(npos+1), rhohail )
  npos = npos + 4

  CALL ibm2flt( mgrib(npos+1), alpharain )
  npos = npos + 4

  CALL ibm2flt( mgrib(npos+1), alphaice )
  npos = npos + 4

  CALL ibm2flt( mgrib(npos+1), alphasnow )
  npos = npos + 4

  CALL ibm2flt( mgrib(npos+1), alphagrpl )
  npos = npos + 4

  CALL ibm2flt( mgrib(npos+1), alphahail )
  npos = npos + 4

!  CALL ibm2flt( mgrib(npos+1), rdummy )
  npos = npos + 4

!  CALL ibm2flt( mgrib(npos+1), rdummy )
  npos = npos + 4

  CALL ibm2flt( mgrib(npos+1), tstop )
  npos = npos + 4

  CALL ibm2flt( mgrib(npos+1), thisdmp )
  npos = npos + 4

  CALL ibm2flt( mgrib(npos+1), latitud )
  npos = npos + 4

  CALL ibm2flt( mgrib(npos+1), ctrlat )
  npos = npos + 4

  CALL ibm2flt( mgrib(npos+1), ctrlon )
  npos = npos + 4

  IF ( totin /= 0 ) THEN
    nstyp = char2i(mgrib(npos+1))
    IF ( nstyp < 0 ) THEN
      nstyp = 1
    END IF
    npos  = npos + 1
!
!-----------------------------------------------------------------------
!
!  Reserve 20 integers for future compatibility.
!
!-----------------------------------------------------------------------
!
    prcin = char2i(mgrib(npos+1))
    npos  = npos + 1

    radin = char2i(mgrib(npos+1))
    npos  = npos + 1

    flxin = char2i(mgrib(npos+1))
    npos  = npos + 1

    snowcin = char2i(mgrib(npos+1))
    npos  = npos + 1

    snowin = char2i(mgrib(npos+1))
    npos  = npos + 1

    nscalarin = char2i(mgrib(npos+1))
    npos  = npos + 1

    p_qcin = char2i(mgrib(npos+1))
    npos  = npos + 1

    p_qrin = char2i(mgrib(npos+1))
    npos  = npos + 1

    p_qiin = char2i(mgrib(npos+1))
    npos  = npos + 1

    p_qsin = char2i(mgrib(npos+1))
    npos  = npos + 1

    p_qhin = char2i(mgrib(npos+1))
    npos  = npos + 1

    DO i=12,20
!      idummy = char2i(mgrib(npos+1))
      npos = npos + 1
    END DO

    DO i=0,16,4
!      CALL ibm2flt( mgrib(npos+i+1), rdummy )
      npos = npos + 4
    END DO
  END IF
!
!-----------------------------------------------------------------------
!
!  If grdout=1 or grdbas=1, write out grid variables
!
!-----------------------------------------------------------------------
!
  IF( grdin == 1 .OR. grdbas == 1 ) THEN

    DO i=1,nx
      CALL ibm2flt( mgrib(npos+1), x(i) )
      npos = npos + 4
    END DO

    DO j=1,ny
      CALL ibm2flt( mgrib(npos+1), y(j) )
      npos = npos + 4
    END DO

    DO k=1,nz
      CALL ibm2flt( mgrib(npos+1), z(k) )
      npos = npos + 4
    END DO

  END IF

  IF ( npos /= hhdlen ) THEN
    WRITE (6,'(1x,a,i6,a,i6/a)')                                        &
        'Length of header was incorrect. npos = ',npos,                 &
        ', hhdlen = ',hhdlen,                                           &
        'Program stopped in GTHISHD'
    CALL arpsstop('arpsstop called from ghishd with npos ',1)
  END IF

  RETURN
END SUBROUTINE rdhishd
!
!##################################################################
!##################################################################
!######                                                      ######
!######               SUBROUTINE IBM2FLT                     ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE ibm2flt( cibm, flt )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Convert GRIB IBM floating representation to machine floating
!  number.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Yuhe Liu
!  11/29/1995
!  MODIFICATIONS:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    cibm      4 bytes character string containing GRIB IBM format
!              floating representation
!
!  OUTPUT:
!
!    flt       Machine dependent floating number
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  CHARACTER (LEN=1) :: cibm(4)

  REAL :: flt
!
!-----------------------------------------------------------------------
!
!  Local variables
!
!-----------------------------------------------------------------------
!
  REAL :: tema

  INTEGER :: ISIGN
  INTEGER :: iexp,imant
!
!-----------------------------------------------------------------------
!
!  External function
!
!-----------------------------------------------------------------------
!
  INTEGER :: char2i
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  iexp = char2i(cibm(1))
  IF ( iexp > 127 ) THEN
    iexp = iexp - 128
    ISIGN = -1
  ELSE
    ISIGN = 1
  END IF

  imant = char2i(cibm(2))*65536                                          &
        + char2i(cibm(3))*256                                            &
        + char2i(cibm(4))

  tema = 16.0**(iexp-70)

  flt = FLOAT(ISIGN*imant)*tema

  RETURN

END SUBROUTINE ibm2flt
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GTIPDS                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE gtipds( pds,ipds )
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
!  11/30/1995
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    pds      PDS (GRIB Section 1).
!
!  OUTPUT:
!
!    ipds     Integer array of PDS
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  CHARACTER (LEN=1) :: pds(28)     ! PDS

  INTEGER :: ipds(25)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i
!
!-----------------------------------------------------------------------
!
!  External function
!
!-----------------------------------------------------------------------
!
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
!  Initialize array ipds
!
!-----------------------------------------------------------------------
!
  DO i=1,25
    ipds(i) = 0
  END DO

!
!-----------------------------------------------------------------------
!
!  Put pds into ipds
!
!-----------------------------------------------------------------------
!
  ipds(1) = char2i(pds(1))*65536                                         &
          + char2i(pds(2))*256                                           &
          + char2i(pds(3))

  ipds(2) = char2i(pds(4))
  ipds(3) = char2i(pds(5))
  ipds(4) = char2i(pds(6))
  ipds(5) = char2i(pds(7))

  i = char2i(pds(8))
  ipds(6) = ishft(IAND(i,128),-7)
  ipds(7) = ishft(IAND(i, 64),-6)

  ipds(8) = char2i(pds(9))
  ipds(9) = char2i(pds(10))
!
!-----------------------------------------------------------------------
!
!  check if level is in two octets or one
!
!-----------------------------------------------------------------------
!
  i = ipds(9)
  IF ( (i >= 1.AND.i <= 100).OR.i == 102.OR.                            &
          i == 103.OR.i == 105.OR.i == 107.OR.                          &
          i == 109.OR.i == 111.OR.i == 113.OR.                          &
          i == 115.OR.i == 125.OR.i == 160.OR.                          &
          i == 200.OR.i == 201 ) THEN
    ipds(11) = char2i(pds(11))*256                                       &
             + char2i(pds(12))        ! one level stored two bytes
    IF ( ipds(11) > 32768 ) THEN
      ipds(11) = 32768 - ipds(11)    ! positive integer for the abs(level)
    END IF
  ELSE
    ipds(10) = char2i(pds(11))
    ipds(11) = char2i(pds(12))
  END IF

  ipds(12) = char2i(pds(13))
  ipds(13) = char2i(pds(14))
  ipds(14) = char2i(pds(15))
  ipds(15) = char2i(pds(16))
  ipds(16) = char2i(pds(17))
  ipds(17) = char2i(pds(18))
!
!-----------------------------------------------------------------------
!
!  Check if time P1 is in two octets or one
!
!-----------------------------------------------------------------------
!
  ipds(20) = char2i(pds(21))

  IF (ipds(20) == 10) THEN
    ipds(18) = char2i(pds(19))*256                                       &
             + char2i(pds(20))
  ELSE
    ipds(18) = char2i(pds(19))
    ipds(19) = char2i(pds(20))
  END IF

  ipds(21) = char2i(pds(22))*256                                         &
           + char2i(pds(23))

  ipds(22) = char2i(pds(24))
  ipds(23) = char2i(pds(25))
  ipds(24) = char2i(pds(26))

  ipds(25) = char2i(pds(27))*256                                         &
           + char2i(pds(28))
  IF ( ipds(25) > 32768 ) THEN
    ipds(25) = 32768 - ipds(25)    ! negative value
  END IF

  RETURN
END SUBROUTINE gtipds
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE GTIGDS                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE gtigds(gds,igds)
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
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    gds      GDS
!
!  OUTPUT:
!
!    igds     Integer array of GDS
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  CHARACTER (LEN=1) :: gds(42)     ! GDS

  INTEGER :: igds(25)
!
!-----------------------------------------------------------------------
!
!  Local GRIB parameters
!
!-----------------------------------------------------------------------
!
  INTEGER :: i
!
!-----------------------------------------------------------------------
!
!  External function
!
!-----------------------------------------------------------------------
!
  INTEGER :: char2i
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Initialize array igds
!
!-----------------------------------------------------------------------
!
  DO i=1,25
    igds(i) = 0
  END DO

!
!-----------------------------------------------------------------------
!
!  Put GDS length in bytes 1,2,3
!
!-----------------------------------------------------------------------
!
  igds(1) = char2i(gds(4))    ! nz use one byte
  igds(2) = char2i(gds(5))    ! = 255, N/A for PV & PL
  igds(3) = char2i(gds(6))    ! data representation
  igds(8) = char2i(gds(17))   ! resolution and u&v component flag
!
!-----------------------------------------------------------------------
!
!  ARPS only support 3 types of projection
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Mercator projection
!
!-----------------------------------------------------------------------
!
  IF ( igds(3) == 1 ) THEN
    igds(4) = char2i(gds(7))*256 + char2i(gds(8))  ! nx
    igds(5) = char2i(gds(9))*256 + char2i(gds(10)) ! ny

    igds(6) = char2i(gds(11))*65536   & ! latitude of 1st point
            + char2i(gds(12))*256                                      &
            + char2i(gds(13))
    IF ( igds(6) > 8388608 ) THEN
      igds(6) = 8388608 - igds(6)    ! negative value
    END IF

    igds(7) = char2i(gds(14))*65536   & ! longitude of 1st point
            + char2i(gds(15))*256                                      &
            + char2i(gds(16))
    IF ( igds(7) > 8388608 ) THEN
      igds(7) = 8388608 - igds(7)    ! negative value
    END IF

    igds(9) = char2i(gds(18))*65536   & ! latitude of last point
            + char2i(gds(19))*256                                      &
            + char2i(gds(20))
    IF ( igds(9) > 8388608 ) THEN
      igds(9) = 8388608 - igds(9)    ! negative value
    END IF

    igds(10) = char2i(gds(21))*65536  & ! longitude of last point
         + char2i(gds(22))*256                                           &
             + char2i(gds(23))
    IF ( igds(10) > 8388608 ) THEN
      igds(10) = 8388608 - igds(10)  ! negative value
    END IF

    igds(11) = char2i(gds(32))*65536  & ! d(lat) in milidegree
         + char2i(gds(33))*256                                           &
             + char2i(gds(34))

    igds(12) = char2i(gds(29))*65536  & ! d(lon) in milidegree
         + char2i(gds(30))*256                                           &
             + char2i(gds(31))

    igds(13) = char2i(gds(24))*65536  & ! true latitude
         + char2i(gds(25))*256                                           &
             + char2i(gds(26))

    igds(14) = char2i(gds(28))
!
!-----------------------------------------------------------------------
!
!  Lambert Conformal projection
!
!-----------------------------------------------------------------------
!
  ELSE IF ( igds(3) == 3 ) THEN
    igds(4) = char2i(gds(7))*256      & ! nx
        + char2i(gds(8))
    igds(5) = char2i(gds(9))*256      & ! ny
        + char2i(gds(10))

    igds(6) = char2i(gds(11))*65536                                      &
            + char2i(gds(12))*256                                        &
            + char2i(gds(13))
    IF ( igds(6) > 8388608 ) THEN
      igds(6) = 8388608 - igds(6)    ! negative value
    END IF

    igds(7) = char2i(gds(14))*65536                                      &
            + char2i(gds(15))*256                                        &
            + char2i(gds(16))
    IF ( igds(7) > 8388608 ) THEN
      igds(7) = 8388608 - igds(7)    ! negative value
    END IF

    igds(9) = char2i(gds(18))*65536                                      &
            + char2i(gds(19))*256                                        &
            + char2i(gds(20))
    IF ( igds(9) > 8388608 ) THEN
      igds(9) = 8388608 - igds(9)    ! negative value
    END IF

    igds(10) = char2i(gds(21))*65536  & ! dx
         + char2i(gds(22))*256                                           &
             + char2i(gds(23))
    igds(11) = char2i(gds(24))*65536  & ! dy
         + char2i(gds(25))*256                                           &
             + char2i(gds(26))

    igds(12) = char2i(gds(27))
    igds(13) = char2i(gds(28))

    igds(15) = char2i(gds(29))*65536  & ! 1st true latitude
         + char2i(gds(30))*256                                           &
             + char2i(gds(31))
    igds(16) = char2i(gds(32))*65536  & ! 2nd true latitude
         + char2i(gds(33))*256                                           &
             + char2i(gds(34))

    igds(17) = char2i(gds(35))*65536  & ! latitude of south pole
         + char2i(gds(36))*256                                           &
             + char2i(gds(37))
    IF ( igds(17) > 8388608 ) THEN
      igds(17) = 8388608 - igds(17)    ! negative value
    END IF

    igds(18) = char2i(gds(38))*65536  & ! longitude of south pole
         + char2i(gds(39))*256                                           &
             + char2i(gds(40))
    IF ( igds(18) > 8388608 ) THEN
      igds(18) = 8388608 - igds(18)    ! negative value
    END IF
!
!-----------------------------------------------------------------------
!
!  Polar stereographic projection
!
!-----------------------------------------------------------------------
!
  ELSE IF ( igds(3) == 5 ) THEN
    igds(4) = char2i(gds(7))*256      & ! nx
        + char2i(gds(8))
    igds(5) = char2i(gds(9))*256      & ! ny
        + char2i(gds(10))

    igds(6) = char2i(gds(11))*65536   & ! latitude of 1st point
        + char2i(gds(12))*256                                            &
            + char2i(gds(13))
    IF ( igds(6) > 8388608 ) THEN
      igds(6) = 8388608 - igds(6)    ! negative value
    END IF

    igds(7) = char2i(gds(14))*65536   & ! longitude of 1st point
        + char2i(gds(15))*256                                            &
            + char2i(gds(16))
    IF ( igds(7) > 8388608 ) THEN
      igds(7) = 8388608 - igds(7)    ! negative value
    END IF

    igds(9) = char2i(gds(18))*65536   & ! true longitude
        + char2i(gds(19))*256                                            &
            + char2i(gds(20))
    IF ( igds(9) > 8388608 ) THEN
      igds(9) = 8388608 - igds(9)    ! negative value
    END IF

    igds(10) = char2i(gds(21))*65536  & ! dx
         + char2i(gds(22))*256                                           &
             + char2i(gds(23))
    igds(11) = char2i(gds(24))*65536  & ! dy
         + char2i(gds(25))*256                                           &
             + char2i(gds(26))

    igds(12) = char2i(gds(27))
    igds(13) = char2i(gds(28))

  END IF

  RETURN
END SUBROUTINE gtigds
