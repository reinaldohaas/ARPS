










!
!  This file contains Linux system specific routines
!  and functions.
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ASNCTL                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE asnctl (string, i , ierr)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Dummy subroutine to substitute Cray function ASNCTL.
!
!-----------------------------------------------------------------------
!
!  Author: Ming Xue
!  Date: 11/30/93
!
!  MODIFICATION HISTORY:
!
!  9/1/94 (Y. Lu)
!  Cleaned up documentation.
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
  INTEGER :: i,ierr
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  RETURN
END SUBROUTINE asnctl
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ASNUNIT                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE asnunit(nchan, string, ierr )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Dummy subroutine to substitute Cray function ASNUNIT.
!
!-----------------------------------------------------------------------
!
!  Author: Ming Xue
!  Date: 11/30/93
!
!  MODIFICATION HISTORY:
!
!  9/1/94 (Y. Lu)
!  Cleaned up documentation.
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

  CHARACTER (LEN=*) :: string
  INTEGER :: nchan, ierr
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  RETURN
END SUBROUTINE asnunit
!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE ASNFILE                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE asnfile(FILE, string, ierr )
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Dummy subroutine to substitute Cray function ASNFILE.
!
!-----------------------------------------------------------------------
!
!  Author: Ming Xue
!  Date: 11/30/93
!
!  MODIFICATION HISTORY:
!
!  9/1/94 (Y. Lu)
!  Cleaned up documentation.
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

  CHARACTER (LEN=*) :: FILE
  CHARACTER (LEN=*) :: string
  INTEGER :: ierr
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  RETURN
END SUBROUTINE asnfile
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE UNIXCMD                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE unixcmd(cmd)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  To executable a system command by making a system call.
!
!-----------------------------------------------------------------------
!
!  Author: Ming Xue
!  Date: 4/15/95
!
!  MODIFICATION HISTORY:
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

  CHARACTER (LEN=*) :: cmd

  INTEGER :: ireturn, system
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!


   !CALL system( cmd )
   ireturn = system ( cmd )


  RETURN
END SUBROUTINE unixcmd

!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CMPRS                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cmprs(filename)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Make a system call to compress the file 'filename' with compress
!  format developed by Free Software Foundation.
!
!-----------------------------------------------------------------------
!
!  Author: Ming Xue
!  Date: 11/30/93
!
!  MODIFICATION HISTORY:
!
!  9/1/94 (Y. Lu)
!  Cleaned up documentation.
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
  CHARACTER (LEN=*)   :: filename
  CHARACTER (LEN=256) :: CHAR
  INTEGER :: lenstr
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  lenstr = LEN( filename )
  CALL strlnth( filename, lenstr)

  !CHAR = 'compress '
  CHAR =  'gzip -f  '
  CHAR(10:9+lenstr) = filename

  IF( 9+lenstr > 256) THEN
    WRITE(6,'(1x,a)')                                                   &
        'Work character char too small in CMPRS, call returned.'
    RETURN
  END IF

  WRITE(6,'(1x,a,a,a)') 'Compressing file ',filename,' ...'
  CALL unixcmd ( CHAR(1:9+lenstr) )

  RETURN
END SUBROUTINE cmprs
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE CMPRSGZ                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE cmprsgz(filename)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Make a system call to compress the file 'filename' with gzip
!  compress format developed by Free Software Foundation.
!
!-----------------------------------------------------------------------
!
!  Author: Ming Xue
!  Date: 11/30/93
!
!  MODIFICATION HISTORY:
!
!  9/1/94 (Y. Lu)
!  Cleaned up documentation.
!
!  07/26/2005 (K. Brewster)
!  Created cmprsgz from a copy of cmprs
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

  CHARACTER (LEN=*)   :: filename
  CHARACTER (LEN=256) :: CHAR
  INTEGER :: lenstr
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  lenstr = LEN( filename )
  CALL strlnth( filename, lenstr)

  IF( 5+lenstr > 256) THEN
    WRITE(6,'(1x,a)')                                                   &
        'Work character char too small in CMPRS, call returned.'
  END IF

  CHAR = 'gzip '
  CHAR(6:5+lenstr) = filename

  WRITE(6,'(1x,a,a,a)') 'Compressing file ',filename,' ...'
  CALL unixcmd ( CHAR(1:5+lenstr) )

  RETURN
END SUBROUTINE cmprsgz
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE UNCMPRS                    ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE uncmprs(filename)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Make a system call to uncompress the file 'filename' with gunzip
!  format developed by Free Software Foundation.
!
!-----------------------------------------------------------------------
!
!  Author: Ming Xue
!  Date: 11/30/93
!
!  MODIFICATION HISTORY:
!
!  9/1/94 (Y. Lu)
!  Cleaned up documentation.
!
!  11/27/1995 (M. Xue)
!  Decompression for both .Z and .gz files.
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

  CHARACTER (LEN=*)   :: filename
  CHARACTER (LEN=256) :: CHAR
  INTEGER :: lenstr
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  lenstr = LEN( filename )
  CALL strlnth( filename, lenstr)

  IF(filename(lenstr-1:lenstr) == '.Z') THEN
    CHAR = 'uncompress '
    CHAR(12:12+lenstr) = filename
    lenstr=lenstr+12
  ELSE IF(filename(lenstr-2:lenstr) == '.gz') THEN
    CHAR = 'gunzip '
    CHAR(8:7+lenstr) = filename
    lenstr = lenstr+8
  ELSE
    WRITE(6,'(1x,a,/1x,a)')                                             &
        'File name does not have the right affix.',                     &
        'No decompressing was done on file ',filename
    RETURN
  END IF

  IF( lenstr > 256) THEN
    WRITE(6,'(1x,a)')                                                   &
        'Work character char too small in UNCMPRS, call returned.'
  END IF

  WRITE(6,'(1x,a,a,a)') 'Decompressing file ',filename,' ...'
  CALL unixcmd( CHAR(1:lenstr) )

  RETURN
END SUBROUTINE uncmprs
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE GRBSBYTE                  ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE grbsbyte(iout,in,iskip,nbits)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Take rightmost nbits bit fields of IN and insert them in
!  bitstrings of IOUT.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dr. Robert C. Gammill, Consultant, NCAR
!  July 1972
!
!  MODIFICATIONS:
!
!  12/20/95 (Yuhe Liu)
!  Converted to ARPS standard format and added documents
!  05/06/03 (Kevin W. Thomas)
!  Code did not handle endian issues properly.
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    in      Integer array to be packed
!    iskip
!    nbits   bits number of packing
!
!  OUTPUT:
!
!    iout    Packed stream
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: in

  INTEGER :: nbits
  INTEGER :: iskip

  INTEGER :: iout(*)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: j, ii
  INTEGER :: icon, INDEX
  INTEGER :: mask, msk
  INTEGER :: movel

  INTEGER :: itemp

  INTEGER :: nbitsw
  DATA nbitsw/32/
!
!-----------------------------------------------------------------------
!
!  Masks table put in decimal so it will compile on any 32 bit
!  computer
!
!-----------------------------------------------------------------------
!
  INTEGER :: masks(32)
  DATA masks /          1,          3,          7,         15,          &
                       31,         63,        127,        255,          &
                      511,       1023,       2047,       4095,          &
                     8191,      16383,      32767,      65535,          &
                   131071,     262143,     524287,    1048575,          &
                  2097151,    4194303,    8388607,   16777215,          &
                 33554431,   67108863,  134217727,  268435455,          &
                536870911, 1073741823, 2147483647,         -1/
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
!  Nbits must be less than or equal to nbitsW
!
!-----------------------------------------------------------------------
!
  icon = nbitsw - nbits
  IF (icon < 0) RETURN
  mask = masks(nbits)
!
!-----------------------------------------------------------------------
!
!  INDEX TELLS HOW MANY WORDS INTO IOUT THE NEXT BYTE IS TO BE STORED.
!
!-----------------------------------------------------------------------
!
  INDEX  = iskip/nbitsw
!
!-----------------------------------------------------------------------
!
!  II tells how many bits in from the left side of the word to store
!  it.
!
!-----------------------------------------------------------------------
!
  ii = MOD(iskip,nbitsw)

  j = IAND(mask,in)
  movel = icon - ii
!
!-----------------------------------------------------------------------
!
!  Byte is to be stored in middle of word.  shift left.
!
!-----------------------------------------------------------------------
!
  IF (movel > 0) THEN
    msk = ishft(mask,movel)
    call flipper(iout(index+1))
    iout(INDEX+1) = ior(IAND(NOT(msk),iout(INDEX+1)),                   &
                        ishft(j,movel))
    call flipper(iout(index+1))
!
!-----------------------------------------------------------------------
!
!  The byte is to be split across a word break.
!
!-----------------------------------------------------------------------
!
  ELSE IF (movel < 0) THEN
    msk = masks(nbits+movel)
    call flipper(iout(index+1))
    call flipper(iout(index+2))
    iout(INDEX+1) = ior(IAND(NOT(msk),iout(INDEX+1)),                   &
                        ishft(j,movel))
    itemp = IAND(masks(nbitsw+movel),iout(INDEX+2))
    iout(INDEX+2) = ior(itemp,ishft(j,nbitsw+movel))
    call flipper(iout(index+1))
    call flipper(iout(index+2))
!
!-----------------------------------------------------------------------
!
!  Byte is to be stored right-adjusted.
!
!-----------------------------------------------------------------------
!
  ELSE
    call flipper(iout(index+1))
    iout(INDEX+1) = ior(IAND(NOT(mask),iout(INDEX+1)),j)
    call flipper(iout(index+1))
  END IF

  RETURN
END SUBROUTINE grbsbyte

!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE FLIPPER                   ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE flipper(value)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Handle Linux endian problems by flipping the bytes around.
!
!  Each time bytes are handle, the routines are called twice.  The first
!  time is to undo endian, so that new data can be inserted.  The second
!  time is to redo endian.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Kevin W. Thomas, CAPS
!  May 6, 2003
!
!  MODIFICATIONS:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    value   Integer to have bytes flipped
!
!  OUTPUT:
!
!    value   The answer.
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  INTEGER :: value
  CHARACTER (KIND=1) :: data(4), itmp
  INTEGER :: in
  EQUIVALENCE (data, in )

  in = value

  itmp = data(1)
  data(1) = data(4)
  data(4) = itmp

  itmp = data(2)
  data(2) = data(3)
  data(3) = itmp

  value = in
  return
END SUBROUTINE flipper

!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE GRBSBYTES                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE grbsbytes(iout,in,iskip,nbits,nskip,nwrd)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Take rightmost nbits bit fields from words of IN and insert them
!  consecutively in bitstrings of IOUT.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dr. Robert C. Gammill, Consultant, NCAR
!  July 1972
!
!  MODIFICATIONS:
!
!  12/05/95 (Yuhe Liu)
!  Converted to ARPS standard format and added documents
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nwrd    Number of word to be packed
!    in      Integer array to be packed
!    iskip
!    nbits   bits number of packing
!    nskip
!
!  OUTPUT:
!
!    iout    Packed stream
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nwrd
  INTEGER :: in(nwrd)

  INTEGER :: nbits
  INTEGER :: iskip
  INTEGER :: nskip

  INTEGER :: iout(nwrd)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, j, ii
  INTEGER :: ibits
  INTEGER :: icon, INDEX
  INTEGER :: mask, msk
  INTEGER :: movel, iwords, istep

  INTEGER :: itemp

  INTEGER :: nbitsw
  DATA nbitsw/32/
!
!-----------------------------------------------------------------------
!
!  Masks table put in decimal so it will compile on any 32 bit
!  computer
!
!-----------------------------------------------------------------------
!
  INTEGER :: masks(32)
  DATA masks /          1,          3,          7,         15,          &
                       31,         63,        127,        255,          &
                      511,       1023,       2047,       4095,          &
                     8191,      16383,      32767,      65535,          &
                   131071,     262143,     524287,    1048575,          &
                  2097151,    4194303,    8388607,   16777215,          &
                 33554431,   67108863,  134217727,  268435455,          &
                536870911, 1073741823, 2147483647,         -1/
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
!  Nbits must be less than or equal to nbitsW
!
!-----------------------------------------------------------------------
!
  icon = nbitsw - nbits
  IF (icon < 0) RETURN
  mask = masks(nbits)
!
!-----------------------------------------------------------------------
!
!  INDEX TELLS HOW MANY WORDS INTO IOUT THE NEXT BYTE IS TO BE STORED.
!
!-----------------------------------------------------------------------
!
  INDEX  = iskip/nbitsw
!
!-----------------------------------------------------------------------
!
!  II tells how many bits in from the left side of the word to store
!  it.
!
!-----------------------------------------------------------------------
!
  ii = MOD(iskip,nbitsw)
!
!-----------------------------------------------------------------------
!
!  Istep is the distance in bits from one byte position to the next.
!
!-----------------------------------------------------------------------
!
  istep  = nbits + nskip
!
!-----------------------------------------------------------------------
!
!  Iwords tells how many words to skip from one byte to the next.
!
!-----------------------------------------------------------------------
!
  iwords = istep / nbitsw
!
!-----------------------------------------------------------------------
!
!  Ibits tells how many bits to skip after skipping iwords.
!
!-----------------------------------------------------------------------
!
  ibits  = MOD(istep,nbitsw)

  DO i = 1,nwrd
    j = IAND(mask,in(i))
    movel = icon - ii
!
!-----------------------------------------------------------------------
!
!  Byte is to be stored in middle of word.  shift left.
!
!-----------------------------------------------------------------------
!
    IF (movel > 0) THEN
      msk = ishft(mask,movel)
      iout(INDEX+1) = ior(IAND(NOT(msk),iout(INDEX+1)),                 &
                          ishft(j,movel))
!
!-----------------------------------------------------------------------
!
!  The byte is to be split across a word break.
!
!-----------------------------------------------------------------------
!
    ELSE IF (movel < 0) THEN
      msk = masks(nbits+movel)
      iout(INDEX+1) = ior(IAND(NOT(msk),iout(INDEX+1)),                 &
                          ishft(j,movel))
      itemp = IAND(masks(nbitsw+movel),iout(INDEX+2))
      iout(INDEX+2) = ior(itemp,ishft(j,nbitsw+movel))
!
!-----------------------------------------------------------------------
!
!  Byte is to be stored right-adjusted.
!
!-----------------------------------------------------------------------
!
    ELSE
      iout(INDEX+1) = ior(IAND(NOT(mask),iout(INDEX+1)),j)
    END IF

    ii = ii + ibits
    INDEX = INDEX + iwords
    IF (ii >= nbitsw) THEN
      ii = ii - nbitsw
      INDEX = INDEX + 1
    END IF

  END DO

  CALL swap4byte(iout,(nwrd*nbits)/nbitsw)

  RETURN
END SUBROUTINE grbsbytes
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE GRBGBYTE                  ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE grbgbyte(in,iout,iskip,nbits)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Extract nwrd bitstrings, nbits bits long, and store them right
!  justified 0 fill, into IOUT.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dr. Robert C. Gammill, Consultant, NCAR
!  May 1972
!
!  MODIFICATIONS:
!
!  12/20/95 (Yuhe Liu)
!  Converted to ARPS standard format and added documents
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: in(*)  ! Integer array to be packed

  INTEGER, INTENT(IN) :: nbits  ! bits number of packing
  INTEGER, INTENT(IN) :: iskip

  INTEGER, INTENT(OUT) :: iout  ! Packed stream
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: ii, itmp(2)
  INTEGER :: icon, INDEX
  INTEGER :: mask
  INTEGER :: movel, mover

  INTEGER :: nbitsw
  DATA nbitsw/32/
!
!-----------------------------------------------------------------------
!
!  Masks table put in decimal so it will compile on any 32 bit
!  computer
!
!  DATA MASKS /Z'00000001',Z'00000003',Z'00000007',Z'0000000F',
!    :            Z'0000001F',Z'0000003F',Z'0000007F',Z'000000FF',
!    :            Z'000001FF',Z'000003FF',Z'000007FF',Z'00000FFF',
!    :            Z'00001FFF',Z'00003FFF',Z'00007FFF',Z'0000FFFF',
!    :            Z'0001FFFF',Z'0003FFFF',Z'0007FFFF',Z'000FFFFF',
!    :            Z'001FFFFF',Z'003FFFFF',Z'007FFFFF',Z'00FFFFFF',
!    :            Z'01FFFFFF',Z'03FFFFFF',Z'07FFFFFF',Z'0FFFFFFF',
!    :            Z'1FFFFFFF',Z'3FFFFFFF',Z'7FFFFFFF',Z'FFFFFFFF'/
!
!-----------------------------------------------------------------------
!
  INTEGER :: masks(32)
  DATA masks /        1,          3,          7,         15,            &
                     31,         63,        127,        255,            &
                    511,       1023,       2047,       4095,            &
                   8191,      16383,      32767,      65535,            &
                 131071,     262143,     524287,    1048575,            &
                2097151,    4194303,    8388607,   16777215,            &
               33554431,   67108863,  134217727,  268435455,            &
              536870911, 1073741823, 2147483647,         -1/
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
!  NBYTE must be less than or equal to nbitsw
!
!-----------------------------------------------------------------------
!
  icon = nbitsw - nbits
  IF (icon < 0) RETURN
  mask = masks(nbits)
!
!-----------------------------------------------------------------------
!
!  Index tells how many words into the array 'in' the next byte appears.
!
!-----------------------------------------------------------------------
!
  INDEX = iskip/nbitsw
!
!-----------------------------------------------------------------------
!
!  II tells how many bits the byte is from the left side of the word.
!
!-----------------------------------------------------------------------
!
  ii = MOD(iskip,nbitsw)
!
!-----------------------------------------------------------------------
!
!    MOVER specifies how far to the right a byte must be moved in order
!    to be right adjusted.
!
!-----------------------------------------------------------------------
!
  mover = icon - ii
!
!-----------------------------------------------------------------------
!
!  Right adjust the byte.
!
!-----------------------------------------------------------------------
!
  IF (mover > 0) THEN
    itmp(1) = in(INDEX+1)
    iout = IAND(ishft(itmp(1),-mover),mask)
!
!-----------------------------------------------------------------------
!
!  The byte is split across a word break.
!
!-----------------------------------------------------------------------
!
  ELSE IF (mover < 0) THEN
    movel = - mover
    mover = nbitsw - movel
    itmp(1) = in(INDEX+1)
    itmp(2) = in(INDEX+2)
    iout = IAND( ior(ishft(itmp(1),movel),                              &
                 ishft(itmp(2),-mover)), mask )
!
!-----------------------------------------------------------------------
!
!  The byte is already right adjusted.
!
!-----------------------------------------------------------------------
!
  ELSE
    itmp(1) = in(INDEX+1)
    iout = IAND(itmp(1),mask)
  END IF

  RETURN
END SUBROUTINE grbgbyte
!
!##################################################################
!##################################################################
!######                                                      ######
!######                 SUBROUTINE GRBGBYTES                 ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE grbgbytes(in,iout,iskip,nbits,nskip,nwrd)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Extract nwrd bitstrings, nbits bits long, and store them right
!  justified 0 fill, into successive words of IOUT.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Dr. Robert C. Gammill, Consultant, NCAR
!  May 1972
!
!  MODIFICATIONS:
!
!  12/05/95 (Yuhe Liu)
!  Converted to ARPS standard format and added documents
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!    nwrd    Number of word to be packed
!    in      Integer array to be packed
!    iskip
!    nbits   bits number of packing
!    nskip
!
!  OUTPUT:
!
!    iout    Packed stream
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: nwrd
  INTEGER :: in(nwrd)

  INTEGER :: nbits
  INTEGER :: iskip
  INTEGER :: nskip

  INTEGER :: iout(nwrd)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: i, ii
  INTEGER :: ibits
  INTEGER :: icon, INDEX
  INTEGER :: mask
  INTEGER :: movel, mover, iwords, istep

  INTEGER :: nbitsw
  DATA nbitsw/32/
!
!-----------------------------------------------------------------------
!
!  Masks table put in decimal so it will compile on any 32 bit
!  computer
!
!  DATA MASKS /Z'00000001',Z'00000003',Z'00000007',Z'0000000F',
!    :            Z'0000001F',Z'0000003F',Z'0000007F',Z'000000FF',
!    :            Z'000001FF',Z'000003FF',Z'000007FF',Z'00000FFF',
!    :            Z'00001FFF',Z'00003FFF',Z'00007FFF',Z'0000FFFF',
!    :            Z'0001FFFF',Z'0003FFFF',Z'0007FFFF',Z'000FFFFF',
!    :            Z'001FFFFF',Z'003FFFFF',Z'007FFFFF',Z'00FFFFFF',
!    :            Z'01FFFFFF',Z'03FFFFFF',Z'07FFFFFF',Z'0FFFFFFF',
!    :            Z'1FFFFFFF',Z'3FFFFFFF',Z'7FFFFFFF',Z'FFFFFFFF'/
!
!-----------------------------------------------------------------------
!
  INTEGER :: masks(32)
  DATA masks /        1,          3,          7,         15,            &
                     31,         63,        127,        255,            &
                    511,       1023,       2047,       4095,            &
                   8191,      16383,      32767,      65535,            &
                 131071,     262143,     524287,    1048575,            &
                2097151,    4194303,    8388607,   16777215,            &
               33554431,   67108863,  134217727,  268435455,            &
              536870911, 1073741823, 2147483647,         -1/
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
!  NBYTE must be less than or equal to nbitsw
!
!-----------------------------------------------------------------------
!
  icon = nbitsw - nbits
  IF (icon < 0) RETURN
  mask = masks(nbits)
!
!-----------------------------------------------------------------------
!
!  Index tells how many words into the array 'in' the next byte appears.
!
!-----------------------------------------------------------------------
!
  INDEX = iskip/nbitsw
!
!-----------------------------------------------------------------------
!
!  II tells how many bits the byte is from the left side of the word.
!
!-----------------------------------------------------------------------
!
  ii = MOD(iskip,nbitsw)
!
!-----------------------------------------------------------------------
!
!  ISTEP is the distance in bits from the start of one byte to the next.
!
!-----------------------------------------------------------------------
!
  istep = nbits + nskip
!
!-----------------------------------------------------------------------
!
!  IWORDS tells how many words to skip from one byte to the next.
!
!-----------------------------------------------------------------------
!
  iwords = istep / nbitsw
!
!-----------------------------------------------------------------------
!
!  IBITS tells how many bits to skip after skipping iwords.
!
!-----------------------------------------------------------------------
!
  ibits = MOD(istep,nbitsw)

  DO i = 1,nwrd
!
!-----------------------------------------------------------------------
!
!    MOVER specifies how far to the right a byte must be moved in order
!    to be right adjusted.
!
!-----------------------------------------------------------------------
!
    mover = icon - ii
!
!-----------------------------------------------------------------------
!
!  Right adjust the byte.
!
!-----------------------------------------------------------------------
!
    IF (mover > 0) THEN
      iout(i) = IAND(ishft(in(INDEX+1),-mover),mask)
!
!-----------------------------------------------------------------------
!
!  The byte is split across a word break.
!
!-----------------------------------------------------------------------
!
    ELSE IF (mover < 0) THEN
      movel = - mover
      mover = nbitsw - movel
      iout(i) = IAND( ior(ishft(in(INDEX+1),movel),                     &
                      ishft(in(INDEX+2),-mover)),mask )
!
!-----------------------------------------------------------------------
!
!  The byte is already right adjusted.
!
!-----------------------------------------------------------------------
!
    ! WDT RLC Handle potential error condition
    ELSE
      IF (INDEX >= nwrd) THEN
        WRITE(6,*) 'GRBGBYTES: ERROR: INDEX > nwrd ', INDEX, nwrd
      ELSE
        iout(i) = IAND(in(INDEX+1),mask)
      END IF
    END IF
!
!-----------------------------------------------------------------------
!
!  Increment ii and index.
!
!-----------------------------------------------------------------------
!
    ii = ii + ibits
    INDEX = INDEX + iwords

    IF (ii >= nbitsw) THEN
      ii = ii - nbitsw
      INDEX = INDEX + 1
    END IF

  END DO

  RETURN
END SUBROUTINE grbgbytes

REAL FUNCTION f_cputime()
!
!-----------------------------------------------------------------------
!
!  Function to measure the CPU time usage.
!
!-----------------------------------------------------------------------
!

!
!  Call cpu_time,     Intrinsic subroutine in Fortran 95
!                     tested with IFORT and G95
!
  REAL :: tscalar

  CALL cpu_time(tscalar)    ! A scalar REAL(INTENT(OUT)) as actual parameter

  IF (tscalar < 0)  THEN    ! A meaningful time cannot be returned
    f_cputime = 0.0
  ELSE
    f_cputime = tscalar
  END IF


  RETURN
END FUNCTION f_cputime

!
!##################################################################
!##################################################################
!######                                                      ######
!######                FUNCTION char2i                       ######
!######                                                      ######
!##################################################################
!##################################################################
!
INTEGER FUNCTION char2i(a)
!
!------------------------------------------------------------------
!
! PURPOSE:
!
!   Moves a bit string from a char*1 to integer
!   Intend to replace the intrinsic function ICHAR, which
!   only supports 0 <= ICHAR(a) <= 127 on the
!   IBM SP.  If "a" is greater than 127 in the collating sequence,
!   ICHAR(a) does not return the expected bit value.
!   This function can be used for all values 0 <= ICHAR(a) <= 255.
!
!   For all other platforms, especially those with little endian
!   coding, this function is just a wrapper function of the
!   intrinsic function ICHAR.
!
! MODIFICATION HISTORY:
!
! USAGE:     I = char2i(a)
!
!   INPUT ARGUMENT :
!
!          a - Character*1 variable that holds the bitstring to extract
!
!   RETURN ARGUMENT :
!
!          char2i - Integer value of the bitstring in character a
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  CHARACTER(LEN=1), INTENT(IN) :: a

  char2i = ICHAR(a)

  RETURN
END FUNCTION char2i

SUBROUTINE inquiredir(dirname,fexist)

  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN)  :: dirname
  LOGICAL,          INTENT(OUT) :: fexist

  INQUIRE(DIRECTORY=TRIM(dirname),EXIST=fexist)

  RETURN
END SUBROUTINE inquiredir

SUBROUTINE makedir(dirname,istatus)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN)  :: dirname
  INTEGER,          INTENT(OUT) :: istatus

!----------------------------------------------------------------------
  INTEGER :: lendir

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  lendir = LEN_TRIM(dirname)

  CALL c_mkdir(TRIM(dirname)//char(0), istatus)

  RETURN
END SUBROUTINE makedir
