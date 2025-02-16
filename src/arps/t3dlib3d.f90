!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE DUMMY                      ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE dummy
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Dummy subroutine to substitute Cray function ASNCTL.
!
!-----------------------------------------------------------------------
!

  RETURN
END SUBROUTINE dummy

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
  INTEGER :: ishell,istat
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  istat = ishell( cmd )

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
!  Make a system call to compress the file 'filename' with
!  using gzip or when gzip fails compress.
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
!  11/27/95 (M. Xue)
!  Options to compress using gzip and compress.
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
  INTEGER :: lenstr,ishell,istat
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  lenstr = LEN( filename )
  CALL strlnth( filename, lenstr)

  CHAR = 'gzip --fast '
  CHAR(13:12+lenstr) = filename

  IF( 12+lenstr > 256) THEN
    WRITE(6,'(1x,a)')                                                   &
        'Work character char too small in CMPRS, call returned.'
  END IF

  WRITE(6,'(1x,a,a,a)') 'Compressing file ',filename,' ...'
  istat = ishell( CHAR(1:12+lenstr) )

  IF( istat /= 0) THEN

    CHAR = 'compress '
    CHAR(10:9+lenstr) = filename
    istat = ishell( CHAR(1:9+lenstr) )

    IF( istat /= 0) WRITE(6,'(1x,a,a,a)')                               &
        'Compression of file ',filename(1:lenstr),                      &
        ' was unsucessful in CMPRS.'

  END IF


  RETURN
END SUBROUTINE cmprs


!
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
  INTEGER :: lenstr,lenfn,ishell,istat
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!

  lenfn = LEN( filename )
  CALL strlnth( filename, lenfn)
  lenstr=lenfn

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

  IF( lenstr > 256 ) THEN
    WRITE(6,'(1x,a)')                                                   &
        'Work character char too small in UNCMPRS, call returned.'
  END IF

  WRITE(6,'(1x,a,a,a)') 'Decompressing file ',filename,' ...'
  istat = ishell( CHAR(1:lenstr) )

  IF( istat /= 0) WRITE(6,'(1x,a,a,a)')                                 &
      'Compression of file ',filename(1:lenfn),                         &
      ' was unsucessful in CMPRS.'

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
!  SBYTE - Insert a number of bit fields. Cray routine.
!
!  Reverses the action of gbytes, taking fields from s and
!  inserting them into a bit string in d. see GBYTE.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: D. Robertson
!  Aug. 1981
!  MODIFICATIONS:
!
!  12/05/95 (Yuhe Liu)
!  Converted to ARPS standard format and added documents
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
  INTEGER :: nbpw
  PARAMETER (nbpw=64)

  INTEGER :: iw, id
  INTEGER :: iskip1, ibits
  INTEGER :: icon, INDEX
  INTEGER :: mask, msk
  INTEGER :: movel, iwords

  INTEGER :: itemp

  INTEGER :: sh1,sh2,sh3
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  id = 1 + iskip/nbpw
  iskip1 = MOD(iskip,nbpw)

  sh1 = iskip1 + nbits
  IF (sh1 > nbpw) GO TO 50

  sh2 = nbpw - sh1
  IF (sh2 < 0) sh2 = nbpw-sh2
!
!-----------------------------------------------------------------------
!
!  Byte goes into 1 word of iout.
!
!-----------------------------------------------------------------------
!
  iout(1) = shift( OR( AND(shift(iout(1),sh1),mask(nbpw-nbits)),        &
                       AND(in,shift(mask(nbits),nbits))),sh2)
  RETURN

  50    CONTINUE
!
!-----------------------------------------------------------------------
!
!  Byte goes into 2 words of iout.
!
!-----------------------------------------------------------------------
!
  sh3 = 2*nbpw-sh1
  iout(1) = OR( AND(iout(id),mask(iskip1)),                             &
                AND(shift(in,sh3),compl(mask(iskip1))) )

  iout(2) = OR( AND(iout(id+1),shift(compl(mask(sh1-nbpw)),nbpw)),      &
                  shift(AND(in,compl(mask(sh3))),sh3) )

  RETURN
END SUBROUTINE grbsbyte
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
!  SBYTES - Insert a number of bit fields. Cray routine.
!
!  Reverses the action of gbytes, taking fields from s and
!  inserting them into a bit string in d. see GBYTES.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: D. Robertson
!  Aug. 1981
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
  INTEGER :: nbpw
  PARAMETER (nbpw=64)

  INTEGER :: iw, id
  INTEGER :: iskip1,ibits
  INTEGER :: icon, INDEX
  INTEGER :: mask, msk
  INTEGER :: movel, iwords, istep

  INTEGER :: itemp

  INTEGER :: sh1,sh2,sh3
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  id = 1 + iskip/nbpw
  iskip1 = MOD(iskip,nbpw)
  istep = nskip + nbits

  DO iw=1,nwrd

    sh1 = iskip1 + nbits
    IF (sh1 > nbpw) GO TO 50

    sh2 = nbpw - sh1
    IF (sh2 < 0) sh2 = nbpw-sh2
!
!-----------------------------------------------------------------------
!
!  Byte goes into 1 word of iout.
!
!-----------------------------------------------------------------------
!
    iout(id) = shift(OR(AND(shift(iout(id),sh1),mask(nbpw-nbits)),      &
                     AND(in(iw),shift(mask(nbits),nbits))),sh2)
    GO TO 65
    50   CONTINUE
!
!-----------------------------------------------------------------------
!
!  Byte goes into 2 words of iout.
!
!-----------------------------------------------------------------------
!
    sh3 = 2*nbpw-sh1
    iout(id) = OR( AND(iout(id),mask(iskip1)),                          &
                   AND(shift(in(iw),sh3),compl(mask(iskip1))) )

    iout(id+1) = OR( AND(iout(id+1),shift(compl(mask(sh1-nbpw)),        &
                                          nbpw)),                       &
                     shift(AND(in(iw),compl(mask(sh3))),sh3) )
    65   CONTINUE
!
!-----------------------------------------------------------------------
!
!  Update starting word and bit position
!
!-----------------------------------------------------------------------
!
    iskip1 = iskip1 + istep
    IF (iskip1 < nbpw) CYCLE

    iskip1 = iskip1 - nbpw
    id = id + 1 + iskip1/nbpw
    iskip1 = MOD(iskip1,nbpw)
  END DO

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
!  GBYTE - Extract a number of bit fields. Cray routine.
!
!  IN is bit string of indefinite length. Gbytes will
!  extract a bitstrings, nbits bits long, and store them
!  right justified 0 fill, into successive words of d. The
!  successive bitstrings start at bit positions
!
!      iskip+1+(iw-1)*(nbits+nskip)
!
!  In the bit string s. i.e. skip iskip bits at the start,
!  and nskip bits between the extracted strings.
!  Bit iskp+1 in a string is found in word is=1+iskip/nbpw in IN,
!  where nbpw is the number of bits per word. the starting bit
!  is found by skipping mod(iskp,nbpw) bits in that word.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: D. Robertson
!  Aug. 1981
!  MODIFICATIONS:
!
!  12/05/95 (Yuhe Liu)
!  Converted to ARPS standard format and added documents
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
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: in(*)

  INTEGER :: nbits
  INTEGER :: iskip

  INTEGER :: iout
!
!-----------------------------------------------------------------------
!
!  Misc. local variables:
!
!-----------------------------------------------------------------------
!
  INTEGER :: nbpw
  PARAMETER (nbpw=64)

  INTEGER :: iw, id
  INTEGER :: iskip1, ibits
  INTEGER :: icon, INDEX
  INTEGER :: mask, msk
  INTEGER :: movel, iwords

  INTEGER :: itemp

  INTEGER :: sh1,sh2,sh3
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  id = 1 + iskip/nbpw
  iskip1 = MOD(iskip,nbpw)

  sh1 = iskip1 + nbits
  IF(sh1 > nbpw) GO TO 50
!
!-----------------------------------------------------------------------
!
!  Byte comes fromm 1 word of IN
!
!-----------------------------------------------------------------------
!
  iout = AND( shift(in(id),sh1), shift(mask(nbits),nbits) )

  RETURN

  50    CONTINUE

  sh1 = sh1 - nbpw
!
!-----------------------------------------------------------------------
!
!  Byte comes from 2 words of IN.
!
!-----------------------------------------------------------------------
!
  iout = OR( shift(AND(in(id),compl(mask(iskip1))),sh1),                &
             shift(AND(in(id+1),mask(sh1)),sh1) )

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
!  GBYTES - Extract a number of bit fields. Cray routine.
!
!  IN is bit string of indefinite length. Gbytes will
!  extract nwrd bitstrings, nbits bits long, and store them
!  right justified 0 fill, into successive words of d. The
!  successive bitstrings start at bit positions
!
!      iskip+1+(iw-1)*(nbits+nskip)
!
!  In the bit string s. i.e. skip iskip bits at the start,
!  and nskip bits between the extracted strings.
!  Bit iskp+1 in a string is found in word is=1+iskip/nbpw in IN,
!  where nbpw is the number of bits per word. the starting bit
!  is found by skipping mod(iskp,nbpw) bits in that word.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: D. Robertson
!  Aug. 1981
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
  INTEGER :: nbpw
  PARAMETER (nbpw=64)

  INTEGER :: iw, id
  INTEGER :: iskip1,ibits
  INTEGER :: icon, INDEX
  INTEGER :: mask, msk
  INTEGER :: movel, iwords, istep

  INTEGER :: itemp

  INTEGER :: sh1,sh2,sh3
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  id = 1 + iskip/nbpw
  iskip1 = MOD(iskip,nbpw)
  istep = nskip + nbits

  DO iw=1,nwrd

    sh1 = iskip1 + nbits
    IF(sh1 > nbpw) GO TO 50
!
!-----------------------------------------------------------------------
!
!  Byte comes fromm 1 word of IN
!
!-----------------------------------------------------------------------
!
    iout(iw) = AND( shift(in(id),sh1), shift(mask(nbits),nbits) )
    GO TO 65
    50   CONTINUE
    sh1 = sh1 - nbpw
!
!-----------------------------------------------------------------------
!
!  Byte comes from 2 words of IN.
!
!-----------------------------------------------------------------------
!
    iout(iw) = OR( shift(AND(in(id),compl(mask(iskip1))),sh1),          &
                   shift(AND(in(id+1),mask(sh1)),sh1) )
    65   CONTINUE
!
!-----------------------------------------------------------------------
!
!  Update starting word and bit position
!
!-----------------------------------------------------------------------
!
    iskip1 = iskip1 + istep
    IF(iskip1 < nbpw) CYCLE

    iskip1 = iskip1 - nbpw
    id = id + 1 + iskip1/nbpw
    iskip1 = MOD(iskip1,nbpw)

  END DO

  RETURN
END SUBROUTINE grbgbytes

REAL FUNCTION f_cputime()

!------------------------------------------------------------------------------
! CRAY DEFINITION FOR TIMING
!------------------------------------------------------------------------------

  f_cputime = 0.0

  RETURN
END FUNCTION f_cputime


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
!   which only supports 0 <= ICHAR(a) <= 127 on the
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

  INTEGER ::  mold                                                      
  CHARACTER(LEN=1) :: a  

  char2i = ICHAR(a)

  RETURN
END FUNCTION char2i

SUBROUTINE inquiredir(dirname,fexist)

  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN)  :: dirname
  LOGICAL,          INTENT(OUT) :: fexist

  INQUIRE(FILE=TRIM(dirname),EXIST=fexist)

  RETURN
END SUBROUTINE inquiredir
