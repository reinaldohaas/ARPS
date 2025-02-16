!
!
!##################################################################
!##################################################################
!######                                                      ######
!######                  F77toF90 convertor                  ######
!######                                                      ######
!##################################################################
!##################################################################
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Convert F77 source codes into F90 source codes
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Alan Miller
!
!  MODIFICATION HISTORY:
!
!    05/01/1999(Yanming Li, CAPS/OU)
!
!    Fixed runtime segmentation fault problem on origin(IRIS)/Cray
!    Modified to output f90 code conforming to ARPS f90 standards
!    Fixed problem with comment/continuation.
!    Added code to remove include statements and incert USE statements
!    before IMPLICIT NONE. This part of code is disable in this version.
!
!-----------------------------------------------------------------------
!
!  Usage:
!    compile: f90 convert_src_to_f90.f90
!    to convert a single f77 source code:
!       1)type a.out and carriage return
!       2)type filename after the prompt, "Enter name of Fortran source file:"
!       3)type one more carriage return to abort the program.
!     to convert a group of f77 source codes:
!       1)create a file, called enlist,  to store the filenames of f77 source
!         code to be converted
!       2)insert a empty line at the end of enlist.
!       3)a.out < enlist
!
!-----------------------------------------------------------------------
!
! Original comments:
!
! Takes Fortran 77 code in standard format and makes some changes to produce
! free-format Fortran 90 code.
! Changes included are:
!     C or c in column 1 replaced with !
!     Continuation denoted by a character in column 6 replaced with & at the
!         end of the previous line.
!     Indenting of code for DO-loops and IF blocks.
!     END of program unit replaced by END SUBROUTINE (/PROGRAM/FUNCTION) name
!     Fortran `keywords' are in upper case, all other words other than those
!         in character strings are converted to lower case.
!     .LT., .EQ., etc. replaced with <, ==, etc.
!     Labels removed from DO loops; all of which will end with END DO.
!         If labels are not referenced, they are removed.
!     Short continued lines are adjoined to the previous line.
!     ENDIF, ELSEIF & GOTO split into separate words.
!     3-way arithmetic IF constructs are converted to IF .. ELSE IF form.
!     Embedded blanks are removed from numbers in DATA statements.
!     INTENT declarations are added for dummy arguments.
!     Some GO TOs are converted to CYCLE or EXIT.
!     Converts CHARACTER * to CHARACTER (LEN=xx) ::.

! To be done:
!     DATA statements to be replaced by assignments on the declaration line.
!     IMPLICIT NONE statements to be included.
!     Declaration of types of unlisted variables.
!     Functions to be converted to ELF90 form, i.e. REAL FUNCTION XYZ(arg)
!         converted to FUNCTION xyz(arg) RESULT(fn_val).

! Known problems
!     Cannot handle character strings or names broken at the end of lines.
!     No attempt to convert BLOCKDATA, COMMON or EQUIVALENCE.
!     Does not convert Hollerith strings, e.g. 31HTHIS IS A COMMENT ...
!     Does not convert computed GOTOs.
!     May do the wrong thing if variable names start with IF or end with DO.
!     INTENTs are sometimes wrong.  In particular, INTENT(IN) arguments are
!         often shown as INTENT(IN OUT).
!     Computed GOTOs are not converted, but if they occur inside DO loops,
!         and one or more of the labels is to the end of a DO loop, then
!         a label may be removed.

! The default extension for the name of the input file is `for'; this can be
! over-ruled by giving the full name (e.g. myprog.f77).   The output file name
! will be the input name (and directory) with extension `.f90'.

! Added conversion of `enddo' to END DO - 13 March 1997
! Corrected bug which occurred when an arithmetic IF within a DO-loop involved
!     a jump to the end of the DO-loop - 17 August 1997.
! ELSEIF, ENDIF & ELSEIF were being split into 2 separate words, and then the
!     last letter converted back to lower case - corrected 17 August 1997.
! Corrected bug which occurred when .LT. (or other comparison) had a blank
!     before and/or after, followed on the same line by a text string, followed
!     by a Fortran word such as THEN or GO TO - 8 December 1997.
! Added (LEN=1) after CHARACTER if length not specified - 9 December 1997.
! Embedded blanks are removed from numerical constants in DATA statements.
!     Added 9 December 1997.
! Added INTENTs and TYPE declarations for dummy arguments - 23 December 1997.
! Corrected problem when DO statement contains a comma immediately after DO,
!     and improved the detection of INTENTs when a dummy argument appears in an
!     IF-expression.  Added extra indentation on continuation lines to improve
!     readability - 13 January 1998
! Corrected a bug which could occur when the last type declaration was matched
!     to a dummy variable and the line deleted - 5 June 1998
! Corrected jumps out of inner nested DO loops, and replaced GO TOs, out of
!     DO loops to the next executable line, with EXIT - 8 June 1998
! Added conversion of CHARACTER * to CHARACTER (LEN=xx) ::
!     including CHARACTER*10 a, d, c*50, d   - 21 June 1998.
! Corrected for case of final command of a DO loop which is not CONTINUE and
!     which flows onto the next line - 29 June 1998.

! Latest revision - 29 June 1998
! Author - Alan Miller @ vic.cmis.csiro.au
! WWW-page:  www.ozemail.com.au/~milleraj


MODULE implicit
! Module to set and reset implicit variable types for use by to_f90.

IMPLICIT NONE
INTEGER, SAVE :: var_type(26) = (/  &
                 1,1,1,1,1,1,1,1,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1 /)
!                a b c d e f g h i j k l m n o p q r s t u v w x y z
CHARACTER (LEN=24), SAVE :: vt(0:7) = (/ 'NO TYPE                 ', &
                                         'REAL                    ', &
                                         'INTEGER                 ', &
                                         'DOUBLE PRECISION        ', &
                                         'LOGICAL                 ', &
                                         'COMPLEX                 ', &
                                         'CHARACTER               ', &
                                         'OTHER TYPE              ' /)

CONTAINS


SUBROUTINE reset_defaults()

var_type(1:8) = 1            ! REAL (A-H)
var_type(9:14) = 2           ! INTEGER (I-N)
var_type(15:26) = 1          ! REAL (O-Z)

RETURN
END SUBROUTINE reset_defaults



SUBROUTINE set_implicit_types(text)
! Read in implicit statement and interpret.

CHARACTER (LEN=*), INTENT(IN OUT) :: text

! Local variables
INTEGER :: ivt, length, start, i, j, pos, left, right
LOGICAL :: first

i = INDEX(text, 'IMPLICIT')
IF (i > 0) text = text(i+8:)
text = ADJUSTL(text)

DO
  IF (text(1:4) == 'NONE') THEN
    var_type = 0
    RETURN
  ELSE IF (text(1:4) == 'REAL') THEN
    ivt = 1
  ELSE IF (text(1:7) == 'INTEGER') THEN
    ivt = 2
  ELSE IF (text(1:24) == 'DOUBLE PRECISION COMPLEX') THEN
    ivt = 7
    vt(7) = 'DOUBLE PRECISION COMPLEX'
  ELSE IF (text(1:16) == 'DOUBLE PRECISION') THEN
    ivt = 3
  ELSE IF (text(1:7) == 'LOGICAL') THEN
    ivt = 4
  ELSE IF (text(1:7) == 'COMPLEX') THEN
    ivt = 5
  ELSE IF (text(1:9) == 'CHARACTER') THEN
    ivt = 6
  ELSE
    ivt = 7
    i = INDEX(text, ' ')
    vt(7) = text(1:i-1)
  END IF

! Interpret the part in brackets, e.g. (a - h, o - z)

  length = LEN_TRIM(text)
  start = 5
  left = INDEX(text(start:length), '(') + start - 1
  IF (left < start) RETURN
  right = INDEX(text(start:length), ')') + start - 1
  IF (right < left) RETURN
                                       ! Interpret text(left+1:right-1)
  first = .TRUE.
  DO pos = left+1, right
    SELECT CASE (text(pos:pos))
      CASE (' ')
        CYCLE
      CASE ('-')
        first = .FALSE.
      CASE (',', ')')
        IF (first) THEN
          var_type(i) = ivt
        ELSE
          var_type(i:j) = ivt
          first = .TRUE.
        END IF
      CASE DEFAULT
        IF (first) THEN
          i = ICHAR(text(pos:pos)) - ICHAR('a') + 1
          IF (i < 1) THEN
            i = ICHAR(text(pos:pos)) - ICHAR('A') + 1
          END IF
        ELSE
          j = ICHAR(text(pos:pos)) - ICHAR('a') + 1
          IF (j < 1) THEN
            j = ICHAR(text(pos:pos)) - ICHAR('A') + 1
          END IF
        END IF
    END SELECT
  END DO

  start = right + 1
  IF (start >= length) RETURN
  text = text(start:length)
  DO
    IF (text(1:1) == ',' .OR. text(1:1) == ' ') THEN
      text = text(2:)
    ELSE
      EXIT
    END IF
  END DO
END DO

RETURN
END SUBROUTINE set_implicit_types



FUNCTION implicit_type(ch) RESULT(vtype)
! Return the variable type given the first character of its name.
! The first character is expected to be lower case, but just in case ..

CHARACTER (LEN=1), INTENT(IN) :: ch
CHARACTER (LEN=24)            :: vtype

! Local variable
INTEGER  :: i, j

i = ICHAR(ch) - ICHAR('a') + 1
IF (i >= 1 .AND. i <= 26) THEN
  j = var_type(i)
  vtype = vt(j)
ELSE
  i = ICHAR(ch) - ICHAR('A') + 1
  IF (i >= 1 .AND. i <= 26) THEN
    j = var_type(i)
    vtype = vt(j)
  ELSE
    vtype = ' '
  END IF
END IF

RETURN
END FUNCTION implicit_type

END MODULE implicit



PROGRAM to_f90
USE implicit
IMPLICIT NONE

TYPE :: code
  CHARACTER (LEN=140)  :: text
  CHARACTER (LEN=  5)  :: label
  TYPE (code), POINTER :: next
END TYPE code

TYPE :: argument
  CHARACTER (LEN=10)       :: name
  INTEGER                  :: intention    ! IN = 1, OUT = 2, IN OUT = 3
  CHARACTER (LEN=24)       :: var_type     ! Room for DOUBLE PRECISION COMPLEX
  INTEGER                  :: dim          ! DIM = 0 for scalars
  CHARACTER (LEN=16)       :: dimensions   ! Not used if DIM = 0
  TYPE (argument), POINTER :: next
END TYPE argument

CHARACTER (LEN=30)       :: f77_name, f90_name
CHARACTER (LEN= 1)       :: tab = CHAR(9)
CHARACTER (LEN=50)       :: prog_unit_name = ' ', blank = ' '
CHARACTER (LEN= 9)       :: delimiters = ' =+-*/,()'
CHARACTER (LEN=10)       :: numbers = '1234567890'
CHARACTER (LEN= 5)       :: lab
CHARACTER (LEN=140)       :: text, vtype
CHARACTER (LEN=140)      :: statement, statement1, statement2
INTEGER                  :: ii,jj
INTEGER (KIND(123))      :: iostatus, pos,posa, count, last, n_marks, pos1(20),  &
                            pos2(20), lab_length, indent, i, i1, i2, length, &
                            left_brac, right_brac, numb_arg
TYPE (code), POINTER     :: head, current, current1,tail, last_line, next_line,  &
                            first_decl, last_decl, start_prog_unit, before_implicit, &
                            end_prog_unit, begin_include, end_include
LOGICAL                  :: asterisk, OK, data_stmnt, first_arg, continuation
TYPE (argument), POINTER :: arg_start, arg, last_arg

integer istatus

INTERFACE
  SUBROUTINE mark_text(text, n_marks, pos1, pos2)
    IMPLICIT NONE
    CHARACTER (LEN = *), INTENT(IN)  :: text
    INTEGER (KIND(123)), INTENT(OUT) :: n_marks, pos1(:), pos2(:)
  END SUBROUTINE mark_text

  SUBROUTINE convert_text(text, n_marks, pos1, pos2)
    IMPLICIT NONE
    CHARACTER (LEN = *), INTENT(IN OUT) :: text
    INTEGER (KIND(123)), INTENT(IN)     :: n_marks
    INTEGER (KIND(123)), INTENT(IN OUT) :: pos1(:), pos2(:)
  END SUBROUTINE convert_text

  SUBROUTINE remove_data_blanks(text)
    IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(IN OUT) :: text
  END SUBROUTINE remove_data_blanks

  FUNCTION last_char( text ) RESULT(ch)
    IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(IN) :: text
    CHARACTER (LEN=1)             :: ch
  END FUNCTION last_char

  FUNCTION find_delimited_name (text, name) RESULT(pos)
    IMPLICIT NONE
    CHARACTER (LEN=*), INTENT(IN) :: text, name
    INTEGER                       :: pos
  END FUNCTION find_delimited_name
END INTERFACE

DO
  WRITE(*, '(a)', ADVANCE='NO')' Enter name of Fortran source file: '
  READ(5, '(a)',end=111, err=111) f77_name

goto 112
111 continue
  print*,'end of file list. '
  Stop
112 continue

  IF (INDEX(f77_name, '.') == 0) THEN
    last = LEN_TRIM(f77_name)
    f77_name(last+1:last+4) = '.for'
  END IF
  OPEN(8, file=f77_name, status='old',IOSTAT=istatus)

  IF( istatus.ne.0 ) then
     print*,'opening of f77_name failed.'
     stop
  endif

  pos = INDEX(f77_name, '.')
  f90_name = f77_name(1:pos) // 'f90'
  OPEN(9, file=f90_name)

!     Set up a linked list containing the lines of code

  NULLIFY(head, tail)
  ALLOCATE(head)
  tail => head
  READ(8, '(a)') head % text
  IF (head % text(1:1) == 'C' .OR. head % text(1:1) == 'c' .OR.   &
      head % text(1:1) == '*') THEN
    head % text(1:1) = '!'
  ELSE IF (head % text(1:1) == tab) THEN
    head % text = '      ' // head % text(2:)
  END IF
  head % label = ' '
  count = 1

  DO
    NULLIFY(current)
    ALLOCATE(current)
    READ(8, '(a)', IOSTAT=iostatus) current % text
    IF (iostatus /= 0) EXIT

!     Change C, c or * in column 1 to !
    IF (current % text(1:1) == 'C' .OR. current % text(1:1) == 'c' .OR.   &
        current % text(1:1) == '*'  ) THEN
        pos = LEN_TRIM(current % text)
        do ii = 2, pos
          text(ii:ii)='-'     ! make a line of '-'s
        enddo
      IF (LEN_TRIM(current % text) >= 1) THEN  ! modified and fixed the 
        current % text(1:1) = '!'              ! missing '!' problem
        if (current % text(2:2)=='#'.and. current % text(20:20) =='#'       &
             .and. current % text(40:40)=='#') then
         ! do ii = 2, pos
         ! text(ii:ii)='-'
         ! enddo
        current % text(2:pos)=text(2:pos)
        endif
      ELSE
        current % text = ' '           ! Leave blank if nothing else on line
      END IF

      posa = INDEX(current % text, '#')
      if( posa > 0 .AND. current%text(2:posa-1) == ' ' ) THEN
         current%text = current%text(1:1) // current%text(posa:)
      endif
      if(current%text(2:6) == ' ' ) THEN
         current%text = current%text(1:1) // current%text(5:)
      endif
 
      current % label = ' '
    ELSE
      current % label = current % text(1:5)
    END IF

    count = count + 1
    IF (current % label(1:1) == tab) THEN          ! Expand tabs
      current % label = ' '
      current % text = '      ' // current % text(2:)
    ELSE IF (current % label(1:1) == '!') THEN
      current % label = ' '   
    ELSE
      current % label= ADJUSTL(current % label)
    END IF

    pos = LEN_TRIM(current % text)  ! remove any tab in the code
    if( pos >= 3 ) THEN
      do ii = 2, pos 
        if ( current % text (ii:ii) == tab ) THEN
          current % text = current % text (1 : ii-1) // '     '//current % text(ii+1:)
        endif
      enddo
    endif
      
    NULLIFY(current % next)
    tail % next => current
    tail => current
  END DO

  WRITE(*, *)'No. of lines read =', count


! Replace continuation line symbol with '&', and put it at end of last line

  current => head
  NULLIFY(last_line)

  DO
!     Look for blanks in columns 1-5 followed by non-blank in column 6.
!     If found, add an ampersand at the end of the previous line.

!     Replace tabs with single spaces
    DO
      pos = INDEX(current % text, tab)
      IF (pos == 0) EXIT
      current % text(pos:pos) = ' '
    END DO

    IF (current % label == '     ' .AND. current % text(6:6) /= ' ' .AND.  &
        current % text(1:1) /= '!') THEN
        pos = INDEX(last_line%text(2:), '!')   !be careful the value of pos, it is 
                                               !different when we use --(1:)
        if( pos > 0)then

           last_line%text =   &
                 last_line%text(1:pos)//'& '//last_line%text(pos+1:) 
        else
           last = len_trim( last_line%text)
           last_line % text(last+2:last+2) = '&'
        endif

        current % text(6:6) = ' '
    endif


    last_line => current
    IF (ASSOCIATED(current, tail)) EXIT
    IF (.NOT. ASSOCIATED(current)) EXIT
    current => current % next
  END DO

 
  current => head
  NULLIFY(last_line)

  DO
!     Look for blanks in columns 1-5 followed by non-blank in column 6.
!     If found, add an ampersand at the end of the previous line.


    pos = index(current%text(2:), '!')
    posa = index(current%text, '&')


!    if( index( last_line%text, '&') > 0 ) then
!        current%text =  current%text(indent:)
!    endif 

if ( index( current%text, '!' ) /= 1 ) then

    if(  posa > 0 .AND. &
        index(last_line%text, '&') == 0 ) then

        statement = adjustl(current%text)
        
        indent = len_trim(current%text) -  &
                len_trim( statement )
        current%text = ADJUSTL( current%text )
    endif

    if ( pos /= 0 .and. current%text(1:pos-1)==' ' ) then
      jj = INDEX(last_line%text(2:), '!')
      if( jj > 0 ) THEN
      do ii=1,jj
       text(ii:ii) = ' '
      enddo
      current%text = text( :jj-1)//ADJUSTL(current % text)
      else
        current%text(1:1) = '!'
      endif

      else if( index( last_line%text, '&') > 0 ) then

       if( current%text(1:indent+5) == ' ' ) then
        current%text =  current%text( indent+5: )
       else 
        current%text = ADJUSTL( current%text )
       endif
        

      else
         current % text = ADJUSTL(current % text)
      endif
endif
    last_line => current
    IF (ASSOCIATED(current, tail)) EXIT
    IF (.NOT. ASSOCIATED(current)) EXIT
    current => current % next
  END DO


!---------------------------------------------------------------------------

  current => head
  NULLIFY(last_line)
  data_stmnt = .FALSE.

  DO
!     Look for blanks in columns 1-5 followed by non-blank in column 6.
!     If found, add an ampersand at the end of the previous line.

!     Mark regions of text which must not have their case changed.
    CALL mark_text(current % text, n_marks, pos1, pos2)

!     Convert cases of regions which are not protected.
    CALL convert_text(current % text, n_marks, pos1, pos2)

!     If line is start of a program unit, record its name
    IF (current % text(1:7) == 'PROGRAM') THEN
      prog_unit_name = current % text(1:50)
    ELSE IF (current % text(1:10) == 'SUBROUTINE') THEN
      pos = INDEX(current % text, '(') - 1
      IF (pos < 0) pos = LEN_TRIM(current % text)
      prog_unit_name = current % text(1:pos)
    ELSE IF (current % text(1:9) == 'BLOCKDATA') THEN
      prog_unit_name = current % text(1:50)
    ELSE
                               ! N.B. 'FUNCTION' could be part of a comment
      pos = INDEX(current % text, 'FUNCTION')
      IF (pos > 0 .AND. INDEX(current % text, '!') == 0 .AND.        &
                           INDEX(current % text, "'") == 0) THEN
        last = INDEX(current % text, '(') - 1
        IF (last < 0) last = LEN_TRIM(current % text)
        prog_unit_name = current % text(pos:last)
      END IF
    END IF

!     If first word is one of INTEGER, REAL, DOUBLE PRECISION, CHARACTER ,
!     LOGICAL or COMPLEX, add :: unless FUNCTION appears on the same line
!     or next non-blank character is '*' as in REAL*8.
!
!mx  - and unless :: appears in the line
!
    IF (INDEX(current % text, 'FUNCTION') == 0 .and.  &
        INDEX(current % text, '::') == 0) THEN
      pos = 0
      IF (INDEX(current % text, 'INTEGER') == 1) THEN
        pos = 9
      ELSE IF (INDEX(current % text, 'REAL') == 1) THEN
        pos = 6
      ELSE IF (INDEX(current % text, 'DOUBLE PRECISION') == 1) THEN
        pos = 18
      ELSE IF (INDEX(current % text, 'CHARACTER') == 1) THEN
        pos = 11
      ELSE IF (INDEX(current % text, 'COMPLEX') == 1) THEN
        pos = 9
      ELSE IF (INDEX(current % text, 'LOGICAL') == 1) THEN
        pos = 9
      END IF

      IF (pos > 0) THEN
        asterisk = INDEX(current % text(pos-1:pos), '*') > 0
        IF (.NOT. asterisk) THEN
          IF (pos /= 11) THEN
            current % text = current % text(1:pos-1) // ':: ' //  &
                             ADJUSTL( current % text(pos:) )
          ELSE                         ! CHARACTER type, default length = 1

!mx added check for 'LEN='

!           if(INDEX(current % text, 'LEN=') == 0 ) then
              current % text = 'CHARACTER (LEN=1) :: ' //   &
                               ADJUSTL( current % text(pos:) )
!           endif
          END IF
        ELSE
          IF (pos == 11) THEN          ! CHARACTER * found
            i1 = INDEX(current % text, '*') + 1
            length = LEN_TRIM(current % text)
                                       ! Get length, could be (*)
            DO
              IF (current % text(i1:i1) /= ' ') EXIT
              IF (i1 >= length) EXIT
              i1 = i1 + 1
            END DO
            IF (current % text(i1:i1) == '(') THEN
              i1 = i1 + 1
              i2 = INDEX(current % text, ')') - 1
            ELSE
              i2 = INDEX(current % text(i1:), ' ') + i1 - 2
            END IF
            current % text = 'CHARACTER (LEN=' // current % text(i1:i2) //  &
                             ') :: ' // ADJUSTL( current % text(i2+2:) )
          END IF
        END IF

                   ! Check for 2 or more lengths in CHARACTER declaration.
                   ! e.g. CHARACTER a, b, c*10, d
                   ! Put 2nd (& later) declarations on separate lines:
                   ! CHARACTER*10 c
        IF (pos == 11) THEN
          pos = INDEX(current % text, '::') + 2
          DO
            i = INDEX(current % text(pos:), '*')
            IF (i == 0) EXIT
            i = i + pos - 1
            length = LEN_TRIM(current % text)
            i1 = INDEX(current % text(:i-1), ',', BACK=.TRUE.)
            i1 = MAX(pos, i1)
            i2 = INDEX(current % text(i+1:), ',')
            IF (i2 == 0) THEN
              i2 = length + 1
            ELSE
              i2 = i2 + i
            END IF
                   ! i1, i2 mark commas at beginning & end of `, name*xx,'
            IF (i1 == pos .AND. i2 == length + 1) THEN
                   ! Only one declaration left on line, e.g.
                   ! CHARACTER :: name*50
              current % text = 'CHARACTER (LEN=' // current % text(i+1:length) &
                               // ') :: ' // ADJUSTL( current % text(i1:i-1) )
              EXIT
            END IF

            ALLOCATE( next_line )
            next_line % next => current % next
            current % next => next_line
            next_line % text = 'CHARACTER' // current % text(i:i2-1) //  &
                               ' ' // current % text(i1+1:i-1)
            IF (i2 < length) THEN
              current % text = current % text(:i1) // current % text(i2+1:length)
            ELSE
              current % text = current % text(:i1-1)
            END IF
          END DO
        END IF
      END IF
    END IF

!     If this is in a DATA statement, eliminate any blanks within numbers
    IF ( current % text(1:4) == 'DATA') THEN
       print*, current%text

!      CALL remove_data_blanks(current % text)
!      last = LEN_TRIM(current % text)
!      data_stmnt = .TRUE.
    END IF

!     If line only contains 'END', add the program unit name
    IF (LEN_TRIM(current % text) == 3 .AND. current % text(1:3) == 'END') THEN
      current % text = current % text(1:3) // ' ' // prog_unit_name
      prog_unit_name = ' '

!     Convert `enddo' to 'END DO'
    ELSE IF (current % text(1:5) == 'enddo') THEN
      current % text = 'END DO' // current % text(6:)
    END IF

    last_line => current
    IF (ASSOCIATED(current, tail)) EXIT
    IF (.NOT. ASSOCIATED(current)) EXIT
    current => current % next
  END DO

!-------------------------------------------------------------------------

!     Now convert Do-loops

  current => head
  WRITE(*, *) '      Converting DO-loops & 3-way IFs'
  DO
    IF (current % text(1:1) /= '!' .AND. current % text(1:1) /= ' ') THEN
      pos = INDEX(current % text, 'DO')
      IF ( pos > 0 .AND. (current % text(pos+2:pos+2) == ' ' .OR.  &
                          current % text(pos+2:pos+2) == ',' ) ) THEN
        IF ( current % text(pos+2:pos+2) == ',' )  &
                          current % text(pos+2:pos+2) = ' '
        IF (pos == 1) THEN
          OK = .TRUE.
        ELSE IF (SCAN(current % text(pos-1:pos-1), delimiters) > 0) THEN
          OK = INDEX(current % text(:pos-1), 'END ') == 0
        ELSE
          OK = .FALSE.
        END IF
        IF (OK) THEN
          text = ADJUSTL( current % text(pos+3:) )
          last = INDEX( text, ' ')
          lab = text(:last-1)
          IF (SCAN(lab(1:1), numbers) == 0) lab = ' '
          lab_length = LEN_TRIM(lab)
          IF (lab_length > 0) THEN
            pos = INDEX(lab, ',')      ! Check for a comma after label
            IF (pos > 0) THEN
               lab(pos:) = ' '
               i = INDEX(current % text, ',')
               current % text(i:i) = ' '
               lab_length = pos - 1
            END IF
            CALL do_loop_fixup(current)
          END IF
        END IF
! Look for IF, then a number as last non-blank character
      ELSE
        pos = INDEX(current % text, 'IF')
        IF (pos > 0) THEN
          last = LEN_TRIM(current % text)
          IF (SCAN(current % text(last:last), numbers) > 0) THEN
            CALL fix_3way_IF(current)
          END IF
        END IF
      END IF
    END IF

    IF (ASSOCIATED(current, tail)) EXIT
    IF (.NOT. ASSOCIATED(current)) EXIT
    current => current % next
  END DO

!-------------------------------------------------------------------------

!     Determine INTENTs for dummy arguments

  WRITE(*, *) '      Determining INTENTs of dummy arguments'

!     Search for either FUNCTION or SUBROUTINE.
!     Extract name of program unit.

  current => head
  NULLIFY(last_line)
  outer_loop: DO
    DO
      IF (current % text(1:1) /= '!' .AND. current % text(1:1) /= ' ') THEN
        IF (current % text(1:7)  == 'PROGRAM' ) THEN
           prog_unit_name = current % text
           EXIT

        else IF (current % text(1:10) == 'SUBROUTINE') THEN
          pos = INDEX(current % text, '(') - 1
          IF (pos < 0) pos = LEN_TRIM(current % text)
          prog_unit_name = current % text(1:pos)
          EXIT
        ELSE
          pos = INDEX(current % text, 'FUNCTION')
          IF (pos > 0) THEN
            last = INDEX(current % text, '(') - 1
            IF (last < 0) last = LEN_TRIM(current % text)
            prog_unit_name = current % text(pos:last)
            EXIT
          END IF
        END IF
      END IF

      last_line => current
      current => current % next
      IF (ASSOCIATED(current, tail)) EXIT outer_loop
    END DO

!     If there is no blank line between this program unit and the previous
!     one, then insert one.

    IF (ASSOCIATED(last_line)) THEN
      IF (LEN_TRIM(last_line % text) > 0) THEN
        CALL insert_and_moveto_newline(last_line)
        last_line % text = ' '
      END IF
    END IF

    ALLOCATE( start_prog_unit )
    start_prog_unit => current

!     Find end of program unit

    DO
      current => current % next
      IF (current % text(1:1) /= '!' .AND. current % text(1:1) /= ' ') THEN
        IF (current % text(1:3) == 'END') THEN
          IF (INDEX(current % text(5:), prog_unit_name) > 0) THEN
            ALLOCATE( end_prog_unit )
            end_prog_unit => current
            EXIT
          END IF
        END IF
      END IF
      IF (ASSOCIATED(current, tail)) EXIT outer_loop
    END DO

    IF( .false. ) then ! commented out by ming

!   Find 'include *.inc' and move them  right before implicit none

    allocate(current1)
    allocate(before_implicit)
    current1 => start_prog_unit 
    do

     if( current1%text(1:8) == 'IMPLICIT') THEN
       before_implicit => last_line
       exit
     endif
     if( current1%text(1:7) == 'INTEGER') THEN
       before_implicit => last_line
       exit
     endif
     if( current1%text(1:4) == 'REAL') THEN
       before_implicit => last_line
       exit
     endif
     if( current1%text(1:6) == 'COMMON') THEN
       before_implicit => last_line
       exit
     endif
     if( current1%text(1:7) == 'INCLUDE') THEN
       before_implicit => last_line
       exit
     endif
     if( current1%text(1:9) == 'PARAMETER') THEN
       before_implicit => last_line
       exit
     endif
     if( current1%text(1:9) == 'CHARACTER') THEN
       before_implicit => last_line
       exit
     endif

      last_line => current1
     if(.NOT. ASSOCIATED(current1))exit
     if(ASSOCIATED(current1, end_prog_unit) ) exit
     current1 => current1%next
    enddo


    allocate(current1)
    current1 => start_prog_unit
    do
     if(current1%text(1:1)/= '!'.AND.INDEX(current1%text,'INCLUDE')>0)then 
       pos = INDEX(current1%text(2:), "'")
       posa = INDEX(current1%text, '.')

       CALL insert_and_moveto_newline(before_implicit) 
       before_implicit%text='USE '//current1%text(pos+2:posa-1)

       CALL insert_and_moveto_newline(before_implicit)
       before_implicit%text='  '

      endif
      last_line => current1
     if(.NOT. ASSOCIATED(current1))exit
     if(ASSOCIATED(current1, end_prog_unit) ) exit
     current1 => current1%next
    enddo

    current1 => start_prog_unit
    do
     next_line=>current1
     if(current1%text(1:1)/= '!'.AND.INDEX(current1%text,'INCLUDE')>0)then
      current1%next => next_line%next
      last_line%next => next_line%next
      deallocate(current1)
      current1 => last_line
      endif
      last_line => current1
     if(.NOT. ASSOCIATED(current1))exit
     if(ASSOCIATED(current1, end_prog_unit) ) exit
     current1 => current1%next
    enddo

    endif ! if (.false. ) block

!     Find first & last declarations

    ALLOCATE( first_decl, last_decl )
    CALL find_declarations( start_prog_unit, end_prog_unit, first_decl, &
                            last_decl )
    IF (.NOT. ASSOCIATED(last_decl)) GO TO 100

!     Extract list of dummy arguments

!    CALL get_arg_list()
      numb_arg = 0
    IF (numb_arg == 0) GO TO 100

!     See if the declarations contain any IMPLICIT statements

    CALL reset_defaults()
    current => first_decl
    DO
      IF( current % text(1:8) == 'IMPLICIT' ) THEN
        statement = current % text(10:)
        CALL set_implicit_types(statement)
      END IF
      IF (ASSOCIATED(current, last_decl)) EXIT
      current => current % next
    END DO

!     Search through the declarations for variable types & dimensions

    CALL get_var_types()

!     Search through rest of code to try to determine the INTENTs

    CALL get_intents()

!     Insert INTENT statements

    arg => arg_start
    DO

      SELECT CASE (arg % intention)
        CASE (0, 3)
          current => first_decl
          do
            pos = INDEX(current%text, arg%name)
            if ( pos > 0) then
              statement1= current%text(:pos-1)
              posa= INDEX(statement1, '::')
              if(posa > 0) then
                statement1 = statement1(:posa-1)
                else
                 posa=pos
              endif
              statement2= current%text(pos:)
              current % text = statement1(:posa-1) // ', INTENT(INOUT)' &
                               //':: '//arg%name// statement2(pos:)
              print*, current%text
              exit
             endif
             IF (ASSOCIATED(current, last_decl)) EXIT
             IF (.NOT. ASSOCIATED(current)) EXIT
             current => current % next
            END DO
        CASE (1)
          current => first_decl
          do
            pos = INDEX(current%text, arg%name)
            if ( pos > 0) then
              statement1= current%text(:pos-1)
              posa= INDEX(statement1, '::')
              if(posa > 0) then
                statement1 = statement1(:posa-1)
                else
                 posa=pos
              endif
              statement2= current%text(pos:)
              current % text = statement1(:posa-1) // ', INTENT(IN   )' &
                               //':: '//arg%name// statement2(pos:)
              print*, current%text
              exit
             endif
             IF (ASSOCIATED(current, last_decl)) EXIT
             IF (.NOT. ASSOCIATED(current)) EXIT
             current => current % next
           END DO
        CASE (2)
          current => first_decl
          do
            pos = INDEX(current%text, arg%name)
            if ( pos > 0) then
              statement1= current%text(:pos-1)
              posa= INDEX(statement1, '::')
              if(posa > 0) then
                statement1 = statement1(:posa-1)
                else
                 posa=pos
              endif
              statement2= current%text(pos:)
              current % text = statement1(:posa-1) // ', INTENT(OUT  )' &
                               //':: '//arg%name// statement2(pos:)
              print*, current%text
              exit
             endif
             IF (ASSOCIATED(current, last_decl)) EXIT
             IF (.NOT. ASSOCIATED(current)) EXIT
             current => current % next
            END DO
      END SELECT

      IF (ASSOCIATED(arg, last_arg)) EXIT
      arg => arg % next
    END DO


!     Insert a blank line after the last declaration if there is not one
!     there already, or a comment.

    next_line => last_decl % next
    IF (next_line % text(1:1) /= ' ' .AND. next_line % text(1:1) /= '!') THEN
      CALL insert_and_moveto_newline(last_decl)
      last_decl % text = ' '
    END IF

!     Move onto the next SUBROUTINE or FUNCTION

    100 current => end_prog_unit
    IF (ASSOCIATED(current, tail)) EXIT
    last_line => current
    current => current % next
    IF (ASSOCIATED(current, tail)) EXIT
  END DO outer_loop

!-------------------------------------------------------------------------

!     Indenting and writing output file

  current => head
  indent = 0
  continuation = .FALSE.
  WRITE(*, *) '      Writing file: ', f90_name
  DO
    
    IF (current % text(1:1) /= '!') THEN 
      IF( INDEX(current % text, 'END ') > 0 )THEN
        indent = MAX(indent-2, 0)
      
      !  WRITE(9, '(a)') blank(:indent) // TRIM(current % text)
         current%text = blank(:indent) // TRIM(current % text)

        continuation = (last_char(current % text) == '&')
      ELSE IF (INDEX(current % text, 'DO ') > 0) THEN

      ! WRITE(9, '(a)') blank(:indent) // TRIM(current % text)
        current%text = blank(:indent) // TRIM(current % text)

        continuation = (last_char(current % text) == '&')
        indent = indent + 2
                                                   ! Temporary reduction in
                                                   ! indentation for `ELSE'
      ELSE IF (INDEX(current % text, 'ELSE') > 0) THEN
        last = MAX(0, indent-2)

      !  WRITE(9, '(a)') blank(:last) // TRIM(current % text)
        current%text = blank(:last) // TRIM(current % text)

        continuation = (last_char(current % text) == '&')
                                                   ! Indent increased if `IF'
                                                   ! is followed by `THEN'
      ELSE IF (INDEX(current % text, 'IF ') > 0 .OR.           &
               INDEX(current % text, 'IF(') > 0) THEN
        current % text =  blank(:indent) // TRIM(current % text)
                                             ! If IF statement runs onto
                                             ! next line, try joining
        last = LEN_TRIM(current % text)
        IF (current % text(last:last) == '&') THEN
          next_line => current % next
          IF (last + LEN_TRIM(next_line % text) < 80) THEN
            current % text(last:last) = ' '
            current % text = TRIM(current % text) // ' ' //  &
                           TRIM(next_line % text)
            current % next => next_line % next
          END IF
        END IF

       ! WRITE(9, '(a)') TRIM(current % text)
         current%text = TRIM(current % text)
        continuation = (last_char(current % text) == '&')
        next_line => current
        DO
          IF (INDEX(next_line % text, ' THEN') > 0 .OR.  &
              INDEX(next_line % text, ')THEN') > 0) THEN
            indent = indent + 2
            EXIT
          ELSE
            IF ( last_char(next_line % text) /= '&') EXIT
          END IF
          next_line => next_line % next
        END DO

      ELSE

!     If line ends with '&', attempt to join on the next line if it is short.

!        last = LEN_TRIM(current % text)
!        IF (last > 0) THEN
!          IF (current % text(last:last) == '&') THEN
!            last = LEN_TRIM(current % text(:last-1))
!            next_line => current % next
!            IF (last + indent + LEN_TRIM(next_line % text) < 78) THEN
!              current % text = current % text(:last) // ' ' // &
!                               TRIM(next_line % text)
!              current % next => next_line % next
!              DEALLOCATE(next_line)
!            END IF
!          END IF
!        END IF

        last = LEN_TRIM(current % text)
        IF (last > 0) THEN
          IF (current % text(last:last) == '&') THEN
              current % text = current % text(:last) 
          END IF
        END IF


        IF (continuation) THEN
         ! WRITE(9, '(a)') blank(:indent+4) // TRIM(current % text)
          current%text=blank(:indent+4) // TRIM(current % text)
        ELSE
         ! WRITE(9, '(a)') blank(:indent) // TRIM(current % text)
            current%text = blank(:indent) // TRIM(current % text)
        END IF

        continuation = (last_char(current % text) == '&')
      END IF
!     Comment line (unchanged)
    ELSE
      !WRITE(9, '(a)') TRIM(current % text)
      current%text = TRIM(current % text)

      continuation = .FALSE.
    END IF
    IF (ASSOCIATED(current, tail)) EXIT
    IF (.NOT. ASSOCIATED(current)) EXIT
    current => current % next
  END DO

  current=>head
  do
    pos = INDEX(current%text(2:), '!')
    posa = INDEX(current%text, ')')
    if(INDEX(current%text,'CHARACTER') > 0 .AND. pos>0 .AND. pos < posa ) then
      length = len_trim(current%text)
      current%text = current%text(1:pos-2)//current%text(posa:length) &
                     //' ' //current%text(pos:posa-1)
    endif
    IF (ASSOCIATED(current, tail)) EXIT
    IF (.NOT. ASSOCIATED(current)) EXIT
    current => current % next
  END DO

  current=>head
  do
    pos = INDEX(current%text, 'LEN=(*)')
    if (pos > 0) THEN
      current%text = current%text(1:pos-1)//'LEN=*'// &
                     current%text(pos+7:)
    endif
    IF (ASSOCIATED(current, tail)) EXIT
    IF (.NOT. ASSOCIATED(current)) EXIT
    current => current % next
  END DO

  current=>head
  do
    pos = INDEX(current%text, '!')
    posa = INDEX(current%text, ')')
    if(INDEX(current%text,'CHARACTER') > 0 .AND. pos>0 .AND. pos < posa ) then
      length = len_trim(current%text)
      current%text = current%text(1:pos-2)//current%text(posa:length) &
                     //' ' //current%text(pos:posa-1)
    endif
    IF (ASSOCIATED(current, tail)) EXIT
    IF (.NOT. ASSOCIATED(current)) EXIT
    current => current % next
  END DO


  current=>head
  do
    if (current%text(1:1) == '!' .OR. INDEX(current%text, 'SUBROUTINE') &
        > 0  .OR. INDEX(current%text, 'PROGRAM')>0 ) THEN
     current%text = TRIM(current % text)

    else
     current%text = '  '//TRIM(current % text)

    endif
    IF (ASSOCIATED(current, tail)) EXIT
    IF (.NOT. ASSOCIATED(current)) EXIT
    current => current % next
  END DO

  current=>head
  indent = 0
  do
   posa= INDEX( current%text(2:), '!')
    if ( posa /= 0 .and. current%text(1:posa-1)==' ' ) then  !modify pos by posa
      jj = INDEX(last_line%text, '!')
      do ii=1,jj
       text(ii:ii) = ' '
      enddo
      current%text= text(:jj-1)//adjustl(current%text)
    endif

    next_line => current%next
    if (INDEX(current%text, 'SUBROUTINE')>0 .AND. INDEX(current%text, '&') &
        > 0 ) THEN
      pos = INDEX(current%text, ' ')
      indent = pos 
           current%text = TRIM(current % text)
!    else if(INDEX(current%text, 'CALL') > 0 .AND. INDEX(current%text, '&') &
!                 >0) then
!      pos = INDEX(current%text, '(')
!      indent = pos 
!           current%text = TRIM(current % text)
!    else if(INDEX(current%text, '=')>0 .AND. INDEX(current%text  &
!                 , '==')==0 .AND. INDEX(current%text, '&') &
!                 >0 ) then
!       pos = index(current%text, '=')
!       indent = pos 
!            current%text = TRIM(current % text)
    else if(indent > 0) then
     current%text = blank(:indent)//ADJUSTL(current % text)
    endif
     if( INDEX ( last_line%text, '&')>0 .AND. index(current%text, '&') &
         == 0) then
         indent =0
      endif 
    last_line => current
    IF (ASSOCIATED(current, tail)) EXIT
    IF (.NOT. ASSOCIATED(current)) EXIT
    current => current % next
  END DO

  current=>head
  do

    posa = INDEX( current%text(2:), '!')
    if( current%text(1:1) == '!' .AND. posa > 0 .AND. current%text(2:posa-1) == ' ' ) THEN
       current%text(1:1) = ' '
    endif

    pos = INDEX(current%text, '&')
    last = LEN_TRIM(current % text)
!    if( pos > 0 ) then
!       print*, 'pos =', pos, '  last= ', last
!    endif

    if ( pos > 0 .AND. pos == last .AND. current%text(1:1) /= '!' ) THEN
     if( last < 73 ) then
     current%text(pos:pos) = ' '
     current%text(73:73) = '&'
     endif

     WRITE(9, '(a)') TRIM(current % text)

    else
    WRITE(9, '(a)') TRIM(current % text)
    endif
    IF (ASSOCIATED(current, tail)) EXIT
    IF (.NOT. ASSOCIATED(current)) EXIT
    current => current % next
  END DO

  CLOSE(8)
  CLOSE(9)
END DO

STOP


CONTAINS


SUBROUTINE do_loop_fixup(start)
!     Convert DO-loops from:    DO xxx i=1,n    To:   DO i=1,n
!                           xxx CONTINUE              END DO
!     start points to the first line of the DO loop
!     lab is the label, and lab_length is its length

TYPE (code), POINTER :: start

!     Local variables

TYPE (code), POINTER :: current, end_loop, new_line
INTEGER              :: pos, len2
LOGICAL              :: continued, referenced, nest_different
CHARACTER (LEN=5)    :: label2

current => start % next                          ! Find end of loop
DO
  IF (current % label == lab(:lab_length)) THEN
    continued = (INDEX(current % text, 'CONTINUE') > 0)
    EXIT
  END IF
  IF (.NOT. ASSOCIATED(current)) RETURN
  current => current % next
END DO

end_loop => current
IF (continued) THEN
  end_loop % text = 'END DO'
ELSE
                                       ! Final line of DO loop is not CONTINUE
                                       ! Guard against long instruction which
                                       ! continues onto next line.
  DO
    IF (last_char(current % text) == '&') THEN
      current => current % next
      end_loop => current
    ELSE
      EXIT
    END IF
  END DO
  ALLOCATE(new_line)
  new_line % text = 'END DO'
  new_line % next => current % next
  new_line % label = ' '
  current % next => new_line
  end_loop => new_line
END IF

!     If possible, replace references to the label with CYCLE.
!     Set referenced = .TRUE. if the label is on a line with a command
!     other than CONTINUE, and there is a jump to that line.

current => start % next
referenced = .FALSE.
nest_different = .FALSE.
DO
  IF (current % text == 'END DO') EXIT
  IF (current % text(1:1) /= '!') THEN
                                                 ! Tested for nested DO's
    IF (INDEX(current % text, 'DO ') > 0 .OR.     &
        INDEX(current % text, 'DO,') > 0) THEN
      pos = INDEX(current % text, lab(:lab_length))
      IF (pos > 0) THEN                          ! Remove label
        current % text = current % text(:pos-1) //    &
                         current % text(pos+lab_length+1:)
                                                 ! Add extra END DO
        ALLOCATE(new_line)
        new_line % text = 'END DO'
        new_line % next => end_loop % next
        end_loop % next => new_line
        end_loop => new_line
      ELSE                             ! Nested but with different label
        IF (.NOT. nest_different) THEN
          pos = INDEX(current % text, 'DO')
          text = ADJUSTL( current % text(pos+3:) )
          last = INDEX(text, ' ')
          label2 = text(:last-1)
          len2 = LEN_TRIM(label2)
          nest_different = .TRUE.
        END IF
      END IF
    ELSE
      pos = INDEX(current % text, 'GO TO ')
      IF (pos > 0 .AND. INDEX(current % text, lab(:lab_length)) > 0) THEN
                             ! Avoid nested loops & computed GOTOs
        IF (continued .AND. .NOT. nest_different .AND.  &
            INDEX( current % text(pos+5:), ')') == 0) THEN
          current % text(pos:) = 'CYCLE'
        ELSE
          referenced = .TRUE.
        END IF
      END IF
    END IF

    IF (INDEX(current % text, 'IF') > 0) THEN
      last = LEN_TRIM(current % text)
      IF (SCAN(current % text(last:last), numbers) > 0) THEN
        CALL fix_3way_IF(current)
      END IF
    END IF
    IF (current % label == lab) THEN
                                               ! Remove label from line
      IF (.NOT. referenced) THEN
        current % label = ' '
        pos = INDEX(current % text, ' ')
        current % text = ADJUSTL(current % text(pos+1:))
      END IF
    ELSE IF (nest_different) THEN
      nest_different = (current % label == label2)
    END IF
  END IF

  IF (.NOT. ASSOCIATED(current)) RETURN
  IF (ASSOCIATED(current, end_loop)) EXIT
  current => current % next
END DO

!     Remove label from the original DO instruction.

pos = INDEX(start % text, lab(:lab_length))
IF (pos > 0) start % text = start % text(:pos-1) //    &
                            start % text(pos+lab_length+1:)

!     Find if next executable line after the end of the DO loop is labelled.
!     If so, it may be possible to replace a GO TO xxx with EXIT.

DO
  current => current % next
  IF (current % text(1:1) /= '!' .AND. current % text(1:1) /= ' ') EXIT
END DO

label2 = current % label
IF (label2 == ' ') RETURN
len2 = LEN_TRIM(label2)

current => start % next
DO
  pos = INDEX(current % text, label2(:len2))
  IF (pos > 0) THEN
    pos = INDEX(current % text, 'GO TO ')
    IF (pos > 0) current % text = current % text(:pos-1) // 'EXIT'
  END IF

  IF (ASSOCIATED(current, end_loop)) EXIT
  current => current % next
END DO

RETURN
END SUBROUTINE do_loop_fixup



SUBROUTINE fix_3way_IF(start)
!     Convert 3-way IFs to IF () THEN .. ELSE IF () THEN .. ELSE

TYPE (code), POINTER :: start

!     Local variables

TYPE (code), POINTER :: current
INTEGER              :: pos1, count, length, pos2, i, lab1, lab2, lab3, lenq, &
                        next_label, lenz
CHARACTER (LEN=1)    :: ch
CHARACTER (LEN=128)  :: quantity
CHARACTER (LEN=3)    :: zero_txt

current => start
length = LEN_TRIM(current % text)

!     Find closing bracket to match the opening bracket.
!     Only cases with the closing bracket on the same line are converted.

pos1 = INDEX(current % text, ' IF') + 3
IF (pos1 == 0) RETURN
pos1 = INDEX(current % text(pos1:), '(') + pos1 - 1
count = 1
pos2 = pos1 + 1
DO
  i = SCAN(current % text(pos2:length), '()')
  IF (i == 0) RETURN
  pos2 = i + pos2 - 1
  IF (current % text(pos2:pos2) == '(') THEN
    count = count + 1
  ELSE
    count = count - 1
  END IF
  IF (count == 0) EXIT
  pos2 = pos2 + 1
END DO

!     See if there are 3 labels after the closing bracket.

READ(current % text(pos2+1:length), *, ERR=100) lab1, lab2, lab3

!     As it is probably very old code, the first alphabetic character in the
!     expression should tell us whether the quantity is REAL or INTEGER.

DO i = pos1+1, pos2-1
  ch = current % text(i:i)
  IF (ch >= 'i' .AND. ch <= 'n') THEN
    zero_txt = '0'
    lenz = 1
    EXIT
  ELSE IF (ch >= 'a' .AND. ch <= 'z') THEN
    zero_txt = '0.0'
    lenz = 3
    EXIT
  ELSE IF (i == pos2-1) THEN
    RETURN
  END IF
END DO

quantity = current % text(pos1:pos2)
lenq = LEN_TRIM(quantity)

!     Find the next executable line to see if it is labelled.
next_label = 0
DO
  IF (.NOT. ASSOCIATED(current)) EXIT
  current => current % next
  IF (current % text(1:1) == '!' .OR. LEN_TRIM(current % text) == 0) CYCLE
  IF (LEN_TRIM(current % label) > 0) READ(current % label, *) next_label
  EXIT
END DO
current => start

IF (lab1 == lab2) THEN
  current % text = current % text(:pos2-1) // ' > ' // zero_txt(:lenz) //  &
                   ') THEN'
  CALL insert_and_moveto_newline(current)
  current % text = ' '
  WRITE(current % text, '(a, i5)') 'GO TO ', lab3
  IF (lab1 /= next_label) THEN
    CALL insert_and_moveto_newline(current)
    current % text = 'ELSE'
    CALL insert_and_moveto_newline(current)
    current % text = ' '
    WRITE(current % text, '(a, i5)') 'GO TO ', lab1
  END IF
  CALL insert_and_moveto_newline(current)
  current % text = 'END IF'

ELSE IF (lab2 == lab3) THEN
  current % text = current % text(:pos2-1) // ' < ' // zero_txt(:lenz) //  &
                   ') THEN'
  CALL insert_and_moveto_newline(current)
  current % text = ' '
  WRITE(current % text, '(a, i5)') 'GO TO ', lab1
  IF (lab2 /= next_label) THEN
    CALL insert_and_moveto_newline(current)
    current % text = 'ELSE'
    CALL insert_and_moveto_newline(current)
    current % text = ' '
    WRITE(current % text, '(a, i5)') 'GO TO ', lab2
  END IF
  CALL insert_and_moveto_newline(current)
  current % text = 'END IF'

ELSE IF (lab1 == lab3) THEN
  current % text = current % text(:pos2-1) // ' == ' // zero_txt(:lenz) //  &
                   ') THEN'
  CALL insert_and_moveto_newline(current)
  current % text = ' '
  WRITE(current % text, '(a, i5)') 'GO TO ', lab2
  IF (lab1 /= next_label) THEN
    CALL insert_and_moveto_newline(current)
    current % text = 'ELSE'
    CALL insert_and_moveto_newline(current)
    current % text = ' '
    WRITE(current % text, '(a, i5)') 'GO TO ', lab1
  END IF
  CALL insert_and_moveto_newline(current)
  current % text = 'END IF'

ELSE
  current % text = current % text(:pos2-1) // ' < ' // zero_txt(:lenz) //  &
                   ') THEN'
  CALL insert_and_moveto_newline(current)
  current % text = ' '
  WRITE(current % text, '(a, i5)') 'GO TO ', lab1
  CALL insert_and_moveto_newline(current)
  current % text = 'ELSE IF ' // quantity(1:lenq-1) // ' == ' // &
                   zero_txt(:lenz) // ') THEN'
  CALL insert_and_moveto_newline(current)
  current % text = ' '
  WRITE(current % text, '(a, i5)') 'GO TO ', lab2
  IF (lab3 /= next_label) THEN
    CALL insert_and_moveto_newline(current)
    current % text = 'ELSE'
    CALL insert_and_moveto_newline(current)
    current % text = ' '
    WRITE(current % text, '(a, i5)') 'GO TO ', lab3
  END IF
  CALL insert_and_moveto_newline(current)
  current % text = 'END IF'

END IF

100 RETURN
END SUBROUTINE fix_3way_IF



SUBROUTINE insert_and_moveto_newline(current)
TYPE (code), POINTER :: current

!     Local variable
TYPE (code), POINTER :: new_line

ALLOCATE(new_line)
new_line % next => current % next
current % next => new_line
current => new_line

RETURN
END SUBROUTINE insert_and_moveto_newline



SUBROUTINE find_declarations( start, tail, first_decl, last_decl )
! Find the first & last declaration lines in a program unit.

TYPE (code), POINTER :: start, tail
TYPE (code), POINTER :: first_decl, last_decl

! Local variables
CHARACTER (LEN=9), PARAMETER :: declaration(13) = (/ 'IMPLICIT ', 'INTEGER  ', &
                                'REAL     ', 'DOUBLE   ', 'LOGICAL  ', &
                                'COMPLEX  ', 'DIMENSION', 'EXTERNAL ', &
                                'DATA     ', 'COMMON   ', 'PARAMETER', &
                                'SAVE     ', 'CHARACTER' /)
TYPE (code), POINTER         :: current
INTEGER                      :: pos, length, i

NULLIFY( first_decl, last_decl )

! Search for first declaration
current => start % next
search1: DO
  IF ( current % text(1:1) /= '!' .AND.  current % text(1:1) /= ' ' ) THEN
    pos = SCAN( current % text(1:13), delimiters )
    IF (pos > 0) THEN
      length = MIN(9, pos - 1)
      IF (length >= 4) THEN
        DO i = 1, 13
          IF ( current % text(:length) == declaration(i)(:length) ) THEN
            first_decl => current
            EXIT search1
          END IF
        END DO
      END IF
    END IF
  END IF

  current => current % next
  IF ( ASSOCIATED( current, tail ) ) RETURN
END DO search1

! Search for last declaration

last_decl => first_decl
DO
  IF ( current % text(1:1) /= '!' .AND.  current % text(1:1) /= ' ' ) THEN
    pos = INDEX( current % text, '=' )
    IF (pos > 0) THEN
      IF (pos < 12) RETURN
      IF (current % text(1:9) /= 'PARAMETER' .AND.  &
          current % text(1:9) /= 'CHARACTER') RETURN
    END IF

    IF ( current % text(1:4) == 'CALL' ) RETURN

    IF ( current % text(1:2) == 'IF' ) THEN
      IF ( current % text(3:3) == ' ' ) RETURN
      IF ( current % text(3:3) == '(' ) RETURN
    END IF

    IF ( current % text(1:3) == 'DO ' ) RETURN

! Skip continuation lines

    DO
      IF ( last_char(current % text) /= '&' ) EXIT
      current => current % next
    END DO

    last_decl => current
  END IF
  current => current % next
  IF ( ASSOCIATED( current, tail ) ) RETURN
END DO

RETURN
END SUBROUTINE find_declarations


SUBROUTINE get_arg_list()
! Extract list of dummy arguments

! Local variables
INTEGER :: pos, last

current => start_prog_unit
numb_arg = 0
DO                                 ! Find '(' if there are any arguments
  pos = INDEX( current % text, '(')
  IF (pos == 0) THEN
    IF ( last_char( current % text ) /= '&' ) RETURN
    current => current % next
  ELSE
    EXIT
  END IF
END DO
pos = pos + 1

NULLIFY( arg_start )
ALLOCATE( arg_start )
first_arg = .TRUE.
DO                                 ! Loop through lines of arguments
  last = SCAN(current % text(pos:), ',)')
  IF (last == 0) THEN
    IF (last_char( current % text ) /= '&' ) EXIT
    current => current % next
    pos = 1
  ELSE
    last = last + pos - 1
    NULLIFY( arg )
    ALLOCATE( arg )
    IF (first_arg) THEN
      IF (LEN_TRIM(current % text(pos:last-1)) == 0) EXIT
      arg_start => arg
      first_arg = .FALSE.
      NULLIFY( last_arg )
      ALLOCATE( last_arg )
    ELSE
      last_arg % next => arg
    END IF
    numb_arg = numb_arg + 1
    last_arg => arg

    arg % name = ADJUSTL( current % text(pos:last-1) )
    arg % intention = 0
    arg % var_type = ' '
    arg % dim = 0
    pos = last + 1
  END IF
END DO

RETURN
END SUBROUTINE get_arg_list



SUBROUTINE get_var_types()
! Search thru the declarations for the types of dummy arguments

current => first_decl
DO
  text = current % text(:30)
  IF (text(:4) == 'REAL' .OR. text(:7) == 'INTEGER' .OR.     &
      text(:6) == 'DOUBLE' .OR. text(:9) == 'CHARACTER' .OR. &
      text(:7) == 'LOGICAL' .OR. text(:7) == 'COMPLEX') THEN
                                   ! Copy the variable type to vtype
    last = INDEX(text, ' ::') - 1
    IF (last < 0) THEN
      last = INDEX(text, '*')
      IF (last == 0) THEN
        last = 24
      ELSE
        last = INDEX(text(last+2:), ' ') + last
      END IF
    END IF
    vtype = text(:last)

    DO                             ! Loop thru any continuation statements
      arg => arg_start
      DO                           ! Loop thru arguments
        pos = find_delimited_name(current % text, arg % name)
        IF (pos > 0) THEN
          arg % var_type = vtype
                                   ! Check for dimensions in brackets
          left_brac = INDEX(current % text(pos+1:), '(')
          right_brac = 0
          IF (left_brac > 0) THEN
            left_brac = left_brac + pos
            i = INDEX(current % text(pos+1:left_brac), ',')
            IF (i == 0) THEN
              right_brac = INDEX(current % text(pos+1:), ')') + pos
              arg % dimensions = current % text(left_brac:right_brac)
              arg % dim = 1
            END IF
          END IF

                                   ! Remove name from the declaration
                                   ! Test whether
                                   ! declaration is left without names
        END IF
        IF (ASSOCIATED(arg, last_arg)) EXIT
        arg => arg % next
      END DO

      IF (last_char(current % text) /= '&') EXIT
      current => current % next
    END DO

  ELSE IF (text(:9) == 'DIMENSION') THEN
    DO                             ! Loop thru any continuation statements
      arg => arg_start
      DO                           ! Loop thru arguments
        pos = find_delimited_name(current % text(10:), arg % name)
        IF (pos > 0) THEN
          pos = pos + 9
          left_brac = INDEX(current % text(pos:), '(') + pos - 1
          right_brac = INDEX(current % text(pos:), ')') + pos - 1
          arg % dimensions = current % text(left_brac:right_brac)
          arg % dim = 1
                                   ! Remove name from the declaration
                                   ! Now remove the comma
                                   ! Test if declaration is left without names
                                 ! Either join onto next line
                                 ! ... this could remove the line to which
                                 !     last_decl points
                                 ! or blank out the declaration
        END IF
        IF (ASSOCIATED(arg, last_arg)) EXIT
        arg => arg % next
      END DO

      IF (last_char(current % text) /= '&') EXIT
      current => current % next
    END DO

  END IF
  IF (ASSOCIATED(current, last_decl)) EXIT
  current => current % next
END DO

!     If there are any arguments for which the type has not been determined,
!     use the implicit types

arg => arg_start
DO
  IF (arg % var_type == ' ')   &
      arg % var_type = implicit_type(arg % name(1:1))
  IF (ASSOCIATED(arg, last_arg)) EXIT
  arg => arg % next
END DO

RETURN
END SUBROUTINE get_var_types


SUBROUTINE get_intents()
! Search thru the body of the current program unit to try to determine
! the intents of dummy arguments.

CHARACTER (LEN=80) :: last_part
INTEGER            :: j, nbrac

DO
  IF (current % text(1:1) /= '!' .AND. current % text(1:1) /= ' ') THEN
    statement = current % text
    IF (statement(1:3) == 'IF ' .OR. statement(1:3) == 'IF(') THEN
                                       ! Split line into two parts
                                       ! IF (condition) | last_part
      i = INDEX(statement, '(')
      length = LEN_TRIM(statement)
      nbrac = 1
      DO j = i+1, length-1
        IF (statement(j:j) == ')') THEN
          nbrac = nbrac - 1
          IF (nbrac == 0) EXIT
        ELSE IF (statement(j:j) == '(') THEN
          nbrac = nbrac + 1
        END IF
      END DO
      IF (j < length) THEN
        last_part = statement(j+1:)
      ELSE
        last_part = ' '
      END IF
      statement = statement(:j)
                                       ! It is assumed that a variable inside
                                       ! an IF-expression cannot be altered
      arg => arg_start
      DO
        i = find_delimited_name(statement, arg % name)
        IF (i > 0) THEN
          IF (arg % intention == 0) arg % intention = 1
        END IF
        IF (ASSOCIATED(arg, last_arg)) EXIT
        arg => arg % next
      END DO
      statement = last_part
    END IF

    pos = INDEX(statement, '=', BACK=.TRUE.)
    IF (pos > 0) THEN
      IF (statement(pos-1:pos-1) /= '=' .AND.  &
          statement(pos-1:pos-1) /= '/' .AND.  &
          statement(pos-1:pos-1) /= '<' .AND.  &
          statement(pos-1:pos-1) /= '>') THEN

                                       ! Look for each argument name;
                                       ! is it before or after '='?
        arg => arg_start
        DO
          i = find_delimited_name(statement, arg % name)
          IF (i > 0) THEN
            IF (i < pos) THEN
              arg % intention = IOR(arg % intention, 2)
            ELSE
              IF (arg % intention == 0) arg % intention = 1
            END IF
          END IF
          IF (ASSOCIATED(arg, last_arg)) EXIT
          arg => arg % next
        END DO
      END IF
    END IF
  END IF

  IF (ASSOCIATED(current, end_prog_unit)) EXIT
  current => current % next
END DO

RETURN
END SUBROUTINE get_intents

END PROGRAM to_f90



SUBROUTINE mark_text(text, n_marks, pos1, pos2)

!     Look for exclamation marks or quotes to find any text which must be
!     protected from case changes.
!     It is assumed that strings are NOT continued from one line to the next.
IMPLICIT NONE

CHARACTER (LEN = *), INTENT(IN)  :: text
INTEGER (KIND(123)), INTENT(OUT) :: n_marks, pos1(:), pos2(:)

!     Local variables
INTEGER (KIND(123))  :: mark, start, pos_exclaim, pos_sngl_quote, pos_dbl_quote

  mark = 1
  start = 1
  DO
    pos_exclaim = INDEX(text(start:80), '!')
    pos_sngl_quote = INDEX(text(start:80), "'")
    pos_dbl_quote = INDEX(text(start:80), '"')
    IF (pos_exclaim == 0) pos_exclaim = 81
    IF (pos_sngl_quote == 0) pos_sngl_quote = 81
    IF (pos_dbl_quote == 0) pos_dbl_quote = 81
    pos1(mark) = MIN(pos_exclaim, pos_sngl_quote, pos_dbl_quote)
    IF (pos1(mark) == 81) THEN                 ! No more protected regions
      n_marks = mark - 1
      EXIT
    ELSE IF (pos_exclaim == pos1(mark)) THEN   ! Rest of line is a comment
      pos1(mark) = pos1(mark) + start - 1
      pos2(mark) = 80
      n_marks = mark
      EXIT
    ELSE IF (pos_sngl_quote == pos1(mark)) THEN     ! Look for matching quote
      pos1(mark) = pos1(mark) + start - 1
      pos2(mark) = INDEX(text(pos1(mark)+1:80), "'")
      IF (pos2(mark) > 0) THEN
        pos2(mark) = pos1(mark) + pos2(mark)
        start = pos2(mark) + 1
        mark = mark + 1
        CYCLE
      ELSE
        pos2(mark) = 80
        n_marks = mark
        EXIT
      END IF
    ELSE IF (pos_dbl_quote == pos1(mark)) THEN      ! Look for matching quote
      pos1(mark) = pos1(mark) + start - 1
      pos2(mark) = INDEX(text(pos1(mark)+1:80), '"')
      IF (pos2(mark) > 0) THEN
        pos2(mark) = pos1(mark) + pos2(mark)
        start = pos2(mark) + 1
        mark = mark + 1
        CYCLE
      ELSE
        pos2(mark) = 80
        n_marks = mark
        EXIT
      END IF
    END IF
  END DO

RETURN
END SUBROUTINE mark_text


SUBROUTINE convert_text(text, n_marks, pos1, pos2)

!     Convert unprotected text to upper case if it is a FORTRAN word,
!     otherwise convert to lower case.
IMPLICIT NONE

CHARACTER (LEN = *), INTENT(IN OUT) :: text
INTEGER (KIND(123)), INTENT(IN)     :: n_marks
INTEGER (KIND(123)), INTENT(IN OUT) :: pos1(:), pos2(:)

!     Local variables

INTEGER (KIND(123))   :: length, inc = ICHAR('A') - ICHAR('a'), indx(27),     &
                         pos, mark, i, i1, j, j1, j2, ptr
LOGICAL               :: matched
CHARACTER (LEN = 11)  :: fortran_word(185)
CHARACTER (LEN = 4)   :: compare(6) = (/ ".LT.", ".LE.", ".EQ.", ".GE.",      &
                                         ".GT.", ".NE." /)
CHARACTER (LEN = 2)   :: replacement(6) = (/ "< ", "<=", "==", ">=", "> ",    &
                                             "/=" /)

DATA fortran_word/                                                            &
   "ABS", "ACCESS", "ACOS", "AIMAG", "AINT", "ALOG", "ALOG10",                &
   "AMAX0", "AMAX1", "AMIN0", "AMIN1", "AMOD", "AND", "ANINT", "APPEND",      &
   "ASIN", "ASSIGN", "ATAN", "ATAN2", "BACKSPACE", "BLANK", "BLOCK",          &
   "BLOCKDATA", "BLOCKSIZE", "CALL", "CCOS", "CDABS", "CDCOS", "CDEXP",       &
   "CDLOG", "CDSIN", "CDSQRT", "CEXP", "CHAR", "CHARACTER", "CLOG",           &
   "CLOSE", "CMPLX", "COMMON", "COMPLEX", "CONJG", "CONTINUE", "COS", "COSH", &
   "CSIN", "CSQRT", "DABS", "DACOS", "DASIN", "DATA", "DATAN", "DATAN2",      &
   "DBLE", "DCMPLX", "DCONJG", "DCOS", "DCOSH", "DELETE", "DEXP",             &
   "DIMAG", "DIMENSION", "DINT", "DIRECT", "DLOG", "DLOG10", "DMAX1",         &
   "DMIN1", "DMOD", "DNINT", "DO", "DOUBLE", "DSIGN", "DSIN", "DSINH",        &
   "DSQRT", "DTAN", "DTANH", "ELSE", "ELSEIF", "END", "ENDFILE", "ENDIF",     &
   "ENTRY", "EQ", "EQUIVALENCE", "ERR", "EXIST", "EXIT", "EXP", "EXTERNAL",   &
   "FILE", "FLOAT", "FMT", "FORM", "FORMAT", "FORMATTED", "FUNCTION",         &
   "GE", "GOTO", "GO", "GT", "IABS", "IAND", "ICHAR", "IDINT", "IDNINT",      &
   "IEOR", "IF", "IFIX", "IMPLICIT", "INCLUDE", "INDEX", "INPUT",             &
   "INQUIRE", "INT", "INTEGER", "INTRINSIC", "IOSTAT", "ISIGN", "KEEP",       &
   "LE", "LEN", "LGE", "LGT", "LLE", "LLT", "LOG", "LOG10", "LOGICAL", "LT",  &
   "MAX", "MAX0", "MAX1", "MIN", "MIN0", "MIN1", "MOD", "NAME", "NAMELIST",   &
   "NAMED", "NE", "NEQV", "NEW", "NEXTREC", "NONE", "NOT", "NUMBER", "OLD",   &
   "OPEN", "OPENED", "OR", "PARAMETER", "PAUSE", "POSITION", "PRECISION",     &
   "PRINT", "PROGRAM", "READ", "REAL", "REC", "RECL", "RETURN",               &
   "REWIND", "SAVE", "SCRATCH", "SEQUENTIAL", "SIGN", "SIN", "SINH", "SNGL",  &
   "SPACE", "SQRT", "STATUS", "STOP", "SUBROUTINE", "TAN", "TANH",            &
   "THEN", "TO", "TYPE", "UNFORMATTED", "UNIT", "UNKNOWN", "WHILE", "WRITE"   &
/
!          A   B   C   D   E   F   G    H    I    J    K    L    M    N    O
!          P    Q    R    S    T    U    V    W    X    Y    Z
DATA indx/ 1, 20, 25, 47, 78, 91, 98, 102, 102, 120, 120, 121, 131, 138, 148,   &
         152, 158, 158, 164, 176, 181, 184, 184, 186, 186, 186, 186/

IF (pos1(1) == 1 .AND. pos2(1) == 80) RETURN      ! Entire line protected

pos = 1
mark = 1
length = LEN_TRIM(text)
DO                                     ! Convert to upper case
  IF (n_marks >= mark .AND. pos == pos1(mark)) THEN
    pos = pos2(mark) + 1
    mark = mark + 1
    IF (pos >= length) EXIT
  END IF
  IF (text(pos:pos) >= 'a' .AND. text(pos:pos) <= 'z')           &
              text(pos:pos) = CHAR ( ICHAR(text(pos:pos)) + inc )
  pos = pos + 1
  IF (pos > length) EXIT
END DO

!     Search for `words' in text.
!     Convert to lower case if they are not FORTRAN words.
i1 = 1
pos = 1
mark = 1
DO
  IF (pos > length) EXIT
  IF (n_marks >= mark .AND. pos >= pos1(mark)) THEN
    pos = pos2(mark) + 1
    i1 = pos
    mark = mark + 1
    IF (pos >= length) EXIT
  END IF

  DO
    IF ((text(pos:pos) >= 'A' .AND. text(pos:pos) <= 'Z')        &
        .OR. (text(pos:pos) >= '0' .AND. text(pos:pos) <= '9')   &
        .OR. text(pos:pos) == '_') THEN
      pos = pos + 1
      CYCLE
    ELSE
      EXIT
    END IF
  END DO

  pos = pos - 1
! Now i1 & pos = positions of 1st & last characters of current string

  IF (pos < i1) THEN                ! Single non-alphanumeric character
    pos = i1 + 1
    i1 = pos
    CYCLE
  END IF

  ptr = ICHAR(text(i1:i1)) - ICHAR('A') + 1
  IF (ptr < 1 .OR. ptr > 26) THEN
    pos = pos + 1
    IF (pos > length) EXIT
    i1 = pos
    CYCLE
  END IF

  matched = .FALSE.
  IF (pos > i1) THEN
    j1 = indx(ptr)
    j2 = indx(ptr+1) - 1
    DO j = j1, j2
      IF (text(i1:pos) == fortran_word(j)) THEN
        matched = .TRUE.
        EXIT
      END IF
    END DO ! j = j1, j2
  END IF

! Replace .LT. with <, etc.
  IF (matched .AND. i1 > 1) THEN
    IF(text(i1-1:i1-1) == '.') THEN
      DO j = 1, 6
        IF (text(i1-1:pos+1) == compare(j)) THEN
          text(i1-1:pos+1) = ' ' // replacement(j) // ' '
          EXIT
        END IF
      END DO
      DO                                 ! Remove excess blanks
        i1 = MAX(i1, 3)
        j1 = INDEX(text(i1-2:pos+2), '  ')
        IF (j1 == 0) EXIT
        j1 = j1 + i1 - 3
        text(j1:) = text(j1+1:)
        pos2(mark) = pos2(mark) - 1      ! Adjust mark positions
        DO i = mark+1, n_marks
          pos1(i) = pos1(i) - 1
          pos2(i) = pos2(i) - 1
        END DO
        pos = pos - 1
      END DO
    END IF
  END IF

! Output line of text to screen if it contains SUBROUTINE or FUNCTION.
! Convert ENDIF to END IF, ELSEIF to ELSE IF, and GOTO to GO TO.
  IF (matched) THEN
    IF (text(i1:pos) == 'SUBROUTINE' .OR. text(i1:pos) == 'FUNCTION') THEN
      WRITE(*, '(1x, a)') text(1:length)
    ELSE IF (text(i1:pos) == 'ENDIF') THEN
      text(i1:) = 'END IF' // text(pos+1:)
      pos = pos + 1
    ELSE IF (text(i1:pos) == 'ELSEIF') THEN
      text(i1:) = 'ELSE IF' // text(pos+1:)
      pos = pos + 1
    ELSE IF (text(i1:pos) == 'GOTO') THEN
      text(i1:) = 'GO TO' // text(pos+1:)
      pos = pos + 1
    END IF
  END IF

! If text is not matched, convert to lower case, if necessary.
  IF (.NOT. matched) THEN
    DO j = i1, pos
      IF (text(j:j) >= 'A' .AND. text(j:j) <= 'Z')              &
              text(j:j) = CHAR ( ICHAR(text(j:j)) - inc )
    END DO
  END IF

  pos = pos + 1
  IF (pos > length) EXIT
  i1 = pos
END DO

END SUBROUTINE convert_text



SUBROUTINE remove_data_blanks(text)
! Remove any blanks embedded between numerical digits in DATA statements

IMPLICIT NONE
CHARACTER (LEN=*), INTENT(IN OUT) :: text

! Local variables
INTEGER            :: length, pos, i1
CHARACTER (LEN=10) :: numbers = '1234567890'

length = LEN_TRIM(text)
i1 = 2
DO
  pos = INDEX(text(i1:length), ' ')
  IF (pos == 0) EXIT
  i1 = i1 + pos - 1
  IF (SCAN(text(i1-1:i1-1), numbers) > 0 .AND.  &
      SCAN(text(i1+1:i1+1), numbers) > 0) THEN
    text = text(:i1-1) // text(i1+1:length)
    length = length - 1
  END IF
  i1 = i1 + 2
  IF (i1 > length) EXIT
END DO

RETURN
END SUBROUTINE remove_data_blanks


FUNCTION last_char( text ) RESULT(ch)
! Return the last character on a line
IMPLICIT NONE

CHARACTER (LEN=*), INTENT(IN) :: text
CHARACTER (LEN=1)             :: ch

! Local variable
INTEGER :: last

last = LEN_TRIM( text )
IF (last == 0) THEN
  ch = ' '
ELSE
  ch = text(last:last)
END IF

RETURN
END FUNCTION last_char


FUNCTION find_delimited_name (text, name) RESULT(pos)
! Find a name in a character string with delimiters either side of it,
! or after it if it starts at position 1.
! An extended version of the intrinsic INDEX.
! pos = the position of the first character of name in text (= 0 if not found).
! N.B. When the name is short (e.g. i or n) it could occur as part of some
!      other name.

IMPLICIT NONE
CHARACTER (LEN=*), INTENT(IN) :: text, name
INTEGER                       :: pos

! Local variables
INTEGER :: i1, ltext, lname

i1 = 1
ltext = LEN_TRIM(text)
lname = LEN_TRIM(name)
DO
  pos = INDEX(text(i1:ltext), TRIM(name))
  IF (pos == 0) RETURN
  pos = pos + i1 - 1
  IF (pos > 1) THEN
    IF ( SCAN(text(pos-1:pos-1), ' <=+-/*,') > 0 ) THEN
      IF ( SCAN(text(pos+lname:pos+lname), ' >=(+-/*,') > 0 ) RETURN
    END IF
  ELSE
    IF ( SCAN(text(pos+lname:pos+lname), ' >=(+-/*,') > 0 ) RETURN
  END IF
  i1 = pos + lname
  IF (i1 + lname > ltext) EXIT
END DO

pos = 0

RETURN
END FUNCTION find_delimited_name

