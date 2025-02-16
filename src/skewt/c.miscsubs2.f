c
c
c                MISCELLANEOUS SUBROUTINES (CHARACTER PROCESSING)
c
c                       Richard Carpenter, Univ. of Okla.
c
c                                  August 1991
c
c
c
c     SUBROUTINE  DESCRIPTION
c     ----------  -----------
c     Getnum1        : Convert a string of numbers into integers
c     Getposn        : Find the location of a string of numbers
c     Getdate        : Get the date and time (Cray)
c     Chardate        : Convert a numerical date obtained with Getdate into character
c     Chrdate2        : Convert a numerical date obtained with Getdate into character
c   F Fnblnk        : Location of first non-blank
c   F Lnblnk        : Location of last non-blank
c     Lwrc        : Convert a string to lowercase
c     Uprc        : Convert a string to uppercase
c   F Rindex        : Finds the right-most location of a string
c     Parse        : Find the words in a character string
c     Ljust        : Left-justify a character string
c     Pslab        : Print an array 
c     Chkexist        : See whether a file exists
c     Csqueeze        : Remove all blanks from a string
c     Csqueeze2        : Remove all consecutive blanks from a string
c
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########        GETNUM1         ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
c
      Subroutine Getnum1 (name, ipos1, ipos2, number)
      Implicit None
      Integer ipos1, ipos2, number, iten, i
      Character*(*) name
c
c     Convert the  Name(Ipos1:Ipos2)  into  Number
c     Note that '0' is ASCII 48.
c
      number         = 0
      iten        = 1
      Do 1000        i=ipos2,ipos1,-1
      number        = number + ( Ichar (name(i:i)) - 48 ) * iten
      iten        = iten * 10
 1000 Continue
      End
c
c
c
c
c                    ########################################
c                    ########################################
c                    ########################################
c                    ########                        ########
c                    ########         GETPOSN        ########
c                    ########                        ########
c                    ########################################
c                    ########################################
c                    ########################################
c
c
      Subroutine Getposn (dumpfile, ipos1, ipos2)
      Implicit None
      Integer        ipos1, ipos2, jchar, i
      Character*(*) dumpfile
c
c     Find the location of the number of the file, Dumpfile(Ipos1:Ipos2)
c
c     Note that 'index' is an m4 command
c     Ipos2     = old Lendmpfl = length of filename and location of last char
c     Ipos1     = leading char of the digit of the filename (after the H or R)
c
      ipos2     = Index (dumpfile,' ') - 1
c      
      Do 100    i=ipos2-1,1,-1
c     jchar     = Ichar (dumpfile(i:i))
c     If (jchar.LT.48 .OR. jchar.GT.57) Go To 101
      If (dumpfile(i:i).LT.'0' .OR. dumpfile(i:i).GT.'9') Goto 101
 100  Continue
 101  Continue
      ipos1     = i + 1
c
      End
c
c
c
c
c                    ########################################
c                    ########################################
c                    ########################################
c                    ########                        ########
c                    ########         GETDATE        ########
c                    ########                        ########
c                    ########################################
c                    ########################################
c                    ########################################
c
c
      Subroutine Getdate (idate) 
      Implicit None
      Integer idate(7)
!      Character*9 char
      CHARACTER*10  datestr, timestr, zonestr
c     External Date, Clock
c
c     Get date and time from Cray
c
c     idate(1) : Day of the week
c     idate(2) : Day of the month
c     idate(3) : Month
c     idate(4) : Year
c     idate(5) : Hour
c     idate(6) : Minute
c     idate(7) : Second
c
!      Call Date (char)                ! returns date in the form "10/26/90"
!      Read (char, 9901) idate(3), idate(2), idate(4)
!      !Call Clock (char)                ! returns time in the form "10:26:00"
!      Call Time (char)                ! returns time in the form "10:26:00"
!      Read (char, 9901) idate(5), idate(6), idate(7)
!      idate(1)        = 0
!      idate(4)        = idate(4) + 1900
! 9901 Format (3(i2,1x))
!  
      idate(1) = 0
      CALL DATE_AND_TIME(datestr,timestr,zonestr)
      READ(datestr,9901) idate(4),idate(2),idate(3)
      READ(timestr,9902) idate(5),idate(6),idate(7)
 9901 FORMAT(i4,1x,2(i2,1x))
 9902 FORMAT(3(i2,1x))

      END
c
c
c
c
c                    ########################################
c                    ########################################
c                    ########################################
c                    ########                        ########
c                    ########        CHARDATE        ########
c                    ########                        ########
c                    ########################################
c                    ########################################
c                    ########################################
c
c
      Subroutine Chardate (idate, cdate) 
      Implicit None
      Integer idate(7), i, j, k
      Character*24 cdate, wkday, month*39
      Data wkday /'   SunMonTueWedThuFriSat'/, 
     >     month /'   JanFebMarAprMayJunJulAugSepOctNovDec'/
c
c     Convert integer date and time into character form.
c        Mon 29-Oct-1990 08:36:00
c        123456789012345678901234
c
      j                = 3 * idate(1)
      k                = 3 * idate(3)
      Write (cdate, 9901) wkday(j+1:j+3), idate(2), month(k+1:k+3), 
     >                           (idate(i),i=4,7)
 9901 Format (a3, i3.2, '-', a3, '-', i4, i3.2, ':', i2.2, ':', i2.2)
c
      End
c
c
c
c
c                    ########################################
c                    ########################################
c                    ########################################
c                    ########                            ########
c                    ########        CHRDATE2        ########
c                    ########                            ########
c                    ########################################
c                    ########################################
c                    ########################################
c
c
      Subroutine Chrdate2 (idate, cdate) 
      Implicit None
      Integer idate(7), i, j, k
      Character*(*) cdate
c
c     Convert integer date and time into character form.
c        10/29/90 08:36:00
c        123456789012345678901234
c
      j                = 3 * idate(1)
      k                = 3 * idate(3)
      Write(cdate, 9901)idate(3),idate(2),idate(4)-1900,(idate(i),i=5,7)
 9901 Format (i2.2,2('/',i2.2),1x, i2.2,2(':',i2.2))
c
      End
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########         LNBLNK         ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Integer Function Lnblnk (string)
      Implicit None
c
c     This routine is identical to the Alliant routine having the same name.
c     Returns 0 if strings is all blanks or uninitialized (all 0's).
c
      Character*(*) string
      Integer        i
c
      Lnblnk        = 0
      If (Ichar(string(1:1)).EQ.0) Return
c
      Do i=Len(string),1,-1
        If (string(i:i).NE.' ') Then
          Lnblnk    = i
          Return
        End If
      End Do
c
c
      End
c
c
c
c
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########         FNBLNK         ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Integer Function Fnblnk (string)
      Implicit None
c
c     Finds the first non-blank character.
c     Returns 0 if strings is all blanks or uninitialized (all 0's).
c
      Character*(*) string
      Integer        i
c
      Fnblnk        = 0
      If (Ichar(string(1:1)).EQ.0) Return
c
      Do i=1,Len(string)
        If (string(i:i).NE.' ') Then
          Fnblnk    = i
          Return
        End If
      End Do
c
c
      End
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########          LWRC          ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Subroutine Lwrc (line)
      Implicit None
c
c     Convert a string to lowercase.
c
      Character*(*) line
      Integer        i, ichr
c
      Do  i=1,Len(line)
      ichr        = Ichar (line(i:i))
      If (ichr.GE.65 .AND. ichr.LE.90)  line(i:i) = Char (ichr+32)
      End Do
      End
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########          UPRC          ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Subroutine Uprc (line)
      Implicit None
c
c     Convert a string to uppercase.
c
      Character*(*) line
      Integer        i, ichr
c
      Do  i=1,Len(line)
      ichr        = Ichar (line(i:i))
      If (ichr.GE.97 .AND. ichr.LE.122)  line(i:i) = Char (ichr-32)
      End Do
      End
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########         RINDEX         ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Integer Function IndexR (s1, s2)
      Implicit None
c
c     Note: Apparently RINDEX is the name of some mysterious AIX function.
c
c     This function works much like the Index function, except it finds the
c     _last_ occurance of string s2 within string s1. Returns 0 if substring
c     not found. 
c
      Character*(*) s1, s2
      Integer        j
c
      IndexR        = 0
c
      Do j=Len(s1),1,-1
      If (Index(s1(j:),s2) .NE. 0) Then
        IndexR        = j
        Return
      End If
      End Do
c
      End
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########         PARSE          ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Subroutine Parse (strng, words, nwords, nwmax)
      Implicit None
      Integer nwords, nwmax, j, j0, nw, indx
      Character*(*) strng, words(nwmax)
      Character*256 string
c
c
c     strng        = Input string
c     words        = Returns an array of of left-justified words
c     nwords        = Number of words detected
c     nwmax        = Dimension of words.
c
c
c
      If (Len(strng).GT.Len(string)) Then
        Print '(a,i4,a)', '% Parse: only the first', Len(string), 
     >                    ' characters of input string will be checked.'
      End If
c
c  Store string in temp var, left-justify, squeeze out multiple blanks.
c
      string        = strng
      Call Ljust (string)
      Call Csqueeze2 (string)
c
c
c     indx        = first occurance of whitespace.
c
c
c
c     Print '(1x,2a)', 'Parse: string:<', string
      indx        = Index (string,' ')
      If (indx.EQ.0) Then        ! all one long word
        nwords        = 1
        words(nwords) = string
        Return
      End If
c
c     Print '(1x,2a)', 'Parse: string:<', string
c
c
c     Generate left-justified words.
c
c
c     j0        = First letter of next word
c
      j0        = 1
      Do nw=1,nwmax
      indx        = Index (string,' ')
      If (indx.LE.1) Then        !THIS DOESN'T ALLOW FOR MULTIPLE BLANKS
        nwords        = nw - 1
c       Print *, 'nwords = ', nwords
        Return
      End If
      words(nw)        = string(:indx-1)
c     Print *, 'Parse: nw,j0,indx-1,words(nw): ', nw,j0,indx-1,words(nw)
      If (indx.LE.1) Then
        nwords        = nw - 1
c       Print *, 'nwords = ', nwords
        Return
      End If
c
c
c     Get ready for next time through
c
      string        = string(indx+1:)
      Call Ljust (string)
      End Do

c
      nwords        = nw - 1
c
      End
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########         LJUST          ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
c        This routine returns a left-justified string.  All whitespace up to
c        the first letter is removed.  If the string is blank, a null string
c        will be returned.  
c        Note: Tabs not handled properly yet.
c
      Subroutine Ljust (string)
      Implicit None
      Character*(*) string
      Integer        Fnblnk        
c
      string        = string(Max(1,Fnblnk(string)):)
c
      End
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########       CHECKEXIST       ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Subroutine CheckExist (file)
      Implicit None
      Character*(*) file
      Logical   exist
c
      Inquire (File=file, exist=exist)
      If (.NOT. exist) Then
        Print *, '*** File does not exist: ', file
        Stop
      End If
      End
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########        CHKEXIST        ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Subroutine Chkexist (file,string)
      Implicit None
      Character*(*) file
      Logical   exist
      Character*(*) string
c
      Inquire (File=file, exist=exist)
      If (.NOT. exist) Then
      Print *, string, 'File does not exist: ', file
      Stop
      End If
      End
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########        CSQUEEZE        ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Subroutine Csqueeze (string)
      Implicit None
      Integer        i, Lnblnk
      Character*(*) string
c
c  Removes all blanks from a string
c
      i                = 1
c
      Do While (i.LT.Lnblnk(string))
      If (string(i:i).EQ.' ') Then
        string(i:) = string(i+1:)
      Else
        i        = i + 1
      End If
      End Do
c
c
      End
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########        CSQUEEZE2       ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Subroutine Csqueeze2 (string)
      Implicit None
      Integer        i, Lnblnk
      Character*(*) string
c
c  Removes consequtive blanks from a string.
c
      i                = 1
c
      Do While (i.LT.Lnblnk(string))
      If (string(i:i+1).EQ.'  ') Then
        string(i:) = string(i+1:)
      Else
        i        = i + 1
      End If
      End Do
c
c
      End
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########          PSLAB         ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c        Print a 2-d array, or XZ slabs of a 3-d array, to an output file
c        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Subroutine Pslab (a,file,string,n132,
     >   nx1,nx2, ny1,ny2, nz1,nz2, i1,i2, j1,j2, k1,k2)
      Implicit None
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-
c
c  INPUT:
c     a                = Array
c     file        = Output filename ('-' for output to screen)
c     string        = String written to output file
c     n132        = +-132 for 132-column output (12 or 10 cols of data);
c                  +-80 otherwise (7 or 6 cols of data);
c                  printed format is f10.5 if n132>=0; 1p,e12.4 otherwise.
c     nx1...nz2        = Upper & lower dimensions of A
c     i1...k2        = Portion of the array to print
c
c  OUTPUT: Formatted output is sent to file.
c
c  NOTES:
c     - To print a 2-d array, set k1=k2; a single XY slab will be printed.
c     - The parameter n132 controls both the width of the output (80/132)
c        as well as the format (f10.5 if n132>=0, 1p,e12.4 otherwise).
c     - 6 or 10 columns of data will be printed, depending on n132.
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-
c
      Integer nx1,nx2,ny1,ny2,nz1,nz2,i1, i2, j1, j2, k1, k2, lu,
     >  i,j,k,n132,ncols,i3,i4,nformat
      Real a(nx1:nx2,ny1:ny2,nz1:nz2)
      Character*(*) file,string
c
c----------------------------------------------------------------------=
c
      If (file.EQ.'-' .OR. file.EQ.' ') Then
        Print *, '% Pslab: printing to screen'
        lu        = 6
      Else
        lu        = 11
        Print *, '% Pslab: ', file
        Open (Unit=lu, File=file, Status='Unknown')
      End If
      Write(lu,*) string
      Write(lu,9100) 'Array dimensions: (', nx1,nx2,ny1,ny2,nz1,nz2,file
      Write(lu,9100) 'Printed region:   (', i1, i2, j1, j2, k1, k2
 9100 Format (1x,a,2(i3,':',i3,','),i3,':',i3,')', :,'   File: ',a)
      Write(lu,*)
c
      nformat        = 10        !width of format field
      If (n132.LT.0) nformat = 12
      ncols        = 76/nformat
      If (Abs(n132).EQ.132) Then
      nformat        = 12
      ncols        = 128/nformat
      End If
c
c
c     XZ slabs of 3-d array
c
      If (k1.NE.k2) Then
      Do i3=i1,i2,ncols
      i4        = Min(i2,i3+ncols-1)
      Do j=j1,j2
      If (n132.GE.0) Then
        Write(lu,1000) 'j=', j, (i,i=i3,i4)
      Else
        Write(lu,1012) 'j=', j, (i,i=i3,i4)
      End If
      Do k=k2,k1,-1
      If (n132.GE.0) Then
        Write(lu,3000) k, (a(i,j,k),i=i3,i4)
      Else
        Write(lu,4000) k, (a(i,j,k),i=i3,i4)
      End If
      End Do
      Write(lu,*)
      End Do
      End Do
c
c
c     XY slab of 2-d array
c
      Else
      Do i3=i1,i2,ncols
      i4        = Min(i2,i3+ncols-1)
      Do k=k1,k2
      If (n132.GE.0) Then
        Write(lu,2000) (i,i=i3,i4)
      Else
        Write(lu,2012) (i,i=i3,i4)
      End If
      Do j=j2,j1,-1
      If (n132.GE.0) Then 
        Write(lu,3000) j, (a(i,j,k),i=i3,i4)
      Else
        Write(lu,4000) j, (a(i,j,k),i=i3,i4)
      End If
      End Do
      Write(lu,*)
      End Do
      End Do
      End If
c
      If (lu.NE.6) Close (lu)
 1000 Format (a,i3,i8,20i10)                !top row of integers: 3-d, 10x format
 1012 Format (a,i3,i8,20i12)                !top row of integers: 3-d, 12x format
 2000 Format (2x,20i10)                        !top row of integers: 2-d, 10x format
 2012 Format (2x,20i12)                        !top row of integers: 2-d, 12x format
 3000 Format (i3,2x,20f10.5)                !f10.5
 4000 Format (i3,2x,1p,20e12.4)                !e12.4
c
      Return
      End
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########         PSLABF         ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c        Print a 2-d array, or XZ slabs of a 3-d array, to an output file
c        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      SUBROUTINE Pslabf (a,file,string,nwid,format,
     >   nx1,nx2, ny1,ny2, nz1,nz2, i1,i2, j1,j2, k1,k2)
      IMPLICIT NONE
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-
c
c  INPUT:
c     a                = Array
c     file        = Output filename ('-' for output to screen)
c     string        = String written to output file
c     nwid        = max number of cols for output, usu 80 or 132.
c     format        = format string
c     nx1...nz2        = Upper & lower dimensions of A
c     i1...k2        = Portion of the array to print
c
c  OUTPUT: FORMATted output is sent to file.
c
c  NOTES:
c     - To print a 2-d array, set k1=k2; a single XY slab will be printed.
c     - The parameter n132 controls both the width of the output (80/132)
c        as well as the format (f10.5 if n132>=0, 1p,e12.4 otherwise).
c     - 6 or 10 columns of data will be printed, depending on n132.
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%-
c
      INTEGER nx1,nx2,ny1,ny2,nz1,nz2,i1, i2, j1, j2, k1, k2, lu,
     >  i,j,k,n132,ncols,j3,j4,i3,i4,nformat, Lnblnk, nwid, jf1,jf2
      REAL a(nx1:nx2,ny1:ny2,nz1:nz2)
      CHARACTER*(*) file,string,format, fmt*40, hfmt*40
c
c----------------------------------------------------------------------=
c
c
      IF (file.EQ.'-' .OR. file.EQ.' ') THEN
c       PRINT *, '% Pslabf: printing to screen'
        lu        = 6
      ELSE
        lu        = 11
        PRINT *, '% Pslabf: ', file
        OPEN (UNIT=lu, FILE=file, STATUS='UNKNOWN')
      END IF
      WRITE(lu,*) string
      WRITE(lu,9100) 'Array dimensions: (', nx1,nx2,ny1,ny2,nz1,nz2,file
      WRITE(lu,9100) 'Printed region:   (', i1, i2, j1, j2, k1, k2
 9100 FORMAT (1x,a,2(i3,':',i3,','),i3,':',i3,')', :,'   File: ',a)
      WRITE(lu,*)
c
c
c  The format MUST have the form  (F8.2), (E8.2), or (G8.2).  IF a G format,
c  then the 1P scale factor will be used.
c  
c
      IF (nwid.LE.0) nwid = 80
c
c     format        = '(F8.2)'
      jf1        = 3        !Index(format,',')
      jf2        = Index(format,'.') - 1
      jf2        = Min(jf2,jf1+1)
      IF (jf2-jf1.LT.0) THEN
        PRINT *, '% Pslabf: unable to parse format'
        PRINT *, '% Pslabf: ', jf1,jf2, format
      END IF
!     print *, 'reading from <', format(jf1:jf2), '> ', jf1, jf2
      IF (jf2.EQ.jf1) THEN
        READ (format(jf1:jf2), Fmt='(I1)') nformat
      ELSE
        READ (format(jf1:jf2), Fmt='(I2)') nformat
      END IF
      IF (nformat.LT.5 .AND. nformat.GT.18) THEN
        PRINT *, '% Pslabf: unable to parse format'
        PRINT *, '% Pslabf: ', nformat, format
        STOP
      END IF
c
      ncols        = (nwid-4)/nformat                ! number of cols of data
      fmt        = '(I3,2X,40' // format(:Lnblnk(format)) // ')'
      IF (Index(format,'G').GT.0)
     >  fmt        = '(I3,2X,1P,40' // format(:Lnblnk(format)) // ')'
      hfmt        = '(A,I3,40I' // format(jf1:jf2) // ')'
      IF (k1.EQ.k2)
     >hfmt        = '(2X,40I' // format(jf1:jf2) // ')'
      PRINT *, '% Pslabf: Fmt=  ', fmt
      PRINT *, '% Pslabf: Fmt=  ', hfmt
c
c
c     XZ slabs of 3-d array
c
      IF (k1.NE.k2) THEN
c
c...XZ slabs of 3-d array
c
      IF (j1.EQ.j2) THEN
      DO i3=i1,i2,ncols
      i4        = Min(i2,i3+ncols-1)
      DO j=j1,j2
        WRITE(lu,Fmt=hfmt) 'j=', j, (i,i=i3,i4)
      DO k=k2,k1,-1
        WRITE(lu,Fmt=fmt) k, (a(i,j,k),i=i3,i4)
      END DO
      WRITE(lu,*)
      END DO
      END DO
c
c
c...YZ slabs of 3-d array
c
      ELSE
      DO j3=j1,j2,ncols
      j4        = Min(j2,j3+ncols-1)
      DO i=i1,i2
        WRITE(lu,Fmt=hfmt) 'i=', i, (j,j=j3,j4)
      DO k=k2,k1,-1
        WRITE(lu,Fmt=fmt) k, (a(i,j,k),j=j3,j4)
      END DO
      WRITE(lu,*)
      END DO
      END DO
      END IF
c
c
c     XY slab of 2-d array
c
      ELSE
      DO i3=i1,i2,ncols
      i4        = Min(i2,i3+ncols-1)
      DO k=k1,k2
        WRITE(lu,Fmt=hfmt) (i,i=i3,i4)
      DO j=j2,j1,-1
        WRITE(lu,Fmt=fmt) j, (a(i,j,k),i=i3,i4)
      END DO
      WRITE(lu,*)
      END DO
      END DO
      END IF
c
      IF (lu.NE.6) CLOSE (lu)
c
c
      END
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########        PARSENUM        ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Subroutine ParseNum (string, numtraj, ntraj)
      Implicit None
      Character*(*) string
      Character*8 format
      Integer        ntraj,numtraj(ntraj)
      Integer        lnb,jj,kk,itraj,itraj1,itraj2,Lnblnk, ipunct,jpunct
c
c...Setup
c
      ntraj        = 0
      lnb        = Lnblnk(string)
      jj        = 1
c
c...Check for comma(44), dash(45)
c
 1000 Continue
c
      Call Punct (string,jj,jpunct,ipunct)
      kk        = jpunct - 1
c     print *, 'jj,kk,lnb:', jj,kk,lnb
c
      If (lnb.LT.kk) Return
c
c...Read single traj (n,) or first element of list (n-)
c
      Write (format,'(a,i1,a)') '(I', kk-jj+1, ')'
      Read (string(jj:kk),Fmt=format) itraj1
cprint *, itraj1
      ntraj        = ntraj + 1
      numtraj(ntraj) = itraj1
c
c
c...Parse to dash
c
      If (ipunct.EQ.45) Then
c     print *, 'dash'
        jj        = kk + 2
        Call Punct (string,jj,jpunct,ipunct)
        kk        = jpunct - 1
c     print *, 'jj,kk,lnb:', jj,kk,lnb
        Write (format,'(a,i1,a)') '(I', kk-jj+1, ')'
        Read (string(jj:kk),Fmt=format) itraj2
        Do itraj=itraj1+1,itraj2
c  print *, itraj
          ntraj        = ntraj + 1
          numtraj(ntraj) = itraj
        End Do
      End If
c
c
      jj        = kk + 2
      If (jj.LE.lnb) Go To 1000
c
      End
c
c
c
c
c                   ########################################
c                   ########################################
c                   ########################################
c                   ########                        ########
c                   ########          PUNCT         ########
c                   ########                        ########
c                   ########################################
c                   ########################################
c                   ########################################
c
c
      Subroutine Punct (string,j1,jpunct,ipunct)
      Implicit None
      Character*(*) string
      Character*1 TAB
      Integer        j1,jpunct,ipunct,jblank,jcomma,jdash,jdot,jtab,Lnblnk
c
c...Find the location of the next punctuation
c   This doesn't account for the last character being non-punctuation
c   (non-blank).
c
c  string        = string to be parsed (input)
c  j1                = index of first character to be scanned (input)
c  j2                = index of first punctuation (output)
c  ipunct        = ASCII code of punctuation (output)
c
c
      TAB        = Char(9)
      jpunct        = Len(string) + 1
c
      jblank        = Index (string(j1:),' ') + j1 - 1
      jcomma        = Index (string(j1:),',') + j1 - 1
      jdash        = Index (string(j1:),'-') + j1 - 1
      jdot        = Index (string(j1:),'.') + j1 - 1
      jtab        = Index (string(j1:),TAB) + j1 - 1
c
      If (jblank.LT.j1) jblank = 9999
      If (jcomma.LT.j1) jcomma = 9999
      If (jdash.LT.j1)        jdash = 9999
      If (jdot.LT.j1)        jdot = 9999
      If (jtab.LT.j1)        jtab = 9999
c
      jpunct        = Min(jblank,jpunct)
      jpunct        = Min(jcomma,jpunct)
      jpunct        = Min(jdash,jpunct)
      jpunct        = Min(jdot,jpunct)
      jpunct        = Min(jtab,jpunct)
c
      ipunct        = Ichar (string(jpunct:jpunct))
c
c     Print *, 'jpunct,ipunct: ', jpunct,ipunct
c
      End
