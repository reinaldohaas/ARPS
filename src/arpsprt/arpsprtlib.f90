!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE STRCAP                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE strcap(inline,outline,LEN)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Capitalize any letters found in input character string "inline"
!  and put result in character string "outline".
!
!  Make use of the fact that ASCII lower-case "a" is 32 greater than
!  upper-case "A" (and similarly for b-thru-z).
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster   OU School of Meteorology  April, 1992
!
!  MODIFICATION HISTORY:
!    8/02/92  Modified from arps2.5 to arps3.0
!
!-----------------------------------------------------------------------
!
!  INPUT :
!    inline     line to be capitalized
!    len        length of inline
!
!  OUTPUT :
!
!    outline    capitalized inline
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE

  INTEGER :: LEN
  CHARACTER (LEN=1) :: inline(LEN)
  CHARACTER (LEN=1) :: outline(LEN)
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  INTEGER :: letasc
  INTEGER :: k
!
!-----------------------------------------------------------------------
!
!  Fortran intrinsic functions
!
!-----------------------------------------------------------------------
!
!  LOGICAL :: LGE,LLE
!  INTEGER :: ICHAR
!  CHARACTER (LEN=1) :: CHAR
!
!-----------------------------------------------------------------------
!
!  Search line for valid characters
!
!-----------------------------------------------------------------------
!
  DO k=1,LEN
    IF(LGE(inline(k),'a') .AND. LLE(inline(k),'z')) THEN
      letasc=ICHAR(inline(k))
      outline(k)=CHAR(letasc-32)
    ELSE
      outline(k)=inline(k)
    END IF
  END DO
  RETURN
END SUBROUTINE strcap
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE PARSLN                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE parsln(inline,reqprm,LEN,maxpar,nparms)
!
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Parses user input line and pulls off character
!  strings in ASCII range A-Z or 1-9.
!
!  Input line should be capitalized BEFORE calling this routine.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster   OU School of Meteorology
!
!  MODIFICATION HISTORY:
!    8/02/92  Modified from arps2.5 to arps3.0
!
!-----------------------------------------------------------------------
!
!  INPUT :
!    inline     character string to be aprsed
!    len        length of inline
!    maxpar     maximum number of substrings to store
!
!  OUTPUT :
!    reqprm     substrings
!    nparms     number of substrings found
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  INTEGER :: LEN,maxpar
  CHARACTER (LEN=1) :: inline(LEN)
  CHARACTER (LEN=2) :: reqprm(maxpar)
  INTEGER :: nparms
!
!-----------------------------------------------------------------------
!
!  Misc. local variables
!
!-----------------------------------------------------------------------
!
  CHARACTER (LEN=2) :: star
  INTEGER :: k,l
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!-----------------------------------------------------------------------
!
!  Set all output parms to missing value, "*"
!  which has ASCII decimal number = 42
!
!-----------------------------------------------------------------------
!
  star='**'
  DO k=1,maxpar
    reqprm(k)=star
  END DO
!
!-----------------------------------------------------------------------
!
!  Search line for valid characters
!
!-----------------------------------------------------------------------
!
  k=1
  l=1
  20   CONTINUE
  IF((LGE(inline(l),'A') .AND. LLE(inline(l),'Z')) .AND.                &
        (LGE(inline(l+1),'A') .AND. LLE(inline(l+1),'Z')) ) THEN
    reqprm(k)(1:1)=inline(l)
    reqprm(k)(2:2)=inline(l+1)
    k=k+1
    IF(k > maxpar) GO TO 100
    l=l+2
  ELSE
    l=l+1
  END IF
  IF(l <= LEN) GO TO 20
  100  CONTINUE
!
!-----------------------------------------------------------------------
!
!  Record number of substrings saved in reqprm
!
!-----------------------------------------------------------------------
!
  nparms=k-1
  RETURN
END SUBROUTINE parsln
!
!##################################################################
!##################################################################
!######                                                      ######
!######                SUBROUTINE KNTARY                     ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

SUBROUTINE kntary(fmtver,                                               &
           grdout,basout,varout,mstout,iceout,trbout,                   &
           n3dary,n3dsclr)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Count the number of 3-D arrays and 3-D scalars that should be
!  written into a history file depending on the switches that
!  have been set and the version.
!
!  For now, we have one version so fmtver is not consulted.
!
!-----------------------------------------------------------------------
!
!  AUTHOR: Keith Brewster
!  7/14/92
!
!  MODIFICATION HISTORY:
!
!-----------------------------------------------------------------------
!
!  INPUT:
!
!  fmtver  Format/Version number a character string describing the
!          format and version number of the grid writing routine
!
!  grdout =0 or 1. If grdout=0, grid variables are not written.
!  basout =0 or 1. If basout=0, base state variables are not dumped.
!  varout =0 or 1. If varout=0, model perturbation variables are not written.
!  mstout =0 or 1. If mstout=0, water variables are not written.
!  iceout =0 or 1. If iceout=0, writing of qi, qs and qh is skipped.
!  trbout =0 or 1. If trbout=0, writing of turbulence param km is skipped.
!
!  OUTPUT:
!
!  n3dary     Total number of 3-D arrays written into history file
!  n3dsclr    Total number of 3-D scalars written into history file
!
!-----------------------------------------------------------------------
!
!  Variable Declarations.
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
!
  CHARACTER (LEN=40) :: fmtver
  INTEGER :: grdout,basout,varout,mstout,iceout,trbout
  INTEGER :: n3dary,n3dsclr
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  n3dary=0
  n3dsclr=0
  IF( grdout == 1 ) n3dary=n3dary+3
  IF( basout == 1 ) THEN
    n3dary=n3dary+7
    n3dsclr=n3dsclr+4
  END IF
  IF( varout == 1 ) THEN
    n3dary=n3dary+5
    n3dsclr=n3dsclr+2
  END IF
  IF( mstout == 1 ) THEN
    n3dary=n3dary+3
    n3dsclr=n3dsclr+3
    IF( iceout == 1 ) THEN
      n3dary=n3dary+3
      n3dsclr=n3dsclr+4
    END IF
  END IF
  IF( trbout == 1 ) THEN
    n3dary=n3dary+1
    n3dsclr=n3dsclr+1
  END IF

  RETURN
END SUBROUTINE kntary
