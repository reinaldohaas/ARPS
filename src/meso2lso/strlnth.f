c
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######                SUBROUTINE STRLNTH                    ######
c     ######                                                      ######
c     ######                Copyright (c) 1992-1994               ######
c     ######    Center for Analysis and Prediction of Storms      ######
c     ######    University of Oklahoma.  All rights reserved.     ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
      SUBROUTINE STRLNTH( string, length )
c
c#######################################################################
c      
c     PURPOSE:
c
c     Return the length of the non-blank part of a character string.
c
c#######################################################################
c
c     AUTHOR: Ming Xue
c     11/20/91
c
c     MODIFICATION HISTORY:
c
c     5/05/92 (M. Xue)
c     Added full documentation.
c
c     9/10/94 (Weygandt & Y. Lu)
c     Cleaned up documentation.
c
c#######################################################################
c
c     INPUT:
c
c       string   A character string
c       length   The declared length of the character string 'string'.
c
c     OUTPUT:
c
c       length   The length of the non-blank part of the string.
c
c#######################################################################
c
c
c#######################################################################
c
c     Variable Declarations.
c
c#######################################################################
c
      implicit none

      character string*(*)      ! A character string for the name of
                                ! this run.
      integer length            ! The length of the non-blank part
                                ! of a string.

      integer i
c
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C
C     Beginning of executable code...
C
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
c

      DO 100 i = length,1,-1

        IF(string(i:i) .ne. ' ') GOTO 200

100   CONTINUE

200   CONTINUE

      length = max(1,i) 

      RETURN
      END
