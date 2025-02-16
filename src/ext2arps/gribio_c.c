/*
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######               GRIBIO (C Functions)                   ######
c     ######                                                      ######
c     ######                     Developed by                     ######
c     ######     Center for Analysis and Prediction of Storms     ######
c     ######                University of Oklahoma                ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
c#######################################################################
c
c     PURPOSE:
c
c     C programs will be called by Fortran program, GRIBREAD and
c     GRIBDUMP.
c
c#######################################################################
c
c     AUTHOR: Yuhe Liu
c     11/01/1995
c
c     MODIFICATION HISTORY:
c     08/23/2004  Add GSIZE function.  Kevin W. Thomas
C
c#######################################################################
c
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/file.h>

#ifdef CRAY
#include <fortran.h>
#elif UNDERSCORE
#define GSIZE  gsize_
#define GOPEN  gopen_
#define GREAD  gread_
#define GWRITE gwrite_
#define GSEEK  gseek_
#define GCLOSE gclose_
#define CHKGRB chkgrb_
#else
#define GSIZE  gsize
#define GOPEN  gopen
#define GREAD  gread
#define GWRITE gwrite
#define GSEEK  gseek
#define GCLOSE gclose
#define CHKGRB chkgrb
#endif

/*
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######               SUBROUTINE GSIZE                       ######
c     ######                                                      ######
c     ######                     Developed by                     ######
c     ######     Center for Analysis and Prediction of Storms     ######
c     ######                University of Oklahoma                ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
*/

void GSIZE( int *size )

/*
c
c#######################################################################
c
c     PURPOSE:
c
c     Return the size of STDIO file pointer.  This will allow the caller
c     to allocate enough space for the "grbunit" variable used by the
c     GOPEN and other routines.
c
c#######################################################################
c
c     AUTHOR: Kevin W. Thomas
c     08/23/2004
c
c     MODIFICATION HISTORY:
C
c#######################################################################
c
c     ARGUMENTS:
c
c     size      Size of the stdio FILE structure.
c
c#######################################################################
c
*/

{
    
    *size = sizeof( FILE );
    return;
}

/*
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######               SUBROUTINE GOPEN                       ######
c     ######                                                      ######
c     ######                     Developed by                     ######
c     ######     Center for Analysis and Prediction of Storms     ######
c     ######                University of Oklahoma                ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
*/

void GOPEN( FILE **funit, char *fname, int *len, char *modes, int *ireturn )

/*
c
c#######################################################################
c
c     PURPOSE:
c
c     Open GRIB file.
c
c#######################################################################
c
c     AUTHOR: Yuhe Liu
c     11/01/1995
c
c     MODIFICATION HISTORY:
C
c#######################################################################
c
c     ARGUMENTS:
c
c     funit    I/O unit of GRIB file opened by gopen
c     fname    File name being opened
c     modes    Open modes
c     ireturn  Return flag
c             =  0, no error
c             = -1, error: invalid open mode
c             = -2, error: can't open the file
c
c#######################################################################
c
*/

{
    char mflag[4];

    int lennam = *len;
    int lenmod = strlen( modes );

    *funit = NULL;
    *ireturn = 0;

    if ( lennam )
        while( fname[--lennam] == ' ' );
    fname[++lennam] = '\0';

/*    if ( lenmod )
 *        while( modes[--lenmod] == ' ' );
 *    modes[++lenmod] = '\0'; */

    mflag[0] = '\0';

    while( *modes && ( strlen(mflag) < 3 ) )
    {
        switch( *modes )
        {
            case 'a':
            case 'A':
                strcat( mflag, "a" );
                break;
            case 'r':
            case 'R':
                strcat( mflag, "r" );
                break;
            case 'w':
            case 'W':
                strcat( mflag, "w" );
                break;
            case '+':
                strcat( mflag, "+" );
                break;
        }
        modes++;
    }

    if ( !strcmp( mflag, "r" ) && !strcmp( mflag, "w" ) )
        strcpy( mflag, "r+w" );

    *funit = fopen( fname, mflag ); 
    if ( *funit == NULL ) {
        perror( fname );
        *ireturn = -2;
        return;
    }

    return;
}

/*
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######               SUBROUTINE GREAD                       ######
c     ######                                                      ######
c     ######                     Developed by                     ######
c     ######     Center for Analysis and Prediction of Storms     ######
c     ######                University of Oklahoma                ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
*/

void GREAD( FILE **funit, void *message, int *nbytes, int *ireturn )

/*
c
c#######################################################################
c
c     PURPOSE:
c
c     Read a block of bytes, nbytes, starting from the current position 
c     in GRIB file.
c
c#######################################################################
c
c     AUTHOR: Yuhe Liu
c     11/01/1995
c
c     MODIFICATION HISTORY:
C
c#######################################################################
c
c     ARGUMENTS:
c
c     funit    I/O unit of GRIB file opened by gopen
c     message  Buffer of the message read in by this subroutine.
c     nbytes   Number of bytes of the message
c     ireturn  Return flag
c             =  0, no error
c             = -1, read error occurs
c
c#######################################################################
c
*/

{
    int iret;

    *ireturn = 0;

    iret = fread( message, 1, *nbytes, *funit );

    if ( iret != *nbytes ) {
        *ireturn = -1;
        perror( "gread" );
        return;
    }
    else *ireturn = iret;

    return;
}

/*
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######               SUBROUTINE GWRITE                      ######
c     ######                                                      ######
c     ######                     Developed by                     ######
c     ######     Center for Analysis and Prediction of Storms     ######
c     ######                University of Oklahoma                ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
*/

void GWRITE( FILE **funit, void *message, int *nbytes, int *ireturn )

/*
c
c#######################################################################
c
c     PURPOSE:
c
c     Write a block of bytes, nbytes, beginning at the position at
c     nbstrt into GRIB file.
c
c#######################################################################
c
c     AUTHOR: Yuhe Liu
c     11/01/1995
c
c     MODIFICATION HISTORY:
C
c#######################################################################
c
c     ARGUMENTS:
c
c     funit    I/O unit of GRIB file opened by gopen
c     message  Buffer of the message write in by this subroutine.
c     nbytes   Number of bytes of the message
c     ireturn  Return flag
c             =  0, no error
c             = -1, write error occurs
c
c#######################################################################
c
*/

{
    int iret;

    *ireturn = 0;

    iret = fwrite( message, 1, *nbytes, *funit );
    if ( iret != *nbytes ) {
        *ireturn = -1;
        perror( "gwrite" );
        return;
    }
    else
        *ireturn = iret;

    return;
}

/*
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######               SUBROUTINE GSEEK                       ######
c     ######                                                      ######
c     ######                     Developed by                     ######
c     ######     Center for Analysis and Prediction of Storms     ######
c     ######                University of Oklahoma                ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
*/

void GSEEK (FILE **funit, int *offset, int *whence, int *ireturn )

/*
c
c#######################################################################
c
c     PURPOSE:
c
c     Seek the position of the GRIB file pointer at offset
c
c#######################################################################
c
c     AUTHOR: Yuhe Liu
c     11/01/1995
c
c     MODIFICATION HISTORY:
C
c#######################################################################
c
c     ARGUMENTS:
c
c     funit    I/O unit of GRIB file opened by gopen
c     offset   The number of offset bytes to be set to GRIB file
c     whence   Flag to set the offset
c             = 0, pointer set to offset
c             = 1, pointer set to current location plus offset
c             = 2, pointer set to the size of the file
c     ireturn  Return flag
c             =  0, no error
c             = -1, error: Invalid seek operation
c             = -2, error in fseek
c
c#######################################################################
c
*/

{
    int iret;

    *ireturn = 0;

    if ( *whence < 0 || *whence > 2 ) {
        perror( "Invalid seek operation in gseek" );
        *ireturn = -1;
        return;
    }
    else if ( *whence < 2 )
        iret = fseek( *funit, *offset, *whence );
    else {
        iret    = fseek( *funit, 0, 2 );
        *offset = ftell( *funit );
    }

    if ( iret != 0 ) {
        perror( "fseek" );
        *ireturn = -2;
        return;
    }
    else *ireturn = iret;

    return;
}

/*
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######               SUBROUTINE GCLOSE                      ######
c     ######                                                      ######
c     ######                     Developed by                     ######
c     ######     Center for Analysis and Prediction of Storms     ######
c     ######                University of Oklahoma                ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
*/

void GCLOSE( FILE **funit, int *ireturn )

/*
c
c#######################################################################
c
c     PURPOSE:
c
c     Close GRIB file.
c
c#######################################################################
c
c     AUTHOR: Yuhe Liu
c     11/01/1995
c
c     MODIFICATION HISTORY:
C
c#######################################################################
c
c     ARGUMENTS:
c
c     funit    I/O unit of GRIB file opened by gopen
c     ireturn  Return flag
c             =  0, no error
c             = -1, error in file close
c
c#######################################################################
c
*/

{
    *ireturn = fclose( *funit );
}

/*
c
c     ##################################################################
c     ##################################################################
c     ######                                                      ######
c     ######               SUBROUTINE CHKBRB                      ######
c     ######                                                      ######
c     ######                     Developed by                     ######
c     ######     Center for Analysis and Prediction of Storms     ######
c     ######                University of Oklahoma                ######
c     ######                                                      ######
c     ##################################################################
c     ##################################################################
c
*/

void CHKGRB( char *fname, int *len, int *modes, int *ireturn )

/*
c
c#######################################################################
c
c PURPOSE:
c
c     Check whether a given GRIB file is indeed in GRIB format and 
C     determine whether it is in grib 1 or grib2.
c
c#######################################################################
c
c     AUTHOR: Yuheng Wang
c     10/20/2007
c
c     MODIFICATION HISTORY:
C
c#######################################################################
c
c     ARGUMENTS:
c
c     fname    File name being opened
C     len      Length of the file name string
c     modes    GRIB mode, 1 - GRIB; 2-GRIB2
c     ireturn  Return flag
c             =  0, no error
c             = -2, error: can't open the file
c
c#######################################################################
c
*/

{
    int lennam = *len;

    FILE *funit;
    char letter;
    char gmsg[5];
    int  i;

    *ireturn = 0;

    if ( lennam )
        while( fname[--lennam] == ' ' );
    fname[++lennam] = '\0';

    funit = fopen( fname, "r" ); 
    if ( funit == NULL ) {
        perror( fname );
        *ireturn = -2;
        return;
    }

    for (i=0;i<4;i++) {
      letter = fgetc( funit );
      if (letter != EOF ) {
	gmsg[i] = letter;
      } else {
	*ireturn = -1;
	perror( "ERROR: not a grib file.\n");
	return;
      }
    }
    gmsg[i] = '\0';

    /* printf("**** file=%s, %s ****\n",fname,gmsg); */
    if (strcmp(gmsg,"GRIB") != 0) {
      printf("ERROR: File %s is not a GRIB file.\n",fname);
      *ireturn = -1;
      return;
    }

    fseek( funit, 7, 0 );
    *modes = fgetc( funit );
        
    /* printf("**** version = %d ****\n",*modes);  */

    *ireturn = fclose( funit );

    return;
}
