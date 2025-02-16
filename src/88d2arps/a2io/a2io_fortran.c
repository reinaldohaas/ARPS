/*
 *	A2io_fortran.c
 *
 *	This is a Fortran interface version.
 *
 *	Do *NOT* include the "len" argument in any Fortran calls!
 *	It is supplied automatically.
 *
 */
 
#pragma ident	"@(#)a2io_fortran.c	5.3	03/04/05	CAPS"

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>

radar_open_f( s, irc, len ) 
char *s;
int *irc;
int len;
{

/*
 *	We need to extract only the valid part of "s".
 */

	char *t;
	char *ptr;
	int i;

	t = (char *) malloc( len );

	if ( t == (char *) NULL ) {
		perror("radar_open_f:  malloc");
		exit(1);
	}

	strncpy(t, s, len );

/*
 *	Change spaces to ascii nulls
 */

	ptr = t;

	for(i=0; i<len; i++)
	{
		if ( *ptr == ' ' ) *ptr = (char) NULL;
		ptr++;
	}

	*irc = radar_open( t );

	free( t );

	return;

}

radar_open_f_( s, irc, len ) 
char *s;
int *irc;
int len;
{
	radar_open_f( s, irc, len );
	return;
}

/*
 *	Note that we call a Fortran interface here so that the blank spaces
 *	are cleaned up.
 */

radar_init_f( s, irc, len )
char *s;
int *irc;
int len;
{
	radar_open_f( s, irc, len );
	return;
}

radar_init_f_( s, irc, len )
char *s;
int *irc;
int len;
{
	radar_init_f( s, irc, len );
	return;
}

set_radar98_f( n )
int *n;
{	
	set_radar98(*n);
	return;
}

set_radar98_f_( n )
int *n;
{
	set_radar98_f(n);
	return;
}

get_radar98_f( n )
int *n;
{	
	*n = get_radar98();
	return;
}

get_radar98_f_( n )
int *n;
{
	get_radar98_f(n);
	return;
}
