/*
 *	find_key...		Read a circular buffer key from a data file
 *				created by the program "make_wdss_key".  The
 *				environmental variable A2IO_KEYFILE points to
 *				this file location.
 *
 *				The default action if A2IO_KEYFILE isn't set is
 *				to do everything that used to occur in the past.
 *
 *				The code assumes the radar name has been set.
 *
 *	User calls:
 *		find_key()	Return the realtime key value.
 */

#pragma ident "@(#)find_key.c	4.1	01/31/01 NSSL"

#include <stdio.h>
#include <stdlib.h>

/*
 *	Circular buffer types, included here so users don't need to include
 *	a2io.h.
 */

#define	CIRC_RAW		0x88d2
#define	CIRC_CREMS_DEALIAS	0x88e0

#define	RT_OFFSET		12
static FILE *fp;

static int eof;

find_key( type, my_radar )
char *type, *my_radar;
{
	char key_buf[80];
	char key_radar[4];
	char *myptr;
	static int first_time = 1;
	int key1, key2, key3;

	if ( first_time )
	{
		myptr = getenv( "A2IO_KEYFILE" );
		if ( myptr == NULL )
			fp = NULL;
		else
		{
			printf("A2IO_KEYFILE:	%s\n",myptr);
			fp = fopen( myptr, "r" );
			if ( fp == NULL )
				perror( myptr );
		}
		first_time = 0;
	}

	else
	{
		if ( fp != NULL )
			rewind( fp );
	}

	eof = 0;

/*
 *	Default behavior.
 */

	if ( fp == NULL )
		return( handle_keys( type, CIRC_RAW, CIRC_CREMS_DEALIAS) );

	while( key_readline( key_buf ) != -1 )
	{
		sscanf( key_buf, "%s %x %x %d",
			key_radar, &key1, &key2, &key3 );
		if ( strncmp( key_radar, my_radar, 4 ) == 0 )
			return ( handle_keys( type, key1, key2 ) );

	}

	error_msg( my_radar, "not found in key file" );

	return ( handle_keys( type, CIRC_RAW, CIRC_CREMS_DEALIAS) );
}

handle_keys( type, raw_key, edit_key )
char *type;
int raw_key, edit_key;
{
	if ( type[0] == 'R' )
		return( raw_key );
	else if ( type[0] == 'E' )
		return( edit_key );
	else
	{	
		error_msg( type, "bad circular buffer type");
		return( raw_key );
	}
}

error_msg( s, t )
char *s, *t;
{
	printf("\n%s:  %s.  Default (which may be WRONG!) is returned.\n",
		s, t );
	return;
}

key_readline( key_buf )
char *key_buf;
{

/*
 *	No buffer overflow problems, because the code that writes this file
 *	doesn't have a problem.
 */

	if ( eof )
		return( -1 );


	for(;;)
	{
		if ( fgets( key_buf, 80, fp ) == NULL )
		{
			eof = 1;
			return( -1 );
		}

/*
 *	Check for comment line.
 */

		if ( key_buf[0] != '#' )
			break;
	}

	return( 0 );
}
