/*
 *	a2io...		Read 88D level 2 (Archive 2) data.
 *
 *	radar_open()	Open radar connection.
 *	radar_init()	Backwards compatability call.
 *	get_rt_mode()	Let the user know whether the data is live/realtime.
 *
 *	Internal routines
 *	bad_mode()	Fatal error message routine.
 */

#pragma ident "@(#)a2io.c	5.5	03/10/04 CAPS"

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include "a2io.h"
#include "extern.h"

#include "proto.h"

radar_open( s )
char *s;
{
	static int have_been_called = 0;
	char my_radar[5];
	int rc;

/*
 *	Figure out what our data source.  After we do that, call our
 *	initialization routine.
 */

/*
 *	Initialize some things first, in the event the user ignores a fail
 *	return value (recently added), and tries to read data anyways.
 */

	read_record_init();

/*
 *	Configure A2IO functions.  They will get changed later if the format
 *	isn't A2IO data.
 */

	funct_init_a2io();

/*
 *	Everything but disk support is removed from this version.
 */

	input_format = DISK;
	rc = radar_open_raw_disk( s );
	if ( rc < 0 )
		return( rc );
	set_console_messages( 0 );
	parse_user_settings();

/*
 *	Multiple calls to "radar_open()" aren't supported, but are legal,
 *	assuming other functions above trap and handle the duplicate event.
 *
 *	Eventually, they will be legal as we may have multiple data streams.
 *	Until this occurs, none of the initialization routines should be called
 *	more than once.
 */

	if ( have_been_called )
		return;

	have_been_called = 1;

/*
 *	Initialization routines.
 */

	get_record_init();

/*
 *	Retrieve radar info.
 */

	get_radar_info();

	return( 0 );
}

radar_init( s )
char *s;
{
	int rc;

	rc = radar_open( s );
	return( rc );
}

/*
 *	Obsolete.
 */

get_rt_mode()
{
	return( 0 );
}

/*
 *	Set the China WSR98-D radar format variable.
 */

set_radar98( n )
int n;
{

	if (n == 0)
		radar98 = 0;
	else
		radar98 = 1;
	  printf(" Variable 'radar98' is now %d\n",radar98);
	return;
}

/*
 *	Get the China WSR98-D radar format variable.
 */

get_radar98()
{
	return(radar98);
}

/*
 *	Nssl_radar_mode is checked in various places.  If it is bad, we
 *	complain and exit.
 */


bad_mode( s )
char *s;
{
	printf("%s:  input_format is %d\n",s,input_format);
	printf("%s:  either 'radar_init' wasn't called\n",s);
	printf("%s:  or internal data has been corrupted.\n",s);
	printf("%s:  fatal error.\n", s);
	exit( 1 );

/* NOTREACHED */
}
