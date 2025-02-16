/*
 *	get_record.c...		Handle data conversion for the user.
 *
 *	There are two formats possible, "tape" and "realtime" ("live").  Each
 *	has different valid data types.  Each has a different length for the
 *	data record.
 *
 *	The default it to use the same format.
 *
 *	get_record()			Return a record.
 *	set_record_format()		Set the record format.
 *
 *	Internal routines...
 *	get_record_init()		Init routine.
 *	get_record_format()		Returns "output_format".
 *	get_record_tape_format()	Convert live format to tape.
 */
 
#ifndef	LINT
static char Get_Record[] = "@(#)get_record.c	3.14	09/14/00";
#endif	/* LINT */

#include <stdio.h>
#include <string.h>
#include <time.h>
#include "a2io.h"
#include "extern.h"
#include "proto.h"

get_record_init()
{
	if ( input_format == TAPE )
		output_format =		TAPE;
	else if ( input_format == DISK )
		output_format =		TAPE;

	else
		bad_mode( "get_record_init" );

	return;
}

get_record( mybuf )
char *mybuf;
{
	int rc = 0;

/*
 *	Sanity check, so bcopy() doesn't core dump.  Note that zero isn't a
 *	currently used record length.
 */

	if ( message_length <= 0 )
		return( 0 );

	if ( output_format == TAPE )
		rc = get_record_tape_format( mybuf );
	else
		bad_mode( "get_record" );

	return( rc );

}

/*
 *	Set the format for get_record() calls.  The default is the same
 *	format as the data.
 */

set_record_format( n )
int n;
{
	if ( n == TAPE )
	{
		output_format = n;
		return(0);
	}

	printf("set_record_format:  bad format %d\n",n);
	return(-1);
}

/*
 *	Return the value.
 */

get_record_format()
{
	return( output_format );
}

/*
 *	User wants TAPE format.
 */

get_record_tape_format( mybuf )
char *mybuf;
{
	struct a2tape24 fake_header;

	int len;

/*
 *	Volume scan (201) record.
 */
	if ( message_type == A2_VOLSCAN_TYPE )
	{
		sprintf( fake_header.archive2, "ARCHIVE2.");
		sprintf( fake_header.vs, "%03d", nssl_scan_num);
/*
 *	Precip alg insists on using this as part of a 4-byte julian date.
 */
		fake_header.not_used = 0;
		fake_header.julian_date = nssl_a2header->julian_date;
		fake_header.milsec = nssl_a2header->milsec;
		bcopy( (char *) &fake_header, mybuf, 24 );
		return( 24 );
	}

/*
 *	Can't convert 202 records, as that has no equivalent tape format.
 */
	if ( message_type == A2_SITE_TYPE )
		return( A2_CANT_CONVERT );

/*
 *	Everything else.
 */
	else
	{
		len = REC_SIZE_TAPE;
/*
 *	Handle special types.
 */
		if ( message_type == A2_ARCHIVE2 ||
			message_type == A2_ARCHIVE2_NCDC ||
			message_type == A2_VS_24_BYTE )
				len = message_length;
		bzero( mybuf, len );
		bcopy( nssl_ptr, mybuf, len );
		return( len );
	}
}
