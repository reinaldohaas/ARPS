/*
 *	read_record_raw.c	Perform i/o.  The i/o may contain more than
 *				one data data record.
 *
 *	Internal
 *	read_record_raw()	Read one data record.
 */
 
#pragma ident "@(#)read_record_raw.c	4.1	01/29/01 CAPS"

#include <stdio.h>
#include <time.h>
#include "a2io.h"
#include "extern.h"
#include "proto.h"

/*
 *	Call the appropriate i/o routine.
 */

read_record_raw( mybuf )
char *mybuf;
{
	int rc = 0;

	if ( input_format == DISK )
		rc = read_record_raw_disk( mybuf, MAX_RECORD_SIZE_TAPE );

	else
		bad_mode( "read_record_raw" );

	return( rc );
}
