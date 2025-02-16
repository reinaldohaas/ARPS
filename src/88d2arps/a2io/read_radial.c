/*
 *	read_radial()			Guarantee a radial is ready or that
 *					end of data has occurred.
 *
 *	User Routines
 *	read_radial()			Return a radial or an EOF.
 */
 
#pragma ident "@(#)read_radial.c	5.1	10/15/01 CAPS"

#include <stdio.h>

#include "config.h"

#include "a2defs.h"

#ifdef	ALLOW_TDWR
#include "tdwrdefs.h"
#endif

#include "proto.h"

#define	RR_DATA		0
#define	RR_EOD		1

/*
 *	Make the next radial of data to the user.  The i/o itself is done in
 *	one of the "read_record" routines.
 */

read_radial()
{
	int rc;

	for(;;)
	{
		rc = read_record();

		if ( rc == A2_DATA_TYPE )
			return( RR_DATA );

		if ( rc == A2_DIGITAL_TYPE )
			return( RR_DATA );

#ifdef	ALLOW_TDWR
		if ( rc == TDWR_2B00 )
			return( RR_DATA );

		if ( rc == TDWR_2B01 )
			return( RR_DATA );
#endif

		if ( rc == RADAR_END_OF_DATA )
			return( RR_EOD );
	}

/* NOTREACHED */

}
