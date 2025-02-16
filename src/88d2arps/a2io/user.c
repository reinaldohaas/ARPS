/*
 *	User routines.
 *
 *	Internal routines
 *	parse_user_settings()		User settings for disk/files/tape.
 *	parse_user_settings_realtime()	User settings for realtime.
 *	update_a2_start_scan()		Part of tape changing process.
 */
 
#pragma ident "@(#)user.c	5.6	03/12/03 CAPS"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "extern.h"

/*
 *	Only retrieve the information for playback.
 */

parse_user_settings()
{
	char *myptr;
	int flag = 0;
	int error;

	error = 0;

/*
 *	Number of user tapes.
 */

	myptr = getenv("A2_NUM_TAPES" );

	if ( myptr != NULL )
	{
		a2_num_tapes = atoi( myptr );
		printf("A2_NUM_TAPES:  %d",a2_num_tapes);
		if ( a2_start_scan < 0 )
		{
			printf("	BAD VALUE!");
			error++;
		}
		printf("\n");
		flag++;
	}

/*
 *	Starting volume scan.
 */

	myptr = getenv("A2_START_SCAN");

	if ( myptr != NULL )
	{
		a2_start_scan = atoi( myptr );
		printf("A2_START_SCAN:	%d",a2_start_scan);
		if ( a2_start_scan < 0 )
		{
			printf("	BAD VALUE!");
			error++;
		}
		printf("\n");
		flag++;
	}

/*
 *	Ending volume scan.
 */

	myptr = getenv("A2_END_SCAN");

	if ( myptr != NULL )
	{
		a2_end_scan = atoi( myptr );
		printf("A2_END_SCAN:	%d",a2_end_scan);
		if ( a2_end_scan < 1 )
		{
			printf("	BAD VALUE!");
			error++;
		}
		printf("\n");
		flag++;
	}

	if ( a2_end_scan < a2_start_scan )
	{
		printf("ERROR:  Ending scan (%d) is less than starting scan (%d).\n",
			a2_end_scan, a2_start_scan );
		error++;
	}

/*
 *	Number of volume scans to process.
 */

	myptr = getenv("A2_NUM_SCAN");

	if ( myptr != NULL )
	{
		a2_num_scan = atoi( myptr );
		printf("A2_NUM_SCAN:	%d",a2_num_scan);
		if ( a2_num_scan < 1 )
		{
			printf("	BAD VALUE!");
			error++;
		}
		printf("\n");
		flag++;
	}

/*
 *	User thresholds.
 */

	myptr = getenv("A2_SPW_THRESHOLD");

	if ( myptr != NULL )
	{
		allow_thresholds = 1;
		spw_threshold = atoi( myptr );
		printf("A2_SPW_THRESHOLD:	%d",spw_threshold );
		if ( spw_threshold < 1 )
		{
			printf("	BAD VALUE!");
			error++;
		}
		printf("\n");
	}

/*
 *	Override the VCP tables for the fixed elevation angle.
 *
 *	This is a request from Steve Vasiloff @ SLC, where the first elevation
 *	angle is below the standard one used.
 */

	myptr = getenv("A2_NO_FIXED_ANGLE");

	if ( myptr != NULL )
	{
		a2_no_fixed_angle = atoi( myptr );
		printf("A2_NO_FIXED_ANGLE:	%d",a2_no_fixed_angle );
		if ( a2_no_fixed_angle < 1 )
		{
			printf("	BAD VALUE!");
			error++;
		}
		printf("\n");
	}

/*
 *	Allow the user to rotate the data if it is TDWR format.
 *	Input is in degrees, and can be a integer or real.  Data is
 *	transformed in the degrees x 100, which the rest of the software
 *	users.
 */

	myptr = getenv("TDWR_ROTATE" );

	if ( myptr != NULL )
	{
		tdwr_rotate = (int) ( atof( myptr ) * 100 + 0.5 );
/*
 *	Believe whatever the user gives us.  Negative values are ok here.
 */
		printf("TDWR_ROTATE [user value x 100]:	%d\n", tdwr_rotate);
	}

/*
 *	Workaround for Australia tape_II bug.
 */

	myptr = getenv("AUS_VS_BUG" );

	if ( myptr != NULL )
	{
		aus_vs_bug = atoi( myptr );
		printf("AUS_VS_BUG:	%d",aus_vs_bug );
		if ( aus_vs_bug < 1 )
		{
			printf("	BAD VALUE!");
			error++;
		}
		printf("\n");
	}

	if ( error )
	{
		printf("FATAL ERROR:  ");
		printf("One or more of the above values are bad.\n");
		exit( 1 );
	}

	if ( flag )
		printf("\n");

	return;
}

/*
 *	User controllable settings for realtime.
 *
 *	Ignore error settings.
 */

parse_user_settings_realtime()
{
	char *myptr;

/*
 *	Override the VCP tables for the fixed elevation angle.
 *
 *	This is a request from Steve Vasiloff @ SLC, where the first elevation
 *	angle is below the standard one used.
 */

	myptr = getenv("A2_NO_FIXED_ANGLE");

	if ( myptr != NULL )
	{
		a2_no_fixed_angle = atoi( myptr );
		printf("A2_NO_FIXED_ANGLE:	%d",a2_no_fixed_angle );
		if ( a2_no_fixed_angle < 1 )
		{
			printf("	BAD VALUE!  Setting ignored.");
		}
		printf("\n");
	}

/*
 *	Allow the user to rotate the data if it is TDWR format.
 *	Input is in degrees, and can be a integer or real.  Data is
 *	transformed in the degrees x 100, which the rest of the software
 *	users.
 */

	myptr = getenv("TDWR_ROTATE" );

	if ( myptr != NULL )
	{
		tdwr_rotate = (int) ( atof( myptr ) * 100 + 0.5 );
/*
 *	Believe whatever the user gives us.  Negative values are ok here.
 */
		printf("TDWR_ROTATE [user value x 100]:	%d\n", tdwr_rotate);
	}

	return;
}

/*
 *	Update the "a2_start_scan" set of variables.
 */

update_a2_start_scan( n )
int n;
{
	char tmp_ptr[20];
	char *myptr;
	int flag = 0;
	int error;

	if ( n == 1 )
	{
		printf("update_a2_start_scan:  logic fault, n is 1\n");
		return;
	}

	a2_start_scan = 0;
	a2_end_scan = 9999;
	a2_num_scan = 999999;

	error = 0;

/*
 *	Starting volume scan.
 */

	sprintf( tmp_ptr, "A2_START_SCAN_%d", n );
	myptr = getenv( tmp_ptr );

	if ( myptr != NULL )
	{
		a2_start_scan = atoi( myptr );
		printf("%s:	%d",tmp_ptr,a2_start_scan);
		if ( a2_start_scan < 0 )
		{
			printf("	BAD VALUE!");
			error++;
		}
		printf("\n");
		flag++;
	}

/*
 *	Ending volume scan.
 */

	sprintf( tmp_ptr, "A2_END_SCAN_%d", n );
	myptr = getenv( tmp_ptr );

	if ( myptr != NULL )
	{
		a2_end_scan = atoi( myptr );
		printf("%s:	%d",tmp_ptr,a2_end_scan);
		if ( a2_end_scan < 1 )
		{
			printf("	BAD VALUE!");
			error++;
		}
		printf("\n");
		flag++;
	}

	if ( a2_end_scan < a2_start_scan )
	{
		printf("ERROR:  Ending scan (%d) is less than starting scan (%d).\n",
			a2_end_scan, a2_start_scan );
		error++;
	}

/*
 *	Number of volume scans to process.
 */

	sprintf( tmp_ptr, "A2_NUM_SCAN_%d", n );
	myptr = getenv( tmp_ptr );

	if ( myptr != NULL )
	{
		a2_num_scan = atoi( myptr );
		printf("%s:	%d",tmp_ptr,a2_num_scan);
		if ( a2_num_scan < 1 )
		{
			printf("	BAD VALUE!");
			error++;
		}
		printf("\n");
		flag++;
	}

	if ( error )
	{
		printf("FATAL ERROR:  ");
		printf("One or more of the above values are bad.\n");
		exit( 1 );
	}

	if ( flag )
		printf("\n");

	return;
}
