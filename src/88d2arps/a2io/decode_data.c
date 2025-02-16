/*
 *	User routines
 *	None
 *
 *	Internal routines
 *	decode_data()		Decode the data, message type 1.
 *	decode_data31()		Decode the data, message type 31.
 */

#pragma ident "@(#)decode_data.c	5.2	10/17/05	CAPS"

#include <stdio.h>
#include <strings.h>
#include <time.h>
#include "a2io.h"
#include "a2io_message_31.h"
#include "const.h"
#include "extern.h"
#include "vcp.h"

static int old_scan_num = -1;
static int old_tilt_num = -1;

/*
 *
 *	Decode the 88D data.
 */

decode_data()
{
	static int old_vcp = -999;

/*
 *	If we aren't a radial data message, we are done.
 */

	if ( message_type != A2_DATA_TYPE ) return;

/*
 *	Make sure the number of gates info is right.  It *might* be wrong
 *	on some of the older Archive II tapes.
 *
 *	Do NOT use nssl_a2data->dbz/vel/spw_ptr!
 */

#ifdef	VERY_OLD_FLAWED_DATA
/*
 *	WARNING...basic level 2 data is modified here!
 */

	if (nssl_a2data->arc_dbz_ptr == 0)
		nssl_a2data->num_gates_dbz = 0;
	if (nssl_a2data->arc_vel_ptr == 0)
		nssl_a2data->num_gates_vel = 0;
#endif

/*
 *	Set the "status_scan" and "status_tilt" flags.
 */

	if ( old_scan_num != nssl_scan_num )
	{
		nssl_scan_status = 1;
		nssl_tilt_status = 1;
	}
	else if ( old_tilt_num != nssl_a2data->elevation_number )
	{
		nssl_scan_status = 0;
		nssl_tilt_status = 1;
	}
	else
	{
		nssl_scan_status = 0;
		nssl_tilt_status = 0;
	}

/*
 *	In the realtime data stream, sweep 2 has a radial status of 3 (new
 *	volume scan) which is wrong.  This causes some algorithms (at least
 *	the 88D PRECIP algorithm) to fail.  We will correct the value.
 */

/*
 *	WARNING...basic level 2 data is modified here!
 */

	if ( nssl_a2data->radial_status == 3 && 
			nssl_a2data->elevation_number != 1 )
		nssl_a2data->radial_status = 0;

/*
 *	Check and see if we have a data restart.  We do not set tilt_status.
 */

	nssl_restart_flag = 0;
	if ( nssl_a2data->radial_status == 0 && !nssl_tilt_status )
	{
		printf("read_record:  data restart detected\n");
		nssl_restart_flag++;
	}

/*
 *	Check for an error state, which seems to occur on some Build 8 tapes,
 *	perhaps only NCDC versions.  It also can occur in the realtime system
 *	in rare cases when XCDedit core dumps at the end of a volume scan.
 *
 *	The error is that new volume scans start without the 24-byte header
 *	record, and for tapes, without an eof.
 */

	if ( nssl_tilt_status && nssl_a2data->elevation_number == 1 )
	{
		if ( !nssl_scan_status )
			printf("read_record:  new volume scan without new header.\n");
	}

	old_scan_num = nssl_scan_num;
	old_tilt_num = nssl_a2data->elevation_number;

/*
 *	Find VCP info.
 */

	switch (nssl_a2data->vcp)
	{
/*
 *	If we're VCP 0, then we might be TDWR data that doesn't use it.
 */
		case 0:
			nssl_vcp = NULL;
			break;

		case 11:
			nssl_vcp = vcpat11;
			break;
		case 21:
			nssl_vcp = vcpat21;
			break;
/*
 *	VCP 22 is being used by Janelle Janish of OSF work.
 */
		case 22:
			nssl_vcp = NULL;
			break;
		case 31:
			nssl_vcp = vcpat31;
			break;
		case 32:
			nssl_vcp = vcpat32;
			break;
/*
 *	VCP 42 thru 49 are used in "research" data.
 *	Ditto VCP 51 thru 54.
 *	Also 55-59 and 61.
 */
		case 42:
		case 43:
		case 44:
		case 45:
		case 46:
		case 47:
		case 48:
		case 49:
		case 51:
		case 52:
		case 53:
		case 54:
		case 55:
		case 56:
		case 57:
		case 58:
		case 59:
		case 61:
/*
 *	VCP's 144-148 are assigned values for Australia.
 */
		case 144:
		case 145:
		case 146:
		case 147:
		case 148:
			nssl_vcp = NULL;
			break;

		case 300:
			nssl_vcp = vcpat300;
			break;
		default:
/*
 *	Only warn once, when the VCP number changes.  This is at the request
 *	of Greg Stumpf.
 */

/*
 *	Turn off this code, as new VCP's have been added.
 */
#ifdef	VCP_WARN
			if ( old_vcp != nssl_a2data->vcp )
				printf("read_record:  %d:  unknown vcp number\n",
					nssl_a2data->vcp);
#endif
			nssl_vcp = NULL;
			break;	
	}

	old_vcp = nssl_a2data->vcp;

/*
 *	Convert milli-seconds to seconds.
 */
 
	zulu_time = nssl_a2data->zulu_time;
	zulu_time /= 1000;

/*
 *	Add the number of days.
 */

	zulu_time += nssl_a2data->mod_julian_date * SECONDS_PER_DAY;

/*
 *	For UNIX, time=0 is December 31, 1969, while it is January 1, 1970
 *	for NEXRAD systems, so we need to subtract off one day.
 */

	zulu_time -= SECONDS_PER_DAY;

/*
 *	Set a pointer to the data record.
 */
	nssl_gmt = gmtime(&zulu_time);

/*
 *	Copy it to dedicated space so that we won't have problems if "a2data"
 *	changes and we want the old data.
 */
	bcopy( (char *) nssl_gmt, &nssl_gmt_data, sizeof(struct tm));
	nssl_gmt = &nssl_gmt_data;

/*
 *	Check for a bogus elevation angle.  This is an indication of a radar
 *	hardware problem.
 *
 *	Don't check if in VCP 54, as large elevation angles are valid.
 *	Don't check if in VCP 144, as large elevation angles are also valid.
 */

	if ( nssl_a2data->elevation > elev_thresh
		&& nssl_a2data->vcp != 0
		&& nssl_a2data->vcp != 54
		&& nssl_a2data->vcp != 57
		&& nssl_a2data->vcp != 144
		&& nssl_a2data->vcp != 145
		&& nssl_a2data->vcp != 146
		&& nssl_a2data->vcp != 147
		&& nssl_a2data->vcp != 148
			)
	{
		printf("read_record:  bogus elevation angle of %.2f degrees.\n",
			(float)get_elev()/100.);

/*
 *	WARNING...basic level 2 data is modified here!
 */

		nssl_a2data->elevation = elev_save;
		printf("read_record:  setting elevation angle to %.2f degrees.\n",
			(float)get_elev()/100.);
	}

	elev_save = nssl_a2data->elevation;
	
	return;

}

/*
 *
 *	Decode the 88D data.
 */

decode_data31()
{
	static int old_vcp = -999;

/*
 *	If we aren't a radial data message, we are done.
 */

	if ( message_type != A2_DIGITAL_TYPE ) return;

/*
 *	Make sure the number of gates info is right.  It *might* be wrong
 *	on some of the older Archive II tapes.
 *
 *	Do NOT use nssl_a2data->dbz/vel/spw_ptr!
 */

#ifdef	VERY_OLD_FLAWED_DATA
/*
 *	WARNING...basic level 2 data is modified here!
 */

	if (nssl_a2data->arc_dbz_ptr == 0)
		nssl_a2data->num_gates_dbz = 0;
	if (nssl_a2data->arc_vel_ptr == 0)
		nssl_a2data->num_gates_vel = 0;
#endif

/*
 *	Set the "status_scan" and "status_tilt" flags.
 */

	if ( old_scan_num != nssl_scan_num )
	{
		nssl_scan_status = 1;
		nssl_tilt_status = 1;
	}
	else if ( old_tilt_num != nssl_a2data31->elevation_number )
	{
		nssl_scan_status = 0;
		nssl_tilt_status = 1;
	}
	else
	{
		nssl_scan_status = 0;
		nssl_tilt_status = 0;
	}

/*
 *	In the realtime data stream, sweep 2 has a radial status of 3 (new
 *	volume scan) which is wrong.  This causes some algorithms (at least
 *	the 88D PRECIP algorithm) to fail.  We will correct the value.
 */

/*
 *	WARNING...basic level 2 data is modified here!
 */

	if ( nssl_a2data31->radial_status == 3 && 
			nssl_a2data31->elevation_number != 1 )
		nssl_a2data31->radial_status = 0;

/*
 *	Check and see if we have a data restart.  We do not set tilt_status.
 */

	nssl_restart_flag = 0;
	if ( nssl_a2data31->radial_status == 0 && !nssl_tilt_status )
	{
		printf("read_record:  data restart detected\n");
		nssl_restart_flag++;
	}

/*
 *	Check for an error state, which seems to occur on some Build 8 tapes,
 *	perhaps only NCDC versions.  It also can occur in the realtime system
 *	in rare cases when XCDedit core dumps at the end of a volume scan.
 *
 *	The error is that new volume scans start without the 24-byte header
 *	record, and for tapes, without an eof.
 */

	if ( nssl_tilt_status && nssl_a2data31->elevation_number == 1 )
	{
		if ( !nssl_scan_status )
			printf("read_record:  new volume scan without new header.\n");
	}

	old_scan_num = nssl_scan_num;
	old_tilt_num = nssl_a2data31->elevation_number;

#define	NOT_YET
#ifndef	NOT_YET
/*
 *	Find VCP info.
 */

	switch (nssl_a2data31->vcp)
	{
/*
 *	If we're VCP 0, then we might be TDWR data that doesn't use it.
 */
		case 0:
			nssl_vcp = NULL;
			break;

		case 11:
			nssl_vcp = vcpat11;
			break;
		case 21:
			nssl_vcp = vcpat21;
			break;
/*
 *	VCP 22 is being used by Janelle Janish of OSF work.
 */
		case 22:
			nssl_vcp = NULL;
			break;
		case 31:
			nssl_vcp = vcpat31;
			break;
		case 32:
			nssl_vcp = vcpat32;
			break;
/*
 *	VCP 42 thru 49 are used in "research" data.
 *	Ditto VCP 51 thru 54.
 *	Also 55-59 and 61.
 */
		case 42:
		case 43:
		case 44:
		case 45:
		case 46:
		case 47:
		case 48:
		case 49:
		case 51:
		case 52:
		case 53:
		case 54:
		case 55:
		case 56:
		case 57:
		case 58:
		case 59:
		case 61:
/*
 *	VCP's 144-148 are assigned values for Australia.
 */
		case 144:
		case 145:
		case 146:
		case 147:
		case 148:
			nssl_vcp = NULL;
			break;

		case 300:
			nssl_vcp = vcpat300;
			break;
		default:
/*
 *	Only warn once, when the VCP number changes.  This is at the request
 *	of Greg Stumpf.
 */

/*
 *	Turn off this code, as new VCP's have been added.
 */
#ifdef	VCP_WARN
			if ( old_vcp != nssl_a2data31->vcp )
				printf("read_record:  %d:  unknown vcp number\n",
					nssl_a2data31->vcp);
#endif
			nssl_vcp = NULL;
			break;	
	}

	old_vcp = nssl_a2data31->vcp;
#endif

/*
 *	Convert milli-seconds to seconds.
 */
 
	zulu_time = nssl_a2data31->zulu_time;
	zulu_time /= 1000;

/*
 *	Add the number of days.
 */

	zulu_time += nssl_a2data31->mod_julian_date * SECONDS_PER_DAY;

/*
 *	For UNIX, time=0 is December 31, 1969, while it is January 1, 1970
 *	for NEXRAD systems, so we need to subtract off one day.
 */

	zulu_time -= SECONDS_PER_DAY;

/*
 *	Set a pointer to the data record.
 */
	nssl_gmt = gmtime(&zulu_time);

/*
 *	Copy it to dedicated space so that we won't have problems if "a2data"
 *	changes and we want the old data.
 */
	bcopy( (char *) nssl_gmt, &nssl_gmt_data, sizeof(struct tm));
	nssl_gmt = &nssl_gmt_data;

#ifndef	NOT_YET
/*
 *	Check for a bogus elevation angle.  This is an indication of a radar
 *	hardware problem.
 *
 *	Don't check if in VCP 54, as large elevation angles are valid.
 *	Don't check if in VCP 144, as large elevation angles are also valid.
 */

	if ( nssl_a2data31->elevation > elev_thresh
		&& nssl_a2data31->vcp != 0
		&& nssl_a2data31->vcp != 54
		&& nssl_a2data31->vcp != 57
		&& nssl_a2data31->vcp != 144
		&& nssl_a2data31->vcp != 145
		&& nssl_a2data31->vcp != 146
		&& nssl_a2data31->vcp != 147
		&& nssl_a2data31->vcp != 148
			)
	{
		printf("read_record:  bogus elevation angle of %.2f degrees.\n",
			(float)get_elev()/100.);

/*
 *	WARNING...basic level 2 data is modified here!
 */

		nssl_a2data31->elevation = elev_save;
		printf("read_record:  setting elevation angle to %.2f degrees.\n",
			(float)get_elev()/100.);
	}
#endif

	elev_save = nssl_a2data31->elevation;
	
	return;

}
