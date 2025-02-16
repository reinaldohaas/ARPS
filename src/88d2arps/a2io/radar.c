/*
 *	Radar...	Retrieval routines for everything but dbz/vel/spw.
 *			Everything is for the most recent record with a
 *			message type of A2_DATA_TYPE.  If this hasn't occurred,
 *			then garbage will be returned!
 *
 *	get_azi()	Returns azimuth * 100.
 *	get_beamwidth()	Returns 95.
 *	get_data_format() Return the format (type) of data.
 *	get_day()	Returns day of month.
 *	get_elev()	Returns the elevation angle * 100.  In the event of a
 *			bogus value, the last good value is returned.  (Some
 *			88Ds have a hardware problem that causes them to
 *			produce ridiculous values.)
 *	get_field_name() Map data numbers to data types.  For example number 0
 *			will return "DBZ".
 *	get_first_gate() Range to first gate in meters, and this *can* be
 *			negative.  If the requested data is not in the radial,
 *			then MISSING will be returned.  The value may or may
 *			not be modified, depending on how "read_record.c" is
 *			compiled.  The OSF requests the values not be modified.
 *	get_fixed_angle() Return the fixed elevation angle *100.  The value
 *			assumes unmodified VCPs.
 *     *get_frequency()	Returns MISSING.
 *	get_gate_spacing() Gate spacing in meters.
 *	get_hour()	Returns the hour.
 *	get_min()	Returns the minute.
 *	get_month()	Returns the month.
 *	get_nyquist()	Returns the nyquist value * 100 (m/s) if velocity is
 *			present for the current radial, otherwise, MISSING.
 *	get_number_of_gates()  Returns the number of gates.  If the field isn't
 *			present, it returns 0.
 *     *get_polarization()  Returns 1.
 *     *get_power()	Returns MISSING.
 *     *get_prf()	Returns 1000.
 *     *get_prt()	Returns MISSING.
 *	get_proj_name()	Returns "REALTIME NEXRAD".
 *     *get_pulse_width()  Returns MISSING.
 *	get_radial_status()  Returns "radial status".
 *	get_restart_flag()  When set, the 88D is restarting the current elevation
 *			scan.  Algorithms should reset themselves when set.
 *	get_scan_status() Returns 1 for a new volume scan, otherwise, 0.
 *     *get_scan_dir()	Don't depend on this code.
 *	get_sec()	Returns the seconds.
 *	get_seq()	Return message sequence number.  This is non-radar
 *			info, and is available as a debugging tool.
 *	get_status()	Returns 1 if the data requested is present, otherwise,
 *			it returns 0.
 *	get_tilt()	Returns the tilt (sweep) number.
 *	get_tilt_status() Returns 1 for a new tilt, otherwise, 0.
 *     *get_tilt_type()	Returns 1.
 *	get_timestamp()	Returns date/time in the same format as time(3).
 *	get_unix_timestamp()
 *			Returns date/time in the same format as time(3).
 *	get_vcp()	Returns VCP number.
 *	get_vel_res()	Return "vel res" number.
 *	get_vs_timestamp()	Returns date/time in the same format as time(3).
 *	get_year()	Returns the year in 2 digit format.
 *
 *	NOTE:  this assumes the current message type is A2_DATA_TYPE.  If not,
 *	       it is an error.  This error is not flagged!
 *
 *	*-  These are old "UF" routines, and don't return anything of value.
 */
 
#pragma ident "@(#)radar.c	4.1	09/15/00 NSSL"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "a2defs.h"	/* NEED NEW NAME */
#include "const.h"
#include "extern.h"
#include "extern_funct.h"
#include "proto.h"

get_azi()
{
	if ( funct_get_azi == NULL )
		radar_oops( "get_azi" );

	return( funct_get_azi() );
}

get_beamwidth()
{
	if ( funct_get_beamwidth == NULL )
		radar_oops( "get_beamwidth" );

	return( funct_get_beamwidth() );
}

get_day()
{
	return( nssl_gmt->tm_mday );
}

get_data_format()
{
	return( data_format );
}

get_elev()
{
	if ( funct_get_elev == NULL )
		radar_oops( "get_elev" );

	return( funct_get_elev() );
}

/*
 *	Provide users with some friendly conversions.
 */

get_field_name( n, s  )
int n;
char *s;
{
	if ( n == DBZ )
		strcpy(s,"DBZ");

	else if ( n == VEL )
		strcpy(s,"VEL");

	else if ( n == SPW )
		strcpy(s,"SPW");

	else if ( n == RVEL )
		strcpy(s, "RVEL" );

	else
		strcpy(s,"???");

	return;
}

get_field_num( s )
char *s;
{
	if ( !strcmp(s, "DBZ") )
		return( DBZ );

	else if ( !strcmp(s, "VEL" ) )
		return( VEL );

	else if ( !strcmp(s, "SPW" ) )
		return( SPW );

	else if ( !strcmp(s, "SW" ) )
		return( SW );

	else if ( !strcmp(s, "RVEL" ) )
		return( RVEL );

	else if ( !strcmp(s, "ZDR" ) )
		return( ZDR );

	else if ( !strcmp(s, "PHI" ) )
		return( PHI );

	else if ( !strcmp(s, "RHO" ) )
		return( RHO );

	else
		return(-1);
}

get_first_gate( n )
int n;
{
	if ( funct_get_first_gate == NULL )
		radar_oops( "get_first_gate" );

	return( funct_get_first_gate( n ) );
}

get_fixed_angle()
{
	if ( funct_get_fixed_angle == NULL )
		radar_oops( "get_fixed_angle" );

	return( funct_get_fixed_angle() );
}

/*
 *	I don't know how to retrieve this.
 */

get_frequency()
{
	return( MISSING );
}

get_gate_spacing(n)
{
	if ( funct_get_gate_spacing == NULL )
		radar_oops( "get_gate_spacing" );

	return( funct_get_gate_spacing( n ) );
}

get_hour()
{
	return( nssl_gmt->tm_hour );
}

get_min()
{
	return( nssl_gmt->tm_min );
}

get_month()
{
	return( nssl_gmt->tm_mon + 1 );
}

get_nyquist()
{
	if ( funct_get_nyquist == NULL )
		radar_oops( "get_nyquist" );

	return( funct_get_nyquist() );
}

get_number_of_gates(n)
int n;
{
	if ( funct_get_number_of_gates == NULL )
		radar_oops( "get_number_of_gates" );

	return( funct_get_number_of_gates( n ) );
}

/*
 *	This value is totally meaningless.
 */

get_polarization()
{
	return( 0 );
}

/*
 *	I'm not sure what this should be...
 */

get_power()
{
	return( MISSING );
}

/*
 *	Ditto.  I seem to recall there is some equation in one of the
 *	packages that wants a number, so we'll give it one.
 *
 *	This value is also totally meaningless.
 */

get_prf()
{
	return( 1000 );
}

/*
 *	Fill it later...
 */

get_prt()
{
/*
	if (get_status(VEL) == 0)
		return( MISSING );
 */
	return( MISSING );
}
	
/*
 *	Yes, I know this name no longers has any real values.
 */

get_proj_name(s)
char *s;
{
	bzero(s,16);
	strncpy(s,"REALTIME NEXRAD",15);
	return;
}

/*
 *	I'm not sure what this is...
 */

get_pulse_width()
{
	return( MISSING );
}

get_radial_status()
{
	if ( funct_get_radial_status == NULL )
		radar_oops( "get_radial_status" );

	return( funct_get_radial_status() );
}

/*
 *	Decode routines set this flag.
 */

get_restart_flag()
{
	return( nssl_restart_flag );
}

/*
 *	Always use our own internal bookkeeping here.
 */

get_scan()
{
	return( nssl_scan_num + vs_offset );
}

get_scan_status()
{
	return( nssl_scan_status );
}

get_scan_dir()
{
	if ( funct_get_scan_dir == NULL )
		radar_oops( "get_scan_dir" );

	return( funct_get_scan_dir() );
}

get_sec()
{
	return( nssl_gmt->tm_sec );
}

get_seq()
{
	if ( funct_get_seq == NULL )
		radar_oops( "get_seq" );

	return( funct_get_seq() );
}

get_status(n)
int n;
{
	if ( funct_get_status == NULL )
		radar_oops( "get_status" );

	return( funct_get_status( n ) );
}

get_tilt()
{
	if ( funct_get_tilt == NULL )
		radar_oops( "get_tilt" );

	return( funct_get_tilt() );
}

get_tilt_status()
{
	return( nssl_tilt_status );
}

/*
 *	Another meaningless number.
 */

get_tilt_type()
{
	return( 1 );
}

get_timestamp()
{
	if ( funct_get_timestamp == NULL )
		radar_oops( "get_timestamp" );

	return( funct_get_timestamp() );
}

/*
 *	Given "julian_date" offset (short) and a seconds offset (long),
 *	in 88D format, compute a standard Unix timestamp.
 */

get_unix_timestamp( julian_date, milsec )
short julian_date;
int milsec;
{
	int ans;

	ans = milsec;
	ans /= 1000;

/*
 *	Add the number of days.
 */

	ans += julian_date * SECONDS_PER_DAY;

/*
 *	For UNIX, time=0 is December 31, 1969, while it is January 1, 1970
 *	for NEXRAD systems, so we need to subtract off one day.
 */

	ans -= SECONDS_PER_DAY;

	return( ans );
}

get_vcp()
{
	if ( funct_get_vcp == NULL )
		radar_oops( "get_vcp" );

	return( funct_get_vcp() );
}

get_vel_res()
{
	if ( funct_get_vel_res == NULL )
		radar_oops( "get_vel_res" );

	return( funct_get_vel_res() );
}

get_vs_timestamp()
{
/*
 *	Currently only set in the A2IO live code.  Later, this needs to be
 *	add to all possibilities.
 */
	return( nssl_vs_timestamp );
}

/*
 *	tm_year returns "year - 1900", and we want two digits only.
 */

get_year()
{
	int ans;

	ans = nssl_gmt->tm_year;
	if ( ans > 99 )
		ans -= 100;
	return( ans );
}

radar_oops( s )
char *s;
{
	printf("radar_oops:  missing function for generic routine '%s'\n",
		s);
	printf("radar_oops:  FATAL ERROR!\n");

	exit( 1 );
}
