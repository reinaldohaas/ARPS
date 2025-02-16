/*
 *	Radar_fortran.c
 *
 *	This is a Fortran interface version.
 *
 *	Some obsolete functions aren't implemented.
 *
 *	Do *NOT* include the "len" argument in any Fortran calls!
 *	It is supplied automatically.
 */
 
#pragma ident	"@(#)radar_fortran.c	6.2	03/15/06	CAPS"

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

get_azi_f( azi )
int *azi;
{
	*azi = get_azi();
	return;
}

get_azi_f_( azi )
int *azi;
{
	get_azi_f( azi );
	return;
}

get_day_f( day )
int *day;
{
	*day = get_day();
	return;
}

get_day_f_( day )
int *day;
{
	get_day_f( day );
	return;
}

get_data_format_f( data_format )
int *data_format;
{
	*data_format = get_data_format();
	return;
}

get_data_format_f_( data_format )
int *data_format;
{
	get_data_format_f( data_format );
	return;
}

get_elev_f( elev )
int *elev;
{
	*elev = get_elev();
	return;
}

get_elev_f_( elev )
int *elev;
{
	get_elev_f( elev );
	return;
}

get_field_name_f( n, s, len )
int *n;
char *s;
int len;
{
	get_field_name( *n, s );
	return;
}

get_field_name_f_( n, s, len )
int *n;
char *s;
int len;
{
	get_field_name_f( n, s, len );
	return;
}

get_field_num_f( s, num, len )
char *s;
int *num;
int len;
{
	char *t;
	char *ptr;
	int i;

	t = (char *) malloc( len + 1 );

	if ( t == (char *) NULL ) {
		perror("radar_open_f:  malloc");
		exit(1);
	}

	bzero(t, len + 1);

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

	*num = get_field_num( t );

	free(t);
	return;
}

get_field_num_f_( s, num, len )
char *s;
int *num;
int len;
{
	get_field_num_f( s, num, len );
	return;
}

get_first_gate_f( n, first_gate )
int *n;
int *first_gate;
{
	*first_gate = get_first_gate( *n );
	return;
}

get_first_gate_f_( n, first_gate )
int *n;
int *first_gate;
{
	get_first_gate_f( n, first_gate );
	return;
}

get_fixed_angle_f( fixed_angle )
int *fixed_angle;
{
	*fixed_angle = get_fixed_angle();
	return;
}

get_fixed_angle_f_( fixed_angle )
int *fixed_angle;
{
	get_fixed_angle_f( fixed_angle );
	return;
}

get_gate_spacing_f( n, gate_spacing )
int *n;
int *gate_spacing;
{
	*gate_spacing = get_gate_spacing( *n );
	return;
}

get_gate_spacing_f_( n, gate_spacing )
int *n;
int *gate_spacing;
{
	get_gate_spacing_f( n, gate_spacing );
	return;
}

get_hour_f( hour )
int *hour;
{
	*hour = get_hour();
	return;
}

get_hour_f_( hour )
int *hour;
{
	get_hour_f( hour );
	return;
}

get_min_f( min )
int *min;
{
	*min = get_min();
	return;
}

get_min_f_( min )
int *min;
{
	get_min_f( min );
	return;
}

get_month_f( month )
int *month;
{
	*month = get_month();
	return;
}

get_month_f_( month )
int *month;
{
	get_month_f( month );
	return;
}

get_nyquist_f( nyquist )
int *nyquist;
{
	*nyquist = get_nyquist();
	return;
}

get_nyquist_f_( nyquist )
int *nyquist;
{
	get_nyquist_f( nyquist );
	return;
}

get_number_of_gates_f( n, number_of_gates )
int *n;
int *number_of_gates;
{
	*number_of_gates = get_number_of_gates( *n );
	return;
}

get_number_of_gates_f_( n, number_of_gates )
int *n;
int *number_of_gates;
{
	get_number_of_gates_f( n, number_of_gates );
	return;
}

get_radial_status_f( radial_status )
int *radial_status;
{
	*radial_status = get_radial_status();
	return;
}

get_radial_status_f_( radial_status )
int *radial_status;
{
	get_radial_status_f( radial_status );
	return;
}

get_restart_flag_f( restart_flag )
int *restart_flag;
{
	*restart_flag = get_restart_flag();
	return;
}

get_restart_flag_f_( restart_flag )
int *restart_flag;
{
	get_restart_flag_f( restart_flag );
	return;
}

get_scan_f( scan )
int *scan;
{
	*scan = get_scan();
	return;
}

get_scan_f_( scan )
int *scan;
{
	get_scan_f( scan );
	return;
}

get_scan_status_f( scan_status )
int *scan_status;
{
	*scan_status = get_scan_status();
	return;
}

get_scan_status_f_( scan_status )
int *scan_status;
{
	get_scan_status_f( scan_status );
	return;
}

get_sec_f( sec )
int *sec;
{
	*sec = get_sec();
	return;
}

get_sec_f_( sec )
int *sec;
{
	get_sec_f( sec );
	return;
}

get_seq_f( seq )
int *seq;
{
	*seq = get_seq();
	return;
}

get_seq_f_( seq )
int *seq;
{
	get_seq_f( seq );
	return;
}

get_status_f( n, status )
int *n;
int *status;
{
	*status = get_status( *n );
	return;
}

get_status_f_( n, status )
int *n;
int *status;
{
	get_status_f( n, status);
	return;
}

get_tilt_f( tilt )
int *tilt;
{
	*tilt = get_tilt();
	return;
}

get_tilt_f_( tilt )
int *tilt;
{
	get_tilt_f( tilt );
	return;
}

get_tilt_status_f( tilt_status )
int *tilt_status;
{
	*tilt_status = get_tilt_status();
	return;
}

get_tilt_status_f_( tilt_status )
int *tilt_status;
{
	get_tilt_status_f( tilt_status );
	return;
}

get_timestamp_f( timestamp )
int *timestamp;
{
	*timestamp = get_timestamp();
	return;
}

get_timestamp_f_( timestamp )
int *timestamp;
{
	get_timestamp_f( timestamp );
	return;
}

get_vcp_f( vcp )
int *vcp;
{
	*vcp = get_vcp();
	return;
}

get_vcp_f_( vcp )
int *vcp;
{
	get_vcp_f( vcp );
	return;
}

get_vs_timestamp_f( vs_timestamp )
int *vs_timestamp;
{
	*vs_timestamp = get_vs_timestamp();
	return;
}

get_vs_timestamp_f_( vs_timestamp )
int *vs_timestamp;
{
	get_vs_timestamp_f( vs_timestamp );
	return;
}

get_year_f( year )
int *year;
{
	*year = get_year();
	return;
}

get_year_f_( year )
int *year;
{
	get_year_f( year );
	return;
}

