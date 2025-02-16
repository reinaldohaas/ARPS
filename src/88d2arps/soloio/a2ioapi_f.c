/*
 *	A2io_fortran.c
 *
 *	This is a Fortran interface version.
 *
 *	Do *NOT* include the "len" argument in any Fortran calls!
 *	It is supplied automatically.
 *
 */

#ifdef UNDERSCORE
#define set_radar98_f         set_radar98_f_
#define set_radar_name_f      set_radar_name_f_
#define get_year_f            get_year_f_
#define get_month_f           get_month_f_
#define get_day_f             get_day_f_
#define get_hour_f            get_hour_f_
#define get_min_f             get_min_f_
#define get_sec_f             get_sec_f_
#define get_azi_f             get_azi_f_
#define get_vcp_f             get_vcp_f_
#define get_field_num_f       get_field_num_f_
#define get_nyquist_f         get_nyquist_f_
#define get_data_field_f      get_data_field_f_
#define get_gate_spacing_f    get_gate_spacing_f_
#define get_first_gate_f      get_first_gate_f_
#define get_number_of_gates_f get_number_of_gates_f_
#define get_fixed_angle_f     get_fixed_angle_f_
#define get_tilt_f            get_tilt_f_
#define get_scan_f            get_scan_f_
#define get_status_f          get_status_f_
#define get_restart_flag_f    get_restart_flag_f_
#define get_latitude_f        get_latitude_f_
#define get_longitude_f       get_longitude_f_
#define get_altitude_f        get_altitude_f_
#define read_radial_f         read_radial_f_
#define radar_open_f          radar_open_f_
#define get_data_field_f      get_data_field_f_
#define get_radar_info_f      get_radar_info_f_
#endif

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <strings.h>


radar_open_f( s, irc, len )
char *s;
int *irc;
int len;
{
   /*
    *  We need to extract only the valid part of "s".
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
     *  Change spaces to ascii nulls
     */

    ptr = t;

    for(i=0; i<len; i++)
    {
        if ( *ptr == ' ' ) *ptr = (char) NULL;
        ptr++;
    }

	*irc = radar_init( t );

	return;
}

set_radar98_f( n )
int *n;
{
	set_radar98(*n);
	return;
}

read_radial_f( irc )
int *irc;
{
    *irc = read_radial();
    return;
}

get_azi_f( azi )
int *azi;
{
    *azi = get_azi();
    return;
}

get_day_f( day )
int *day;
{
    *day = get_day();
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

get_first_gate_f( n, first_gate )
int *n;
int *first_gate;
{
    *first_gate = get_first_gate( *n );
    return;
}

get_fixed_angle_f( fixed_angle )
int *fixed_angle;
{
    *fixed_angle = get_fixed_angle();
    return;
}

get_gate_spacing_f( n, gate_spacing )
int *n;
int *gate_spacing;
{
    *gate_spacing = get_gate_spacing( *n );
    return;
}

get_hour_f( hour )
int *hour;
{
    *hour = get_hour();
    return;
}

get_min_f( min )
int *min;
{
    *min = get_min();
    return;
}

get_month_f( month )
int *month;
{
    *month = get_month();
    return;
}

get_nyquist_f( nyquist )
int *nyquist;
{
    *nyquist = get_nyquist();
    return;
}

get_number_of_gates_f( n, number_of_gates )
int *n;
int *number_of_gates;
{
    *number_of_gates = get_number_of_gates( *n );
    return;
}

get_restart_flag_f( restart_flag )
int *restart_flag;
{
    *restart_flag = get_restart_flag();
    return;
}

get_scan_f( scan )
int *scan;
{
    *scan = get_scan();
    return;
}

get_sec_f( sec )
int *sec;
{
    *sec = get_sec();
    return;
}

get_status_f( n, status )
int *n;
int *status;
{
    *status = get_status( *n );
    return;
}

get_tilt_f( tilt )
int *tilt;
{
    *tilt = get_tilt();
    return;
}

get_vcp_f( vcp )
int *vcp;
{
    *vcp = get_vcp();
    return;
}

get_year_f( year )
int *year;
{
    *year = get_year();
    return;
}

get_altitude_f( altitude )
int *altitude;
{
    *altitude = get_altitude();
    return;
}

get_latitude_f( latitude )
int *latitude;
{
    *latitude = get_latitude();
    return;
}

get_longitude_f( longitude )
int *longitude;
{
    *longitude = get_longitude();
    return;
}

set_radar_name_f( s, alen )
char *s;
int  *alen;
{
    set_radar_name( s, *alen );
    return;
}

get_data_field_f( num, data, n, irc, len )
int *num;
float data[];
int *n;
int *irc;
{
    *irc = get_data_field( *num, data, *n );
    return;
}

void get_radar_info_f( )
{
  return;
}
