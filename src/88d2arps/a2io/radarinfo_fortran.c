/*
 *	Radarinfo_fortran.c
 *
 *	This is a Fortran interface version.
 *
 *	Some obsolete functions aren't implemented.
 *
 *	Do *NOT* include the "len" argument in any Fortran calls!
 *	It is supplied automatically.
 */
 
#pragma ident	"@(#)radarinfo_fortran.c	5.2	03/03/05	CAPS"

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <strings.h>

get_altitude_f( altitude )
int *altitude;
{
	*altitude = get_altitude();
	return;
}

get_altitude_f_( altitude )
int *altitude;
{
	get_altitude_f( altitude );
	return;
}

get_latitude_f( latitude )
int *latitude;
{
	*latitude = get_latitude();
	return;
}

get_latitude_f_( latitude )
int *latitude;
{
	get_latitude_f( latitude );
	return;
}

get_latitude_dms_f( d, m, s )
int *d, *m, *s;
{
	get_latitude_dms( d, m, s );
	return;
}

get_latitude_dms_f_( d, m, s )
int *d, *m, *s;
{
	get_latitude_dms_f( d, m, s );
	return;
}

get_longitude_f( longitude )
int *longitude;
{
	*longitude = get_longitude();
	return;
}

get_longitude_f_( longitude )
int *longitude;
{
	get_longitude_f( longitude );
	return;
}

get_longitude_dms_f( d, m, s )
int *d, *m, *s;
{
	get_longitude_dms( d, m, s );
	return;
}

get_longitude_dms_f_( d, m, s )
int *d, *m, *s;
{
	get_longitude_dms_f( d, m, s );
	return;
}

get_radar_info_f()
{
	get_radar_info();
}

get_radar_info_f_()
{
	get_radar_info_f();
}

set_radar_name_f( s, len )
char *s;
int len;
{
	set_radar_name( s );
	return;
}

set_radar_name_f_( s, len )
char *s;
int len;
{
	set_radar_name_f( s, len);
	return;
}
