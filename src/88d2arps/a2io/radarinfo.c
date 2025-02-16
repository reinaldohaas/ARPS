/*
 *	radarinfo.c...	A set of routines to return info from "radarinfo.dat".
 *			See rules before get_radar_info().
 *
 *	User Routines
 *	get_altitude()		Station altitude in meters.
 *	get_latitude()		Station latitude * 10000.
 *	get_latitude_dms()	Station latitude in degrees, minutes, seconds.
 *	get_longitude()		Station longitude * 10000.
 *	get_longitude_dms()	Station longitude in degrees, minutes, seconds.
 *	get_radar_name()	Get the radar name.
 *	get_site_name()		Get the site name.
 *	set_radar_name()	Set radar name, this *is* the 4-ltr id.
 *
 *	Internal Routines
 *	dms_df()		Convert degrees/min/sec to a whole number.
 *	get_radar_info()	Read/decode radarinfo.dat.
 */

#pragma ident "@(#)radarinfo.c	4.2	01/30/01 NSSL"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include "a2defs.h"
#include "proto.h"

static int radar_name_set = 0;
static char radar_name[4];
static char site_name[8];
static int lat, lon;
static int lat1,lat2,lat3,lon1,lon2,lon3;
static int alt;

/*
 *	The altitude info (station elevation) may now be in the data file.
 */

get_altitude()
{
	return( alt );
}


get_latitude()
{
	return( lat );
}

get_latitude_dms( d, m, s )
int *d, *m, *s;
{
	*d = lat1;
	*m = lat2;
	*s = lat3;

	return;
}

get_longitude()
{
	return( lon );
}

get_longitude_dms( d, m, s )
int *d, *m, *s;
{
	*d = lon1;
	*m = lon2;
	*s = lon3;

	return;
}

/*
 *	If the user has set a radar name via the environmental variable
 *	RADARNAME, make it immediately available.  If "radar_name_set"
 *	is zero, we don't change it so that the complaints from
 *	get_radar_info() don't come out of sequence.
 *
 *	We need this code because for live data because:
 *
 *		radar_init() calls radar_init_raw_live()
 *			which calls get_radar_name()
 *			then calls find_key()
 *
 *	using the info from get_radar_name().  If the program hasn't called
 *	set_radar_name(), the user gets a NULL string.
 */

get_radar_name(s)
char *s;
{
	char *ptr;

	bzero(s,4);

	if ( radar_name_set )
		strncpy( s, radar_name, 4 );
	else
	{
		strncpy( s, "UNKN", 4 );
		ptr = getenv("RADARNAME");
		if ( ptr != NULL )
			strncpy( s, ptr, 4 );
	}
	return;
}

get_site_name(s)
char *s;
{
	bzero(s,8);
	strncpy(s,site_name,8);
	return;
}

set_radar_name(s)
char *s;
{
	strncpy(radar_name,s,4);
	radar_name_set++;

	return;
}

/*
 *	Radar name, site name, latitude, and longitude are all important
 *	radar parameters.  Unfortunately, none of them are in the actual
 *	data themselves.
 *
 *	The user may set the radar name by:
 *	o  set_radar_name( some_name );	[in code]
 *	o  setenv RADARNAME some_name	[ignored if set_radar_name() used]
 *
 *	The program retrieves the info from the file:
 *	o  setenv RADARFILE some_file	[data file]
 *	o  from "./radarinfo.dat"	[if RADARFILE not used]
 *
 */

#define	SIZE	120

static FILE *fp;

get_radar_info()
{
	int rc;
	int nline;
	char line[SIZE];
	char *default_name = "radarinfo.dat";
	char name[4], site[8];
	char *ptr;

/*
 *	Initialize info to missing.  Do NOT initialize "radar_name" here as
 *	the user may have already called set_radar_name().
 */

	strncpy( site_name, "Unknown ", 8 );
	lat = 0;
	lon = 0;
	alt = MISSING;

/*
 *	Retrieve the radar name.
 *
 *	First, look and see if "set_radar_name()" has been called.  If it
 *	hasn't, we look for an environmental variable.
 */

	if ( radar_name_set == 0 )
	{
		ptr = getenv("RADARNAME");
		if (ptr == NULL)
		{
			printf("WARNING:  no call to set_radar_name, and no environmental variable RADARNAME\n");
			printf("WARNING:  radar name set to 'UNKN'\n");
			strncpy(radar_name,"UNKN",4);
			goto finish;
		}
		else
			strncpy(radar_name,ptr,4);
	}

/*
 *	Open the radar info file.
 *
 *	Look for an environment variable.
 *
 *	If the environmental variable isn't set,  look for the radar file in
 *	the current directory.  If it is set, and the file doesn't exist,
 *	complain and still look for the file in the current directory.
 */

	ptr = getenv("RADARFILE");
	if ( ptr != NULL )
	{
		fp = fopen(ptr,"r");
		if (fp == NULL)
			perror(ptr);
	}

	if (fp == NULL)
	{
		fp = fopen(default_name,"r");
		if (fp == NULL)
		{
			printf("WARNING:  can't find a valid '%s' file.\n",
			       default_name);
			goto finish;
		}
	}

	fgets(line,SIZE,fp);			/* eat header line */

	nline = 0;
	while( fgets(line,SIZE,fp) != NULL )
	{
		nline++;
		rc = sscanf(line,"%4s %8c %d %d %d %d %d %d %d",
			name,site,&lat1,&lat2,&lat3,&lon1,&lon2,&lon3,&alt);
		if ( rc < 8 || rc > 9 )
		{
			printf("get_radar_info:  %d arguments returned.  ", rc);
/*
			printf("Radarinfo file is corrupt.  Exit.\n");
			exit( 1 );
*/
            printf("WARNING:  Radarinfo file is corrupt at line %d, rc=%d\n",
				nline,rc);
		}
		if ( rc == 8 )
			alt = MISSING;
#ifdef	DEBUG_RADARINFO
		printf("name %4s site %8.8s lat1 %d lat2 %d lat3 %d ",
			name, site, lat1, lat2, lat3 );
		printf("lon1 %d lon2 %d lon3 %d alt %d\n",
			lon1, lon2, lon3, alt );
#endif
		if (strncmp( radar_name, name, 4 ) == 0)
		{
			strncpy(site_name, site, 8);
			break;
		}
	}

	fclose(fp);

	fp = NULL;

	lat = dms_df( lat1, lat2, lat3 );
	lon = dms_df( lon1, lon2, lon3 );

	if ( lat == 0 )
		printf("WARNING:  radar '%4.4s' not found in data file\n",
			radar_name);

finish:	;

	if ( lat == 0 )
	{
		printf("WARNING:  latitude and longitude set to 0\n");
		if ( strncmp( radar_name, "UNKN", 4 ) == 0 )
			printf("INFO:     will watch for NCDC header file\n");
	}
#ifdef	DEBUG
	printf("\nRadar:	%-4.4s\nSite:	%-8.8s\nLat.:	%f\nLong.:	%f\n\n",
		radar_name, site_name, lat/100000., lon/100000.);
#endif	/* DEBUG */

	return;
}

/*
 *	Convert degrees, minutes, seconds to an integer (floating point
 *	value times 100000).
 */

dms_df( a, b, c )
int a, b, c;
{
	float tmp;
	int ans;
	int sign;
	int aval;

	if ( a >= 0 )
	{
		sign = 1;
		aval = a;
	}
	else
	{
		sign = -1;
		aval = -a;
	}

	tmp = ( (float) aval + (float) b / 60.0 + (float) c / 3600.0 );

	ans = sign * (int) ( tmp * 100000 );

	return( ans );
}
