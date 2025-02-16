/*
 *
 *    radarinfo.c
 *
 *	This is a subsection of the WSR-88D archive-II I/O library
 *    codes, a2io.c,  written by Kevin Thomas of NSSL, and
 *    shamelessly stolen and adapted by Keith Brewster, CAPS.
 *    It reads the a2io radar information file, radarinfo.dat
 *
 *    Keith Brewster, May, 1996
 *
 */
 
#include <stdio.h>
#include <strings.h>

#define MISSING -999

int radar_name_set =0;
char radar_name[4];
char site_name[8];
int lat;
int lon;
int alt;

int lat1,lat2,lat3;
int lon1,lon2,lon3;

/*
 *	The altitude info (station elevation) may now be in the data file.
 */

int get_altitude()
{
	return( alt );
}

int get_latitude()
{
	return( lat );
}

int get_longitude()
{
	return( lon );
}

void get_radar_name(s)
char *s;
{
	strncpy(s,radar_name,4);
	return;
}

void get_site_name(s)
char *s;
{
	strncpy(s,site_name,8);
	return;
}

void set_radar_name(s)
char *s;
{
	strncpy(radar_name,s,4);
	radar_name_set++;
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
 */

#define	SIZE	120

FILE *fp;

void get_radar_info()
{
	int rc;
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
		ptr = (char *)getenv("RADARNAME");
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

	ptr = (char *)getenv("RADARFILE");
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
			printf("WARNING:  can't find a valid '%s' file\n",
				default_name);
			goto finish;
		}
	}

	fgets(line,SIZE,fp);			/* eat header line */

	while( fgets(line,SIZE,fp) != NULL )
	{
		rc = sscanf(line,"%4s %8c %d %d %d %d %d %d %d",
			name,site,&lat1,&lat2,&lat3,&lon1,&lon2,&lon3,&alt);
		if ( rc < 8 || rc > 9 )
		{
			printf("get_radar_info:  %d arguments returned.  ",rc);
			printf("Radarinfo file is corrupt.  Exit.\n");
			exit( 1 );
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

	lat = dms_df( lat1, lat2, lat3 );
	lon = dms_df( lon1, lon2, lon3 );

	if ( lat == 0 )
		printf("WARNING:  radar '%s' not found in data file\n",
			radar_name);

finish:	;

	if ( lat == 0 )
		printf("WARNING:  latitude and longitude set to 0\n");
#ifdef	DEBUG_A2IO
	printf("\nRadar:	%-4.4s\nSite:	%-8.8s\nLat.:	%f\nLong.:	%f\n\n",
		radar_name, site_name, lat/100000., lon/100000.);
#endif	/* DEBUG_A2IO */

	return;
}
/*
 *    Convert degrees, minutes, seconds to an integer (floating point
 *    value times 100000).
 */

int dms_df(a,b,c)
int a,b,c;
{
      float tmp;
      int ans;

      tmp = ( (float) a + (float) b / 60.0 + (float) c / 3600.0 );

      ans = (int) ( tmp * 100000 );

      return( ans );
}
