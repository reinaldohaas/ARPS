/*
 *	read_record_raw_disk_dir.c	Perform i/o for disk directories.
 *
 *	Internal
 *	radar_open_raw_disk_dir()	Open routine for directories.
 *	radar_close_raw_disk_dir()	Closing routine for directories.
 *	read_record_raw_disk_dir()	I/O routine for disk.
 *	process_disk()			Make an index of the data on disk.
 *	vs_cmp()			Function used by qsort().
 *	get_disk_format()		File format (KWT, LDM)
 *	match()				Pattern matching routine.
 *	decode_vs_number()		Extract VS number from file name.
 *	set_disk_eod()			Set the "end of data" flag.
 *	disk_advance_flag()		Allow another routine to advance to
 *					the next file without performing i/o.
 */
 

#pragma ident "@(#)read_record_raw_disk_dir.c	5.3	03/12/04 CAPS"

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <dirent.h>
#include <malloc.h>
#include <string.h>
#include <strings.h>
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/mman.h>
#include "a2io.h"
#include "config.h"
#include "extern.h"

static char *pdir;
static char *name;

static DIR *dirp = NULL;
static struct dirent *dp;

#define	MAX_NUM	5000
#define	MAX_LEN	30

struct	diskdir {
	char name[MAX_LEN];		/* name of data file */
	long long vs;			/* volume scan #, already decoded */
					/* long long to hold LDM timestamp */
} ;

struct diskdir diskdir[MAX_NUM];
struct diskdir *pdiskdir;

#define	DONE_DATA	-999

static int ninfo = 0;

static int file_current = -1;
static int time_to_open = 1;
static long long old_vs = -1;
static int disk_no_data = 0;

/*
static int disk_format = DISK;
 */

/*
 *	Open/initialization routine.
 */

radar_open_raw_disk_dir( s )
char *s;
{
	if ( disk_no_data )
	{
		disk_no_data = 0;
		file_current = -1;
		time_to_open = 1;
		free( name );
		free( pdir );
	}

	process_disk( s );

	return( 0 );
}

/*
 *	Terminate this data set.
 */

radar_close_raw_disk_dir()
{
	dirp = NULL;
	return;
}

/*
 *	Perform i/o for the user.
 */

read_record_raw_disk_dir( mybuf, mybuf_size )
char *mybuf;
int mybuf_size;
{
	int len;
	int rc;

	if ( disk_no_data )
	{
		printf("User program ignored 'end of data set'.  Exit.\n");
		exit( 1 );
	}

/*
 *	Give the user one chance.
 */

	if ( dirp == NULL )
	{
		printf("read_record_raw_disk_dir:  No directories open.\n");
		printf("read_record_raw_disk_dir:  Returning EOF.\n");

		set_disk_eod();
		return( RADAR_END_OF_DATA );
	}

loop:	;

	if ( time_to_open )
	{
		file_current++;

		pdiskdir = &diskdir[file_current];
		if ( pdiskdir->vs == DONE_DATA )
		{
			set_disk_eod();
			return( RADAR_END_OF_DATA );
		}

/*
 *	Make sure we don't have a duplicate number.  If we do, then we had
 *	both a gzipped and a non-gzipped file present.
 *
 *	Same things applies to bzipped files, though this isn't supposed to
 *	happen.
 */

		if ( pdiskdir->vs == old_vs )
			goto loop;
		old_vs = pdiskdir->vs;

/*
 *	Lots of space.
 */

		len = strlen( pdir ) + MAX_LEN;

		name = malloc( len );
		bzero( name, len );

		sprintf( name, "%s/%s", pdir, pdiskdir->name );

/*
 *	Let the file reader stuff take over.
 */

		rc = radar_open_raw_disk_file( name );

/*
 *	No need to tell the user about an error, as it has already been said.
 */

		if ( rc < 0 )
			goto loop;

		time_to_open = 0;

	}

	rc = read_record_raw_disk_file( mybuf, mybuf_size );

	if ( rc <= 0 )
	{
		time_to_open++;
		goto loop;
	}

	return( rc );

}

/*
 *	Make a table of the volume scans available.
 */

process_disk( s )
char *s;
{
	char rname[4];
	int vs_cmp();
	int format_kwt = 0;
	int format_ldm = 0;
	int format_ncdc = 0;
	int format_caps = 0;
	int format_tdwr = 0;
	long long vnum;
	long long decode_vs_number();
	int flag;

/*
 *	Save the directory name for later use.
 */

	pdir = ( char * ) malloc( strlen( s ) + 4 );

	if ( pdir == NULL )
	{
		perror("process_disk:  malloc");
		exit( 1 );
	}

	strcpy( pdir, s );

	pdiskdir = diskdir;

	dirp = opendir( s );

/*
 *	If we are here, and we aren't a directory, fail, as this should be
 *	impossible.
 */

	if ( dirp == NULL )
	{
		printf("process_data:  opendir %s failed???", s );
		set_disk_eod();
	}

	ninfo = 0;

	format_kwt = 0;

	get_radar_name( rname );

	for(;;)
	{
		dp = readdir( dirp );
		if ( dp == NULL )
			break;

/*
 *	Name isn't malloc'd at this time.
 */
		name = dp->d_name;

		disk_format = get_disk_format( name );

		if ( disk_format == DISK_KWT )
			format_kwt++;

		else if ( disk_format == DISK_LDM )
			format_ldm++;
		else if ( disk_format == DISK_NCDC )
			format_ncdc++;
		else if ( disk_format == DISK_CAPS )
			format_caps++;
		else if ( disk_format == DISK_TDWR )
			format_tdwr++;
		else
			continue;

/*
 *	We can't want to handle multiple formats.
 */

		flag = 0;
		flag += ( format_kwt >= 1 );
		flag += ( format_ldm >= 1 );
		flag += ( format_ncdc >= 1 );
		flag += ( format_caps >= 1 );

		if ( flag > 1 )
		{
			if ( format_kwt )
				printf("process_disk:  format KWT found.\n");
			if ( format_ldm )
				printf("process_disk:  format LDM found.\n");
			if ( format_ncdc )
				printf("process_disk:  format NCDC found.\n");
			if ( format_caps )
				printf("process_disk:  format CAPS found.\n");
			if ( format_tdwr )
				printf("process_disk:  format TDWR found.\n");
			printf("process_disk:  only one format is allowed.\n");
 			exit( 1 );
		}

		vnum = decode_vs_number( name );

/*
 *	Sanity.
 */
		if ( vnum < 0 )
			continue;

		strcpy( pdiskdir->name, name );

		pdiskdir->vs = vnum;

		pdiskdir++;
		ninfo++;

		if ( ninfo >= MAX_NUM )
		{
			printf("process_disk:  MAX_NUM is too small.\n");
			exit( 1 );
		}

	}

	pdiskdir->vs = DONE_DATA;

	closedir( dirp );

	qsort( diskdir, ninfo, sizeof( struct diskdir), vs_cmp );

	return;
}

long long
decode_vs_number( s )
char *s;
{
	char tmp[20];
	char *ptr;

#ifndef	NO_ATOLL

	if ( disk_format == DISK_KWT )
	{
		strcpy( tmp, s );

		ptr = strchr( tmp, '.' );
		if ( ptr == NULL )
		{
			printf("decode_vs_number:  %s:  can't find a period\n",
					tmp );
			return( -999 );
		}
		*ptr = (char) NULL;

		return( atoll( &tmp[5] ) );
	}

/*
 *	"vs" number is really a time stamp here.
 */

	else if ( disk_format == DISK_LDM )
	{
		bzero( tmp, 20 );
		strncpy( tmp, s, 14 );
		return( atoll( tmp ) );
	}

/*
 *	Here too, though it is a little more complicated to load.
 */

	else if ( disk_format == DISK_NCDC )
	{
		bzero( tmp, 20 );
		strncpy( tmp, &s[4], 8 );
		strncpy( &tmp[8], &s[13], 6 );
		return( atoll( tmp ) );
	}

/*
 *	Here three.
 */

	else if ( disk_format == DISK_CAPS )
	{
		bzero( tmp, 20 );
		strncpy( tmp, &s[5], 8 );
		strncpy( &tmp[8], &s[14], 6 );
		return( atoll( tmp ) );
	}

	else if ( disk_format == DISK_TDWR )
	{
		bzero( tmp, 20 );
		strncpy( tmp, &s[17], 2 );
/*
		strncpy( &tmp[8], &s[14], 6 );
 */
		printf("tmp=%s\n",tmp);
		return( atoll( tmp ) );
	}
#endif

/*
 *	Sanity.
 */

	return( -999 );

}

vs_cmp( s, t )
struct diskdir *s;
struct diskdir *t;
{
	return( s->vs - t->vs );
}

get_disk_format( s )
char *s;
{
	char rname[5];
	char str[20];

	if ( match( "NX88d#.0", s ) )
		return( DISK_KWT );
	if ( match( "NX88d##.0", s ) )
		return( DISK_KWT );
	if ( match( "NX88d###.0", s ) )
		return( DISK_KWT );
	if ( match( "##############.ldm", s ) )
		return( DISK_LDM);
	if ( match( "????########_######", s ) )
		return( DISK_NCDC );

/*
 *	The backslash is necessary to avoid a nasty Sun f77 "feature".  When
 *	? ? - (spaces added so the bug isn't tripped) appears, it gets changed
 *	to a tilde character.
 */

	if ( match( "???\?_########_######.ridds", s ) )
		return( DISK_CAPS );

  	get_radar_name(rname);
	sprintf(str,"%4.4s???????????\?-##-*", rname);
	if ( match(str,s) )
		return( DISK_TDWR );
/*
	if ( match( "T??????????????\?-##-*", s ) )
		return (DISK_TDWR);
 */

	return( DISK );
}

/*
 *	Try to match file types.
 *
 *	Symbols:
 *		#	number
 *		?	any character
 *		*	the rest of the strings matches if one char is present.
 *			NOTHING WILL BE CHECKED AFTER "*"!
 */


#define	MATCH	1
#define	NOMATCH	0

match( format, file )
char *format, *file;
{
	char *ptr1, *ptr2;
	char *ptr;

	ptr1 = format;

/*
 *	Remove path information.
 */

	ptr = rindex( file, '/' );

	if ( ptr == NULL )
		ptr2 = file;
	else
	{
		ptr++;
		ptr2 = ptr;
	}

	while( *ptr1 )
	{
		if ( *ptr2 == (char) NULL )
			return( NOMATCH );

		if ( *ptr1 == '#' )
		{
			if ( !isdigit( (int) *ptr2 ) )
				return( NOMATCH );
		}
		else if ( *ptr1 == '?' )
			;
		else if ( *ptr1 == '*' )
			return( MATCH );
		else
		{
			if ( *ptr1 != *ptr2 )
				return( NOMATCH );
		}
		ptr1++;
		ptr2++;
		
	}

	if ( *ptr2 == (char) NULL )
		return( MATCH );

/*
 *	Allow ".gz" extensions.
 */

	if ( strlen( ptr2 ) == 3 && ( strcmp( ptr2, ".gz" ) == 0 ) )
		return( MATCH );

/*
 *	And ".bz2".
 */

	if ( strlen( ptr2 ) == 4 && ( strcmp( ptr2, ".bz2" ) == 0 ) )
		return( MATCH );

/*
 *	And ".Z".
 */

	if ( strlen( ptr2 ) == 2 && ( strcmp( ptr2, ".Z" ) == 0 ) )
		return( MATCH );

	return( NOMATCH );
}

/*
 *	Set "end of data".
 */

set_disk_eod()
{
/*
 *	Sanity.
 */

	disk_no_data++;

	return;
}

/*
 *	Allow an external routine ("read_record_tape()") to tell us to advance
 *	to the next file (volume scan search).
 *
 *	NEW CODE REQUIRED!!!
 */

disk_advance_file()
{
/*
 *	If the input is a single file, it may be captured live data, so it
 *	may not be at the beginning of a volume scan.  Be polite and don't
 *	close it.
 *
 *	If the input is a directory, don't make this assumption, even if
 *	there is only one file out there.  We might change this in the future.
 */

	time_to_open++;

	return;
}
