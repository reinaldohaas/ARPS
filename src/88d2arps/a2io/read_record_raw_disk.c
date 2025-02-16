/*
 *	read_record_raw_disk.c	Control i/o for disk data.
 *
 *	Internal
 *	radar_open_raw_disk()		Open routine for disk.
 *	radar_close_raw_disk()		Close routine for disk.
 *	read_record_raw_disk()		I/O routine for disk.
 *
 *	The original read_record_raw_disk.c is now split into:
 *
 *		read_record_raw_disk_dir.c
 *		read_record_raw_disk_file.c
 *
 *	This routine just acts as a controller.
 */
 
#pragma ident "@(#)read_record_raw_disk.c	4.2	02/23/01 NSSL"

#include <stdio.h>
#include <stdlib.h>
#include <dirent.h>

static DIR *dirp;

static int dup_call = 0;

int (*data_io_open)();
int (*data_io)();

int radar_open_raw_disk_dir(), read_record_raw_disk_dir();
int radar_open_raw_disk_file(), read_record_raw_disk_file();

/*
 *	Open/initialization routine.
 */

radar_open_raw_disk( s )
char *s;
{
	int fd, rc;

/*
 *	Close all open data sets.
 */

	radar_close_raw_disk();

/*
 *	Give the user one chance.
 */

	dup_call = 0;

/*
 *	Opens on directories will be successful.
 */

	fd = open( s, 0 );

/*
 *	Go no further if we have an error.
 */

	if ( fd < 0 )
	{
		perror( s );
		return( -1 );
	}

	close( fd );

/*
 *	Are we a directory?
 */

	dirp = opendir( s );

	if ( dirp == ( DIR * ) NULL )
	{
		data_io_open = radar_open_raw_disk_file;
		data_io = read_record_raw_disk_file;
	}
	else
	{
		closedir( dirp );

		data_io_open = radar_open_raw_disk_dir;
		data_io = read_record_raw_disk_dir;
	}


	if ( data_io_open == NULL )
	{
		printf("radar_open_raw_disk:  data_io_open is NULL???\n");
		return( -1 );
	}

	rc = data_io_open( s );

	return( rc );
}

/*
 *	Close open data sets.  Nothing harmful will happen if nothing is open.
 */

radar_close_raw_disk()
{
	radar_close_raw_disk_dir();
	radar_close_raw_disk_file();

	data_io_open = NULL;
	data_io = NULL;

	return( 0 );
}

/*
 *	Perform i/o for the user.  Let the called routine handle errors.
 */

read_record_raw_disk( mybuf, mybuf_size )
char *mybuf;
int mybuf_size;
{
	int rc;

	if ( data_io == NULL )
	{
/*
 *	Just fail.  The user program had its chance.
 */
		if ( dup_call )
		{
			printf("read_record_raw_disk:  FATAL:  no open datasets\n");
			exit( 1 );
		}
		printf("read_record_raw_disk:  no data set is open\n");
		dup_call++;
		return( -1 );
	}

	dup_call = 0;

	rc = data_io( mybuf, mybuf_size );

	return( rc );
}
