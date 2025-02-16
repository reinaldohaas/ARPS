/*
 *	read_record_raw_disk_mmap.c	Handle memory mapping for disk files.
 *
 *	Internal
 *	open_mmap()			Open memory mapping.
 *	close_mmap()			Close memory mapping.
 *	seek_mmap()			Seek to a location in the mmap'd file.
 *	read_mmap()			Perform i/o.
 */
 
#pragma ident "@(#)read_record_raw_disk_mmap.c	6.1	03/30/06 CAPS"

#include <stdio.h>
#include <fcntl.h>
#include <malloc.h>
#include <strings.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>

#include "config.h"

#include "proto.h"

static char *map = ( char * ) -1; 

static int save_size = 0;

/*
 *	Open the file.  Return the number of bytes available.  This is the
 *	actual number of bytes, not in any decoded (such as ungzipped)
 *	format.
 */

open_mmap( name )
char *name;
{
	struct stat stat;
	int rc;
	int fd;

/*
 *	Make sure we don't mess up memory.
 */

	if ( map != ( char * ) -1 )
		munmap( map, save_size );

	fd = open( name, 0 );

	if ( fd < 0 )
	{
		perror( name );
		return( 0 );
	}

	rc = fstat( fd, &stat );

	if ( rc < 0 )
	{
		perror( name );
		close( fd );
		return( 0 );
	}

	if ( stat.st_size == 0 )
	{
		close( fd );
		return( 0 );
	}

	save_size = stat.st_size;

/*
 *	Always map the entire file.  This way, we don't have to worry about
 *	page size logic.  The actual seek location is handled below.
 *
 *	If mmap() fails fall back to traditional read().  LEMIEUX has a
 *	scratch filesystem that doesn't supported mmap().
 */

	map = mmap( 0, save_size, PROT_READ, MAP_SHARED, fd, 0 );

	if ( map == (char *) -1 )
	{
		printf("open_mmap:  apparently not supported on filesystem.\n");
		printf("open_mmap:  falling back to 'read'.\n");

		map = (char *) malloc(save_size);
		if ( map == (char *) NULL )
		{
			perror("open_mmap:  error allocation memory for read");
			return(-1);
		}
		rc = read(fd,map,save_size);
		if (rc != save_size)
		{
			perror("open_mmap:  read failed");
			return(-1);
		}
	}

	close( fd );

/*
 *	Advise the system we are doing sequential i/o.  Not all systems have
 *	this call.
 *
 *	Not all systems have this implemented in the kernel, so we'll just
 *	ignore the return code.
 */

#ifdef	MADV_SEQUENTIAL
	rc = madvise( map, save_size, MADV_SEQUENTIAL );

/*
	if ( rc < 0 )
		perror("open_mmap:  WARNING:  madvise failed");
 */
#endif

	rc = set_memory( (unsigned char *) map, save_size );

	if ( rc != 0 )
		return( rc );

	return( save_size );
}

close_mmap()
{
	if ( map != ( char * ) -1 )
	{
		munmap( map, save_size );
		map = ( char * ) -1;
	}

	return( 0 );
}

seek_mmap( n )
int n;
{
	return( seek_memory( n ) );
}

/*
 *	Perform i/o.
 */

read_mmap( mybuf, size )
unsigned char *mybuf;
int size;
{
	int rc;

	rc = read_from_memory( mybuf, size );

	if ( rc <= 0 )
	{
		if ( map != ( char * ) -1 )
		{
			munmap( map, save_size );
			map = ( char * ) -1;
		}
	}

	return( rc );
}
