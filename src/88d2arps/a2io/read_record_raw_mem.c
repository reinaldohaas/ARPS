/*
 *	read_record_raw_mem.c		Make data from memory (possibly via
 *					mmap()) available.
 *
 *	Internal
 *	set_memory()			Set up memory pointers.
 *	read_from_memory()		"Read" data.
 *	seek_memory()			Seek to location.
 */
 
#pragma ident "@(#)read_record_raw_mem.c	5.3	10/15/01 CAPS"

#include <stdio.h>
#include <fcntl.h>
#include <strings.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "config.h"

#include "proto.h"
#include "mem.h"

static unsigned char *ptr = NULL;		/* data pointer */

static int fd_mem = -1;
static FILE *fp_mem = NULL;

static int mem_start = 0;
static int mem_end = 0;
static int mem_bytes_read = 0;
static int mem_bytes_available = 0;
static int total_bytes_available = 0;

/*
 *	Set up our pointer.  This must be done before any "standard" open
 *	or access calls.
 */

set_memory( myptr, mysize )
unsigned char *myptr;
int mysize;
{
	char tmpname[20];
	int no_null;
/*
 *	Do this only once so we have a valid open file pointer in memory.
 *	The file itself won't be open!
 */
	if ( fp_mem == NULL )
	{
		no_null = 0;
		fp_mem = fopen( "/dev/null", "r" );

		if ( fp_mem == NULL )
		{
			no_null = 1;
/*
 *	Catamount systems don't have a /dev/null!  Note that the file we
 *	do use has to be "w", not "r"!
 */

/*
			perror( "set_memory:  /dev/null");
 */
			sprintf(tmpname,"/tmp/null%d",getpid());
			fp_mem = fopen( tmpname, "w" );
/*
 *	Give up on error.
 */
			if ( fp_mem == NULL )
			{
				perror(tmpname);
				return( -1 );
			}
		}
		fd_mem = fileno( fp_mem );
		fclose( fp_mem );
		if ( no_null ) unlink(tmpname);
	}

	ptr = myptr;

	mem_start = 0;
	mem_end = mysize;

	mem_bytes_read = 0;
	mem_bytes_available = mem_end;
	total_bytes_available = mem_bytes_available;

	return( 0 );
}

/*
 *	Perform i/o.
 */

read_from_memory( mybuf, size )
unsigned char *mybuf;
int size;
{
	int request;

	if ( fd_mem == -1 )
	{
		printf("read_from_memory:  set_memory needs to be called first\n");
		return( 0 );
	}

	if ( ptr == NULL )
	{
		printf("read_from_memory:  no data in memory!\n");
		return( 0 );
	}

/*
 *	Sanity checks.
 */

	if ( size <= 0 )
		return( -1 );

	if ( mem_bytes_available == 0 )
		return( 0 );


/*
 *	Return if there isn't any more data to read.
 */

	if ( mem_bytes_read == mem_bytes_available )
		return( 0 );

	request = size;

	if ( mem_bytes_read + size > mem_bytes_available )
		request = mem_bytes_available - mem_bytes_read;

/*
 *	Sanity.  We don't want bcopy() to core dump.
 */

	if ( request <= 0 )			/* shouldn't happen */
	{
		printf("read_from_memory:  request is %d???\n", request);
		return( 0 );
	}

	bcopy( &ptr[mem_bytes_read+mem_start], mybuf, (unsigned) request );

	mem_bytes_read += request;

	return( request );
}

/*
 *	Seek to a location in memory.  This makes a handy "rewind" operation.
 */

seek_memory( num )
int num;
{
/*
 *	Maybe I should add a lseek_memory().
 */
	if ( num == MAGIC_SEEK )
		return( mem_bytes_read );

	if ( num < 0 || num > mem_bytes_available )
 
	{
		printf("seek_memory:  %d is not valid\n", num);
		return( -1 );
	}

	mem_bytes_read = num;

	return( mem_bytes_available - mem_bytes_read );
 
}
