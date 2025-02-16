/*
 *	read_record_raw_disk_compress.c	Read disk data that is compressed.
 *
 *	Internal routines
 *	open_compress			Start the compress process.
 *	read_compress()			Perform i/o.
 */
 
#pragma ident "@(#)read_record_raw_disk_compress.c	5.5	01/20/04 CAPS"

#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include <strings.h>
#include <unistd.h>
#include <netinet/in.h>

#include "config.h"

#include "proto.h"

#ifdef	ALLOW_COMPRESS

/*
 *	We need to make a guess about the size, depending on what is in the
 *	data.
 */

#include "a2io.h"

#ifdef	ALLOW_TDWR
#include "tdwrio.h"
#include "tdwr_extern.h"
#endif

static int bytes_available = 0;
static int bytes_read = 0;

#define	INTERNAL_SIZE	4

static unsigned char internal_buf[INTERNAL_SIZE];
static unsigned char *internal_ptr;

static int use_internal = 0;

static int fd_compress = -1;

#define	LARGE_NUMBER_OF_RECORDS	10000

/*
 *	Size of compressed file.
 *
 *	We will read a few bytes, and try to determine whether we are TDWR
 *	or A2IO, and return a reasonable answer for the decoder routines.
 */

/*
 *	Get the size of the file.  If we return 0, then it is an error state,
 *	and we get advanced to the next file.
 */

open_compress( name )
char *name;
{
	int rc;
	int save_size;
	unsigned char buf[2];
#ifdef	ALLOW_TDWR
	short int tdwr_value;
	short int *ptr;
#endif

/*
 *	Make sure we are a compressed file.
 */

	rc = read_from_memory( buf, 2 );

	if ( rc != 2 )
		return( 0 );

	if ( buf[0] != 037 && buf[1] != 0235 )
	{
		printf("open_compress:  '%s' is not in compress format.\n",
			name);
		return( 0 );
	}

	if ( fd_compress != -1 )
	{
		close( fd_compress );
		fd_compress = -1;
	}

	rc = seek_memory( 0 );

	if ( rc <= 0 )
		return( 0 );

/*
 *	We still maintain our own bookkeeping, though open/close/read is
 *	really not handled here.
 */

	fd_compress = open( name, 0 );

	if ( fd_compress < 0 )
	{
		printf("open_compress:  '%s':  unable to open file\n", name );
		perror( name );
		return( -1 );
	}

	rc = compress_handle_open( name );

	if ( rc <= 0 )
		return( rc );

	bytes_read = 0;
	bytes_available = INTERNAL_SIZE;

/*
 *	We can't use read_from_memory(), we we need to decode the data.
 */

	rc = read_compress( (char *) internal_buf, INTERNAL_SIZE );

	if ( rc != INTERNAL_SIZE )
	{
		bytes_available = 0;
		return( 0 );
	}

#ifdef	ALLOW_TDWR
	ptr = ( short int * ) &internal_buf;
	tdwr_value = *ptr;

/*
 *	We may have to swap bytes here.
 */

	tdwr_value = htons ( tdwr_value );

	ptr = &tdwr_value;
#endif


/*
 *	We might be A2.  Check to see if we have a 24-byte header record.
 */

	if ( strncmp( (char *) internal_buf, "ARCH", 4 ) == 0 ||
		strncmp( (char *) internal_buf, "AR2V", 4 ) == 0 )
	{
		save_size = A2_VS_24_BYTE_SIZE +
			REC_SIZE_TAPE * LARGE_NUMBER_OF_RECORDS;
	}

#ifdef	ALLOW_TDWR
/*
 *	Are we TDWR?
 */

	else if ( *ptr >= 0x2b00 && *ptr <= 0x2c02 )
		save_size = TDWR_NORMAL_SIZE * LARGE_NUMBER_OF_RECORDS;
#endif

/*
 *	No 24-byte header.
 */

	else
		save_size = REC_SIZE_TAPE * LARGE_NUMBER_OF_RECORDS;


	use_internal = INTERNAL_SIZE;
	internal_ptr = internal_buf;

/*
 *	Really make things available.
 */

	bytes_read = 0;
	bytes_available = save_size;

	return( bytes_available );

}

/*
 *	Perform i/o.
 */

read_compress( mybuf, size )
char *mybuf;
int size;
{
	int request;
	int user_bytes_read;
	int copy_size;

/*
 *	Sanity checks.
 */

	if ( size <= 0 )
		return( -1 );

	if ( bytes_available == 0 )
		return( 0 );
	
	request = size;

	if ( bytes_read + size > bytes_available )
		request = bytes_available - bytes_read;

/*
 *	Sanity.
 */

	if ( request <= 0 )				/* shouldn't happen */
	{
		printf("read_compress:  request is %d???\n", request);
		close( fd_compress );
		fd_compress = -1;
		return( -1 );
	}

/*
 *	Do we have some data in memory?
 */

	if ( use_internal )
	{
		if ( request < use_internal )
			copy_size = request;
		else
			copy_size = use_internal;
		bcopy( (char *) internal_ptr, mybuf, copy_size );
		request -= copy_size;
		user_bytes_read = compress_handle_read( &mybuf[copy_size],
			request );
		user_bytes_read += copy_size;
		use_internal -= copy_size;
		if ( use_internal > 0 )
			internal_ptr = &internal_ptr[copy_size];
		else
			use_internal = 0;
	}
	else
		user_bytes_read = compress_handle_read( mybuf, request );

	if ( user_bytes_read == 0 )				/* EOF */
	{
		close( fd_compress );
		fd_compress = -1;
	}

	if ( bytes_read > bytes_available )
	{
		printf("read_compress:  bytes_read > bytes_available???\n");
		return( -1 );
	}

	bytes_read += user_bytes_read;

	return( user_bytes_read );

}

#else

#define	NO_CAN_DO	\
		printf("\nCOMPRESS access is not built in this binary.\n"); \
		exit( 1 );

open_compress( name )
char *name;
{
	NO_CAN_DO;
}

read_compress( mybuf, size )
char *mybuf;
int size;
{
	NO_CAN_DO;
}

#endif
