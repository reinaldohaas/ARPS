/*
 *	read_record_raw_disk_bzip2.c	Read disk data that is bzipped.
 *
 *	Internal routines
 *	open_bzip2			Start the BZIP2 process.
 *	read_bzip2			Perform the I/O.
 *
 */
 
#pragma ident "@(#)read_record_raw_disk_bzip2.c	5.6	01/20/04 CAPS"

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <unistd.h>
#include <netinet/in.h>

#include "config.h"

#include "proto.h"

#ifdef	ALLOW_BZIP2

#include "bzlib.h"

/*
 *	Drivers shouldn't depend on the data format, but we do, due to a
 *	BZIP2 design flaw.  It seems that the uncompressed record size isn't
 *	recorded in the data.  The program author has confirmed this.
 */

#include "a2io.h"

#ifdef	ALLOW_TDWR
#include "tdwrio.h"
#include "tdwr_extern.h"
#endif

static int bytes_available = 0;
static int bytes_read = 0;

static BZFILE *fbz;

#define	INTERNAL_SIZE	4

static unsigned char internal_buf[INTERNAL_SIZE];
static unsigned char *internal_ptr;

static int use_internal = 0;

#define	LARGE_NUMBER_OF_RECORDS	10000

/*
 *	Size of bzipped file.
 *
 *	We will read a few bytes, and try to determine whether we are TDWR
 *	or A2IO, and return a reasonable answer for the decoder routines.
 */

/*
 *	Get the size of the file.  If we return 0, then it is an error state,
 *	and we get advanced to the next file.
 */

open_bzip2( name )
char *name;
{
	int rc;
	int save_size;
	unsigned char buf[4];
#ifdef	ALLOW_TDWR
	short int tdwr_value;
	short int *ptr;
#endif

/*
 *	Make sure we are a bzipped file.
 */

	rc = read_from_memory( buf, 3 );

	if ( rc != 3 )
		return( 0 );

/*
 *	We'll ignore the fourth byte, which should be a number.
 */

	if ( buf[0] != 'B' && buf[1] != 'Z' && buf[2] != 'h' )
	{
		printf("open_bzip2:  '%s': is not in bzip2 format.\n",
			name);
		return( 0 );
	}

/*
 *	Make sure we are closed.
 */

	if ( fbz != NULL )
	{
#ifndef	BZIP2_OLD
		BZ2_bzclose( fbz );
#else
		bzclose( fbz );
#endif
		fbz = NULL ;
	}

	rc = seek_memory( 0 );

	if ( rc <= 0 )
		return( 0 );

/*
 *	Even if we're going thru the mmap() routines, we still open the
 *	original file name.  This will make the reader happy.
 */

#ifndef	BZIP2_OLD
	fbz = BZ2_bzopen( name, "r" );
#else
	fbz = bzopen( name, "r" );
#endif

	if ( fbz == NULL )
	{
		printf("open_bzip2:  '%s':  unable to open file\n", name);
		perror( name );
		return( -1 );
	}

	bytes_read = 0;
	bytes_available = INTERNAL_SIZE;

/*
 *	Can't use read_from_memory() here, as we need the data decoded.
 */

	rc = read_bzip2( (char *) internal_buf, INTERNAL_SIZE );

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

	tdwr_value = htons( tdwr_value );

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

	bytes_read = 0;
	bytes_available = save_size;

	return( bytes_available );
}

/*
 *	Perform i/o.
 */

read_bzip2( mybuf, size )
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

	if ( bytes_read == bytes_available )
		return( 0 );

	request = size;

	if ( bytes_read + size > bytes_available )
		request = bytes_available - bytes_read;

/*
 *	Sanity.
 */

	if ( request <= 0 )				/* shouldn't happen */
	{
		printf("read_bzip2:  request is %d???\n", request);

#ifndef	BZIP2_OLD
		BZ2_bzclose( fbz );
#else
		bzclose( fbz );
#endif
		fbz = NULL;

		return( -1 );
	}

/*
 *	When BZ2_bzread is called, a lot of data is read in and buffered for
 *	later use.  When "use_internal" is set, it means we read in data to
 *	process a couple of bytes.  Since there isn't a "rewind" feature, we
 *	either close and reopen the file, which is horribly inefficient, or
 *	buffer up the four bytes we read, then return then now.
 */

	if ( use_internal )
	{
		if ( request < use_internal )
			copy_size = request;
		else
			copy_size = use_internal;
		bcopy( (char *) internal_ptr, mybuf, copy_size );
		request -= copy_size;
#ifndef	BZIP2_OLD
		user_bytes_read = BZ2_bzread( fbz, &mybuf[copy_size], request );
#else
		user_bytes_read = bzread( fbz, &mybuf[copy_size], request );
#endif
		user_bytes_read += copy_size;
		use_internal -= copy_size;
		if ( use_internal > 0 )
			internal_ptr = &internal_ptr[copy_size];
		else
			use_internal = 0;
	}
	else
#ifndef	BZIP2_OLD
		user_bytes_read = BZ2_bzread( fbz, mybuf, request );
#else
		user_bytes_read = bzread( fbz, mybuf, request );
#endif

	if ( user_bytes_read == 0 )				/* EOF */
	{
#ifndef	BZIP2_OLD
		BZ2_bzclose( fbz );
#else
		bzclose( fbz );
#endif
		fbz = NULL;
	}

	if ( bytes_read > bytes_available )
	{
		printf("read_bzip2:  bytes_read > bytes_available???\n");
		return( -1 );
	}

	bytes_read += user_bytes_read;

	return( user_bytes_read );
}

#else

#define	NO_CAN_DO	\
		printf("\nBZIP2 access is not built in this binary.\n"); \
		exit( 1 );

open_bzip2( name )
char *name;
{
	NO_CAN_DO;
}

read_bzip2( mybuf, size )
char *mybuf;
int size;
{
	NO_CAN_DO;
}

#endif
