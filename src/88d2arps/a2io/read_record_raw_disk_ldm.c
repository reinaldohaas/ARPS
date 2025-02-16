/*
 *	read_record_raw_disk_ldm.c	Read LDM data files.
 *
 *	Internal routines
 *	open_ldm()			Start the LDM process.
 *	read_ldm()			Perform i/o.
 *
 *	DANGER WILL ROBINSON.
 *	This set of routines assumes Level 2 data.  No other format is valid.
 *	Non-A2 data will generate errors, and maybe a core dump!
 */
 
#pragma ident "@(#)read_record_raw_disk_ldm.c	5.6	05/13/04 CAPS"

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <unistd.h>
#include <netinet/in.h>
#include "a2io.h"			/* Need some constants */

#include "config.h"

#include "proto.h"

#ifdef	ALLOW_BZIP2

#include "bzlib.h"

static int bytes_available = 0;
static int bytes_read = 0;

static int save_size = 0;

static int have_24_byte_hdr = 0;

#define	DO_READ	-99999

/*
 *	Size of actual data.
 *
 *	We have no way of knowing the size of the actual data without reading
 *	and decoding it.  We *can* determine if we are the right format, and
 *	if we have a 24-byte header record.  That is the best that we can do.
 */

/*
 *	Get the size of the file.  If we return 0, then it is an error state,
 *	and we get advanced to the next file.
 */

open_ldm( name )
char *name;
{
	int rc;
	int offset;
	unsigned char buf[31];

/*
 *	Read some data.
 */

	rc = read_from_memory( buf, 31 );

	if ( rc != 31 )
		return( 0 );

	if ( strncmp( (char *) buf, "ARCH", 4 ) == 0 ||
		strncmp( (char *) buf, "AR2V", 4 ) == 0 )
	{
		have_24_byte_hdr = 1;
		offset = 28;
	}
	else
		offset = 4;

/*
 *	Make sure we are a bzipped data.
 *
 *	We'll ignore the fourth byte, which should be a number.
 */

	if ( buf[offset] != 'B' || buf[offset + 1] != 'Z' ||
		buf[offset + 2] != 'h' )
	{
		printf("get_size_disk_ldm:  '%s' is not in LDM raw format.\n",
			name);
		return( 0 );
	}

/*
 *	Now let memory routines take over.
 */

	rc = seek_memory( 0 );

	if ( rc < 0 )
		return( rc );

/*
 *	24 byte header + 10000 records
 */

	save_size = (offset - 4 ) + REC_SIZE_TAPE * 10000;

	bytes_available = DO_READ;

	return( save_size );

}

/*
 *	Perform i/o.
 */

read_ldm( mybuf, size )
char *mybuf;
int size;
{
#define	BLOCK_MAX	3000000
	static char block[BLOCK_MAX], oblock[BLOCK_MAX];
	int length;
	unsigned int olength;
	int rc;

/*
 *	Sanity checks.
 */

	if ( size <= 0 )
		return( -1 );

	if ( bytes_available == 0 )
		return( 0 );

/*
 *	Handle some "impossible" events.
 */

	if ( have_24_byte_hdr )
	{
		if ( size >= A2_VS_24_BYTE_SIZE )
		{
/*
 *	Turn the flag off so we don't try to process it again.
 */
			have_24_byte_hdr = 0;
			rc = read_from_memory( (unsigned char * ) mybuf, 
				A2_VS_24_BYTE_SIZE );
			return( rc );
		}
		else
		{
	printf("read_disk_ldm:  %d bytes requested, should have been %d\n",
			size, A2_VS_24_BYTE_SIZE );
			printf("read_disk_ldm:  returning EOF\n");
			return( 0 );
		}
	}
	if ( size == A2_VS_24_BYTE_SIZE )
	{
	printf("read_disk_ldm:  %d bytes requested, should have been %d\n",
			size, REC_SIZE_TAPE );
			printf("read_disk_ldm:  returning EOF\n");
			return( 0 );
	}

/*
 *	Do we have any data in memory?
 */

	if ( bytes_available == DO_READ )
	{

		rc = read_from_memory( ( unsigned char *) &length, 4 );

		if ( rc != 4 )
		{
			bytes_available = 0;
			return( 0 );
		}

		length = htonl( length );

/*
 *	A negative length indicates the last block of the file.  We'll just
 *	flip the sign and let the rest of the code handle EOF issues.
 */

		if ( length < 0 )
			length = -length;

/*
 *	Check for data corruption.  If the data trips this, then there is a
 *	problem at the data origination site.
 */

		if ( length > BLOCK_MAX ) {
printf("read_disk_ldm:  Requesting %d bytes, though arrays are %d bytes.\n",
				length,BLOCK_MAX);
printf("read_disk_ldm:  The input data file is likely corrupt.  EOF forced.\n");
			bytes_available = 0;
			return( 0 );
		}

		rc = read_from_memory( (unsigned char *) block, length );

		if ( rc != length )
		{
			bytes_available = 0;
			return( 0 );
		}

		olength = sizeof( oblock );
#ifndef	BZIP2_OLD
		rc = BZ2_bzBuffToBuffDecompress( oblock, &olength, 
			block, length, 0, 0 );
#else
		rc = bzBuffToBuffDecompress( oblock, &olength, 
			block, length, 0, 0 );
#endif
		if ( rc != 0 )
		{
		printf("BZIP decompression problem, need %d, array size %d\n",
			olength,BLOCK_MAX);
			bytes_available = 0;
			return( 0 );
		}

		bytes_available = olength;
	}

	if ( bytes_read + size > bytes_available )
	{
	printf("read_disk_ldm:  bytes_read %d + size %d > bytes_available %d\n",
		bytes_read, size, bytes_available );
		bytes_available = 0;
		return( bytes_available );
	}

	bcopy( &oblock[bytes_read], mybuf, size );

	bytes_read += size;

	if ( bytes_read == bytes_available )
	{
		bytes_read = 0;
		bytes_available = DO_READ;
	}

	return( size );
}

#else

#define	NO_CAN_DO	\
		printf("\nBZIP2 access is not built in this binary.\n"); \
		exit( 1 );

open_ldm( name )
char *name;
{
	NO_CAN_DO;
}

read_ldm( mybuf, size )
char *mybuf;
int size;
{
	NO_CAN_DO;
}

#endif
