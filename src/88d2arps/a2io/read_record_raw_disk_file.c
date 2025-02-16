/*
 *	read_record_raw_disk_file.c	Perform i/o for disk data.
 *
 *	Internal
 *	radar_open_raw_disk_file()	Open routine for files.
 *	radar_close_raw_disk_file()	Close routine for files.
 *	read_record_raw_disk_file()	I/O routine for disk.
 *	file_seek_data()		LDM hook for file continuation.
 *	get_file_format()		Return file format.
 */
 
#pragma ident "@(#)read_record_raw_disk_file.c	5.5	01/20/04 CAPS"

#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>

#include "a2io.h"
#include "a2io_message_31.h"
#include "extern.h"
#include "ncdc.h"

#include "config.h"

#ifdef	ALLOW_TDWR
#include "tdwrio.h"
#include "tdwr_extern.h"
#endif

#include "proto.h"

static char *name = NULL;

static int time_to_open = 1;
static int rec_len;			/* length of record */
static int save_rec_len;
static int file_no_data = 1;

static int data_size = 0;
static int seek_offset = 0;

static int data_continues = 0;

static int ok_data;

static int size;

static int no_header = 0;

/*
 *	We want this, even if TDWR is disabled.  In this case, it is always
 *	zero.
 */

static int is_tdwr = 0;

/*
 *	IO readers, to simplify the reader.
 */

static int (*file_io_reader_init)();
static int (*file_io_reader)();

/*
 *	The function names exist, even if the features aren't compiled in.
 *	When called, the user will get a friendly message saying that the
 *	feature isn't compiled in.
 */

int open_bzip2(), read_bzip2();
int open_compress(), read_compress();
int open_gzip(), read_gzip();
int open_ldm(), read_ldm();
int open_mmap(), read_mmap();

/*
 *	Initialization routine.  The actual open() is postponed due to the
 *	needs of LDM, which doesn't want to re-read everything multiple
 *	times.
 */

radar_open_raw_disk_file( s )
char *s;
{
#define	MAX_BUF_LEN	512		/* large buffer, really too large */
	unsigned char buf[MAX_BUF_LEN];
	int n;
	int rc, rc1;
	int mid, mlen;

	file_no_data = 0;
	time_to_open = 1;
	data_continues = 0;

	if ( name != NULL )
		free( name );

	ok_data = 0;

	size = open_mmap( s );

/*
 *	User gets EOF on read_record_raw_disk_file() call.
 */

	if ( size <= 0 )
		return( -1 );

	ok_data = 1;

/*
 *	Save the name for later use.
 */

	n = strlen( s ) + 1;
	name = malloc( n );

	if ( name == NULL )
	{
		perror("radar_open_raw_disk_file:  malloc failure");
		exit( 1 );
	}

	strncpy( name, s , n );

	rc = read_mmap( buf, MAX_BUF_LEN );

	if ( rc <= 0 )
	{
		printf("radar_open_raw_disk_file:  no data???\n");
		return( -1 );
	}

	file_format = get_file_format( name, buf, rc );

	is_tdwr = 0;

#ifdef	ALLOW_TDWR
/*
 *	Since we have nibbled at the data, let's make the decision right now
 *	if we're TDWR format.
 */

	tdwr_pmsg  = (struct tdwr_msg_hdr *) buf;
	endian_tdwr_msg_hdr();
	mid = tdwr_pmsg->message_id;
	mlen = tdwr_pmsg->message_length;

	if (TDWR_DATA_RADIAL && mlen == TDWR_NORMAL_SIZE)
		is_tdwr = 1;
	if (mid == 0x2b02 && mlen == TDWR_2B02_SIZE)
		is_tdwr = 1;
	if (mid == 0x2c02 && mlen == TDWR_2C02_SIZE)
		is_tdwr = 1;
	if (mid == 0x4206 && mlen == TDWR_4206_SIZE)
		is_tdwr = 1;
#endif

	rc = seek_mmap( 0 );

/*
 *	Sanity.
 */

	if ( rc != size )
	{
		printf("read_open_raw_disk_file:  rc %d size %d\n", rc, size);
		return( -1 );
	}

	return( 0 );
}

/*
 *	Close routine.
 */

radar_close_raw_disk_file()
{
	close_mmap();

	return( 0 );
}

/*
 *	Perform i/o for the user.
 */

read_record_raw_disk_file( mybuf, mybuf_size )
char *mybuf;
int mybuf_size;
{
	int rc, rc1;
	int missing_24_byte = 0;
	int req;
	int mid;
	int alen = sizeof(struct a2header) + 12;

	if ( file_no_data )
	{
		printf("User program ignored 'end of data set'.  Exit.\n");
		exit( 1 );
	}


loop:

/*
 *	The actual open has to occur here, as the user may have done a
 *	seek_data() call after doing the radar_init() call (LDM only).
 */

	if ( time_to_open )
	{
/*
 *	Return EOF on first read if the open failed.
 */

		if ( !ok_data )
		{
			file_no_data++;
			return( RADAR_END_OF_DATA );
		}

/*
 *	Probably bug:  LDM format file that uses BZIP2 compression just
 *	won't work with this logic!
 */

		if ( file_format != FILE_MMAP )
			size = file_io_reader_init( name );

		else
		{
/*
 *	Format FILE_MMAP/FILE_LDM, we've already opened the file.  Just seek.
 *
 *	The user may get strange failures if file_seek_data() is called with
 *	non-LDM data.
 */
			if ( data_continues )
			{
				seek_offset = data_size;
				size = seek_mmap( seek_offset );
				data_size += size;
				data_continues = 0;
			}
			else
			{
				seek_offset = 0;
				data_size = size;
			}
		}

/*
 *	If the file is size zero, just return now.
 */

		if ( !size )
		{
			file_no_data++;
			return( RADAR_END_OF_DATA );
		}

/*
 *	Check size of file.  We'll allow for "out of sequence" NCDC header
 *	files and normal 88D header files.
 *
 *	Some compression schemes can't return the real size, so they fake it
 *	by returning a valid number after checking to see if a 24-byte header
 *	record is present.
 */

		rec_len = 0;

/*
 *	Make sure we're not TDWR, so we don't get in the wrong code!
 */

		if ( size == A2_ARCHIVE2_SIZE && ! is_tdwr )
			rec_len = A2_ARCHIVE2_SIZE;

		else if ( size == NCDC_SIZE && ! is_tdwr )
			rec_len = NCDC_SIZE;

		else if ( ( ( size - A2_VS_24_BYTE_SIZE ) % REC_SIZE_TAPE ) 
				== 0 && ! is_tdwr )
			rec_len = REC_SIZE_TAPE;

/*
 *	Check for LDM fault...no 24-byte header record.
 *	Only complain once per volume scan.
 */

		else if ( ( size % REC_SIZE_TAPE ) == 0 && ! is_tdwr )
		{
			if ( !seek_offset )
		printf("\nWARNING:  24-byte header record is missing!\n");
			missing_24_byte = 1;
			rec_len = REC_SIZE_TAPE;
		}

#ifdef	ALLOW_TDWR
		else if (is_tdwr)
		{
			rc1 = file_io_reader( mybuf, 4 );

			if ( rc1 <= 0 )
			{
				file_no_data++;
				return(RADAR_END_OF_DATA);
			}

			tdwr_pmsg  = (struct tdwr_msg_hdr *) mybuf;
			endian_tdwr_msg_hdr();
			mid = tdwr_pmsg->message_id;
			rec_len = tdwr_pmsg->message_length;

/*
 *	Only these have to be altered.
 */

			if (mid == 0x2b00)
				rec_len = TDWR_2B00_SIZE;
			if (mid == 0x2b01)
				rec_len = TDWR_2B01_SIZE;
	
/*
 *	Undo, as it will be redone later.
 */
			endian_tdwr_msg_hdr();
			rc = file_io_reader( &mybuf[4], rec_len-4 );

			if ( rc <= 0 )
			{
				return(RADAR_END_OF_DATA);
			}
			return(rec_len);
		}
#endif
/*
 *	Checking the size of a file will get us into trouble if not all
 *	tilts use the same message id.  Do this right.
 */
#ifdef	BAD_METHOD
/*
 *	TDWR Data
 */
		else if ( ( size % TDWR_NORMAL_SIZE ) == 0 )
			rec_len = TDWR_NORMAL_SIZE;
		else if ( ( size % TDWR_2B00_SIZE ) == 0 )
			rec_len = TDWR_2B00_SIZE;
		else if ( ( size % TDWR_2B01_SIZE ) == 0 )
			rec_len = TDWR_2B01_SIZE;
/*
 *	TDWR Tape Header File
 */
		else if ( size == TDWR_2B02_SIZE )
			rec_len = TDWR_2B02_SIZE;
#endif

		else
		{
/*
printf("ELSE!\n");
 */
/*
 *	Message type 31 uses variable length records, so we can no longer
 *	determine if the file size is correct.
 */
#ifdef	PRE_BUILD10
			printf("\nWARNING:  can't determine true size of file '%s'\n",name);
			printf("WARNING:  data may be corrupt, incomplete, or non-A2 data\n");
			printf("WARNING:  it will be assumed to be tape format data.\n\n");
#endif	/* PRE_BUILD10 */
			rec_len = REC_SIZE_TAPE;
		}

/*
 *	Sanity.
 */

		if ( rec_len > mybuf_size )
		{
			printf("read_record_raw_disk_file:  need %d bytes for buffer, but only %d bytes allow\n",
			rec_len, mybuf_size );
			printf("read_record_raw_disk_file:  fatal error\n");
			exit( 1 );
		}

/*
 *	If we are a header file, we only have one record, so we read it and
 *	return.
 */

		if ( rec_len == A2_ARCHIVE2_SIZE || rec_len == NCDC_SIZE )
		{
			rc = file_io_reader( mybuf, rec_len );

			if ( rc <= 0 )		/* EOF */
			{
				file_no_data++;
				return( RADAR_END_OF_DATA );
			}
		}

		time_to_open = 0;
		save_rec_len = rec_len;
		req = A2_VS_24_BYTE_SIZE;


/*
 *	We might be missing the 24-byte header record, or perhaps we've
 *	been told to skip over it (LDM data files).
 */

		if ( missing_24_byte )
			req = save_rec_len;

#ifdef	ALLOW_TDWR
	if (is_tdwr){
/*
 *	Nibble.
 */
		rc1 = file_io_reader( mybuf, 4);

			if ( rc1 <= 0 )
			{
				file_no_data++;
				return( RADAR_END_OF_DATA );
			}

			req = message_length;
		}
#endif

		rc = file_io_reader( mybuf, req );

		if ( rc <= 0 )
		{
			file_no_data++;
			return( RADAR_END_OF_DATA );
		}

/*
 *	Make sure we're not missing the header record.
 */
		if ( strncmp( mybuf, "ARCHIVE2", 8 ) == 0 ||
			strncmp( mybuf, "AR2V", 4 ) == 0 )
		{
			no_header = 0;
			return( rc );
		}
		else
			no_header = 1;
	}

/*
 *	Sanity.
 */

	data_continues = 0;

	if (! is_tdwr )
	{
		if (!no_header) {
			rc1 = file_io_reader( mybuf, alen );
		} else {
			no_header = 0;
			rc1 = file_io_reader( &mybuf[A2_VS_24_BYTE_SIZE],
				alen-A2_VS_24_BYTE_SIZE);
			rc1 += A2_VS_24_BYTE_SIZE;
		}
	} else {
		rc1 = file_io_reader( mybuf, rec_len);
	}

	if ( rc1 > 0 && ! is_tdwr )
	{
		nssl_a2header = (struct a2header *) mybuf;
		endian_header();
		if ( nssl_a2header->message_type == A2_DIGITAL_TYPE )
		{
			rec_len = 2 * nssl_a2header->message_size + 12;
			rc = file_io_reader( &mybuf[alen], rec_len - alen );
#ifdef	DEBUG
nssl_a2data31 = (struct a2data31 *) &mybuf[sizeof(struct a2header)];
endian_data31();
printf("%4.4s %d %f %d %d %d %d %d %f\n",
nssl_a2data31->radar_id,
nssl_a2data31->azimuth_number,
nssl_a2data31->azimuth,
nssl_a2data31->compression_indicator,
nssl_a2data31->radial_length,
nssl_a2data31->azimuth_resolution_spacing,
nssl_a2data31->radial_status,
nssl_a2data31->elevation_number,
nssl_a2data31->elevation);
/*
printf("%4.4s %u %d %d %f\n",
nssl_a2data31->radar_id,
nssl_a2data31->zulu_time,
nssl_a2data31->mod_julian_date,
nssl_a2data31->azimuth_number,
nssl_a2data31->azimuth);
 */
printf("data block count [4-9] %d ",nssl_a2data31->data_block_count);
printf("pointers: %d %d %d %d %d %d %d %d %d\n",
	nssl_a2data31->vol_ptr,
	nssl_a2data31->elv_ptr,
	nssl_a2data31->rad_ptr,
	nssl_a2data31->ref_ptr,
	nssl_a2data31->vel_ptr,
	nssl_a2data31->spw_ptr,
	nssl_a2data31->zdr_ptr,
	nssl_a2data31->phi_ptr,
	nssl_a2data31->rho_ptr);
endian_data31();
#endif

		} else {
			rc = file_io_reader( &mybuf[alen],
				REC_SIZE_TAPE - alen );
		}
		rc += alen;
/*
 *	Reverse the action, because we'll do it again later.
 */
		endian_header();
	} else {
		rc = rc1;
	}

	if ( rc <= 0 )
	{
		file_no_data++;
		return( RADAR_END_OF_DATA );
	}
	return( rc );
}

/*
 *	We are an LDM file that has been partially processed.  Set a flag
 *	so that when the actual i/o occurs, we won't get lost looking for
 *	a 24-byte header record.
 *
 *	We can't check to see if the data is gzipped yet, as the data format
 *	doesn't get checked until the first read occurs, though be definition,
 *	no one should be trying this.
 *
 *	This function is designed for the "ldm_ridds" program, so it isn't
 *	considered a user function.
 *
 */

file_seek_data()
{
/*
 *	Sanity.
 *
 *	We can't do this for compressed files.
 */

	if ( file_format != FILE_MMAP && file_format != FILE_LDM )
	{
		printf("seek_data:  it is not possible to use this call with compressed data\n");
		exit( 1 );
	}
	data_continues = 1;

	return( 0 );
}

/*
 *	Recognize the disk format, and return the information.  This format
 *	is how the data is written to disk, as in GZIP, BZIP2, etc, not as
 *	in A2 vs TDWR.
 *
 *	We make no determination of the validity of the data, we just identify
 *	what kind it is.
 */

/*
 *	We'll ignore the fourth byte, which should be a number.
 */

#define	bzip2_ok( offset_val ) \
		buf[ offset_val ] == 'B' && \
		buf[ offset_val + 1 ] == 'Z' && \
		buf[ offset_val + 2 ] == 'h'

get_file_format( my_name, buf, len )
char *my_name;
unsigned char *buf;
int len;
{
	int offset;
	int n;

/*
 *	Look for GZIP and Compressed format.
 */

	if ( len >= 2 )
	{
		if ( buf[0] == 037 && buf[1] == 0213 )
		{
			file_io_reader_init = open_gzip;
			file_io_reader = read_gzip;
			return( FILE_GZIP );
		}

		if ( buf[0] == 037 && buf[1] == 0235 )
		{
			file_io_reader_init = open_compress;
			file_io_reader = read_compress;
			return( FILE_COMPRESS );
		}
	}

/*
 *	Look for BZIP2 format.
 */

	if ( len >= 3 )
	{
		if ( bzip2_ok( 0 ) )
		{
			file_io_reader_init = open_bzip2;
			file_io_reader = read_bzip2;
			return( FILE_BZIP2 );
		}
	}

/*
 *	Look for LDM, which contains BZIP2 records.
 */

	if ( len >= 31 )
	{
		if ( strncmp( (char *) buf, "ARCH", 4 ) == 0 ||
			strncmp( (char *) buf, "AR2V", 4 ) == 0 )
			offset = 28;
		else
			offset = 4;

		if ( bzip2_ok( offset ) )
		{
			file_io_reader_init = open_ldm;
			file_io_reader = read_ldm;
			return( FILE_LDM );
		}
	}

/*
 *	Put in some sanity.
 */

	n = strlen( my_name );

	if ( n >= 3 )
	{
		if ( strcmp( &my_name[n-3], ".gz" ) == 0 )
			printf("WARNING:  '%s' is not a gzip file!\n", my_name );
	}

	if ( n >= 4 )
	{
		if ( strcmp( &my_name[n-4], ".bz2" ) == 0 )
			printf("WARNING:  '%s' is not a bzip2 file!\n", my_name );
		if ( strcmp( &my_name[n-4], ".ldm" ) == 0 )
			printf("WARNING:  '%s' is not an LDM file!\n", my_name );
	}

	if ( n >= 6 )
	{
		if ( strcmp ( &my_name[n-6], ".ridds" ) == 0 )
			printf("WARNING:  '%s' is not an LDM file!\n", my_name );
	}

	if ( n >= 2 )
	{
		if ( strcmp( &my_name[n-2], ".Z" ) == 0 )
			printf("WARNING:  '%s' is not a compressed file!\n",
				my_name );
	}

	file_io_reader_init = open_mmap;
	file_io_reader = read_mmap;

	return( FILE_MMAP );
}
