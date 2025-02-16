/*
 *	read_record_raw_disk_gzip.c	Read disk data that is gzipped.
 *
 *	Internal routines
 *	open_gzip()			Start the GZIP process.
 *	read_gzip()			Perform i/o.
 */
 
#pragma ident "@(#)read_record_raw_disk_gzip.c	4.4	03/07/01 NSSL"

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "config.h"

#include "proto.h"

#ifdef	ALLOW_GZIP

#include "zlib.h"


static int bytes_available = 0;
static int bytes_read = 0;

static int save_size = 0;

static gzFile fgz;

/*
 *	Retrieve the size of the original data.
 *
 *	The code to find the file size when not gzip'd came from the source
 *	code for version 1.2.4.  The documentation files are not included in
 *	this source.  Reference files can be found in the GNU distribution
 *	available on "prep.ai.mit.edu".  Code is in subdirectories of /pub/gnu.
 *
 *	Some info from the above location is included below.
 */

/* gzip (GNU zip) -- compress files with zip algorithm and 'compress' interface
 * Copyright (C) 1992-1993 Jean-loup Gailly
 * The unzip code was written and put in the public domain by Mark Adler.
 * Portions of the lzw code are derived from the public domain 'compress'
 * written by Spencer Thomas, Joe Orost, James Woods, Jim McKie, Steve Davies,
 * Ken Turkowski, Dave Mack and Peter Jannesen.
 *
 * See the license_msg below and the file COPYING for the software license.
 * See the file algorithm.doc for the compression algorithms and file formats.
 */

/*
static char  *license_msg[] = {
"   Copyright (C) 1992-1993 Jean-loup Gailly",
"   This program is free software; you can redistribute it and/or modify",
"   it under the terms of the GNU General Public License as published by",
"   the Free Software Foundation; either version 2, or (at your option)",
"   any later version.",
"",
"   This program is distributed in the hope that it will be useful,",
"   but WITHOUT ANY WARRANTY; without even the implied warranty of",
"   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the",
"   GNU General Public License for more details.",
"",
"   You should have received a copy of the GNU General Public License",
"   along with this program; if not, write to the Free Software",
"   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.",
0};
 */
typedef unsigned char  uch;
typedef unsigned short ush;
typedef unsigned long  ulg;

/* Macros for getting two-byte and four-byte header values */
#define SH(p) ((ush)(uch)((p)[0]) | ((ush)(uch)((p)[1]) << 8))
#define LG(p) ((ulg)(SH(p)) | ((ulg)(SH((p)+2)) << 16))

/*
 *	Initialize things for Gzip format.  The data size is returned.
 *	If zero, then it is an error, and we get advanced to the next file.
 */

open_gzip( name )
char *name;
{
	int size;
	int rc;

	unsigned char buf[4];

/*
 *	Make sure we are a gzipped file.
 */

	rc = read_from_memory( buf, 2 );

	if ( rc != 2 )
		return( 0 );

	if ( buf[0] != 037 && buf[1] != 0213 )
	{
		printf("open_gzip:  '%s':  is not in gzipped format.\n", name);
		return( 0 );
	}

	size = seek_memory( 0 );

	if ( size < 0 )
		return( 0 );

	rc = seek_memory( size - 4 );

	if ( rc < 4 )
		return( 0 );

	rc = read_from_memory( buf, 4 );

	if ( rc != 4 )
		return( 0 );

	save_size = LG( buf );


/*
 *	Make sure we are closed.  This is necessary if the user does something
 *	out of sequence, and we are left with unused/abused memory.
 */

	if ( fgz != NULL )
	{
		gzclose( fgz );
		fgz = NULL;
	}

	rc = seek_memory( 0 );

	if ( rc <= 0 )
		return( 0 );

/*
 *	Even if we're going thru the mmap() routines, we still open the
 *	original file name.
 */

	fgz = gzopen( name, "r" );

	if ( fgz == NULL )
	{
		printf("open_gzip:  '%s':   unable to open file\n", name);
		perror( name );
		return( -1 );
	}

	bytes_read = 0;
	bytes_available = save_size;

	return( bytes_available );
}

/*
 *	Perform i/o.
 */

read_gzip( mybuf, size )
char *mybuf;
int size;
{
	int request;
	int user_bytes_read;

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
		printf("read_gzip:  request is %d???\n", request);
		gzclose( fgz );
		fgz = NULL;

		return( -1 );
	}

	user_bytes_read = gzread( fgz, mybuf, request );

	if ( user_bytes_read == 0 )				/* EOF */
	{
		gzclose( fgz );
		fgz = NULL;
	}

	if ( user_bytes_read > bytes_available )
	{
		printf("read_gzip:  bytes_read > bytes_available???\n");
		return( -1 );
	}

	bytes_read += user_bytes_read;

	return( user_bytes_read );
}

#else	/* ALLOW_GZIP */

#define	NO_CAN_DO	\
	printf("\nGZIP access is not built into this binary.\n"); \
	exit( 1 );

open_gzip( name )
char *name;
{
	NO_CAN_DO;
}

read_gzip( mybuf, size )
char *mybuf;
int size;
{
	NO_CAN_DO;
}

#endif	/* ALLOW_GZIP */
