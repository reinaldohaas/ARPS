/*
 *	compress.c			Allow the reading of compressed files.
 *
 *	Internal routines
 *	compress_handle_open()		Open data file.
 *	compress_handle_read()		Perform I/O.
 *
 *	Unlike "gzip" and "bzip2", there are no public domain packages to
 *	read data directory.  Most of this code, and routines listed below,
 *	are taken from the version 4.1 of the Public Domain implementation
 *	of compress.  The rest of the code is the hooks into A2IO.
 *
 *	At the current time, only NCDC ftp files are compressed.
 *
 *	Compress internal routines
 *	decompress()		Decompression routine.
 *	getcode()		Perform I/O.
 */

#pragma ident "@(#)compress.c	4.6	03/07/01 NSSL"

#  define BITS 16

# define HSIZE	69001		/* 95% occupancy */

typedef long int	code_int;

typedef long int	  count_int;

 typedef	unsigned char	char_type;
static char_type magic_header[] = { "\037\235" };	/* 1F 9D */

/* Defines for third byte of header */
#define BIT_MASK	0x1f
#define BLOCK_MASK	0x80
/* Masks 0x40 and 0x20 are free.  I think 0x20 should mean that there is
   a fourth header byte (for expansion).
*/
#define INIT_BITS 9			/* initial number of bits/code */

#include <stdio.h>
#include <ctype.h>
#include <fcntl.h>
#include <strings.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "config.h"

/*
 *	Get prototypes unique to this file.
 */

#define	COMPRESS_UNIQUE

#include "proto.h"

static int n_bits;			/* number of bits/code */
static int maxbits = BITS;		/* user settable max # bits/code */
static code_int maxcode;		/* maximum code, given n_bits */
static code_int maxmaxcode = 1 << BITS;	/* should NEVER generate this code */
# define MAXCODE(n_bits)	((1 << (n_bits)) - 1)

static count_int htab [HSIZE];
static unsigned short codetab [HSIZE];
#define codetabof(i)	codetab[i]

/*
 * To save much memory, we overlay the table used by compress() with those
 * used by decompress().  The tab_prefix table is the same size and type
 * as the codetab.  The tab_suffix table needs 2**BITS characters.  We
 * get this from the beginning of htab.  The output stack uses the rest
 * of htab, and contains characters.  There is plenty of room for any
 * possible stack (stack used to be 8000 characters).
 */

#define tab_prefixof(i)	codetabof(i)
# define tab_suffixof(i)	((char_type *)(htab))[i]
# define de_stack		((char_type *)&tab_suffixof(1<<BITS))

static code_int free_ent = 0;		/* first unused entry */

static code_int getcode();

/*
 * block compression parameters -- after all codes are used up,
 * and compression rate changes, start over.
 */
static int block_compress = BLOCK_MASK;
static int clear_flg = 0;
/*
 * the next two codes should not be changed lightly, as they must not
 * lie within the contiguous general code space.
 */ 
#define FIRST	257	/* first free entry */
#define	CLEAR	256	/* table clear output code */

static int fd = -1;

#define	M	2000

static unsigned char compress_buf[M];

static int left = 0;
static int loc = 0;

#undef	DEBUG
#ifdef	DEBUG
main()
{
	int rc;

	rc = compress_handle_open( "q.Z" );

	for(;;)
	{
		rc = decompress( 1, compress_buf, M );
		if ( rc <= 0 )
			break;
		write(1,compress_buf,rc);
	}
}
#endif

compress_handle_open( name )
char *name;
{
	unsigned char dummy;
	unsigned char magic[2];
	int rc;

/*
 *	Tell the software to "reset" itself so we don't get old data from the
 *	internal buffer.
 */

	clear_flg = 1;

/*
 *	In case the user didn't process the entire file.
 */

	if ( fd > 0 )
	{
		close( fd );
		fd = -1;
	}

	fd = open( name, 0 );

	rc = read_from_memory( magic, 2 );

	if ( rc != 2 )
	{
		printf("compress_handle_open:  %s:  unexpected EOF\n", name);
		close( fd );
		fd = -1;
		return( -1 );
	}

	if ( magic[0] != magic_header[0] || magic[1] != magic_header[1] )
	{
		printf("compress_handle_open:  %s:  not in compressed format\n",
			 name);
		close( fd );
		fd = -1;
		return( -1 );
	}

	rc = read_from_memory( &dummy, 1 );

	if ( rc != 1 )
	{
		printf("compress_handle_open:  %s:  unexpected EOF\n", name );
		close( fd );
		fd = -1;
		return( -1 );
	}

	clear_flg = 1;

	maxbits = dummy;

	block_compress = maxbits & BLOCK_MASK;
	maxbits &= BIT_MASK;
	maxmaxcode = 1 << maxbits;

	if ( maxbits > BITS )
	{
		printf("compress_handle_open:  %s:  error:  maxbits=%d BITS=%d\n",
			name, maxbits, BITS);
		close( fd );
		fd = -1;
		return( -1 );
	}

	left = 0;
	rc = decompress( 0, &dummy, 1 );
	return( rc );
}

compress_handle_read( mybuf, mybuf_size )
char *mybuf;
int mybuf_size;
{
	int user_left, user_loc;
	int kount;
	int bytes_to_copy;
	int bytes_read;

	bzero( mybuf, (unsigned) mybuf_size );

	kount = 0;
	user_left = mybuf_size;
	user_loc = 0;
loop:
	if ( left == 0 )
	{
		bytes_read = decompress( 1, compress_buf, M );

		if ( bytes_read == 0 )
			return( kount );
		left = bytes_read;
		loc = 0;
	}
	else
		bytes_read = left;

	bytes_to_copy = user_left;
	if ( bytes_read < bytes_to_copy )
		bytes_to_copy = bytes_read;

	bcopy( &compress_buf[loc], &mybuf[user_loc], (unsigned) bytes_to_copy );
	left -= bytes_to_copy;
	loc += bytes_to_copy;

	user_left -= bytes_to_copy;
	user_loc += bytes_to_copy;

	kount += bytes_to_copy;

	if ( user_left != 0 )
		goto loop;

	return( kount );

}

static char_type rmask[9] = {0x00, 0x01, 0x03, 0x07, 0x0f, 0x1f, 0x3f, 0x7f, 0xff};

/*
 *  Decompress data.  Modified to handle local needs.
 */

decompress( flag, mybuf, mybuf_size )
int flag;
unsigned char *mybuf;
int mybuf_size;
{
    static char_type *stackp;
    static int finchar;
    static code_int code, oldcode, incode;
    static int state = 0;
	unsigned char *ptr;
	int kount;

/*
 *	Sanity.
 */

	if ( fd < 0 )
		return( 0 );

	if ( flag )
		goto loop;

	state = 0;

    /*
     * As above, initialize the first 256 entries in the table.
     */
    maxcode = MAXCODE(n_bits = INIT_BITS);
    for ( code = 255; code >= 0; code-- ) {
	tab_prefixof(code) = 0;
	tab_suffixof(code) = (char_type)code;
    }
    free_ent = ((block_compress) ? FIRST : 256 );

    finchar = oldcode = getcode();
    if(oldcode == -1)	/* EOF already? */
	return( 0 );			/* Get out of here */
/*
 *	We aren't returning any data, but we do want to set a "good" return
 *	code.
 */
	return( 1 );

loop:

	state++;

	if ( state == 1 )
	{
    /* first code must be 8 bits = char */
	mybuf[0] =  (char) finchar;
    stackp = de_stack;
		return( 1 );
	}

	code = getcode();

	if ( code == -1 )
	{
		close( fd );
		fd = -1;
		return( 0 );
	}

	if ( (code == CLEAR) && block_compress ) {
	    for ( code = 255; code >= 0; code-- )
		tab_prefixof(code) = 0;
	    clear_flg = 1;
	    free_ent = FIRST - 1;
	    if ( (code = getcode ()) == -1 )	/* O, untimely death! */
	    {
		close( fd );
		fd = -1;
		return( 0 );
	    }
	}
	incode = code;
	/*
	 * Special case for KwKwK string.
	 */
	if ( code >= free_ent ) {
            *stackp++ = (char_type) finchar;
	    code = oldcode;
	}

	/*
	 * Generate output characters in reverse order
	 */
	while ( code >= 256 ) {
	    *stackp++ = tab_suffixof(code);
	    code = tab_prefixof(code);
	}
	*stackp++ = finchar = tab_suffixof(code);

#ifdef	ORIG
	/*
	 * And put them out in forward order
	 */
	do
	    putchar ( *--stackp );
	while ( stackp > de_stack );
#endif
	kount = stackp - de_stack;
	if ( kount >  mybuf_size )
	{
		printf("compress_handle_read:  error:  need %d bytes, but only have room for %d\n",
			kount, mybuf_size );
		return( 0 );
	}

	ptr = mybuf;
	do
	{
		*ptr = *--stackp;
		ptr++;
	} while (stackp > de_stack );
	

	/*
	 * Generate the new entry.
	 */
	if ( (code=free_ent) < maxmaxcode ) {
	    tab_prefixof(code) = (unsigned short)oldcode;
	    tab_suffixof(code) = (char_type) finchar;
	    free_ent = code+1;
	} 
	/*
	 * Remember previous code.
	 */
	oldcode = incode;
	return( kount );
}

/*****************************************************************
 * TAG( getcode )
 *
 * Read one code from the standard input.  If EOF, return -1.
 * Inputs:
 * 	stdin
 * Outputs:
 * 	code or -1 is returned.
 */

code_int
getcode() 
{
    register code_int code;
    static int offset = 0, size = 0;
    static char_type buf[BITS];
    register int r_off, bits;
    register char_type *bp = buf;

    if ( clear_flg > 0 || offset >= size || free_ent > maxcode ) {
	/*
	 * If the next entry will be too big for the current code
	 * size, then we must increase the size.  This implies reading
	 * a new buffer full, too.
	 */
	if ( free_ent > maxcode ) {
	    n_bits++;
	    if ( n_bits == maxbits )
		maxcode = maxmaxcode;	/* won't get any bigger now */
	    else
		maxcode = MAXCODE(n_bits);
	}
	if ( clear_flg > 0) {
    	    maxcode = MAXCODE (n_bits = INIT_BITS);
	    clear_flg = 0;
	}
	size = read_from_memory( buf, n_bits );
	if ( size <= 0 )
	    return -1;			/* end of file */
	offset = 0;
	/* Round size down to integral number of codes */
	size = (size << 3) - (n_bits - 1);
    }
    r_off = offset;
    bits = n_bits;
	/*
	 * Get to the first byte.
	 */
	bp += (r_off >> 3);
	r_off &= 7;
	/* Get first part (low order bits) */
	code = (*bp++ >> r_off);
	bits -= (8 - r_off);
	r_off = 8 - r_off;		/* now, offset into code word */
	/* Get any 8 bit parts in the middle (<=1 for up to 16 bits). */
	if ( bits >= 8 ) {
	    code |= *bp++ << r_off;
	    r_off += 8;
	    bits -= 8;
	}
	/* high order bits. */
	code |= (*bp & rmask[bits]) << r_off;
    offset += n_bits;

    return code;
}
