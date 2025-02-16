/*
 *	User routines
 *	read_record()		Reads a data record.  Return value is the
 *				message type.  Note that there are many local
 *				additions.
 *	set_console_messages()	Allow or prevent 88D console messages.
 *
 *	Internal routines
 *	adjust_a2_pointer()	Internal routine to do pointer bookkeeping.
 *	check_seq_number()	Report on missing data records.
 *	read_record_init()	Internal initialization routine.
 *	find_message_type()	Internal routine to compute message type.
 *	display_console_message()
 *				Display 88D console message.
 *
 *	NOTE:  the return code from read_record() has changed dramatically
 *	       from version 2 of A2IO.
 */

#pragma ident	"@(#)read_record.c	5.9	08/30/04	CAPS"

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <time.h>
#include <netinet/in.h>

#include "config.h"

#include "a2io.h"
#include "a2io_message_31.h"
#include "const.h"
#include "extern.h"
#include "ncdc.h"
#include "read_record.h"

#ifdef	ALLOW_TDWR
#include "tdwrio.h"
#include "tdwr_struct.h"
#endif

#include "vcp.h"

#include "proto.h"

static int allow_console_messages = 0;

read_record()
{
	int rc = 0;

	if ( input_format == DISK )
		rc = read_record_tape( 0 );
	else
		bad_mode("read_record");

	check_seq_number( rc );

	return( rc );
}


/*
 *	Keep track of where we are in the data, so that the next time that
 *	read_record() is called, the right things happen.
 */

adjust_a2_pointer( msg, offset )
char *msg;
int offset;
{
	left -= offset;
	new += offset;

/*
 *	If "offset" is bad, make it fault.
 */

	if ( offset <= 0 )
		left = -1;

	if ( left < 0 )
	{
		printf("read_record:  logic fault: %s\n",msg);
		left = 0;
	}
}


/*
 *	Check to see if we've dropped any records.  This is a replacement
 *	for the code that used to be handled by read_radial() which didn't
 *	take into account non-data records.
 */

#define	ANY_SEQ		-1

check_seq_number( n )
int n;
{
	static int old_seq = ANY_SEQ;
	int seq;
	int diff;

/*
 *	If we aren't a valid 88D message (this includes local NSSL
 *	messages), just return.
 */

	if ( n < 1 )
		return;
	if ( n > 200 )
		return;

	seq = nssl_a2header->seq;

#define	TYPE_2		2					/* ??? */
#define	TYPE_3		3					/* ??? */
#define	TYPE_7		7					/* ??? */
#define	TYPE_13		13					/* ??? */
#define	TYPE_15		15					/* ??? */
#define	TYPE_18		18					/* ??? */
#define	TYPE_20		20				/* operator message */
#define	TYPE_21		21				/* operator message */
#define	TYPE_MPDA	22


/*
 *	Some message types have sequence numbers totally messed up.  Accept
 *	the sequence number and allow the next sequence number to be anything.
 */

	if ( message_type == TYPE_2 ||
		message_type == TYPE_3 ||
		message_type == TYPE_7 ||
		message_type == TYPE_13 ||
		message_type == TYPE_15 ||
		message_type == TYPE_18 ||
		message_type == TYPE_20 ||
		message_type == TYPE_21 ||
		message_type == TYPE_MPDA )
	{
		old_seq = ANY_SEQ;
		seq = ANY_SEQ;
	}

/*
 *	Message type 2 seems to have a sequence number of 0, at least on older
 *	88D tapes.  We'll correct it.
 */

	if ( message_type == TYPE_2 && seq == 0 )
	{
		old_seq = ANY_SEQ;
		seq = ANY_SEQ;
	}

	if ( old_seq != ANY_SEQ )
	{
		diff = seq - old_seq;
		if ( old_seq > seq )			/* wraparound */
			diff += 0x7fff;

/*
 *	We expect a difference of one, so account for that.
 */

		diff--;

/*
 *	Only complain if more than one is missing.  A missing record after
 *	a type 2 message seems to be common, so this might be an 88D bug.
 */

		if ( diff > 1 )
			printf("read_record:  %d data record(s) missing\n",diff );
	}

	old_seq = seq;
	if ( old_seq == 0x7fff )			/* reset */
		old_seq = 0;

	return;
}

/*
 *	Internal initializations.
 */

read_record_init()
{

/*
 *	Make pointers point to something so we don't core dump on program
 *	abuse (access data without a successful read_record() call).
 *	THE INFO RETURNED WILL BE GARBAGE!
 */

	nssl_ptr = &nssl_buf[0];
	nssl_a2header = ( struct a2header *) nssl_ptr;
	nssl_a2data = ( struct a2data *) nssl_ptr;
	nssl_a2data31 = ( struct a2data31 *) nssl_ptr;

/*
 *	Make gmt point to something useful (actually useless) so the program
 *	doesn't core dump if someone tries to access date info without going
 *	through a success read_record() call.
 */

	bzero( (char*) &nssl_gmt_data, sizeof(struct tm));
	nssl_gmt = &nssl_gmt_data;

/*
 *	Compute maximum elevation threshold in the units the data comes in.
 */

	elev_thresh = (unsigned short) 
		( (float) MAX_ELEV / (float) A_SCALE + 0.5);

/*
 *	Compute the size of some structures for later use.
 */

	a2size 	    = sizeof( struct a2header );
	a2size_scan = a2size + sizeof( struct a2scan );
	a2size_site = a2size + sizeof( struct a2site );

	return;
}


/*
 *	Find the message type and set the size of the data in bytes.
 *
 *	Return 0 if everything is ok, and -1 if there is a size inconsistency.
 *	When -1 is return, processing is terminated and the next i/o read
 *	occurs automatically.
 */

#define	A2_ARCHIVE2_NCDC_SIZE	NCDC_SIZE

/*
 *	n		The amount of unprocessed data.  For tapes, this could
 *			be multiple records.
 *	min_size	Absolute smallest acceptable size for a valid 88D data
 *			record.
 *	rec_size	Value to assign "message_length" for valid 88D messages.
 */

find_message_type( n, min_size, rec_size )
int n, min_size, rec_size;
{

/*
 *	Set possible data pointers.
 */

	nssl_ptr = &nssl_buf[ new ];
	nssl_a2header = (struct a2header *) nssl_ptr;

#ifdef	ALLOW_TDWR
	tdwr_pmsg = ( struct tdwr_msg_hdr *) nssl_ptr;
#endif

/*
 *	Records that begin "ARCHIVE2" are definitely some kind of tape
 *	header record.  Ditto for "AR2V####".
 */
	if ( strncmp( nssl_buf, "ARCHIVE2", 8 ) == 0 ||
		strncmp( nssl_buf, "AR2V", 4 ) == 0 )
	{
		endian_header();
/*
 *	Tape header records.
 */
		if ( n == A2_ARCHIVE2_SIZE )
		{
			message_type 	= A2_ARCHIVE2;
			message_length	= A2_ARCHIVE2_SIZE;
			return;
		}

		else if ( n == A2_ARCHIVE2_NCDC_SIZE )
		{
			message_type	= A2_ARCHIVE2_NCDC;
			message_length	= A2_ARCHIVE2_NCDC_SIZE;
			return;
		}

		else if ( n == A2_VS_24_BYTE_SIZE )
		{
			message_type   = A2_VS_24_BYTE;
			message_length = A2_VS_24_BYTE_SIZE;
			return;
		}

/*
 *	If we are here, something has gone very wrong.  Force a skip to the
 *	next record.
 */
		printf("find_message_type:  unidentifiable record that begins 'ARCHIVE2'\n");
		printf("find_message_type:  or 'AR2V'\n");
		message_type	= A2_UNKNOWN_TYPE;
		message_length	= n;
		return;
	}

/*
 *	Turn off, as A2 data sometimes trips this!
 */

#ifdef	ALLOW_TDWR

	endian_tdwr_msg_hdr();

/*
 *	Look for TDWR data.
 */

	if ( tdwr_pmsg->message_id == 0x2b00 )
	{
		message_type	= TDWR_2B00;
		message_length	= TDWR_2B00_SIZE;
		return;
	}

	if ( tdwr_pmsg->message_id == 0x2b01 )
	{
		message_type	= TDWR_2B01;
		message_length	= TDWR_2B01_SIZE;
		return;
	}

	if ( tdwr_pmsg->message_id == 0x2b02 )
	{
		message_type	= TDWR_2B02;
		message_length	= TDWR_2B02_SIZE;
		return;
	}

	if ( tdwr_pmsg->message_id == 0x2c00 )	/* CHECK THIS!!! */
	{
		message_type	= TDWR_2C00;
		message_length	= TDWR_2C00_SIZE;
		return;
	}

	if ( tdwr_pmsg->message_id == 0x2c01 )
	{
		message_type	= TDWR_2C01;
		message_length	= TDWR_2C01_SIZE;
		return;
	}

	if ( tdwr_pmsg->message_id == 0x2c02 )
	{
		message_type	= TDWR_2C02;
		message_length	= TDWR_2C02_SIZE;
		return;
	}

	if ( tdwr_pmsg->message_id == 0x4206 )
	{
		message_type	= 0x4206;
		message_length	= TDWR_4206_SIZE;
		return;
	}

/*
 *	A second call to endian_tdwr_msg_hdr() undoes what the first one did.
 */

	endian_tdwr_msg_hdr();
#endif

//**************************zkf********************
    if (radar98 == 1)
       swap_header();     /* Swap head for IBM with China radar98 data only  */
	else
       endian_header();

//**************************zkf********************

/*
 *	All we are left with is valid 88D messages.
 */

//***************************zkf***********************
    if(radar98==1)
        message_type = nssl_a2header->rda_channel;
    else
        message_type = nssl_a2header->message_type;

//***************************zkf***********************

	if ( message_type == A2_CONSOLE_MSG && allow_console_messages )
		display_console_message();

/*
 *	Digital message, with the real data.  Size is variables.
 */

	if ( message_type == A2_DIGITAL_TYPE )
	{
		message_length = 2 * nssl_a2header->message_size + 12;
		return;
	}

/*
 *	Check and see if we are exactly the record length of a LIVE or TAPE
 *	message.  If so, handle this as a special case, as we know the answer,
 *	and don't want to get confusing complaints.
 */


	if ( n == REC_SIZE_TAPE )
	{
		message_length = n;
		return;
	}


/*
 *	See if we have enough data, if not, complain to the user.
 */

	if ( n < min_size )
	{
printf("find_message_type:  product type %d, bytes expected %d bytes found %d\n",
		message_type, min_size, n );
/*
		printf("find_message_beep_type:  data is probably corrupt!\n");
 */
		printf("find_message_type:  data is probably corrupt!\n");
		message_length = n;
		return;
	}

/*
 *	If we've gotten this far, then we assign the record size that we were
 *	told to assign.
 */

	message_length = rec_size;
/*
printf("message_type=%d message_length=%d\n",message_type,message_length);
 */

	return;
}

/*
 *	Allow or disallow console messages.  The default (handled by
 *	radar_init() ) is to enable for realtime and disabled otherwise.
 */

set_console_messages( n )
int n;
{
	if ( n )
		allow_console_messages = 1;
	else
		allow_console_messages = 0;
	return;
}

/*
 *	Display 88D console messages.
 */

display_console_message()
{
	struct a2console *nssl_a2console;
	struct tm *tm;
#define	MYLEN	A2_MAX_CONSOLE_LEN + 1
	char mybuf[ MYLEN ];
	char *ptr;
	int len = 0;
	int warn = 0;
	long ztime;

/*
 *	Sanity, in the event a user called this routine themselves.
 */

	if ( message_type != A2_CONSOLE_MSG )
		return;

	nssl_a2console = ( struct a2console * ) &nssl_ptr[ a2size ];

	nssl_a2console->len = htons( nssl_a2console->len );

	len = nssl_a2console->len;

	if ( len > A2_MAX_CONSOLE_LEN )
	{
		len = A2_MAX_CONSOLE_LEN;
		warn = 1;
	}

/*
 *	Sanity.
 */

	if ( len <= 0 )
	{
		printf("display_console_message:  bad message length of %d\n",
			len );
		return;
	}

	bzero( mybuf, MYLEN );
	bcopy( nssl_a2console->msg, mybuf, len );

/*
 *	Clean up white space.
 */

	ptr = &mybuf[ len - 1];

	while( *ptr == ' ')
	{
		*ptr = (char) NULL;
		--ptr;
	}

/*
 *	Figure out the time/date of the message.  We replicate the logic for
 *	type 1 data messages here.
 */

/*
 *	Convert milli-seconds to seconds.
 */
 
	ztime = nssl_a2header->milsec;
	ztime /= 1000;

/*
 *	Add the number of days.
 */

	ztime += nssl_a2header->julian_date * SECONDS_PER_DAY;

/*
 *	For UNIX, time=0 is December 31, 1969, while it is January 1, 1970
 *	for NEXRAD systems, so we need to subtract off one day.
 */

	ztime -= SECONDS_PER_DAY;

/*
 *	Set a pointer to the data record.
 */

	tm = gmtime(&ztime);

	printf("\n%02d/%02d/%02d %02d:%02d:%02d - %s\n\n",
			tm->tm_mon+1, tm->tm_mday, tm->tm_year % 100,
			tm->tm_hour, tm->tm_min, tm->tm_sec,
			mybuf);

	if ( warn )
		printf("display_console_message:  WARNING:  console message truncated.\n");

	return;
}
