/*
 *	read_record_tape()	Reads a radial.  Return value is the message
 *				type.  Note that there are many local
 *				additions.
 *	Internal routines
 *	read_record_tape()	Read a tape or disk record.
 *	handle_tape_change()	Handle the tape change.
 *
 *	Note...although the word "tape" appears in the filename, this really
 *	only handles disk files.  The original version was designed to read
 *	tapes.
 */

#pragma ident "@(#)read_record_tape.c	5.6	03/12/03 CAPS"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "config.h"

#include "a2io.h"
#include "a2io_message_31.h"
#include "const.h"
#include "extern.h"
#include "ncdc.h"

#ifdef	ALLOW_TDWR
#include "tdwrio.h"
#include "tdwr_extern.h"
#endif

#include "vcp.h"

#include "proto.h"

static int force_site_type_record = 0;

static int out_of_seq_201_record = 0;

static unsigned short int save_julian_date = 0;
static int save_milsec = 0;

static int tapes_processed = 0;

static int trigger = 0;

static int eod = 0;

#ifdef ALLOW_PKU
static int aus_force_vs_record = 1;
#else
static int aus_force_vs_record = 0;
#endif

static int dbz_type = -999;
static int vel_type = -999;
static int spw_type = -999;

/*
 *	Although the variable name is "format_is_tape", this release no
 *	longer supports tape data.
 */

read_record_tape( format_is_tape )
int format_is_tape;
{
/*
int fd, rc;
 */
	char radar_name[16];		/* user requested radar name */
	int n = 0;			/* number of bytes read */

	static int scan_flag = 0;	/* error counter looking for scan # */
	struct a2tape24 fake_header;
	static int handled_bug = 0;
	int rec_size;

	int i, k;
	int32_t *pdbc;
/*
	int dbc[10];
 */

	if (nssl_a2data != NULL )
	{
	if ( nssl_a2data->radial_status == 4 && aus_vs_bug && !handled_bug )
		aus_force_vs_record++;
	}
/*
 *	Fix bug on Australia tapes.
 */
	if ( aus_force_vs_record )
	{
		aus_force_vs_record = 0;
		handled_bug = 1;
		nssl_scan_num++;

		sprintf( fake_header.archive2, "ARCHIVE2.");
		sprintf( fake_header.vs, "%03d", nssl_scan_num );
		fake_header.not_used = 0;
		fake_header.julian_date = 0;
		fake_header.milsec = 0;
		bcopy( (char *) &fake_header, nssl_buf, A2_VS_24_BYTE_SIZE );
		message_type = A2_VS_24_BYTE;
		message_length = A2_VS_24_BYTE_SIZE;
		return( message_type );
	}

	handled_bug = 0;

/*
 *	Start loop, as we may have to reread...
 */

reread:	;

/*
 *	88D on tape is usually blocked, meaning that a read will return more
 *	than one record.  If left isn't zero, then that means that there is
 *	more data to unblock.
 *
 *	TDWR data is blocked also, but at a different rate.
 *
 *	For realtime, there is one record per record.
 *
 *	A previous version of this routine had to worry about i/o errors and
 *	EOF's.  This is now all handled in read_record_raw_tape().
 */

	if ( left == 0 )				/* need to do i/o */
	{
		rec_size = REC_SIZE_TAPE;
		if (file_format == FILE_MMAP || disk_format == DISK_TDWR)
			rec_size = REC_SIZE_GENERIC;

		n = read_record_raw_disk( nssl_buf, rec_size);

		if ( n == RADAR_END_OF_DATA )
		{
			eod = 1;

			if ( eod == 0 )
			{
				left = 0;
				goto reread;
			}
			message_type = n;
			message_length = 0;
			return( message_type );
		}

		left = n;
		new = 0;
	}

	nssl_ptr = &nssl_buf[new];

/*
 *	Set pointer to current data location and figure out our message type.
 *
 *	Size inconstancies are also handled here.
 */

	find_message_type( left, REC_SIZE_TAPE, REC_SIZE_TAPE );

/*
 *	Handle tape header and volume scan messages.
 */

	if ( message_type == A2_ARCHIVE2 )
	{
		adjust_a2_pointer( "Tape Header", message_length );
		return( message_type );
	}

	if ( message_type == A2_ARCHIVE2_NCDC )
	{
		get_radar_name( radar_name );
		nssl_dptr = &nssl_ptr[ 0 ];
		ncdc_header = ( struct ncdc_header *) nssl_dptr;

		if ( strncmp( radar_name, "UNKN", 4 ) == 0 )
		{
			printf("INFO:     setting radar name to '%4.4s'\n",
				ncdc_header->id );
			set_radar_name( ncdc_header->id );
			get_radar_info();
		}
		else
		{
			if ( strncmp( radar_name, ncdc_header->id, 4 ) != 0 )
			{
printf("WARNING:  User said radar name is '%4.4s'.  NCDC tape said radar name is '%4.4s'.\n",
				radar_name, ncdc_header->id );
			printf("WARNING:  Name '%s' is believed.\n",radar_name);
			}
		}

		adjust_a2_pointer( "Tape Header NCDC", message_length );
		return( message_type );
	}

	out_of_seq_201_record = 0;

/*
 *	Handle an unexpected event:
 *		A 201 (A2_VOLSCAN_TYPE on an 88D tape.  This can happen on
 *		early RIDDS tapes.
 */

	if ( message_type == A2_VOLSCAN_TYPE )
	{
/*
 *	This is ok if the data is in "realtime" format.
 */
printf("read_record:  found 201 record (realtime) instead of 24-byte header (tape)\n");
/*
 *	Save the previous volume scan number.
 */
		out_of_seq_201_record = nssl_scan_num;
		message_length = A2_VS_24_BYTE_SIZE;
		message_type = A2_VS_24_BYTE;	/* so we fall thru */

/*
 *	When all is done, A2_VS_24_BYTE_SIZE bytes will be skipped before the
 *	routine returns.  The catch is that at this point, we have garbage
 *	data that needs to be skipped.  We'll skip the excess now, with the
 *	last part skipped before the code returns.
 */

		adjust_a2_pointer( "Unexpected 201 record",
			left - A2_VS_24_BYTE_SIZE );
	}

/*
 *	NEED TO FIX THIS SO IT WILL WORK FOR TDWR ALSO!
 *
 *	I must have fixed it, as I don't see any problems.  :-)
 */

	if ( message_type == A2_VS_24_BYTE )
	{
		nssl_dptr	= &nssl_ptr[0];
		nssl_a2tape24	= ( struct a2tape24 *) nssl_dptr;
		nssl_scan_num	= atoi( nssl_a2tape24->vs );

		if ( out_of_seq_201_record )
		{
/*
 *	We don't believe the current volume scan number, because it was the
 *	realtime volume scan number, not the tape version, which aren't
 *	related.
 *
 *	People reading the data itself which have to handle this case.
 */

			nssl_scan_num	= out_of_seq_201_record;

/*
 *	Create the A2_VS_24_BYTE record.
 */

			sprintf( nssl_a2tape24->archive2, "ARCHIVE2." );
			sprintf( nssl_a2tape24->vs, "%03d", nssl_scan_num );

/*
 *	We don't have date/time info, so we'll use the last record as a
 *	"best guess".
 */

			nssl_a2tape24->julian_date = save_julian_date;
			nssl_a2tape24->milsec = save_milsec;
		}

/*
 *	The user may have requested a range of volume scans.  If so, handle
 *	that request.
 */

		if ( nssl_scan_num < a2_start_scan )
		{
			printf("Found volume scan %d\n", nssl_scan_num );
			nssl_scan_num = VS_SEARCH;
		}

		if ( nssl_scan_num != VS_SEARCH )
		{
			adjust_a2_pointer( "File Header", message_length );
			--a2_num_scan;
			eod = 0;
			if ( a2_num_scan < 0 )
				eod++;
			if ( a2_end_scan < nssl_scan_num )
				eod++;

/*
 *	Save the volume scan number if we're not at EOD.
 */
			if ( !eod )
				trigger = nssl_scan_num;

			if ( eod )
			{
				set_disk_eod();
				message_type = A2_END_OF_DATA;
				message_length = 0;
				return( message_type );
			}

/*
 *	We need to return the right "volume scan" record number.
 */

			return( message_type );
		}

/*
 *		Fall thru into search mode.
 */

	}

#ifdef	ALLOW_TDWR
/*
 *	nssl_scan_num gets set in set_tdwr_info().
 */

	if ( message_type >= TDWR_2B00 && message_type <= TDWR_2C02 )
		set_tdwr_info();
#endif

/*
 *	If we don't have a valid volume scan number yet, we need to go into
 *	a search mode.  We won't return data until we do have a number.
 *	The reason we don't want to return data, even if it available, is that
 *	we want algorithms to have a full volume scan of data.
 *
 *	For "realtime" only, if you have a "nexrad.cfg" file in the directory
 *	where the program is running from, and "LOCAL_STATUS_KEY" is set
 *	to 0x88d1, than the getvolid() will return a valid volume scan number,
 *	and data will be *immediately* available, even if it in the middle of
 *	a volume scan.
 *
 *	Volume scan searching really is a leftover from the tape days.  We'll
 *	now accept anything, so the code is now obsolete, though I'll leav it
 *	in for this release.
 */
//***************************zkf****************
#ifdef	OBSOLETE
    if (radar98 != 1){
      if ( nssl_scan_num == VS_SEARCH )
	    {

        /*
        *	If we are in volume scan search mode, don't complain.
        */
		    if ( a2_start_scan > 1 )
		    	scan_flag = 1;
		    if ( scan_flag == 0 )
		    	printf("read_record:  looking for valid scan number\n");
		    scan_flag++;
		    if ( scan_flag > WARN_FREQ )
		    	scan_flag = 0;
		    adjust_a2_pointer( "looking for vs number", message_length );
		    disk_advance_file();
		    goto reread;
	    }
    }
#endif
//******************************zkf*************

/*
 *	At this point, we now have a valid volume scan number, so we can now
 *	process data.
 */

	if ( message_type == A2_DATA_TYPE )
	{
		save_julian_date = nssl_a2header->julian_date;
		save_milsec = nssl_a2header->milsec;
		nssl_dptr = &nssl_buf[ new + a2size ];
		nssl_a2data = ( struct a2data *) nssl_dptr;
//*************************zkf******************
        if (radar98 == 1)
           swap_data();     /* Swap data for IBM with radar98 data only */
		    else
           endian_data1();
//*************************zkf******************


		decode_data();
	}

#ifdef	ALLOW_TDWR
	else if ( TDWR_DATA_RADIAL )
		decode_data_tdwr();
#endif
	else if ( message_type == A2_DIGITAL_TYPE )
	{
#ifdef	DEBUG
printf("\n\n\n");
#endif
		save_julian_date = nssl_a2header->julian_date;
		save_milsec = nssl_a2header->milsec;
		nssl_dptr = &nssl_buf[ new + a2size ];
		nssl_a2data31 = ( struct a2data31 *) nssl_dptr;
		endian_data31();

		decode_data31();

/*
 *	Pointer order isn't guaranteed, however, we've assumed an order.  We
 *	now have to unscramble them, but first, make a copy of them!
 */

		pdbc =  (int32_t *) &nssl_a2data31->vol_ptr;
		for(i=1; i<10; i++)
		{
			map31[i] = 0;
		}

		for(i=1; i<=nssl_a2data31->data_block_count; i++)
		{
#ifdef	DEBUG
printf("TYPE=%d *pdbc=%d",i,*pdbc);
#endif
/*
 *	Arbitrary bad number
 */

			if (*pdbc > 10000)
			{
				adjust_a2_pointer("corrupt dual pol info",-1);
				return(MISSING);
			}
#ifdef	DEBUG
/*
printf("TYPE:  %d/%4.4s\n",i,&nssl_dptr[*pdbc]);
 */
printf("  value=%4.4s\n",&nssl_dptr[*pdbc]);
#endif
			if ( strncmp(&nssl_dptr[*pdbc],"RVOL",4) == 0)
				map31[VOL_31] = *pdbc;
			if ( strncmp(&nssl_dptr[*pdbc],"RELV",4) == 0)
				map31[ELV_31] = *pdbc;
			if ( strncmp(&nssl_dptr[*pdbc],"RRAD",4) == 0)
				map31[RAD_31] = *pdbc;
			if ( strncmp(&nssl_dptr[*pdbc],"DREF",4) == 0)
				map31[REF_31] = *pdbc;
			if ( strncmp(&nssl_dptr[*pdbc],"DVEL",4) == 0)
				map31[VEL_31] = *pdbc;
			if ( strncmp(&nssl_dptr[*pdbc],"DSW",3) == 0)
				map31[SPW_31] = *pdbc;
			if ( strncmp(&nssl_dptr[*pdbc],"DZDR",4) == 0)
				map31[ZDR_31] = *pdbc;
			if ( strncmp(&nssl_dptr[*pdbc],"DPHI",4) == 0)
				map31[PHI_31] = *pdbc;
			if ( strncmp(&nssl_dptr[*pdbc],"DRHO",4) == 0)
				map31[RHO_31] = *pdbc;

			pdbc++;
		}
#ifdef	DEBUG
printf("SLOTTERS:  ");
for(k=1; k<10; k++)
printf("%d/%d  ",k,map31[k]);
printf("\n");
#endif

		if ( map31[VOL_31] )
		{
			nssl_a2data31_vol = (struct a2data31_vol *)
				&nssl_dptr[map31[VOL_31]];

			endian_data31_vol();

/*
printf("VCP:  %d\n",nssl_a2data31_vol->vcp);
printf("LAT/LONG:  %f %f\n",nssl_a2data31_vol->latitude,nssl_a2data31_vol->longitude);
 */
		}

		if ( map31[ELV_31] )
		{
			nssl_a2data31_elv = (struct a2data31_elv *)
				&nssl_dptr[map31[ELV_31]];

			endian_data31_elv();

		}

		if ( map31[RAD_31] )
		{
			nssl_a2data31_rad = (struct a2data31_rad *)
				&nssl_dptr[map31[RAD_31]];

			endian_data31_rad();

/*
printf("NYQ:  %d\n",nssl_a2data31_rad->nyquist);
 */
		}


/*
 *	There may be alignment issues with some of these variables, so we
 *	have to make a dedicated copy.
 */


		if ( map31[REF_31] )
		{
			bcopy((char *) &nssl_dptr[map31[REF_31]],
				(char *) &a2data31_dbz,
				sizeof(struct a2data31_data));
			nssl_a2data31_dbz = &a2data31_dbz;

/*
printf("REF_PTR=%d\n",nssl_a2data31->ref_ptr);
printf("Data Block Type %c\n",nssl_a2data31_dbz->data_block_type);
printf("ng %d\n",nssl_a2data31_dbz->number_of_gates);
printf("Data name %3.3s\n",nssl_a2data31_dbz->data_name);
printf("Scale/Offset %f/%f\n",nssl_a2data31_dbz->scale,nssl_a2data31_dbz->offset);
 */

			endian_data31_data(0);
		}

		if ( map31[VEL_31] )
		{
			bcopy((char *) &nssl_dptr[map31[VEL_31]],
				(char *) &a2data31_vel,
				sizeof(struct a2data31_data));
			nssl_a2data31_vel = &a2data31_vel;

			endian_data31_data(1);
		}

		if ( map31[SPW_31] )
		{
			bcopy((char *) &nssl_dptr[map31[SPW_31]],
				(char *) &a2data31_spw,
				sizeof(struct a2data31_data));

			nssl_a2data31_spw = &a2data31_spw;

			endian_data31_data(2);

		}

		if ( map31[ZDR_31] )
		{
			bcopy((char *) &nssl_dptr[map31[ZDR_31]],
				(char *) &a2data31_zdr,
				sizeof(struct a2data31_data));

			nssl_a2data31_zdr = &a2data31_zdr;

			endian_data31_data(3);
		}

		if ( map31[PHI_31] )
		{
			bcopy((char *) &nssl_dptr[map31[PHI_31]],
				(char *) &a2data31_phi,
				sizeof(struct a2data31_data));

			nssl_a2data31_phi = &a2data31_phi;

			endian_data31_data(4);
		}

		if ( map31[RHO_31] )
		{
			bcopy((char *) &nssl_dptr[map31[RHO_31]],
				(char *) &a2data31_rho,
				sizeof(struct a2data31_data));

			nssl_a2data31_rho = &a2data31_rho;

			endian_data31_data(5);
		}

	}

	adjust_a2_pointer( "done decoding data", message_length );

	return( message_type );

}
