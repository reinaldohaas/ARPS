/*
 *	get_data_tdwr.c...		Data retrieval routines.
 *
 *	get_data_field_tdwr()		Returns floating point values.
 *	get_data_field_raw_tdwr()	Returns raw data values (short ints).
 *	get_data_field_uchar_tdwr()	Return raw data values (unsigned chars).
 *	get_scale_values_tdwr()		Recoding scale info for raw data values.
 *
 *	Internal routines
 *	set_tdwr_info()			TDWR bookkeeping structures
 *	decode_data_tdwr()		Decode TDWR data.
 *	print_header()			Print out the TDWR tape header record.
 *	get_data_field_raw_tdwr()	TDWR equivalent of get_data_field_raw()
 *	get_tdwr_radial_status()	Return "tdwr_radial_status".
 */
 
#pragma ident "@(#)get_data_tdwr.c	5.5	03/13/03 CAPS"

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <malloc.h>
#include <time.h>
#include <inttypes.h>

#include "config.h"

#ifdef	ALLOW_TDWR

#include "a2defs.h"
#include "extern.h"

#include "a2io.h"			/* to get TDWR_FORMAT */
#include "extern.h"
#include "tdwrio.h"
#include "tdwr_extern.h"

#include "proto.h"

static int tdwr_radial_status;

#define	MISSING_88D		0
#define	RANGE_FOLDED_88D	1

set_tdwr_info()
{
	int mloc, bloc;

	data_format = TDWR_FORMAT;

	funct_init_tdwr();

	mloc = sizeof( struct tdwr_msg_hdr );
	bloc = sizeof( struct tdwr_basedata_hdr ) + mloc;

/*
 *	We set pointers for everything, though not all will be used.
 */

	tdwr_pmsg  = (struct tdwr_msg_hdr *) nssl_ptr;
	tdwr_phdr  = (struct tdwr_basedata_hdr * ) &nssl_ptr[mloc];
	tdwr_p2b00 = (struct tdwr_prod_2b00 *) &nssl_ptr[bloc];
	tdwr_p2b01 = (struct tdwr_prod_2b01 *) &nssl_ptr[bloc];
	tdwr_p2b02 = (struct tdwr_prod_2b02 *) &nssl_ptr[mloc];

/*
 *	Handle endian problems for real data to prevent corruption.  Tdwr_pmsg
 *	is already handled, as we had to determine message type.
 */

	if ( TDWR_DATA_RADIAL )
	{
		endian_tdwr_basedata_hdr();
/*
 *	Volume scan info is in every data radial.  Only accept the value if
 *	the start of volume scan is flagged.
 */

		if ( tdwr_phdr->start_of_volume_scan )
			nssl_scan_num = tdwr_phdr->volume_scan_count;

		if ( tdwr_pmsg->message_id == 0x2b00 )
			endian_tdwr_2b00();

		else if ( tdwr_pmsg->message_id == 0x2b01 )
			endian_tdwr_2b01();

		else
	printf("set_tdwr_info:  funny message id of %x\n",tdwr_pmsg->message_id);
	}

	return( 0 );
}

decode_data_tdwr()
{
	static int tdwr_old_scan_num = -1;
	static int tdwr_old_tilt_num = -1;
	static int tdwr_old_vsn = -1;
	static int tdwr_rec = 0;
	long tmp;

	tdwr_rec++;

/*
	if ( tdwr_pmsg->message_id == 0x2b02 )
	{
		print_header();
		return( 0 );
	}
 */
#ifdef	DEBUG
	print_header();
#endif

	if ( !TDWR_DATA_RADIAL )
		return( 0 );

/*
 *	Process time.
 */

	tmp = (long) tdwr_phdr->timestamp;

	nssl_gmt = gmtime( (long*) &tmp);
	bcopy((char*) nssl_gmt, &nssl_gmt_data, sizeof(struct tm));
	nssl_gmt = &nssl_gmt_data;

/*
 *	Set the "status_scan" and "status_tilt" flags.
 */

	if ( tdwr_old_scan_num != tdwr_phdr->volume_scan_count  )
	{
		nssl_scan_status = 1;
		nssl_tilt_status = 1;
	}
	else if ( tdwr_old_tilt_num != tdwr_phdr->tilt_number )
	{
		nssl_scan_status = 0;
		nssl_tilt_status = 1;
	}
	else
	{
		nssl_scan_status = 0;
		nssl_tilt_status = 0;
	}

	if ( tdwr_old_scan_num != tdwr_phdr->volume_scan_count || 
		tdwr_old_tilt_num != tdwr_phdr->tilt_number )
	{
		tdwr_old_scan_num = tdwr_phdr->volume_scan_count;
		tdwr_old_tilt_num = tdwr_phdr->tilt_number;
	}

/*
 *	Set the equivalent of 88D "Radial Status", as some software may
 *	depend on the values being valid.
 */

	if ( tdwr_phdr->start_of_elevation_scan )
	{
		if ( tdwr_phdr->start_of_volume_scan )
			tdwr_radial_status = 3;
		else
			tdwr_radial_status = 0;
	}
	else if ( tdwr_phdr->end_of_elevation_scan )
	{
		if ( tdwr_phdr->end_of_volume_scan )
			tdwr_radial_status = 4;
		else
			tdwr_radial_status = 2;
	}
	else
		tdwr_radial_status = 1;

#ifdef	DEBUG
	printf("rec %d ",tdwr_rec);
	printf("%x ",tdwr_pmsg->message_id);
	printf("vsc %d ",tdwr_phdr->volume_scan_count);
	printf("azi %5.2f elev %5.2f ",tdwr_phdr->current_azimuth,
		tdwr_phdr->current_elevation);
	printf("sov %d eov %d ",tdwr_phdr->start_of_volume_scan,
		tdwr_phdr->end_of_volume_scan);
	printf("tilt # %d ",tdwr_phdr->tilt_number);
	printf("soe %d eoe %d ",tdwr_phdr->start_of_elevation_scan,
		tdwr_phdr->end_of_elevation_scan);
	printf("rs %d ",tdwr_radial_status);
	printf("\n");
#endif	/* DEBUG */


/*
 *	We don't play around with "nssl_scan_num", because we've already
 *	handled it.
 */

	if ( tdwr_old_vsn >= 0 )
	{
		if ( tdwr_phdr->volume_scan_count == tdwr_old_vsn )
			/* EMPTY */ ;
		else if ( tdwr_phdr->volume_scan_count == tdwr_old_vsn + 1 )
			/* EMPTY */ ;
/*
 *	This must have been a problem on an early TDWR data tape.
 */
		else
		printf("\nWARNING:  volume scan inconsistancy.  Previous %d.  Current %d\n",
			tdwr_old_vsn, tdwr_phdr->volume_scan_count );
	}
	tdwr_old_vsn = tdwr_phdr->volume_scan_count;

/*
 *	Check and see if we have a data restart.  There seem to be two
 *	different flags.  We'll check for both.
 *
 *	IMPORTANT:  I've not see either of these in any data.
 *
 */

	nssl_restart_flag = 0;

	if ( tdwr_phdr->volume_scan_restart )
	{
		printf("decode_data_tdwr:  volume scan restart detected\n");
		nssl_restart_flag++;
	}

	if ( tdwr_phdr->elevation_scan_restart )
	{
		printf("decode_data_tdwr:  elevation scan restart detected\n");
		nssl_restart_flag++;
	}

	return( 0 );
}

print_header()
{
/*
	printf("Message Id:  %x Message Length:  %d\n",
		tdwr_pmsg->message_id, tdwr_pmsg->message_length );
 */
	if (tdwr_pmsg->message_id == 0x2b02) {
		printf("Base_Data_Date_Stamp:  %-8.8s\n",
			tdwr_p2b02->base_data_date_stamp);
		printf("Base_Data_Time_Stamp:  %-8.8s\n",
			tdwr_p2b02->base_data_time_stamp);
		printf("Base_Data_Tape_Header_Text:  %-64.64s\n\n",
			tdwr_p2b02->base_data_header_text);
	}

	return( 0 );
}

/*
 *	Retrieve a decode data field.
 *
 *	Return values:
 *		0		requested field returned
 *		2, 4		requested field returns, value is "vel_res"
 *		MISSING		requested field not present this radial,
 *				or not present at all
 *
 *	For MISSING, all data values are set to MISSING_DATA.
 *
 *	Arguments are:
 *		num		data field number from a2defs.h
 *		data		real array for the data
 *		n		the number of data points
 *
 *	All "snr" hooks have been removed.
 *
 *	All of the above remarks apply to 88D data.  TDWR data *does* have
 *	SNR data.
 */

get_data_field_tdwr( num, data, n )
int num;
float data[];
int n;
{
	int i, rc;
	int sf, offset;
	float *fptr;
	float f1;
	unsigned short int *iptr, *idata;

/*
 *	Initialize data.
 */

	f1 = MISSING_DATA;

	fptr = data;
	for( i=0; i<n; i++)
		*fptr++ = f1;

	if ( get_status( num ) == 0 )
		return( MISSING );

	rc = get_scale_values( num, &sf, &offset );

	if ( rc < 0 )
		return( MISSING );

	idata = ( unsigned short int * ) malloc( n * 2 );

	if ( idata == NULL )
	{
		perror("get_data_field_tdwr:  malloc: " );
		return( MISSING );
	}

	rc = get_data_field_raw_tdwr( num, idata, n );

/*
printf("NUM=%d sf=%d offset=%d\n",num,sf,offset);
 */
	iptr = idata;
	fptr = data;

	if ( rc >= 0 )
	{
		for(i=0; i<n; i++)
		{
			if ( *iptr == 0 )
				*fptr = MISSING_DATA;
			else
				*fptr = (float) (*iptr) / (float) (sf) -
						( float ) offset;
/*
printf("i=%d iptr=%d fptr=%f\n",i,*iptr,*fptr);
 */
			iptr++;
			fptr++;
		}
	}

	free( idata );

	return( rc );

}

/*
 *	Retrieve a raw data field, and put it in short ints.
 *
 *	Return values:
 *		0		requested field returned
 *		2, 4		requested field returns, value is "vel_res"
 *		MISSING		requested field not present this radial,
 *				or not present at all
 *
 *	For MISSING, all data values are set to MISSING_DATA.
 *
 *	Arguments are:
 *		num		data field number from defs.h
 *		data		real array for the data
 *		n		the number of data points
 *
 *	Use get_scale_values() to return scaling.
 *
 */

get_data_field_raw_tdwr( num, data, n )
int num;
unsigned short int data[];
int n;
{
	int i, k;
	int num_gates;
	unsigned short int *iptr;

/*
 *	Initialize data to missing.
 */

	iptr = data;

	for(i=0; i<n; i++)
		*iptr++ = 0;

	if ( get_status( num ) == 0 )
		return( MISSING );

	num_gates = get_number_of_gates( num );

	k = n;
	if ( num_gates < n )
		k = num_gates;

	iptr = data;

#define	SNR_BAD		0

/*
 *	Reflectivity data cleanup are non-TDWR based.
 *
 *	For message id 0x2b00:
 *		Require SNR threshold.
 *
 *	For message id 0x2b01:
 *		For gates 1-600:
 *			Require "compressed_valid_flag" to be set *and* the
 *			SNR threshold.
 *		For gates 601-1984:
 *			Require SNR threshold.  Ignore "compressed_valid_flag"
 *			which is always zero.
 *
 *	The "SNR threshold" requires a non-zero value.  Zero values seem to
 *	be "ring" garbage, or causes large grey areas inside RADS.
 *
 *	Gate spacing...
 *
 *	For message id 0x2b00:
 *		Always 150 meters
 *
 *	For message id 0x2b01:
 *		For gates 1-900:
 *			150 meters
 *		For gates 901-1984:
 *			300 meters
 *
 *	Consequently, for message id 0x2b01, everything beyond gate 900
 *	will be marked as BAD!
 */

	if ( tdwr_pmsg->message_id == 0x2b00 )
	{
		for(i=0; i<k; i++)
		{
			if ( num == DBZ )
			{
				*iptr = tdwr_p2b00->tdwr_data[i].dbz;
				if ( tdwr_p2b00->tdwr_data[i].snr == SNR_BAD )
					*iptr = 0;
			}

			if ( num == RVEL )
				*iptr = tdwr_p2b00->tdwr_data[i].vel_raw;

			if ( num == VEL )
			{
				*iptr = tdwr_p2b00->tdwr_data[i].vel_dealias;

/*
 *	Delete the point if it has been declared to be bad.
 */

			if ( !tdwr_p2b00->tdwr_data[i].
				compressed_conditioned_valid_velocity_flag )
				*iptr = 0;

			}

			if ( num == SNR )
				*iptr = tdwr_p2b00->tdwr_data[i].snr;

			if ( num == SPW )
				*iptr = tdwr_p2b00->tdwr_data[i].spw;

			iptr++;
		}
	}
	else if ( tdwr_pmsg->message_id == 0x2b01 )
	{
		for(i=0; i<k; i++)
		{
			if ( num == DBZ )
			{
				*iptr = tdwr_p2b01->dbz[i];

/*
 *	Delete the point if it has been declared to be bad, but only check
 *	the first 600 gates, as it is always zero beyond!
 */

				if ( i < TDWR_MAX_2B00 && 
					!tdwr_p2b01->
					tdwr_flags[i].compressed_valid_flag )
					*iptr = 0;

				if ( tdwr_p2b01->snr[i] == SNR_BAD )
					*iptr = 0;

#define	GATE_LIMIT_150M_GATE_SPACING	900
				if ( i >= GATE_LIMIT_150M_GATE_SPACING )
					*iptr = 0;
			}

			if ( num == SNR )
				*iptr = tdwr_p2b01->snr[i];
			iptr++;
		}
	}
	else 
	{
	printf("get_data_field_raw_tdwr:  logic flaw...message_id is %x\n",
		tdwr_pmsg->message_id);
		exit( 1 );
	}

	return( 0 );
}

/*
 *	Retrieve a raw data field, and put it in unsigned characters.  This
 *	is the way the data arrives.
 *
 *	Return values:
 *		0		requested field returned
 *		2, 4		requested field returns, value is "vel_res"
 *		MISSING		requested field not present this radial,
 *				or not present at all
 *
 *	For MISSING, all data values are set to MISSING_DATA.
 *
 *	Arguments are:
 *		num		data field number from defs.h
 *		data		real array for the data
 *		n		the number of data points
 *
 *	Use get_scale_values() to return scaling.
 *
 */

get_data_field_uchar_tdwr( num, data, n )
int num;
unsigned char data[];
int n;
{
	int i;
	int rc;
	unsigned char *uptr;
	unsigned short int *iptr, *idata;

/*
 *	Initialize data to missing.
 */

	bzero( data, n );

	if ( get_status( num ) == 0 )
		return( MISSING );

	if ( num == VEL )
	{
		printf("get_data_field_uchar_tdwr:  TDWR velocity data uses 16 bits.\n");
		printf("get_data_field_uchar_tdwr:  It can't be copied to an 8 bit variable.\n");
		return( MISSING );
	}

	idata = ( unsigned short int * ) malloc( n * 2 );

	if ( idata == NULL )
	{
		perror("get_data_field_uchar_tdwr:  malloc: " );
		return( MISSING );
	}

	rc = get_data_field_raw_tdwr( num, idata, n );

	iptr = idata;
	uptr = data;

	if ( rc >= 0 )
	{
		for(i=0; i<n; i++)
		{
			*uptr = ( unsigned char ) *iptr;
			iptr++;
			uptr++;
		}
	}

	free( idata );

	return( rc );
}

/*
 *	88D data comes in the range 0-255 (unsigned) where 0 is missing and
 *	1 is range folded.  The data can be converted to real (float) values
 *	by an equation that looks like:
 *
 *		(real) real_val = (real) f1 * (int) data_val - (real) f2
 *
 *	"data_val" comes from the A2 data.  "f1" and "f2" from the 88D
 *	documentation are are in function get_data_field().  Note that there
 *	are two possible sets of numbers for velocity.
 *
 *	Computations in floating point isn't practical in realtime, so we
 *	need to do everything in integer arithmetic.  Fortunately, 1/f1 is
 *	always an integer.
 *
 *	When transformed, the equation looks like:
 *
 *		(int) scaled_real_val = (int) data_val - offset
 *
 *	The scale factor is "sf", which is always 1/f1.  Fortunately, that
 *	is always an integer. (scale_real_val = real_val * sf )
 *
 *	To go from real value to 88D values:
 *
 *		(int) data_val = (int) real_val * sf + offset
 *
 *	This routine returns "sf" (scale factor) and "offset".
 *
 *	The user needs to handle the missing and range folded values before
 *	applying the above equation.
 *
 *	Note:  the user should check to make sure the function does not
 *	       return MISSING before using the scaling information.
 *
 *	The purpose of the recomputations is to remove floating point
 *	computations from this code.  Any software that needs to retrieve
 *	floating values (such as NSSL RADS) has all the information necessary.
 *
 *	TDWR follow similar rules, though data sizes are different.
 *	DBZ, SPW, and SNR are 8 bits.  VEL and RVEL are 16 bits.
 */

get_scale_values_tdwr( num, sf, offset )
int num, *sf, *offset;
{

/*
 *	Initialize data to missing.
 */

	if ( get_status( num ) == 0 )
		return( MISSING );

	if ( num == DBZ )
	{
		*sf = 2;
		*offset = 30;
	}
	else if ( num == VEL || num == RVEL )
	{
		*sf = 4;
		*offset = 80;
	}
	else if ( num == SPW )
	{
		*sf = 4;
		*offset = 0;
	}
	else if ( num == SNR )
	{
		*sf = 2;
		*offset = 0;
	}
	else
	{
	printf("get_scale_values_tdwr:  %d:  'num' is corrupt.  help!\n",num);
		return( MISSING );
	}

	return( 0 );
}

get_tdwr_radial_status()
{
	return( tdwr_radial_status );
}

#endif

/*
 *	No one is using this.
 */

#ifdef	ALLOW_TDWR
#ifdef	DEBUG
dumps()
{
	printf("vsc %d\n",tdwr_phdr->volume_scan_count);
	printf("eovs %d sovs %d dummy %d ss %d\n",
		tdwr_phdr->end_of_volume_scan, tdwr_phdr->start_of_volume_scan, tdwr_phdr->dummy1, tdwr_phdr->scan_strategy);
	printf("ptp %d d2 %d dri %d sold %d sopd %d\n",
		tdwr_phdr->peak_transmitter_power,tdwr_phdr->dummy2,
		tdwr_phdr->dummy_record_indicator,
		tdwr_phdr->start_of_live_data, tdwr_phdr->start_of_playback_data);
	printf("tn %d eoes %d soes %d crmn %d res %d\n",
		tdwr_phdr->tilt_number, tdwr_phdr->end_of_elevation_scan,
		tdwr_phdr->start_of_elevation_scan, tdwr_phdr->clutter_residue_map_number,
		tdwr_phdr->reserved1);
/*
 *	skip
 */

	printf("current elev %f\n",tdwr_phdr->current_elevation);
	printf("asr %f\n",tdwr_phdr->angular_scan_rate);

/*
 *	skip
 */

	printf("current_azimuth %f\n",tdwr_phdr->current_azimuth );
	printf("total_noise_power %f\n",tdwr_phdr->total_noise_power );
	printf("timestamp %d\n",tdwr_phdr->timestamp);

	return;
}
#endif	/* DEBUG */
#endif	/* ALLOW_TDWR */
