/*
 *	get_data.c...		Data retrieval routines.
 *
 *	get_data_field()	Returns floating point values.
 *	get_data_field_raw()	Returns raw data values (short ints).
 *	get_data_field_uchar()	Return raw data values (unsigned chars).
 *	get_scale_values()	Recoding scale info for raw data values.
 */

#pragma ident "@(#)get_data_a2io.c	5.2	10/17/05 CAPS"

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <time.h>
#include "a2io.h"
#include "a2io_message_31.h"
#include "extern.h"
#include "extern_funct.h"
#include "proto.h"

#define	MISSING_88D		0
#define	RANGE_FOLDED_88D	1

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
 *		num		data field number from defs.h
 *		data		real array for the data
 *		n		the number of data points
 *
 *	All "snr" hooks have been removed.
 *
 *	All of the above remarks apply to 88D data.  TDWR data *does* have
 *	SNR data.
 */

get_data_field_a2io( num, data, n )
int num;
float data[];
int n;
{
	int i, k;
	int num_gates;
	int start_pointer;
	int flag;
	float *fptr;
	float f1, f2;
	unsigned char *vptr;

	if ( message_type == A2_DIGITAL_TYPE )
	{
		get_data_field_a2io_31( num, data, n );
		return(0);
	}

	if ( message_type != A2_DATA_TYPE )
		return(0);

/*
 *	Initialize data.
 */

	f1 = MISSING_DATA;

	fptr = data;
	for( i=0; i<n; i++)
		*fptr++ = f1;

	if ( funct_get_status == NULL )
		radar_oops( "get_status" );

	if ( funct_get_status( num ) == 0 )
		return( MISSING );

	flag = nssl_a2data->vel_res;

	if ( num == DBZ )
	{
		flag = 0;
		num_gates = nssl_a2data->num_gates_dbz;
#ifdef	VERY_OLD_FLAWED_DATA
		start_pointer = nssl_a2data->arc_dbz_ptr;
#else
		start_pointer = nssl_a2data->dbz_ptr;
#endif
		f1 = 0.5;
		f2 = 33.0;
	}
	else if ( num == VEL || num == RVEL )
	{
		num_gates = nssl_a2data->num_gates_vel;
#ifdef	VERY_OLD_FLAWED_DATA
		start_pointer = nssl_a2data->arc_vel_ptr;
#else
		start_pointer = nssl_a2data->vel_ptr;
#endif
		if ( nssl_a2data->vel_res == 2 )
		{
			f1 = 0.5;
			f2 = 64.5;
		}
		else if ( nssl_a2data->vel_res == 4 )
		{
			f1 = 1.0;
			f2 = 129.0;
		}
		else
		{
printf("get_data_field:  %d:  can't decode velocity this radial, bad vel_res.\n",
				nssl_a2data->vel_res);
			return( MISSING );
		}
#ifdef	VERY_OLD_FLAWED_DATA
		threshold( nssl_a2data->num_gates_vel, nssl_a2data->arc_vel_ptr,
			nssl_a2data->arc_spw_ptr );
#else
		threshold( nssl_a2data->num_gates_vel, nssl_a2data->vel_ptr,
			nssl_a2data->spw_ptr );
#endif
	}
	else if ( num == SPW )
	{
		num_gates = nssl_a2data->num_gates_vel;
#ifdef	VERY_OLD_FLAWED_DATA
		start_pointer = nssl_a2data->arc_spw_ptr;
#else
		start_pointer = nssl_a2data->spw_ptr;
#endif
		f1 = 0.5;
		f2 = 64.5;
	}
	else
	{
		printf("get_data_field:  %d:  'num' is corrupt.  help!\n",num);
		return( MISSING );
	}

	k = n;
	if ( num_gates < n )
		k = num_gates;
	fptr = data;
	vptr = ( unsigned char * ) &nssl_dptr[start_pointer];

	for(i=0; i<k; i++)
	{
		if ( *vptr == MISSING_88D )
			*fptr = MISSING_DATA;
		else if ( *vptr == RANGE_FOLDED_88D )
			*fptr = RANGE_FOLDED_DATA;
		else
			*fptr = f1 * (float) *vptr - f2;
		fptr++;
		vptr++;
	}
	return( flag );
}


get_data_field_a2io_31( num, data, n )
int num;
float data[];
int n;
{
	int i, k;
	int num_gates;
	float *fptr;
	float f1, f2;
	unsigned char *vptr;
	unsigned short *sptr;
	int size;

/*
 *	Initialize data.
 */

	f1 = MISSING_DATA;

	fptr = data;
	for( i=0; i<n; i++)
		*fptr++ = f1;

	if ( funct_get_status == NULL )
		radar_oops( "get_status" );

	if ( funct_get_status( num ) == 0 )
		return( MISSING );

	if ( num == DBZ )
	{
		f1 = 1.0 / nssl_a2data31_dbz->scale;
		f2 = nssl_a2data31_dbz->offset / nssl_a2data31_dbz->scale;
/*
printf("DBZ:  f1=%f f2=%f (should be 0.5 and 33.0)\n",f1,f2);
 */
		num_gates = nssl_a2data31_dbz->number_of_gates;
		vptr = (unsigned char *) &nssl_dptr[map31[REF_31] + 28];
		size = 1;
	}
	else if ( num == VEL )
	{
		f1 = 1.0 / nssl_a2data31_vel->scale;
		f2 = nssl_a2data31_vel->offset / nssl_a2data31_vel->scale;
/*
printf("VEL:  f1=%f f2=%f (should be 0.5 and 64.5)\n",f1,f2);
 */

		num_gates = nssl_a2data31_vel->number_of_gates;
		vptr = (unsigned char *) &nssl_dptr[map31[VEL_31] + 28];
		size = 1;
	}
	else if ( num == SPW )
	{
		f1 = 1.0 / nssl_a2data31_spw->scale;
		f2 = nssl_a2data31_spw->offset / nssl_a2data31_spw->scale;
/*
printf("SPW:  f1=%f f2=%f (should be 0.5 and 64.5)\n",f1,f2);
 */

		num_gates = nssl_a2data31_spw->number_of_gates;
		vptr = (unsigned char *) &nssl_dptr[map31[SPW_31] + 28];
		size = 1;
	}
/*
 * "should be" values for the rest are wrong.
 */
	else if ( num == ZDR )
	{
		f1 = 1.0 / nssl_a2data31_zdr->scale;
		f2 = nssl_a2data31_zdr->offset / nssl_a2data31_zdr->scale;
/*
printf("ZDR:  f1=%f f2=%f (should be 0.5 and 64.5)\n",f1,f2);
 */

		num_gates = nssl_a2data31_zdr->number_of_gates;
		vptr = (unsigned char *) &nssl_dptr[map31[ZDR_31] + 28];
		size = 1;
	}
	else if ( num == PHI )
	{
		f1 = 1.0 / nssl_a2data31_phi->scale;
		f2 = nssl_a2data31_phi->offset / nssl_a2data31_phi->scale;
/*
printf("PHI:  f1=%f f2=%f (should be 0.5 and 64.5)\n",f1,f2);
 */

		num_gates = nssl_a2data31_phi->number_of_gates;
		vptr = (unsigned char *) &nssl_dptr[map31[PHI_31] + 28];
		sptr = (unsigned short *) &nssl_dptr[map31[PHI_31] + 28];
		size = 2;
	}
	else if ( num == RHO )
	{
		f1 = 1.0 / nssl_a2data31_rho->scale;
		f2 = nssl_a2data31_rho->offset / nssl_a2data31_rho->scale;
/*
printf("RHO:  f1=%f f2=%f (should be 0.5 and 64.5)\n",f1,f2);
 */

		num_gates = nssl_a2data31_rho->number_of_gates;
		vptr = (unsigned char *) &nssl_dptr[map31[RHO_31] + 28];
		size = 1;
	}
	else
	{
		printf("get_data_field:  %d:  'num' is corrupt.  help!\n",num);
		return( MISSING );
	}

	k = n;
	if ( num_gates < n )
		k = num_gates;
	fptr = data;

	for(i=0; i<k; i++)
	{
		if (size == 1) {
			if ( *vptr == MISSING_88D )
				*fptr = MISSING_DATA;
			else if ( *vptr == RANGE_FOLDED_88D )
				*fptr = RANGE_FOLDED_DATA;
			else
				*fptr = f1 * (float) *vptr - f2;
			vptr++;
		}
		if (size == 2) {
			*sptr = htons(*sptr);
/*
printf("*vptr=%d and %o...",*vptr,*vptr);
vptr++;
printf("*vptr=%d and %o...",*vptr,*vptr);
vptr++;
printf("*sptr=%d and %o\n",*sptr,*sptr);
printf("*sptr=%d and %o\n",*sptr,*sptr);
 */
			if ( *sptr == MISSING_88D )
				*fptr = MISSING_DATA;
			else if ( *sptr == RANGE_FOLDED_88D )
				*fptr = RANGE_FOLDED_DATA;
			else
				*fptr = f1 * (float) *sptr - f2;
			sptr++;
		}
/*
printf("%2d = %f\n",i,*fptr);
 */
		fptr++;
	}
	return( 0 );
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

get_data_field_raw_a2io( num, data, n )
int num;
unsigned short int data[];
int n;
{
	int i, k;
	int num_gates;
	int start_pointer;
	int flag;
	unsigned char *vptr;
	unsigned short int *iptr;

	if ( message_type == 31 )
	{
		printf("get_data_field_raw_a2io:  not message31 ready\n");
		exit(1);
	}

/*
 *	Initialize data to missing.
 */

	iptr = data;
	for( i=0; i<n; i++)
		*iptr++ = 0;

	if ( funct_get_status == NULL )
		radar_oops( "get_status" );

	if ( funct_get_status( num ) == 0 )
		return( MISSING );

	flag = nssl_a2data->vel_res;

	if ( num == DBZ )
	{
		flag = 0;				/* just to be safe */
		num_gates = nssl_a2data->num_gates_dbz;
#ifdef	VERY_OLD_FLAWED_DATA
		start_pointer = nssl_a2data->arc_dbz_ptr;
#else
		start_pointer = nssl_a2data->dbz_ptr;
#endif
	}
	else if ( num == VEL || num == RVEL )
	{
		num_gates = nssl_a2data->num_gates_vel;
#ifdef	VERY_OLD_FLAWED_DATA
		start_pointer = nssl_a2data->arc_vel_ptr;
		threshold( nssl_a2data->num_gates_vel, nssl_a2data->arc_vel_ptr,
			nssl_a2data->arc_spw_ptr );
#else
		start_pointer = nssl_a2data->vel_ptr;
		threshold( nssl_a2data->num_gates_vel, nssl_a2data->vel_ptr,
			nssl_a2data->spw_ptr );
#endif
	}
	else if ( num == SPW )
	{
		num_gates = nssl_a2data->num_gates_vel;
#ifdef	VERY_OLD_FLAWED_DATA
		start_pointer = nssl_a2data->arc_spw_ptr;
#else
		start_pointer = nssl_a2data->spw_ptr;
#endif
	}
	else
	{
		printf("get_data_field_raw:  %d:  'num' is corrupt.  help!\n",num);
		return( MISSING );
	}

	k = n;
	if ( num_gates < n )
		k = num_gates;
	iptr = data;
	vptr = ( unsigned char * ) &nssl_dptr[start_pointer];

	for(i=0; i<k; i++)
	{
		*iptr = (unsigned short) *vptr;
		iptr++;
		vptr++;
	}
	return( flag );
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

get_data_field_uchar_a2io( num, data, n )
int num;
unsigned char data[];
int n;
{
	int i, k;
	int num_gates;
	int start_pointer;
	int flag;
	unsigned char *vptr;

	if ( message_type == 31 )
	{
		printf("get_data_field_raw_uchar_a2io:  not message31 ready\n");
		exit(1);
	}

/*
 *	Initialize data to missing.
 */

	bzero( data, n );

	if ( funct_get_status == NULL )
		radar_oops( "get_status" );

	if ( funct_get_status( num ) == 0 )
		return( MISSING );

	flag = nssl_a2data->vel_res;

	if ( num == DBZ )
	{
		flag = 0;				/* just to be safe */
		num_gates = nssl_a2data->num_gates_dbz;
#ifdef	VERY_OLD_FLAWED_DATA
		start_pointer = nssl_a2data->arc_dbz_ptr;
#else
		start_pointer = nssl_a2data->dbz_ptr;
#endif
	}
	else if ( num == VEL || num == RVEL )
	{
		num_gates = nssl_a2data->num_gates_vel;
#ifdef	VERY_OLD_FLAWED_DATA
		start_pointer = nssl_a2data->arc_vel_ptr;
		threshold( nssl_a2data->num_gates_vel, nssl_a2data->arc_vel_ptr,
			nssl_a2data->arc_spw_ptr );
#else
		start_pointer = nssl_a2data->vel_ptr;
		threshold( nssl_a2data->num_gates_vel, nssl_a2data->vel_ptr,
			nssl_a2data->spw_ptr );
#endif
	}
	else if ( num == SPW )
	{
		num_gates = nssl_a2data->num_gates_vel;
#ifdef	VERY_OLD_FLAWED_DATA
		start_pointer = nssl_a2data->arc_spw_ptr;
#else
		start_pointer = nssl_a2data->spw_ptr;
#endif
	}
	else
	{
		printf("get_data_field_uchar:  %d:  'num' is corrupt.  help!\n",num);
		return( MISSING );
	}

	k = n;
	if ( num_gates < n )
		k = num_gates;

	vptr = ( unsigned char * ) &nssl_dptr[start_pointer];

	bcopy( vptr, data, k );

	return( flag );
}

/*
 *	Threshold velocity data.
 *
 *	For now, this function is in the A2IO driver, as we TDWR processing
 *	won't use this, at least for now.
 */

threshold( vel_ng, vel_ptr, spw_ptr )
unsigned short int vel_ng, vel_ptr, spw_ptr;
{
	int sf, offset;
	int rc;
	int i;
	int points_saved, points_deleted;
	unsigned char value;
	unsigned char *vptr_vel, *vptr_spw;

/*
 *	Make sure we are supposed to do this.
 */

	if ( allow_thresholds == 0 )
		return;

/*
 *	This shouldn't happen, as we've already check for the presense of VEL.
 */

	if ( vel_ng == 0 )
	{
		printf("threshold:  vel_ng is zero ???\n");
		return;
	}

/*
 *	Find our conversion factors.  We'll do this once per radial, though
 *	in practice, it probably isn't necessary to do it so much.
 */

	if ( funct_get_scale_values == NULL )
		radar_oops( "get_scale_values" );

	rc = funct_get_scale_values( SPW, &sf, &offset );

/*
 *	More sanity checks.
 */

	if ( rc == MISSING || sf == 0 || offset == 0 )
	{
		printf("threshold:  return from get_scale_values %d %d %d\n",
			rc, sf, offset );
		return;
	}

/*
 *	Convert the threshold value to an 88D value.
 */

	value = (unsigned char) (spw_threshold * sf + offset);

/*
 *	Sanity.  Make sure we don't destroy "RANGE_FOLDED_88D".
 */

	if ( value <= 1 )
		value = 2;

	vptr_vel = ( unsigned char *) &nssl_dptr[ vel_ptr ];
	vptr_spw = ( unsigned char *) &nssl_dptr[ spw_ptr ];

	points_saved = 0;
	points_deleted = 0;

/*
 *	WARNING...basic level 2 data is modified here!
 */

	for(i=0; i<(int)vel_ng; i++)
	{
		if ( *vptr_vel != MISSING_88D && *vptr_vel != RANGE_FOLDED_88D )
		{
			if ( *vptr_spw != MISSING_88D && *vptr_spw != RANGE_FOLDED_88D )
			{
				if ( *vptr_spw >= value )
				{
					points_deleted++;
					*vptr_vel = MISSING_88D;
				}
				else
					points_saved++;
			}
		}
		vptr_vel++;
		vptr_spw++;
	}

/*
#define	THRESHOLD_STATS
 */
#ifdef	THRESHOLD_STATS
	printf("%d/%d points saved/deleted\n",points_saved,points_deleted);
#endif

	return;
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

get_scale_values_a2io( num, sf, offset )
int num, *sf, *offset;
{

	if ( message_type == 31 )
	{
		printf("get_scale_values_a2io not message 31 ready\n");
		exit(1);
	}

/*
 *	Initialize data to missing.
 */

	if ( funct_get_status == NULL )
		radar_oops( "get_status" );

	if ( get_status( num ) == 0 )
		return( MISSING );

	if ( num == DBZ )
	{
		*sf = 2;				/* 1/0.5 */
		*offset = 66;				/* 33 * 2 */
	}
	else if ( num == VEL || num == RVEL )
	{
		if ( nssl_a2data->vel_res == 2 )
		{
			*sf = 2;			/* 1/0.5 */
			*offset = 129;			/* 64.5 * 2 */
		}
		else if ( nssl_a2data->vel_res == 4 )
		{
			*sf = 1;			/* 1/1 */
			*offset = 129;			/* 129 * 1 */
		}
		else
		{
printf("get_scale_values_a2io:  %d:  can't decode velocity this radial, bad vel_res.\n",
				nssl_a2data->vel_res);
			return( MISSING );
		}
	}
	else if ( num == SPW )
	{
		*sf = 2;				/* 1/0.5 */
		*offset = 129;				/* 64.5 * 2 */
	}

	else
	{
	printf("get_scale_values_a2io:  %d:  'num' is corrupt.  help!\n",num);
		return( MISSING );
	}

	return( 0 );
}
