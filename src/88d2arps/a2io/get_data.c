/*
 *	get_data.c...		Data retrieval routines.
 *
 *	get_data_field()	Returns floating point values.
 *	get_data_field_raw()	Returns raw data values (short ints).
 *	get_data_field_uchar()	Return raw data values (unsigned chars).
 *	get_scale_values()	Recoding scale info for raw data values.
 */
 
#pragma ident "@(#)get_data.c	3.12	04/09/01 NSSL"

#include <stdio.h>
#include <time.h>

#include "a2io.h"
#include "extern.h"
#include "extern_funct.h"
#include "proto.h"

/*
 *	Retrieve a decode data field.
 *
 *	Return values:
 *		0		requested field returned
 *		2, 4		requested field returns, value is "vel_res"
 *		MISSING		requested field not present this radial,
 *				or not present at all
 *
 *	NOTE:  TDWR always returns zero or MISSING.
 *
 *	For MISSING, all data values are set to MISSING_DATA.
 *
 *	Arguments are:
 *		num		data field number from defs.h
 *		data		real array for the data
 *		n		the number of data points
 */

get_data_field( num, data, n )
int num;
float data[];
int n;
{
	if ( funct_get_data_field == NULL )
		radar_oops("get_data_field");

	return( funct_get_data_field( num, data, n ));
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
 *	NOTE:  TDWR always returns zero or MISSING.
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

get_data_field_raw( num, data, n )
int num;
unsigned short int data[];
int n;
{
	if ( funct_get_data_field_raw == NULL )
		radar_oops("get_data_field_raw");

	return( funct_get_data_field_raw( num, data, n )) ;
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
 *	NOTE:  TDWR always returns zero or MISSING.
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

get_data_field_uchar( num, data, n )
int num;
unsigned char data[];
int n;
{
	if ( funct_get_data_field_uchar == NULL )
		radar_oops("get_data_field_uchar");

	return( funct_get_data_field_uchar( num, data, n ));
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

get_scale_values( num, sf, offset )
int num, *sf, *offset;
{
	if ( funct_get_scale_values == NULL )
		radar_oops("get_scale_values");

	return( funct_get_scale_values( num, sf, offset ));
}
