/*
 *	Tdwrdefs.h...
 *
 *	TDWR parameters...
 */

#pragma ident "@(#)tdwrdefs.h	4.1	01/30/01 NSSL"

#include "radardefs.h"

/*
 *	IT IS EXTREMELY IMPORTANT THAT NONE OF THESE ASSIGNMENTS CONFLICT
 *	WITH ANYTHING IN "a2defs.h"!!!!!!
 *
 *	Data types.
 *
 *	       210	0x2b00 record, 600 gates dbz, vel, rvel, spw
 *	       211	0x2b01 record, 1984 gates dbz, spw
 *	       212	0x2b02 record, tape header record
 *	       213	0x2c00 record, LLWAS2 sensor record
 *	       214	0x2c01 record, LLWAS3 data record
 *	       215	0x2c02 record, LLWAS2 data record
 */

#define	TDWR_CANT_CONVERT	RADAR_CANT_CONVERT
#define	TDWR_END_OF_DATA	RADAR_END_OF_DATA
#define	TDWR_UNKNOWN_TYPE	RADAR_UNKNOWN_TYPE
#define	TDWR_2B00		210
#define	TDWR_2B01		211
#define	TDWR_2B02		212
#define	TDWR_2C00		213
#define	TDWR_2C01		214
#define	TDWR_2C02		215
