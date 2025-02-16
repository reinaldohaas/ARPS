/*
 *	Radardefs.h...
 *
 *	Misc parameters...
 */

#ifndef	_RADARDEFS_H
#define	_RADARDEFS_H

#pragma ident "@(#)radardefs.h	4.1	01/30/01 NSSL"

/*
 *
 *	Data types.
 *
 *	        -2	Data format can't be converted (such as converting
 *				a tape to realtime record)
 *	        -1	End of data
 *		 0	Unknown data type (data may be corrupt)
 */

#define	RADAR_CANT_CONVERT	-2
#define	RADAR_END_OF_DATA	-1
#define	RADAR_UNKNOWN_TYPE	 0

#endif	/* _RADARDEFS_H */
