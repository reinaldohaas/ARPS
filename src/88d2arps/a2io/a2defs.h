/*
 *	A2defs.h...
 *
 *	Various A2 data types, some are local additions.
 */

#pragma ident "@(#)a2defs.h	5.1	01/20/04 CAPS"

#include "radardefs.h"

/*
 *	Data types.
 *
 *		 1-200	Valid 88D messages.
 *		 1	Data radial.
 *		20	User message entered on system console.
 *	       201	NSSL Volume Scan Message.
 *	       202	NSSL Site/Timestamp Message.
 *	       210-215	See "tdwrdefs.h".
 *	       241	"ARCHIVE2" tape header record, 8-byte version.
 *	       242	"ARCHIVE2" tape header record, NCDC version.
 *	       243	Volume Scan Message (24-byte tape header record).
 *
 *	Starting in 2004, "ARCHIVE2" is replaced with "AR2V####", where "####"
 *	is a version number.  This software is backwards compatible.
 */

#define	A2_CANT_CONVERT		RADAR_CANT_CONVERT
#define	A2_END_OF_DATA		RADAR_END_OF_DATA
#define	A2_UNKNOWN_TYPE		RADAR_UNKNOWN_TYPE
#define	A2_DATA_TYPE		1
#define	A2_CONSOLE_MSG		20
#define	A2_DIGITAL_TYPE		31
#define	A2_VOLSCAN_TYPE		201
#define	A2_SITE_TYPE		202
#define	A2_ARCHIVE2		241
#define	A2_ARCHIVE2_NCDC	242
#define	A2_VS_24_BYTE		243

/*
 *	Product types.  Some are obsolete, though still defined.  Skipped
 *	numbers are some of the deleted (very old WDSS) products.
 */

#define	DBZ			0
#define	VEL			1
#define	SPW			2
#define	SW			2	/* this is what the NEXRAD docs use */
#define	SNR			3
#define	RVEL			9	/* raw velocity (TDWR only) */
#define	ZDR			10	/* differential reflectivity */
#define	PHI			11	/* differential phase */
#define	RHO			12	/* differential correlation coef */

#define	MISSING			-999

#define	MISSING_DATA		-999.0
#define	RANGE_FOLDED_DATA	-888.0
