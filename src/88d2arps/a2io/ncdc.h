/*
 *	ncdc.h...
 *
 *	NCDC format.
 */

#pragma ident "@(#)ncdc.h	4.1	01/30/01 NSSL"

#define	NCDC_SIZE	31616

/*
 *	Note that "ncdc_valid" is an unused slot in the NCDC format, so we
 *	use so that our software can determine if the data is valid.
 */

struct ncdc_header {
	char archive2[8];			/* always ARCHIVE2 */
	char id[4];				/* 4-letter site id */
	char tape_num[6];			/* NCDC tape number */
	char blank19;				/* blank */
	char tape_date[9];			/* date tape was written */
	char blank29;				/* blank */
	char tape_time[8];			/* time tape was written */
	char blank38;				/* blank */
	char tape_source[5];			/* either RDASC or NCDC */
	char wban[5];				/* WBAN number of NEXRAD site */
	char tape_mode[5];			/* tape output mode */
	char vol_num[5];			/* volume number */
	char blank59[6];			/* blank */
	char internal[31552];			/* who knows what is here */
} ;
