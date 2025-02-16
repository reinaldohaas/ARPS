/*
 *	A2io.h...
 *
 *	A2 parameters...
 */

#pragma ident "@(#)a2io.h	5.5	10/17/05 CAPS"

#include <inttypes.h>
#include "a2defs.h"

/*
 *	Data structures.
 *
 *	Note that "ctm" is only in tape format.  Live routines handle this
 *	properly.
 */

struct a2header {
	char ctm[12];					/* not used */
	unsigned short message_size;			/* not used */
   #ifdef ALLOW_PKU
        unsigned short int  message_type;                       /* was short */
   #else
        unsigned char rda_channel;                      /* not used */
        unsigned char message_type;                     /* was short */
   #endif

	short seq;
	unsigned short julian_date;			/* not used */
	uint32_t milsec;				/* not used */
	unsigned short num_mes_seg;			/* not used */
	unsigned short cur_mes_seg;			/* not used */
} ;

/*
 *	Type "1" record.
 */

/*
 *	A little history...
 *
 *	In the early level 2 archive tapes, variables "dbz_ptr", "vel_ptr", and
 *	"spw_ptr" were wrong.  Variables "arc_dbz_ptr", "arc_vel_ptr" and
 *	"arc_spw_ptr" were always right.  Later, both sets of points were
 *	correct, and identical.  This code chose to use the "arc_*" names.
 *
 *	Fast forward to ORDA.  The "arc_*" variables are no longer used.
 *	Except when built for old data set names, these variable names are
 *	changed to guarantee that nothing is usin them.
 *
 *	FYI, the "might be wrong" below was another early bug.
 */

struct a2data {
	uint32_t zulu_time;
	unsigned short mod_julian_date;
	unsigned short unamb_range;
	unsigned short azimuth;
	unsigned short azimuth_number;
	unsigned short radial_status;			/* might be wrong! */
	unsigned short elevation;
	unsigned short elevation_number;
	short range_first_gate_dbz;
	short range_first_gate_vel;
	unsigned short gate_size_dbz;
	unsigned short gate_size_vel;
	unsigned short num_gates_dbz;
	unsigned short num_gates_vel;
	unsigned short cut_sector_number;
	int32_t calibration_constant;
	unsigned short dbz_ptr;
	unsigned short vel_ptr;
	unsigned short spw_ptr;
	unsigned short vel_res;
	unsigned short vcp;
	unsigned short dummy[4];
#ifdef	VERY_OLD_FLAWED_DATA
	unsigned short arc_dbz_ptr;
	unsigned short arc_vel_ptr;
	unsigned short arc_spw_ptr;
#else
	unsigned short arc_dbz_ptr_obsolete;
	unsigned short arc_vel_ptr_obsolete;
	unsigned short arc_spw_ptr_obsolete;
#endif
	unsigned short nyquist;
	short atmos;
	unsigned short tover;
	unsigned short radial_spot_blanking_status;
} ;

/*
 *	Type "20" record.
 */

#define	A2_MAX_CONSOLE_LEN	410

struct a2console {
	short len;
	char msg[A2_MAX_CONSOLE_LEN];
} ;

/*
 *	Type "201" record.
 */

struct a2scan {
	short vsn;
	short pad;
};

/*
 *	Type "202" record.
 */

struct a2site {
	unsigned char st_name[64];
	unsigned char st_site[32];
	uint32_t reset_tm;
};

/*
 *	Type "1003" (our number scheme) record.
 */

struct a2tape24 {
	char archive2[9];
	char vs[3];			/* convert with atoi() */
	short int not_used;
	unsigned short int julian_date;
	uint32_t milsec;
	char filler[4];
};

/*
 *	Note that the work "TAPE" in many variable names is of historical
 *	nature.  Data used to be processed from tape.  Now, everything is
 *	on disk, however, I'm not going to change these names.
 */

#define	MAX_RECORD_SIZE_TAPE	145920

/*
 *	Data input/output format types.  Type "DISK" isn't applicable to
 *	output format.
 *
 *	Note that TAPE below refers to the data format written to tapes, which
 *	for Level 2 data is the same, except of course, there are end of
 *	record marks on the tape.
 */

#define	UNKNOWN			-1		/* no call to radar_init() */
#define	TAPE			0		/* tape format */
#define	DISK			3		/* input is disk, 88D format */
#define	DISK_KWT		4		/* input is dir, KWT format */
#define	DISK_LDM		6		/* LDM (modified bzip2) */
#define	DISK_NCDC		7		/* input is dir, NCDC format */
#define	DISK_CAPS		8		/* LDM, CAPS naming scheme */
#define	DISK_TDWR		9		/* TDWR format */

/*
 *	Identify which reader to use.
 */

#define	FILE_MMAP		0		/* Plain disk file */
#define	FILE_GZIP		1		/* Gzip file */
#define	FILE_BZIP2		2		/* Bzip2 file */
#define	FILE_LDM		3		/* LDM (modified bzip2) file */
#define	FILE_COMPRESS		4		/* Compressed file */
 
/*
 *	Data format (actual data type).  Even though this files is "a2io.h",
 *	non-A2IO formats are or will be supported.
 */

#define	UNKNOWN_FORMAT		-1		/* To be determined format */
#define	A2_FORMAT		0		/* A2 format */
#define	TDWR_FORMAT		2		/* TDWR format */

/*
 *	Data record sizes.
 *
 *	REC_SIZE_GENERIC must be as big as the largest possible record
 *	size in any data set.  It is only valid for FILE_MMAP.
 *
 */

#define	REC_SIZE_TAPE		2432
#define	REC_SIZE_GENERIC	20000

#define	A2_ARCHIVE2_SIZE	8
#define	A2_VS_24_BYTE_SIZE	24

/*
 *	Flag that says we need to go into "search" mode for volume scan.
 *	Zero used to be the choice, however, LDM uses 0 at the volume scan
 *	reset time of 00Z.
 */

#define	VS_SEARCH		-99
