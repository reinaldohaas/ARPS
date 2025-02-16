/*
 *	extern.h...	External version of things in read_radial.h
 */

#pragma ident "@(#)extern.h	5.2	08/30/04 CAPS"

/*******************************************************************************
 *	Data and pointers to data.
 ******************************************************************************/

extern struct a2header	*nssl_a2header;		/* header part of record */
extern struct a2data	*nssl_a2data;		/* ptr for A2_DATA_TYPE */
extern struct a2data31	*nssl_a2data31;		/* ptr for A2_DIGITAL_TYPE */

extern struct a2data31_vol	*nssl_a2data31_vol;
extern struct a2data31_elv	*nssl_a2data31_elv;
extern struct a2data31_rad	*nssl_a2data31_rad;
extern struct a2data31_data	*nssl_a2data31_dbz;
extern struct a2data31_data	*nssl_a2data31_vel;
extern struct a2data31_data	*nssl_a2data31_spw;
extern struct a2data31_data	*nssl_a2data31_zdr;
extern struct a2data31_data	*nssl_a2data31_phi;
extern struct a2data31_data	*nssl_a2data31_rho;

extern struct a2data31_data	a2data31_dbz;
extern struct a2data31_data	a2data31_vel;
extern struct a2data31_data	a2data31_spw;
extern struct a2data31_data	a2data31_zdr;
extern struct a2data31_data	a2data31_phi;
extern struct a2data31_data	a2data31_rho;

extern struct a2scan	*nssl_a2scan;		/* ptr for A2_VOLSCAN_TYPE */
extern struct a2site	*nssl_a2site;		/* ptr for A2_SITE_TYPE */
extern struct a2tape24	*nssl_a2tape24;		/* ptr for A2_VS_24_BYTE */
extern struct ncdc_header *ncdc_header;		/* ptr for A2_ARCHIVE2_NCDC */

extern char nssl_buf[];				/* where the data goes */
extern char *nssl_ptr;				/* beginning of record */
extern char *nssl_dptr;				/* data part of record */
extern int map31[];				/* mapping of data */

extern unsigned short int *nssl_vcp;		/* ptr to vcp structure */

/******************************************************************************
 *	This data gets passed back thru user calls.
 *****************************************************************************/

extern int nssl_scan_status;			/* one when new volume scan */
extern int nssl_tilt_status;			/* one when new tilt */

extern int nssl_restart_flag;			/* tilt restart flag */

extern int nssl_scan_num;			/* current volume scan num */

extern int nssl_vs_timestamp;			/* from A2_SITE_TYPE record */

/******************************************************************************
 *	Input/output/datarequest format.
 *****************************************************************************/

extern int input_format;			/* input format */
extern int output_format;			/* output format */
extern int data_format;				/* actual data format type */
extern disk_format;				/* DISK_* in a2io.h */
extern int file_format;				/* FILE_* in a2io.h */

/******************************************************************************
 *	Bookkeeping information.
 *****************************************************************************/

extern int left;				/* bytes left to process */
extern int new;					/* byte loc in "nssl_buf" */

extern int message_type;			/* message type */
extern int message_length;			/* length of the message */

/******************************************************************************
 *	Time handling.
 *****************************************************************************/

extern struct tm *nssl_gmt;			/* time structure */
extern struct tm nssl_gmt_data;			/* dedicated space time */
extern long zulu_time;				/* time conversion */

/******************************************************************************
 *	Data thresholding.
 *****************************************************************************/

extern int allow_thresholds;			/* allow user thresholds */
extern int spw_threshold;			/* spw threshold */

/******************************************************************************
 *	Elevation error handling.
 *****************************************************************************/

extern unsigned short elev_save;	/* same type as a2data->elevation */
extern unsigned short elev_thresh;	/* ditto */

/******************************************************************************
 *	Tape volume number handling.
 *****************************************************************************/

extern int a2_start_scan;			/* starting volume scan */
extern int a2_end_scan;				/* ending volume scan */
extern int a2_num_scan;				/* number of vs's to process */

/******************************************************************************
 *	Misc radar.
 *****************************************************************************/

extern int a2_no_rads;			/* NX disk files, not RADS */
extern int a2_no_fixed_angle;		/* Use first azi, not fixed */

/******************************************************************************
 *	Multiple tape handling.
 *****************************************************************************/

extern short int a2_num_tapes;			/* number of input tapes */
extern short int vs_offset;			/* volume scan offset */

/******************************************************************************
 *	TDWR handling.
 *****************************************************************************/

extern int tdwr_rotate;				/* Angle to rotate the data */

/******************************************************************************
 *	Misc.
 *****************************************************************************/

extern int a2size;				/* sizeof structure */
extern unsigned short int a2size_scan;		/* sizeof structure */
extern unsigned short int a2size_site;		/* sizeof structure */

/******************************************************************************
 *	Fix for Australia archive bug.
 *****************************************************************************/

extern int aus_vs_bug;				/* enable bug fix */

/******************************************************************************
 *	China WSR-98D
 *****************************************************************************/

extern radar98;					/* Are we this radar format? */
