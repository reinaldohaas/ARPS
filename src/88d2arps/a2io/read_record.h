/*
 *	read_record.h		read_record() parameters
 */

#pragma ident "@(#)read_record.h	5.2	08/30/04 CAPS"

/******************************************************************************
 *	Data and pointers to data.
 *****************************************************************************/

struct a2header	*nssl_a2header;			/* header part of record */
struct a2data	*nssl_a2data;			/* ptr for A2_DATA_TYPE */
struct a2data31	*nssl_a2data31;			/* ptr for A2_DIGITAL_TYPE */

struct a2data31_vol	*nssl_a2data31_vol;
struct a2data31_elv	*nssl_a2data31_elv;
struct a2data31_rad	*nssl_a2data31_rad;
struct a2data31_data	*nssl_a2data31_dbz;
struct a2data31_data	*nssl_a2data31_vel;
struct a2data31_data	*nssl_a2data31_spw;
struct a2data31_data	*nssl_a2data31_zdr;
struct a2data31_data	*nssl_a2data31_phi;
struct a2data31_data	*nssl_a2data31_rho;

struct a2data31_data	a2data31_dbz;
struct a2data31_data	a2data31_vel;
struct a2data31_data	a2data31_spw;
struct a2data31_data	a2data31_zdr;
struct a2data31_data	a2data31_phi;
struct a2data31_data	a2data31_rho;

struct a2scan	*nssl_a2scan;			/* ptr for A2_VOLSCAN_TYPE */
struct a2site	*nssl_a2site;			/* ptr for A2_SITE_TYPE */
struct a2tape24	*nssl_a2tape24;			/* ptr for A2_VS_24_BYTE */
struct ncdc_header *ncdc_header;		/* ptr for A2_ARCHIVE2_NCDC */

char nssl_buf[MAX_RECORD_SIZE_TAPE];		/* where the data goes */
char *nssl_ptr;					/* beginning of record */
char *nssl_dptr;				/* data part of record */
int map31[10];					/* mapping of data */

unsigned short int *nssl_vcp;			/* ptr to vcp structure */

/******************************************************************************
 *	This data gets passed back thru user calls.
 *****************************************************************************/

int nssl_scan_status;				/* one when new volume scan */
int nssl_tilt_status;				/* one when new tilt */

int nssl_restart_flag = 0;			/* tilt restart flag */

int nssl_scan_num = VS_SEARCH;			/* current volume scan num */

int nssl_vs_timestamp = 0;			/* from A2_SITE_TYPE record */

/******************************************************************************
 *	Input/output/data format.
 *****************************************************************************/

int input_format = UNKNOWN;			/* input format */
int output_format;				/* output format */
int data_format = A2_FORMAT;			/* actual data format type */
int disk_format = DISK;				/* DISK_* in a2io.h */
int file_format = UNKNOWN;			/* FILE_* in a2io.h */

/******************************************************************************
 *	Bookkeeping information.
 *****************************************************************************/

int left = 0;					/* bytes left to process */
int new = 0;					/* byte loc in "nssl_buf" */

int message_type;				/* message type */
int message_length;				/* length of the message */

/******************************************************************************
 *	Time handling.
 *****************************************************************************/

struct tm *nssl_gmt;				/* time structure */
struct tm nssl_gmt_data;			/* dedicated space for time */
long zulu_time;					/* time conversion */

/******************************************************************************
 *	Data thresholding.
 *****************************************************************************/

int allow_thresholds = 0;			/* allow user thresholds */
int spw_threshold;				/* spw threshold */

/******************************************************************************
 *	Elevation error handling.
 *****************************************************************************/

unsigned short elev_save = 0;		/* same type as a2data->elevation */
unsigned short elev_thresh;		/* ditto */

/******************************************************************************
 *	Tape volume number handling, LDM starts at 0.
 *****************************************************************************/

int a2_start_scan = 0;				/* starting volume scan */
int a2_end_scan = 9999;				/* ending volume scan */
/*
 *	Make this one a large number, so that LDM doesn't fail after 9999
 *	volume scans (see code in read_record_tape()).
 */
int a2_num_scan = 999999;			/* number of vs's to process */

/******************************************************************************
 *	Misc radar.
 *****************************************************************************/

int a2_no_rads = 0;				/* NX disk files, not RADS */
int a2_no_fixed_angle = 0;			/* Use first azi, not fixed */

/******************************************************************************
 *	Multiple tape handling.
 *****************************************************************************/

short int a2_num_tapes = 1;			/* number of input tapes */
short int vs_offset = 0;			/* volume scan offset */

/******************************************************************************
 *	TDWR handling.
 *****************************************************************************/

int tdwr_rotate = 0;				/* Angle to rotate the data */

/******************************************************************************
 *	Misc.
 *****************************************************************************/

int a2size;					/* sizeof structure */
unsigned short int a2size_scan;			/* sizeof structure */
unsigned short int a2size_site;			/* sizeof structure */

/******************************************************************************
 *	Fix for Australia archive bug.
 *****************************************************************************/

int aus_vs_bug = 0;				/* enable bug fix */

/******************************************************************************
 *	China WSR-98D
 *****************************************************************************/

int radar98 = 0;				/* Are we this radar format? */
