/*
 *	a2io_message_31.h
 *
 *	Message type 31, valid Open Build 10.0 or later.
 */

#pragma ident "@(#)a2io_message_31.h	6.1	10/17/07 CAPS"

#include <inttypes.h>

/*
 *	It's another radial data type! 
 */

struct a2data31 {
	unsigned char radar_id[4];
	uint32_t zulu_time;
	unsigned short mod_julian_date;
	unsigned short azimuth_number;
	float azimuth;
	unsigned char compression_indicator;
	unsigned char spare;				/* not used */
	unsigned short radial_length;
	unsigned char azimuth_resolution_spacing;
	unsigned char radial_status;
	unsigned char elevation_number;
	unsigned char cut_sector_number;
	float elevation;
	unsigned char radial_split_blanking_status;
	unsigned char azimuth_indexing_mode;
	unsigned short data_block_count;
	uint32_t vol_ptr;
	uint32_t elv_ptr;
	uint32_t rad_ptr;
	uint32_t ref_ptr;
	uint32_t vel_ptr;
	uint32_t spw_ptr;
	uint32_t zdr_ptr;
	uint32_t phi_ptr;
	uint32_t rho_ptr;
} ;

struct a2data31_vol {
	unsigned char data_block_type;
	unsigned char data_name[3];
	unsigned short ltrup;
	unsigned char version_number_major;
	unsigned char version_number_minor;
	float latitude;
	float longitude;
	short int site_height;
	unsigned short feedhorn_height;
	float calibration_constant;
	float horiz_shvtx_power;
	float vert_shvtx_power;
	float sys_diff_refl;
	float init_sys_diff_phase;
	unsigned short vcp;
	unsigned short spare;
} ;

struct a2data31_elv {
	unsigned char data_block_type;
	unsigned char data_name[3];
	unsigned short ltrup;
	short int atmos;
	float calibration_constant;
} ;

struct a2data31_rad {
	unsigned char data_block_type;
	unsigned char data_name[3];
	unsigned short ltrup;
	short int unambiguous_range;
	float noise_level_horiz;
	float noise_level_vert;
	unsigned short nyquist;
	unsigned short spare;
} ;

struct a2data31_data {
	unsigned char data_block_type;
	unsigned char data_name[3];
	uint32_t reserved;
	unsigned short number_of_gates;
	short int first_gate;
	unsigned short gate_spacing;
	unsigned short tover;
	unsigned short snr_thresh;
	unsigned char control_flags;
	unsigned char data_word_size;
	float scale;
	float offset;
	unsigned char ptr;
};

/*
 *	Data types, this is the order we assume the data is in.  It may not
 *	be, so we have to map to this order.
 */

#define	VOL_31	1
#define	ELV_31	2
#define	RAD_31	3
#define	REF_31	4
#define	VEL_31	5
#define	SPW_31	6
#define	ZDR_31	7
#define	PHI_31	8
#define	RHO_31	9
