/*
 *	Tdwrio.h
 *
 *	TDWR parameters...
 */

#pragma ident "@(#)tdwrio.h	5.1	10/10/01 CAPS"

/*
 *	TDWR Constants
 */

#include <sys/types.h>
#include <inttypes.h>
#include "tdwrdefs.h"

#define	TDWR_MAX_2B00		600
/*
 *	There are 1984 range gates, however, the gate spacing is 150 meters
 *	with 135km, then 300 meters afterwards.  ARPS can't handle this, so
 *	we'll truncate at 135km.
 *
#define	TDWR_MAX_2B01		1984
 */
#define	TDWR_MAX_2B01		900

#define	TDWR_NORMAL_SIZE	6144
#define	TDWR_2B00_SIZE		4852
#define	TDWR_2B01_SIZE		6004
#define	TDWR_2B02_SIZE		1024
#define	TDWR_2C00_SIZE		6144		/* CHECK THIS!!! */
#define	TDWR_2C01_SIZE		1236

/*
 *	The size of the record is 64 byte.  Tapes records are apparently
 *	padded, while realtime isn't.
 */

#define	TDWR_2C02_SIZE		64
#define	TDWR_2C02_TAPE_SIZE	1280

#define	TDWR_4206_SIZE		124

/*
 *	Check to make sure this isn't a header or LLWAA record.
 */

#define	TDWR_DATA_RADIAL	\
	 ( tdwr_pmsg->message_id == 0x2b00 || tdwr_pmsg->message_id == 0x2b01 )

/*
 *	The data.
 */

struct tdwr_msg_hdr {
	u_short	message_id;
	u_short	message_length;
};

struct tdwr_basedata_hdr {
	u_short	volume_scan_count;
	u_short	end_of_volume_scan:1;
	u_short	start_of_volume_scan:1;
	u_short dummy1:6;
	u_short	scan_strategy:8;

	u_short	peak_transmitter_power;
	u_short	dummy2:13;
	u_short	dummy_record_indicator:1;
	u_short	start_of_live_data:1;
	u_short	start_of_playback_data:1;

	u_short	tilt_number:8;
	u_short	end_of_elevation_scan:1;
	u_short	start_of_elevation_scan:1;
	u_short	clutter_residue_map_number:3;
	u_short	reserved1:3;			/* do in two pieces because */
	u_short	reserved2:3;			/* it crosses a word boundary */
	u_short	obscuration_flagging_indicator:1;
	u_short spike_removal_indicator:1;
	u_short velocity_dealiasing_scan:2;
	u_short sector_scan:1;
	u_short	microburst_aloft_scan:1;
	u_short precip_scan_and_resolution:2;
	u_short wind_shift_scan:1;
	u_short low_elevation_scan:1;
	u_short microburst_surface_scan:1;
	u_short gust_front_scan:1;
	u_short low_prf_scan:1;

	float	current_elevation;

	float	angular_scan_rate;

	u_short	pulse_repetition_interval;
	u_short dwell_id:2;
	u_short	solar_indicator:1;
	u_short dummy3:1;
	uint32_t pulses_per_dwell:12;

	u_short final_range_sample;
	u_short range_samples_per_dwell;

	float	current_azimuth;

	float	total_noise_power;

	uint32_t	timestamp;

	u_short	base_data_type;
	u_short dummy4:11;
	u_short	wind_field_model_initializer:1;
	u_short	incomplete_volume_scan:1;
	u_short	incomplete_elevation_scan:1;
	u_short	volume_scan_restart:1;
	u_short	elevation_scan_restart:1;

	u_short	integer_azimuth;
	u_short	load_shed_final_sample;
};

struct tdwr_prod_2b00 {
	struct	tdwr_data {
		u_char	dbz;
		u_char	snr;
		u_short	vel_raw;
		u_char	spw;
		u_short	dummy:2;
		u_short	compressed_dealias_algorithm_failure_flag:3;
		u_short	compressed_point_target_filter_flag:1;
		u_short	compressed_conditioned_valid_velocity_flag:1;
		u_short	compressed_conditioned_valid_flag:1;
		u_short	vel_dealias;
	} tdwr_data[600];

	char filler[1292];
};

struct tdwr_prod_2b01 {
	u_char	dbz[1984];
	u_char	snr[1984];

	struct tdwr_flags {
		u_char dummy1:1;
		u_char compressed_valid_flag:1;
		u_char dummy2:3;
		u_char compressed_point_target_filter_flag:1;
		u_char dummy3:2;
	} tdwr_flags[1984];
	char filler[140];
};

struct tdwr_prod_2b02 {
	char	base_data_date_stamp[8];
	char	base_data_time_stamp[8];
	char	base_data_header_text[64];
	char	filler[940];
};

/*
 *	The following three definitions aren't used anywhere, and are leftover
 *	from an earlier generation of this software.  They are left in as
 *	handy debugging tools.
 */

#ifdef	TDWR_DEBUG

#define	TDWR_VEL(value)	( (float) value * 0.25 - 80.0 )
#define	TDWR_SNR(value)	( (float) value * 0.50 )
#define	TDWR_SPW(value)	( (float) value * 0.25 )

#endif

struct tdwr_prod_2b02 *tdwr_prod_2b02;
