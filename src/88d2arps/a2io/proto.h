int	open_bzip2( char * );
int	read_bzip2( char *, int );


/*
 *	read_record_raw_disk_compress.c
 */

int	open_compress( char * );
int	read_compress( char *, int );

/*
 *	compress.c
 */

int	compress_handle_open( char * );
int	compress_handle_read( char *, int );

/*
 *	These are internal to "compress.c", so we don't want other to get
 *	complaints from others.
 */

#ifdef	COMPRESS_UNIQUE
static	int decompress( int, unsigned char *, int );
static	long int getcode();
#endif	/* COMPRESS_UNIQUE */


void	endian_header();
void	endian_data();

void	endian_scan();
void	endian_site();
void	endian_tape24();

float	endian_float_fix( float);

void	endian_msg_hdr();
void	endian_tdwr_basedata_hdr();
void	endian_tdwr_2b00();
void	endian_tdwr_2b01();

unsigned short tdwr_get_bits( unsigned short, int, int );


int	radar_open_raw_disk_file( char * );
int	read_record_raw_disk_file( char *s, int );
int	file_seek_data();
int	get_file_format( char *, unsigned char *, int );


int	open_gzip( char * );
int	read_gzip( char *, int );


int	open_ldm( char * );
int	read_ldm( char *, int );


int	set_memory( unsigned char *, int );
int	read_from_memory( unsigned char *, int );
int	seek_memory( int );

/*
 *	The rest of these are replacement calls for standard functions.
 *	Types are directory from the man pages.  The "const" are removed.
 */

int	open_memory( char *, int );
FILE	*fopen_memory( char *, char * );
FILE	*fdopen_memory( int, char * );
int	close_memory( int );
int	fclose_memory ( FILE * );
size_t	fread_memory( void *, size_t, size_t, FILE * );
int	fgetc_memory( FILE * );
int	ungetc_memory( int, FILE * );


int	open_mmap( char * );
int	close_mmap();
int	seek_mmap( int );
int	read_mmap( unsigned char *, int );


int get_azi();
int get_azi_a2io();
int get_azi_tdwr();

int get_beamwidth();
int get_beamwidth_a2io();
int get_beamwidth_tdwr();

int get_day();
int get_data_format();

int get_elev();
int get_elev_a2io();
int get_elev_tdwr();

int get_field_name( int, char * );
int get_field_num( char * );

int get_first_gate( int );
int get_first_gate_a2io( int );
int get_first_gate_tdwr( int );

int get_fixed_angle();
int get_fixed_angle_a2io();
int get_fixed_angle_tdwr();

int get_frequency();

int get_get_gate_spacing( int );
int get_get_gate_spacing_a2io( int );
int get_get_gate_spacing_tdwr( int );

int get_hour();
int get_min();
int get_month();

int get_nyquist();
int get_nyquist_a2io();
int get_nyquist_tdwr();

int get_number_of_gates( int );
int get_number_of_gates_a2io( int );
int get_number_of_gates_tdwr( int );

int get_polarization();
int get_power();
int get_prf();
int get_prt();
int get_proj_name( char * );
int get_pulse_width();

int get_radial_status();
int get_radial_status_a2io();
int get_radial_status_tdwr();

int get_restart_flag();
int get_scan();
int get_scan_status();

int get_scan_dir();
int get_scan_dir_a2io();
int get_scan_dir_tdwr();

int get_sec();

int get_seq();
int get_seq_a2io();
int get_seq_tdwr();

int get_status( int );
int get_status_a2io( int );
int get_status_tdwr( int );

int get_tilt();
int get_tilt_status();
int get_tilt_type();

int get_timestamp();
int get_timestamp_a2io();
int get_timestamp_tdwr();

int get_unix_timestamp();

int get_vcp();
int get_vcp_a2io();
int get_vcp_tdwr();

int get_vel_res();
int get_vel_res_a2io();
int get_vel_res_tdwr();

int get_vs_timestamp();
int get_year();

int radar_oops( char * );


int	set_tdwr_info();
int	decode_data_tdwr();
int	print_header();
int	get_data_field_raw_tdwr( int, unsigned short *, int );
int	get_tdwr_radial_status();


/*
 *	User routines.  These are the only ones that a user program should
 *	call.
 */

/*
 *	a2io.c
 */

int radar_open( char * );
int radar_init( char * );			/* backwards compatible */
int get_rt_mode();

/*
 *	find_key.c
 */

int find_key( char *, char * );

/*
 *	get_data.c
 */

int get_data_field( int, float *, int );
int get_data_field_raw( int, unsigned short *, int );
int get_data_field_uchar( int, unsigned char *, int );
int get_scale_values ( int, int *, int * );

/*
 *	get_record_init.c
 */

int get_record( char * );
int set_record_format( int );
int get_record_format();

/*
 *	radar.c 
 */

int get_azi();
int get_beamwidth();
int get_day();
int get_data_format();
int get_elev();

int get_field_name( int, char * );
int get_field_num( char * );
int get_first_gate( int );
int get_fixed_angle();
int get_frequency();

int get_get_gate_spacing( int );
int get_hour();
int get_min();
int get_month();
int get_nyquist();
int get_number_of_gates( int );

int get_polarization();
int get_power();
int get_prf();
int get_prt();
int get_proj_name( char * );
int get_pulse_width();

int get_radial_status();
int get_restart_flag();
int get_scan();
int get_scan_status();
int get_scan_dir();
int get_sec();
int get_seq();
int get_status( int );

int get_tilt();
int get_tilt_status();
int get_tilt_type();
int get_timestamp();

int get_unix_timestamp();
int get_vcp();
int get_vel_res();
int get_vs_timestamp();
int get_year();

/*
 *	radarinfo.c
 */

int get_altitude();
int get_latitude();
int get_latitude_dms( int *, int *, int * );
int get_longitude();
int get_longitude_dms( int *, int *, int * );
int get_radar_name ( char * );
int get_site_name( char * );

/*
 *	read_radial.c
 */

int read_radial();

/*
 *	read_record.c
 */

int read_record();
int set_console_messages( int );
