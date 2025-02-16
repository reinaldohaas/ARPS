/*
 *	endian.c			Fix Little Endian/Big Endian Problems.
 *
 *	Internal
 *	endian_header()			Fix 88d header message.
 *	endian_data()			Fix type 1 data message.
 *	endian_scan()			Fix 201 messages.
 *	endian_site()			Fix 202 messages.
 *	endian_tape24()			Fix 24-byte header records.
 *	endian_tdwr_msg_hdr()		Fix TDWR msg header.
 *	endian_tdwr_basedata_hdr()	Fix TDWR basedata message.
 *	*endian_tdwr_2b00()		Not implemented.
 *	*endian_tdwr_2b01()		Not implemented.
 *	tdwr_get_bits()			Internal
 *
 *	* = not used or tested, probably will go away
 *
 *	NOTE:  	  TDWR implementation is incomplete.
 *	WARNING:  TDWR implementation hasn't been tested since this file
 *		  was redesigned!!!
 */

#pragma ident	"@(#)endian.c	5.11	10/17/05"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/types.h>
#include <netinet/in.h>

#include "config.h"

#include "a2io.h"
#include "a2io_message_31.h"
#include "extern.h"

#include "proto.h"

#define	H	nssl_a2header
#define	D	nssl_a2data
#define	D31	nssl_a2data31
#define	D31_vol	nssl_a2data31_vol
#define	D31_rad	nssl_a2data31_rad
#define	D31_dbz	nssl_a2data31_dbz
#define	D31_vel	nssl_a2data31_vel
#define	D31_spw	nssl_a2data31_spw
#define	D31_zdr	nssl_a2data31_zdr
#define	D31_phi	nssl_a2data31_phi
#define	D31_rho	nssl_a2data31_rho

#ifdef	ALLOW_TDWR
#include "tdwrio.h"
#include "tdwr_extern.h"

#define	TH	tdwr_pmsg
#define	TB	tdwr_phdr
#define	T0	tdwr_p2b00
#define	T1	tdwr_p2b01
#define	T2	tdwr_p2b02
#endif


//----------------------for IBM/AIX------------------
swap_short(unsigned short tran_data)    //transfer unsighed short variables
{
 unsigned short conx,d1,d2,d3,d4,dout;
 dout=0;
 conx=~(~0<<8);
 d1=conx&tran_data;
 tran_data=tran_data>>8;
 d2=conx&tran_data;
 d1=d1<<8;
 dout=d1|d2;
return(dout);
}
 
swap_long(unsigned long tran_data)     //transfer unsighed long variables
{
 unsigned long conx,d1,d2,d3,d4,dout;
 dout=0;
 conx=~(~0<<8);
 d1=conx&tran_data;
 tran_data=tran_data>>8;
 d2=conx&tran_data;
 tran_data=tran_data>>8;
 d3=conx&tran_data;
 tran_data=tran_data>>8;
 d4=conx&tran_data;
 d1=d1<<24;
 d2=d2<<16;
 d3=d3<<8;
 dout=d1|d2|d3|d4;
return(dout);
}
 
void swap_header()
{
#ifdef IBM4
        H->message_size                 = swap_short( H->message_size );
        H->seq                          = swap_short( H->seq );
        H->julian_date                  = swap_short( H->julian_date );
        H->milsec                       = swap_long( H->milsec );
        H->num_mes_seg                  = swap_short( H->num_mes_seg );
        H->cur_mes_seg                  = swap_short( H->cur_mes_seg );
#endif
 
        return;
}
 
void swap_data()
{
#ifdef IBM4
        D->zulu_time                    = swap_long( D->zulu_time );
        D->mod_julian_date              = swap_short( D->mod_julian_date );
        D->unamb_range                  = swap_short( D->unamb_range );
        D->azimuth                      = swap_short( D->azimuth );
        D->azimuth_number               = swap_short( D->azimuth_number );
        D->radial_status                = swap_short( D->radial_status );
        D->elevation                    = swap_short( D->elevation );
        D->elevation_number             = swap_short( D->elevation_number );
        D->range_first_gate_dbz         = swap_short( D->range_first_gate_dbz );
        D->range_first_gate_vel         = swap_short( D->range_first_gate_vel );
        D->gate_size_dbz                = swap_short( D->gate_size_dbz );
        D->gate_size_vel                = swap_short( D->gate_size_vel );
        D->num_gates_dbz                = swap_short( D->num_gates_dbz );
        D->num_gates_vel                = swap_short( D->num_gates_vel );
        D->cut_sector_number            = swap_short( D->cut_sector_number );
        D->calibration_constant         = swap_long( D->calibration_constant );
        D->dbz_ptr                      = swap_short( D->dbz_ptr );
        D->vel_ptr                      = swap_short( D->vel_ptr );
        D->spw_ptr                      = swap_short( D->spw_ptr );
        D->vel_res                      = swap_short( D->vel_res );
        D->vcp                          = swap_short( D->vcp );
        D->dummy[0]                     = swap_short( D->dummy[0] );
        D->dummy[1]                     = swap_short( D->dummy[2] );
        D->dummy[2]                     = swap_short( D->dummy[2] );
        D->dummy[3]                     = swap_short( D->dummy[3] );
#ifdef  VERY_OLD_FLAWED_DATA
        D->arc_dbz_ptr                  = swap_short( D->arc_dbz_ptr );
        D->arc_vel_ptr                  = swap_short D->arc_vel_ptr );
        D->arc_spw_ptr                  = swap_short( D->arc_spw_ptr );
#endif
        D->nyquist                      = swap_short( D->nyquist );
        D->atmos                        = swap_short( D->atmos );
        D->tover                        = swap_short( D->tover );
#endif

        return;
}
//-------------------------for IBM/AIX ------------------------------


void endian_header()
{
	H->message_size	 		= htons( H->message_size );
	H->seq				= htons( H->seq );
	H->julian_date			= htons( H->julian_date );
	H->milsec			= htonl( H->milsec );
	H->num_mes_seg			= htons( H->num_mes_seg );
	H->cur_mes_seg			= htons( H->cur_mes_seg );

	return;
}

void endian_data1()
{
	D->zulu_time			= htonl( D->zulu_time );
	D->mod_julian_date		= htons( D->mod_julian_date );
	D->unamb_range			= htons( D->unamb_range );
	D->azimuth			= htons( D->azimuth );
	D->azimuth_number		= htons( D->azimuth_number );
	D->radial_status		= htons( D->radial_status );
	D->elevation			= htons( D->elevation );
	D->elevation_number		= htons( D->elevation_number );
	D->range_first_gate_dbz		= htons( D->range_first_gate_dbz );
	D->range_first_gate_vel		= htons( D->range_first_gate_vel );
	D->gate_size_dbz		= htons( D->gate_size_dbz );
	D->gate_size_vel		= htons( D->gate_size_vel );
	D->num_gates_dbz		= htons( D->num_gates_dbz );
	D->num_gates_vel		= htons( D->num_gates_vel );
	D->cut_sector_number		= htons( D->cut_sector_number );
	D->calibration_constant		= htonl( D->calibration_constant );
	D->dbz_ptr			= htons( D->dbz_ptr );
	D->vel_ptr			= htons( D->vel_ptr );
	D->spw_ptr			= htons( D->spw_ptr );
	D->vel_res			= htons( D->vel_res );
	D->vcp	 			= htons( D->vcp );
	D->dummy[0]			= htons( D->dummy[0] );
	D->dummy[1]			= htons( D->dummy[2] );
	D->dummy[2]			= htons( D->dummy[2] );
	D->dummy[3]			= htons( D->dummy[3] );
#ifdef	VERY_OLD_FLAWED_DATA
	D->arc_dbz_ptr			= htons( D->arc_dbz_ptr );
	D->arc_vel_ptr			= htons( D->arc_vel_ptr );
	D->arc_spw_ptr			= htons( D->arc_spw_ptr );
#endif
	D->nyquist			= htons( D->nyquist );
	D->atmos			= htons( D->atmos );
	D->tover			= htons( D->tover );

	return;
}

void endian_data31()
{
	D31->zulu_time			= htonl( D31->zulu_time );
	D31->mod_julian_date		= htons( D31->mod_julian_date );
	D31->azimuth_number		= htons( D31->azimuth_number );

	D31->azimuth			= endian_float_fix(D31->azimuth);
	D31->radial_length		= htons( D31->radial_length );
	D31->elevation			= endian_float_fix(D31->elevation);

	D31->data_block_count		= htons( D31->data_block_count );

	D31->vol_ptr			= htonl( D31->vol_ptr);
	D31->elv_ptr			= htonl( D31->elv_ptr);
	D31->rad_ptr			= htonl( D31->rad_ptr);
	D31->ref_ptr			= htonl( D31->ref_ptr);
	D31->vel_ptr			= htonl( D31->vel_ptr);
	D31->spw_ptr			= htonl( D31->spw_ptr);
	D31->zdr_ptr			= htonl( D31->zdr_ptr);
	D31->phi_ptr			= htonl( D31->phi_ptr);
	D31->rho_ptr			= htonl( D31->rho_ptr);

	return;
}

void endian_data31_vol()
{
	D31_vol->latitude		= endian_float_fix(D31_vol->latitude);
	D31_vol->longitude		= endian_float_fix(D31_vol->longitude);
	D31_vol->vcp			= htons( D31_vol->vcp );

	return;
}

void endian_data31_elv()
{
	return;
}

void endian_data31_rad()
{
	D31_rad->nyquist		= htons( D31_rad->nyquist );
	return;
}

void endian_data31_data(n)
int n;
{
	if ( n == 0 )
	{
		D31_dbz->number_of_gates= htons(D31_dbz->number_of_gates);
		D31_dbz->first_gate	= htons(D31_dbz->first_gate);
		D31_dbz->gate_spacing	= htons(D31_dbz->gate_spacing);
		D31_dbz->tover		= htons(D31_dbz->tover);
		D31_dbz->snr_thresh	= htons(D31_dbz->snr_thresh);
		D31_dbz->scale		= endian_float_fix(D31_dbz->scale);
		D31_dbz->offset		= endian_float_fix(D31_dbz->offset);
#ifdef	DEBUG
printf("DBZ:  %f %f\n",D31_dbz->scale,D31_dbz->offset);
#endif
	}
	if ( n == 1 )
	{
		D31_vel->number_of_gates= htons(D31_vel->number_of_gates);
		D31_vel->first_gate	= htons(D31_vel->first_gate);
		D31_vel->gate_spacing	= htons(D31_vel->gate_spacing);
		D31_vel->tover		= htons(D31_vel->tover);
		D31_vel->snr_thresh	= htons(D31_vel->snr_thresh);
		D31_vel->scale		= endian_float_fix(D31_vel->scale);
		D31_vel->offset		= endian_float_fix(D31_vel->offset);
#ifdef	DEBUG
printf("VEL:  %f %f\n",D31_vel->scale,D31_vel->offset);
#endif
	}
	if ( n == 2 )
	{
		D31_spw->number_of_gates= htons(D31_spw->number_of_gates);
		D31_spw->first_gate	= htons(D31_spw->first_gate);
		D31_spw->gate_spacing	= htons(D31_spw->gate_spacing);
		D31_spw->tover		= htons(D31_spw->tover);
		D31_spw->snr_thresh	= htons(D31_spw->snr_thresh);
		D31_spw->scale		= endian_float_fix(D31_spw->scale);
		D31_spw->offset		= endian_float_fix(D31_spw->offset);
#ifdef	DEBUG
printf("SPW:  %f %f\n",D31_spw->scale,D31_spw->offset);
#endif
	}
	if ( n == 3 )
	{
		D31_zdr->number_of_gates= htons(D31_zdr->number_of_gates);
		D31_zdr->first_gate	= htons(D31_zdr->first_gate);
		D31_zdr->gate_spacing	= htons(D31_zdr->gate_spacing);
		D31_zdr->tover		= htons(D31_zdr->tover);
		D31_zdr->snr_thresh	= htons(D31_zdr->snr_thresh);
		D31_zdr->scale		= endian_float_fix(D31_zdr->scale);
		D31_zdr->offset		= endian_float_fix(D31_zdr->offset);
#ifdef	DEBUG
printf("ZDR:  %f %f\n",D31_zdr->scale,D31_zdr->offset);
#endif
	}
	if ( n == 4 )
	{
		D31_phi->number_of_gates= htons(D31_phi->number_of_gates);
		D31_phi->first_gate	= htons(D31_phi->first_gate);
		D31_phi->gate_spacing	= htons(D31_phi->gate_spacing);
		D31_phi->tover		= htons(D31_phi->tover);
		D31_phi->snr_thresh	= htons(D31_phi->snr_thresh);
		D31_phi->scale		= endian_float_fix(D31_phi->scale);
		D31_phi->offset		= endian_float_fix(D31_phi->offset);
#ifdef	DEBUG
printf("PHI:  %f %f\n",D31_phi->scale,D31_phi->offset);
#endif
	}
	if ( n == 5 )
	{
		D31_rho->number_of_gates= htons(D31_rho->number_of_gates);
		D31_rho->first_gate	= htons(D31_rho->first_gate);
		D31_rho->gate_spacing	= htons(D31_rho->gate_spacing);
		D31_rho->tover		= htons(D31_rho->tover);
		D31_rho->snr_thresh	= htons(D31_rho->snr_thresh);
		D31_rho->scale		= endian_float_fix(D31_rho->scale);
		D31_rho->offset		= endian_float_fix(D31_rho->offset);
#ifdef	DEBUG
printf("RHO:  %f %f\n",D31_rho->scale,D31_rho->offset);
#endif
	}
	return;
}

void endian_scan()
{

	nssl_a2scan->vsn		= htons( nssl_a2scan->vsn );

/*
 *	nssl_a2scan->pad is just a filler.
 */

	return;
}

void endian_site()
{
	nssl_a2site->reset_tm		= htonl( nssl_a2site->reset_tm );

	return;
}

void endian_tape24()
{
	nssl_a2tape24->julian_date	= htons( nssl_a2tape24->julian_date );
	nssl_a2tape24->milsec		= htonl( nssl_a2tape24->milsec );

	return;
}

#ifdef	ALLOW_TDWR

void endian_tdwr_msg_hdr()
{

	TH->message_id			= htons( TH->message_id );
	TH->message_length		= htons( TH->message_length );

	return;
}

void endian_tdwr_basedata_hdr()
{
	unsigned short *ptr;
	unsigned short tmp;

	ptr = (unsigned short *) TB;

	TB->volume_scan_count		= htons( TB->volume_scan_count );
	ptr[1]				= htons( ptr[1] );

	tmp = ptr[1];

	TB->end_of_volume_scan		= tdwr_get_bits( tmp, 0, 1 );
	TB->start_of_volume_scan	= tdwr_get_bits( tmp, 1, 1 );
	TB->dummy1			= tdwr_get_bits( tmp, 2, 6 );
	TB->scan_strategy		= tdwr_get_bits( tmp, 8, 8 );

	TB->peak_transmitter_power	= htons( TB->peak_transmitter_power );

	ptr[3]				= htons( ptr[3] );

	tmp = ptr[3];

	TB->dummy2			= tdwr_get_bits( tmp, 0, 13 );
	TB->dummy_record_indicator	= tdwr_get_bits( tmp, 13, 1 );
	TB->start_of_live_data		= tdwr_get_bits( tmp, 14, 1 );
	TB->start_of_playback_data	= tdwr_get_bits( tmp, 15, 1 );

	ptr[4]				= htons( ptr[4] );

	tmp = ptr[4];

	TB->tilt_number			= tdwr_get_bits( tmp, 0, 8 );
	TB->end_of_elevation_scan	= tdwr_get_bits( tmp, 8, 1 );
	TB->start_of_elevation_scan	= tdwr_get_bits( tmp, 9, 1 );
	TB->clutter_residue_map_number	= tdwr_get_bits( tmp, 10, 3 );
	TB->reserved1			= tdwr_get_bits( tmp, 13, 3 );

	ptr[5]				= htons( ptr[5] );

	tmp = ptr[5];

	TB->reserved2			= tdwr_get_bits( tmp, 0, 3 );
	TB->obscuration_flagging_indicator	= tdwr_get_bits( tmp, 3, 1 );
	TB->spike_removal_indicator	= tdwr_get_bits( tmp, 4, 1 );
	TB->velocity_dealiasing_scan	= tdwr_get_bits( tmp, 5, 2 );
	TB->sector_scan			= tdwr_get_bits( tmp, 7, 1 );
	TB->microburst_aloft_scan	= tdwr_get_bits( tmp, 8, 1 );
	TB->precip_scan_and_resolution	= tdwr_get_bits( tmp, 9, 2 );
	TB->wind_shift_scan		= tdwr_get_bits( tmp, 11, 1 );
	TB->low_elevation_scan		= tdwr_get_bits( tmp, 12, 1 );
	TB->microburst_surface_scan	= tdwr_get_bits( tmp, 13, 1 );
	TB->gust_front_scan		= tdwr_get_bits( tmp, 14, 1 );
	TB->low_prf_scan		= tdwr_get_bits( tmp, 15, 1 );

	TB->current_elevation		= 
				endian_float_fix( TB->current_elevation);
	TB->angular_scan_rate		= 
				endian_float_fix( TB->angular_scan_rate);
	TB->pulse_repetition_interval	= htons
				( TB->pulse_repetition_interval );

	ptr[11]				= htons( ptr[11] );

	tmp = ptr[11];

	TB->dwell_id			= tdwr_get_bits( tmp, 0, 2 );
	TB->solar_indicator		= tdwr_get_bits( tmp, 2, 1 );
	TB->dummy3			= tdwr_get_bits( tmp, 3, 1 );
	TB->pulses_per_dwell		= tdwr_get_bits( tmp, 4, 12 );

	TB->final_range_sample		= htons( TB->final_range_sample );
	TB->range_samples_per_dwell	= htons( TB->range_samples_per_dwell );
	TB->current_azimuth		= 
				endian_float_fix( TB->current_azimuth);
	TB->total_noise_power		= 
				endian_float_fix( TB->total_noise_power);
	TB->timestamp			= htonl( TB->timestamp );
	TB->base_data_type		= htons( TB->base_data_type );

	ptr[21]				= htons( ptr[21] );

	tmp = ptr[21];

	TB->dummy4			= tdwr_get_bits( tmp, 0, 11 );
	TB->wind_field_model_initializer= tdwr_get_bits( tmp, 11, 1 );
	TB->incomplete_volume_scan	= tdwr_get_bits( tmp, 12, 1 );
	TB->incomplete_elevation_scan	= tdwr_get_bits( tmp, 13, 1 );
	TB->volume_scan_restart		= tdwr_get_bits( tmp, 14, 1 );
	TB->elevation_scan_restart	= tdwr_get_bits( tmp, 15, 1 );

	TB->integer_azimuth		= htons( TB->integer_azimuth );
	TB->load_shed_final_sample	= htons( TB->load_shed_final_sample );
	return;
}

/*
 *	Some of the 0x2b00 data is short ints.
 */

void endian_tdwr_2b00()
{
	unsigned short *ptr;
	unsigned short tmp;
	int i;

	for( i=0; i<TDWR_MAX_2B00; i++ )
	{
		T0->tdwr_data[i].vel_raw = htons( T0->tdwr_data[i].vel_raw);

		ptr 			= (unsigned short *) &T0->tdwr_data[i];
		ptr[2] 			= htons( ptr[2]);

		tmp 			= ptr[2];

		T0->tdwr_data[i].spw	= tdwr_get_bits( tmp, 0, 8 );
		T0->tdwr_data[i].dummy	= tdwr_get_bits( tmp, 8, 2 );

		T0->tdwr_data[i].compressed_dealias_algorithm_failure_flag
					= tdwr_get_bits( tmp, 10, 3 );

		T0->tdwr_data[i].compressed_point_target_filter_flag
					= tdwr_get_bits( tmp, 13, 1 );

		T0->tdwr_data[i].compressed_conditioned_valid_velocity_flag
					= tdwr_get_bits( tmp, 14, 1 );

		T0->tdwr_data[i].compressed_conditioned_valid_flag 
					= tdwr_get_bits( tmp, 15, 1 );

		T0->tdwr_data[i].vel_dealias = 
			htons( T0->tdwr_data[i].vel_dealias);
	}
	return;
}

/*
 *	Nothing to do.
 */

void endian_tdwr_2b01()
{
	return;
}

unsigned short tdwr_get_bits( unsigned short num, int start, int kount )
{
	unsigned short tmp;

	tmp = num << start;
	tmp = tmp >> (16 - kount);

	return( tmp );
}

#endif

union e {
	long x;
	float y;
} e;

/*
 *	A htonf() routine doesn't exist, so we write our own.
 */

float endian_float_fix( float f )
{
	e.y = f;
	e.x = htonl( e.x );
	return( e.y );
}
