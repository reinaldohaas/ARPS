/*
 *	Internal routines
 */

#pragma ident "@(#)funct_init.c	5.1	10/10/01 CAPS"

#include <stdio.h>

#include "config.h"

#include "funct_init.h"
#include "proto.h"

funct_init_a2io()
{
	funct_get_azi		= get_azi_a2io;
	funct_get_beamwidth	= get_beamwidth_a2io;
	funct_get_elev		= get_elev_a2io;
	funct_get_first_gate	= get_first_gate_a2io;
	funct_get_fixed_angle	= get_fixed_angle_a2io;
	funct_get_gate_spacing	= get_gate_spacing_a2io;
	funct_get_number_of_gates
				= get_number_of_gates_a2io;
	funct_get_nyquist	= get_nyquist_a2io;
	funct_get_radial_status	= get_radial_status_a2io;
	funct_get_scan_dir	= get_scan_dir_a2io;
	funct_get_seq		= get_seq_a2io;
	funct_get_status	= get_status_a2io;
	funct_get_tilt		= get_tilt_a2io;
	funct_get_timestamp	= get_timestamp_a2io;
	funct_get_vcp		= get_vcp_a2io;
	funct_get_vel_res	= get_vel_res_a2io;

	funct_get_data_field	= get_data_field_a2io;
	funct_get_data_field_raw
				= get_data_field_raw_a2io;
	funct_get_data_field_uchar
				= get_data_field_uchar_a2io;
	funct_get_scale_values	= get_scale_values_a2io;

	return;
}

#ifdef	ALLOW_TDWR

funct_init_tdwr()
{
	funct_get_azi		= get_azi_tdwr;
	funct_get_beamwidth	= get_beamwidth_tdwr;
	funct_get_elev		= get_elev_tdwr;
	funct_get_first_gate	= get_first_gate_tdwr;
	funct_get_fixed_angle	= get_fixed_angle_tdwr;
	funct_get_gate_spacing	= get_gate_spacing_tdwr;
	funct_get_number_of_gates
				= get_number_of_gates_tdwr;
	funct_get_nyquist	= get_nyquist_tdwr;
	funct_get_radial_status	= get_radial_status_tdwr;
	funct_get_scan_dir	= get_scan_dir_tdwr;
	funct_get_seq		= get_seq_tdwr;
	funct_get_status	= get_status_tdwr;
	funct_get_tilt		= get_tilt_tdwr;
	funct_get_timestamp	= get_timestamp_tdwr;
	funct_get_vcp		= get_vcp_tdwr;
	funct_get_vel_res	= get_vel_res_tdwr;

	funct_get_data_field	= get_data_field_tdwr;
	funct_get_data_field_raw
				= get_data_field_raw_tdwr;
	funct_get_data_field_uchar
				= get_data_field_uchar_tdwr;
	funct_get_scale_values	= get_scale_values_tdwr;

	return;
}

#endif
