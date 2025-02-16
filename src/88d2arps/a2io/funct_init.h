/*
 *	Define a bunch of functions that are also variable names.
 */

#ifndef	_FUNCT_INIT_H
#define	_FUNCT_INIT_H

#pragma ident "@(#)funct_init.h	4.3	03/23/01 NSSL"

/*
 *	radar.c
 */

int	(*funct_get_azi)();
int	(*funct_get_beamwidth)();
int	(*funct_get_elev)();
int	(*funct_get_first_gate)();
int	(*funct_get_fixed_angle)();
int	(*funct_get_gate_spacing)();
int	(*funct_get_number_of_gates)();
int	(*funct_get_nyquist)();
int	(*funct_get_radial_status)();
int	(*funct_get_scan_dir)();
int	(*funct_get_seq)();
int	(*funct_get_status)();
int	(*funct_get_tilt)();
int	(*funct_get_timestamp)();
int	(*funct_get_vcp)();
int	(*funct_get_vel_res)();

int	get_azi_a2io();
int	get_beamwidth_a2io();
int	get_elev_a2io();
int	get_first_gate_a2io();
int	get_fixed_angle_a2io();
int	get_gate_spacing_a2io();
int	get_number_of_gates_a2io();
int	get_nyquist_a2io();
int	get_radial_status_a2io();
int	get_scan_dir_a2io();
int	get_seq_a2io();
int	get_status_a2io();
int	get_tilt_a2io();
int	get_timestamp_a2io();
int	get_vcp_a2io();
int	get_vel_res_a2io();

int	get_azi_tdwr();
int	get_beamwidth_tdwr();
int	get_elev_tdwr();
int	get_first_gate_tdwr();
int	get_fixed_angle_tdwr();
int	get_gate_spacing_tdwr();
int	get_number_of_gates_tdwr();
int	get_nyquist_tdwr();
int	get_radial_status_tdwr();
int	get_scan_dir_tdwr();
int	get_seq_tdwr();
int	get_status_tdwr();
int	get_tilt_tdwr();
int	get_timestamp_tdwr();
int	get_vcp_tdwr();
int	get_vel_res_tdwr();

/*
 *	get_data.c
 */

int	(*funct_get_data_field)();
int	(*funct_get_data_field_raw)();
int	(*funct_get_data_field_uchar)();
int	(*funct_get_scale_values)();

int	get_data_field_a2io();
int	get_data_field_raw_a2io();
int	get_data_field_uchar_a2io();
int	get_scale_values_a2io();

int	get_data_field_tdwr();
int	get_data_field_raw_tdwr();
int	get_data_field_uchar_tdwr();
int	get_scale_values_tdwr();

#endif	/* _FUNCT_INIT_H */
