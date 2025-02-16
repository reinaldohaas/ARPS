/*
 *	Define a bunch of functions that are also variable names.
 */

#ifndef	_EXTERN_FUNCT_H
#define	_EXTERN_FUNCT_H

#pragma ident "@(#)extern_funct.h	5.1	10/10/01 CAPS"

#include "config.h"

/*
 *	radar.c
 */

extern int	(*funct_get_azi)();
extern int	(*funct_get_beamwidth)();
extern int	(*funct_get_elev)();
extern int	(*funct_get_first_gate)();
extern int	(*funct_get_fixed_angle)();
extern int	(*funct_get_gate_spacing)();
extern int	(*funct_get_number_of_gates)();
extern int	(*funct_get_nyquist)();
extern int	(*funct_get_radial_status)();
extern int	(*funct_get_scan_dir)();
extern int	(*funct_get_seq)();
extern int	(*funct_get_status)();
extern int	(*funct_get_tilt)();
extern int	(*funct_get_timestamp)();
extern int	(*funct_get_vcp)();
extern int	(*funct_get_vel_res)();

extern int	get_azi_a2io();
extern int	get_beamwidth_a2io();
extern int	get_elev_a2io();
extern int	get_first_gate_a2io();
extern int	get_fixed_angle_a2io();
extern int	get_gate_spacing_a2io();
extern int	get_number_of_gates_a2io();
extern int	get_nyquist_a2io();
extern int	get_radial_status_a2io();
extern int	get_scan_dir_a2io();
extern int	get_seq_a2io();
extern int	get_status_a2io();
extern int	get_tilt_a2io();
extern int	get_timestamp_a2io();
extern int	get_vcp_a2io();
extern int	get_vel_res_a2io();

#ifdef	ALLOW_TDWR
extern int	get_azi_tdwr();
extern int	get_beamwidth_tdwr();
extern int	get_elev_tdwr();
extern int	get_first_gate_tdwr();
extern int	get_fixed_angle_tdwr();
extern int	get_gate_spacing_tdwr();
extern int	get_number_of_gates_tdwr();
extern int	get_nyquist_tdwr();
extern int	get_radial_status_tdwr();
extern int	get_scan_dir_tdwr();
extern int	get_seq_tdwr();
extern int	get_status_tdwr();
extern int	get_tilt_tdwr();
extern int	get_timestamp_tdwr();
extern int	get_vcp_tdwr();
extern int	get_vel_res_tdwr();
#endif

/*
 *	get_data.c
 */

extern int	(*funct_get_data_field)();
extern int	(*funct_get_data_field_raw)();
extern int	(*funct_get_data_field_uchar)();
extern int	(*funct_get_scale_values)();

extern int	get_data_field_a2io();
extern int	get_data_field_raw_a2io();
extern int	get_data_field_uchar_a2io();
extern int	get_scale_values_a2io();

#ifdef	ALLOW_TDWR
extern int	get_data_field_tdwr();
extern int	get_data_field_raw_tdwr();
extern int	get_data_field_uchar_tdwr();
extern int	get_scale_values_tdwr();
#endif

#endif	/* _EXTERN_FUNCT_H */
