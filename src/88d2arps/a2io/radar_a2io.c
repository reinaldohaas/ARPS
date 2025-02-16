/*
 *	A2IO specific routines.
 */
 
#include "config.h"

#pragma ident "@(#)radar_a2io.c	5.1	10/09/01 CAPS - range to first gate values are not modified"

#include <stdio.h>
#include <time.h>
#include "a2io.h"
#include "a2io_message_31.h"
#include "const.h"
#include "extern.h"
#include "proto.h"

get_azi_a2io()
{
	if ( message_type == A2_DATA_TYPE )
	  return( (int) ( (float) (nssl_a2data->azimuth) * A_SCALE + 0.5 ) );
	if ( message_type == A2_DIGITAL_TYPE )
	  return( (int) ( nssl_a2data31->azimuth * 100.0 ) );
	return(MISSING);
}

/*
 *	This is just a constant for 88Ds.
 */

get_beamwidth_a2io()
{
	return( 95 );

}

/*
 *	Current elevation angle.  Use get_fixed_angle() to get the fixed
 *	elevation angle.
 *
 */

get_elev_a2io()
{
	if ( message_type == A2_DATA_TYPE )
	  return( (int) ( (float) (nssl_a2data->elevation) * A_SCALE + 0.5 ) );
	if ( message_type == A2_DIGITAL_TYPE )
	  return( (int) ( nssl_a2data31->elevation * 100.0 ) );
	return(MISSING);
}

get_first_gate_a2io( n )
int n;
{
	if ( message_type == A2_DATA_TYPE )
	{
		if (get_status(n) == 0)
			return( MISSING );
		if ( n == DBZ)
			return( nssl_a2data->range_first_gate_dbz );
		else if ( n == VEL || n == SPW || n == RVEL )
			return( nssl_a2data->range_first_gate_vel );
		else
			return( 0 );
	}
	if ( message_type == A2_DIGITAL_TYPE )
	{
		if ( get_status(n) == 0)
			return(MISSING);
		if ( n == DBZ )
			return(nssl_a2data31_dbz->first_gate);
		else if ( n == VEL )
			return(nssl_a2data31_vel->first_gate);
		else if ( n == SPW )
			return(nssl_a2data31_spw->first_gate);
		else if ( n == ZDR )
			return(nssl_a2data31_zdr->first_gate);
		else if ( n == PHI )
			return(nssl_a2data31_phi->first_gate);
		else if ( n == RHO )
			return(nssl_a2data31_rho->first_gate);
		return(0);
	}
	return(MISSING);
}

/*
 *	Fixed elevation angle.  Used get_elev() to get the elevation angle
 *	of the current radial.
 *
 *	For message type 31, just return the current elevation angle.
 */

get_fixed_angle_a2io()
{
	int magic;
	int afix;

	if (message_type == A2_DATA_TYPE )
	{
/*
 *	If we don't have a valid VCP number, then just send back the current
 *	elevation angle.
 *
 *	Ditto if the user requested it.
 */

		if ( nssl_vcp == NULL || a2_no_fixed_angle )
			return( get_elev() );

		magic = (nssl_a2data->elevation_number - 1) * 17 + 6;
		afix = (float) nssl_vcp[magic] * A_SCALE + 0.5 ;

		return ( (int) afix );
	}

	if ( message_type == A2_DIGITAL_TYPE )
		return( get_elev() );
	return(MISSING);
}

get_gate_spacing_a2io(n)
{
	if (get_status(n) == 0)
		return( MISSING );

	if ( message_type == A2_DATA_TYPE )
	{
		if (n == DBZ)
			return( nssl_a2data->gate_size_dbz );
		else if ( n == VEL || n == SPW || n == RVEL )
			return( nssl_a2data->gate_size_vel );
		else
			return( 0 );
	}
	if ( message_type == A2_DIGITAL_TYPE )
	{
		if ( get_status(n) == 0)
			return(MISSING);
		if ( n == DBZ )
			return(nssl_a2data31_dbz->gate_spacing);
		else if ( n == VEL )
			return(nssl_a2data31_vel->gate_spacing);
		else if ( n == SPW )
			return(nssl_a2data31_spw->gate_spacing);
		else if ( n == ZDR )
			return(nssl_a2data31_zdr->gate_spacing);
		else if ( n == PHI )
			return(nssl_a2data31_phi->gate_spacing);
		else if ( n == RHO )
			return(nssl_a2data31_rho->gate_spacing);
		return(0);
	}
	return(MISSING);
}

get_nyquist_a2io()
{
	if (get_status( VEL ) == 0)
		return( MISSING );
	if ( message_type == A2_DATA_TYPE )
		return( (int) nssl_a2data->nyquist );
	if ( message_type == A2_DIGITAL_TYPE )
	{
/*
 *	Will eventually need to make this location independent.
 */
		if ( nssl_a2data31->rad_ptr > 0 )
			return(nssl_a2data31_rad->nyquist );
		else
			return(MISSING);
	}
	return(MISSING);
}

get_number_of_gates_a2io(n)
int n;
{
	if ( message_type == A2_DATA_TYPE )
	{
		if (n == DBZ)
			return( nssl_a2data->num_gates_dbz );
		else if (n == VEL || n == SPW || n == RVEL )
			return( nssl_a2data->num_gates_vel );
		else
			return( 0 );
	}
	if ( message_type == A2_DIGITAL_TYPE )
	{
		if ( get_status(n) == 0)
			return(MISSING);
		if ( n == DBZ )
			return(nssl_a2data31_dbz->number_of_gates);
		else if ( n == VEL )
			return(nssl_a2data31_vel->number_of_gates);
		else if ( n == SPW )
			return(nssl_a2data31_spw->number_of_gates);
		else if ( n == ZDR )
			return(nssl_a2data31_zdr->number_of_gates);
		else if ( n == PHI )
			return(nssl_a2data31_phi->number_of_gates);
		else if ( n == RHO )
			return(nssl_a2data31_rho->number_of_gates);
		return(0);
	}
	return(MISSING);
}

get_radial_status_a2io()
{
	if ( message_type == A2_DATA_TYPE )
		return( nssl_a2data->radial_status );
	if ( message_type == A2_DIGITAL_TYPE )
		return( nssl_a2data31->radial_status );
	return(MISSING);
}


/*
 *	Obsolete.
 */

get_scan_dir_a2io()
{
	return( 0 );
}

get_seq_a2io()
{
	return( nssl_a2header->seq );
}

get_status_a2io(n)
int n;
{
	int val;

/*
 *	Make sure we are a data message first!
 */

	if ( message_type == A2_DATA_TYPE )
	{
		if (n == DBZ)
			val = nssl_a2data->num_gates_dbz;
		else if (n == VEL || n == SPW || n == RVEL)
			val = nssl_a2data->num_gates_vel;
		else
			val = 0;
	
		if (val > 0)
			return( 1 );
		else
			return( 0 );
	}

	if ( message_type == A2_DIGITAL_TYPE )
	{
		if (n == DBZ)
			val = map31[REF_31];
		else if ( n == VEL )
			val = map31[VEL_31];
		else if ( n == SPW )
			val = map31[SPW_31];
		else if (n == ZDR )
			val = map31[ZDR_31];
		else if (n == PHI )
			val = map31[PHI_31];
		else if (n == RHO )
			val = map31[RHO_31];
		else
			val = 0;
	
		if (val > 0)
			return( 1 );
		else
			return( 0 );
	}

	return( 0 );

}

get_tilt_a2io()
{
	if ( message_type == A2_DATA_TYPE )
		return( nssl_a2data->elevation_number );

	if ( message_type == A2_DIGITAL_TYPE )
		return( nssl_a2data31->elevation_number );
	return(MISSING);
}

get_timestamp_a2io()
{
	return( zulu_time );
}

get_vcp_a2io()
{
	if ( message_type == A2_DATA_TYPE )
		return( nssl_a2data->vcp );
	if ( message_type == A2_DIGITAL_TYPE )
		return( nssl_a2data31_vol->vcp );
	return(MISSING);
}

get_vel_res_a2io()
{
	if ( message_type == A2_DATA_TYPE )
		return( nssl_a2data->vel_res );
	if ( message_type == A2_DIGITAL_TYPE )
		return(-888);
	return(MISSING);
}
