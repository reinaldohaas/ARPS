/*
 *	TDWR specific routines.
 */

#pragma ident "@(#)radar_tdwr.c	5.1	10/10/01 CAPS"

#include <stdio.h>
#include <time.h>

#include "config.h"

#include "a2defs.h"		/* NEED TO CHANGE THIS NAME! */
#include "const.h"
#include "extern.h"
#include "proto.h"

#ifdef	ALLOW_TDWR
#include "tdwrio.h"
#include "tdwr_extern.h"
#include "tdwr_nyq.h"
#endif

get_azi_tdwr()
{
	int tdwr_azi;

#ifdef	ALLOW_TDWR
	if ( TDWR_DATA_RADIAL )
	{
		tdwr_azi = (int)
			( tdwr_phdr->current_azimuth * 100. + 0.5 );
		if ( !tdwr_rotate )
			return( tdwr_azi);
		else
		{
			tdwr_azi += tdwr_rotate;

			while( tdwr_azi >= 36000 )
				tdwr_azi -= 36000;
			while ( tdwr_azi < 0 )
				tdwr_azi += 36000;
			return( tdwr_azi );
		}
	}
	else
#endif
		return( MISSING );
}

/*
 *	This is just a constant for 88Ds.
 */

get_beamwidth_tdwr()
{
	return( MISSING );
}


/*
 *	Current elevation angle.  Use get_fixed_angle() to get the fixed
 *	elevation angle.
 *
 */

get_elev_tdwr()
{
#ifdef	ALLOW_TDWR
	if ( TDWR_DATA_RADIAL )
		return ( (int) ( tdwr_phdr->current_elevation * 100. ) );
	else
#endif
		return( MISSING );
}

/*
 *	Gate spacing is uniform for all TDWR fields.
 */

get_first_gate_tdwr( n )
int n;
{
#ifdef	ALLOW_TDWR
	if ( TDWR_DATA_RADIAL )
		return( 150 );
	else
#endif
		return( MISSING );
}

/*
 *	Fixed elevation angle.  Used get_elev() to get the elevation angle
 *	of the current radial.
 */

get_fixed_angle_tdwr()
{
#ifdef	ALLOW_TDWR
	if ( TDWR_DATA_RADIAL )
		return( get_elev() );
	else
#endif
		return ( MISSING );
}

get_gate_spacing_tdwr(n)
{
#ifdef	ALLOW_TDWR
	if ( TDWR_DATA_RADIAL )
		return( 150 );
	else
#endif
		return( MISSING );
}

get_nyquist_tdwr()
{
#ifdef	ALLOW_TDWR
	if ( get_status( VEL ) == 0)
		return( MISSING );

	ptdwr_nyq = tdwr_nyq;

	while ( ptdwr_nyq->pri != 0 )
	{
		if ( tdwr_phdr->pulse_repetition_interval 
			== ptdwr_nyq->pri )
			return( (int) ptdwr_nyq->nyq );
		ptdwr_nyq++;
	}
#endif
	return( MISSING );
}

get_number_of_gates_tdwr(n)
int n;
{
#ifdef	ALLOW_TDWR
	if ( TDWR_DATA_RADIAL )
	{
		if ( tdwr_pmsg->message_id == 0x2b00 )
			return( TDWR_MAX_2B00 );

		if ( tdwr_pmsg->message_id == 0x2b01 )
		{
			if ( n == VEL || n == RVEL || n == SPW )
				return( 0 );
			return( TDWR_MAX_2B01 );
		}

/*
 *	Not possible, but we'll put it here anyways.
 */

		return( 0 ); 
	}
	else
#endif
		return( 0 );
}

get_radial_status_tdwr()
{
#ifdef	ALLOW_TDWR
	return( get_tdwr_radial_status() );
#else
	return( MISSING );
#endif
}

get_scan_dir_tdwr()
{
	return( MISSING );
}

get_seq_tdwr()
{
	return( MISSING );
}

get_status_tdwr(n)
int n;
{
#ifdef	ALLOW_TDWR
	if ( TDWR_DATA_RADIAL )
	{
		if ( n == VEL || n == SPW || n == RVEL )
		{
			if ( tdwr_pmsg->message_id == 0x2b01 )
				return( 0 );
			else
				return( 1 );
		}

		if ( n == DBZ || n == SNR )
			return( 1 );

		return( 0 );
	}
	else
#endif
		return( 0 );

}

get_tilt_tdwr()
{
#ifdef	ALLOW_TDWR
	if ( TDWR_DATA_RADIAL )
	{
		return( tdwr_phdr->tilt_number );
	}
	else
#endif
		return( MISSING );
}

get_timestamp_tdwr()
{
#ifdef	ALLOW_TDWR
	return( tdwr_phdr->timestamp );
#else
	return( MISSING );
#endif
}

get_vcp_tdwr()
{
	return( 0 );
	
}

get_vel_res_tdwr()
{
	return( 0 );
}
