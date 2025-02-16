/*
 * dorade.h --
 *
 *	This file declares structures and functions that access and manage data
 *	in DORADE sweep files.
 *
 *	The struct's are based on those found in the soloii source distribution
 *	from National Center for Atmospheric Research (NCAR), Earth Observing
 *	Laboratory, Research Data Program (RDP), downloaded from
 *	ftp.atd.ucar.edu.
 *
 * Copyright (c) 2007 Gordon D. Carrie
 *
 * Licensed under the Open Software License version 2.1
 *
 * Please send feedback to user0@tkgeomap.org
 *
 * @(#) $Id: dorade.h,v 4.0 2013/02/08 18:31:44 ywang Exp $
 *
 **********************************************************************
 *
 */

#ifndef DORADE_H_
#define DORADE_H_

#include <stdio.h>
#include "alloc.h"
#include "hash.h"
#include "radar.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Binary formats
 */

# define      DD_8_BITS 1
# define     DD_16_BITS 2
# define     DD_24_BITS 3
# define   DD_32_BIT_FP 4
# define   DD_16_BIT_FP 5

/*
 * Dorade scan modes
 */

enum Dorade_ScanMode {
    DORADE_CALIBRATION, DORADE_PPI, DORADE_COPLANE, DORADE_RHI,
    DORADE_VERTICAL_POINTING, DORADE_TARGET_MANUAL
};

/*
 * Structures of this type store information about a parameter.  The sweep
 * structure declared below has one of these for each parameter in the file.
 * The values mostly come from PARM descriptors.
 */

struct Dorade_ParmDesc {
    char use;			/* If true, this parameter is o.k. to use.
				 * If false, this parameter has been flagged
				 * for deletion. */
    char  name[9];		/* Name of parameter being described. */
    char  description[41];	/* Detailed description of this parameter. */
    char  units[9];		/* Units parameter is written in. */
    short interpulse_time;	/* Inter-pulse periods used.
				 * bits 0-1 = frequencies 1-2. */
    short xmitted_freq;		/* Frequencies used for this parameter. */
    float recvr_bandwidth;	/* Effective receiver bandwidth for this
				 * parameter, MHz.*/
    short pulse_width;		/* Effective pulse width of parameter, m. */
    short polarization;		/* Polarization of the radar beam, in na
				 * 0 Horizontal,
				 * 1 Vertical,
				 * 2 Circular,
				 * 3 Elliptical */
    short num_samples;		/* Number of samples used in estimate */
    short binary_format;	/* Binary format of radar data. */
    char  threshold_field[9];	/* Name of parameter upon which this
				 * parameter is thresholded (ascii
				 * characters NONE if not
				 * thresholded). */
    float threshold_value;	/* Value of threshold in ? */
    float scale;		/* Scale factor for parameter. */
    float bias;			/* Bias factor for parameter. */
    long  bad_data;		/* Bad data flag. */
    float min;			/* Smallest value seen in this file */
};

/*
 * Structures of this type store ray info and platform info for one ray.
 * The sweep structure declared below has one of these for each ray.
 */

struct Dorade_RayHdr {
    char good;			/* If true, ray is usable.  If false,
				 * ray should be skipped. */
    struct GeoTime_Jul time;	/* Ray time */
    Angle azimuth;		/* Azimuth (RYIB) */
    Angle elevation;		/* Elevation (RYIB) */
    float peak_power;		/* Last measured peak transmitted
				 * power in kw.  (RYIB) */
    Angle true_scan_rate;	/* Actual scan rate in angle units /second.
				 * (RYIB) */
    long  ray_status;		/* 0 = normal, 1 = transition, 2 = bad.
				 * (RYIB) */
    Angle longitude;		/* Antenna Longitude (Eastern
				 * Hemisphere is positive, West
				 * negative) (ASIB) */
    Angle latitude;		/* Antenna Latitude (Northern
				 * Hemisphere is positive, South
				 * Negative) (ASIB) */
    float altitude_msl;		/* Antenna Altitude above mean sea
				 * level (MSL) in m  (ASIB) */
    float altitude_agl;		/* Antenna Altitude above ground level
				 * (AGL) in m  (ASIB) */
    float ew_velocity;		/* Antenna east-west ground speed
				 * (towards East is positive) in m/sec
				 * (ASIB) */
    float ns_velocity;		/* Antenna north-south ground speed
				 * (towards North is positive) in m/sec
				 * (ASIB) */
    float vert_velocity;	/* Antenna vertical velocity (Up is
				 * positive) in m/sec  (ASIB) */
    Angle heading;		/* Antenna heading (angle between
				 * rotodome rotational axis and true
				 * North, clockwise (looking down)
				 * positive) (ASIB) */
    Angle roll;			/* Roll angle of aircraft tail section
				 * (Horizontal zero, Positive left wing up)
				 * (ASIB) */
    Angle pitch;		/* Pitch angle of rotodome (Horizontal
				 * is zero positive front up)
				 * (ASIB) */
    Angle drift_angle;		/* Antenna drift Angle. (angle between
				 * platform true velocity and heading,
				 * positive is drift more clockwise
				 * looking down) (ASIB) */
    Angle rotation_angle;	/* Angle of the radar beam with
				 * respect to the airframe (zero is
				 * along vertical stabilizer, positive
				 * is clockwise) (ASIB) */
    Angle tilt;			/* Angle of radar beam and line normal
				 * to longitudinal axis of aircraft,
				 * positive is towards nose of
				 * aircraft) (ASIB) */
    float ew_horiz_wind;	/* east - west wind velocity at the
				 * platform (towards East is positive)
				 * in m/sec  (ASIB) */
    float ns_horiz_wind;	/* North - South wind velocity at the
				 * platform (towards North is
				 * positive) in m/sec  (ASIB) */
    float vert_wind;		/* Vertical wind velocity at the
				 * platform (up is positive) in m/sec.
				 * (ASIB) */
    Angle heading_change;	/* Heading change rate in angle units/second.
				 * (ASIB) */
    Angle pitch_change;		/* Pitch change rate in angle units/second.
				 * (ASIB) */
};

/*
 * Structures of this type store the contents of a DORADE sweep file.
 * It assumes one collection of rays with the same set of parameters
 * from one sensor (radar).
 */

struct Dorade_Sweep {
    char radar_name[9];		/* (SSWB, RADD) */
    char proj_name[21];		/* Project number or name. (VOLD) */
    long n_parms;		/* Number of parameters.  (SSWB) */
    long n_use_parms;		/* Number of parameters in use (i.e. use flag
				 * in Dorade_ParmDesc struct is true) */
    short vol_num;		/* Volume number into current tape. (VOLD) */
    struct GeoTime_Jul time;	/* Time data taken (VOLD) */
    short n_sensors;		/* Total Number of sensor descriptors that
				 * follow. (VOLD) */
    float radar_const;		/* Radar/lidar constant in ??  (RADD) */
    float peak_power;		/* Typical peak power of the sensor in kw.
				 * Pulse energy is really the peak_power
				 * pulse_width (RADD) */
    float noise_power;		/* Typical noise power of the sensor in dBm.
				 * (RADD) */
    float receiver_gain;	/* Gain of the receiver in db. (RADD) */
    float antenna_gain;		/* Gain of the antenna in db.  (RADD) */
    float system_gain;		/* System Gain in db. (Ant G - WG loss)
				 * (RADD) */
    Angle horz_beam_width;	/* Horizontal beam width .
				 * beam divergence in milliradians
				 * is equivalent to beamwidth  (RADD) */
    Angle vert_beam_width;	/* Vertical beam width. (RADD) */
    short radar_type;		/* Radar Type (RADD)
				 * 0 = Ground,
				 * 1 = Airborne Fore,
				 * 2 = Airborne Aft,
				 * 3 = airborne Tail,
				 * 4 = Airborne Lower Fuselage,
				 * 5 = Shipborne. */
    enum Dorade_ScanMode scan_mode;	
    				/* Scan Mode (RADD)
				 * 0 = Calibration,
				 * 1 = PPI (constant Elevation)
				 * 2 = Coplane,
				 * 3 = RHI (Constant Azimuth),
				 * 4 = Vertical Pointing,
				 * 5 = Target (Stationary),
				 * 6 = Manual,
				 * 7 = Idle (Out of Control). */
    float req_rotat_vel;	/* Requested rotational velocity of
				 * the antenna in angle units/sec. (RADD) */
    short compression;		/* Compression flag */
    short data_reduction;	/* Data Reduction algorithm:
				 * 1 = none,
				 * 2 = between 2 angles,
				 * 3 = Between concentric circles,
				 * 4 = Above/below certain altitudes.
				 * (RADD) */
    float data_red_parm0;	/* 1 = smallest positive angle, degrees,
				 * 2 = inner circle diameter, km,
				 * 4 = minimum altitude, km.  (RADD) */
    float data_red_parm1;	/* 1 = largest positve angle, degrees,
				 * 2 = outer circle diameter, km,
				 * 4 = maximum altitude, km.  (RADD) */
    GeoPt radar_location;	/* Radar latitude/longitude (RADD) */
    float radar_altitude;	/* Altitude of radar above msl in km. (RADD) */
    float eff_unamb_vel;	/* Effective unambiguous velocity, m/s.
				 * (RADD) */
    float eff_unamb_range;	/* Effective unambiguous range, km. (RADD) */
    short num_freq_trans;	/* Number of frequencies transmitted. (RADD) */
    short num_ipps_trans;	/* Number of different inter-pulse periods
				 * transmitted.  (RADD) */
    float freq1;		/* Frequency 1.  (RADD) */
    float freq2;		/* Frequency 2.  (RADD) */
    float freq3;		/* Frequency 3.  (RADD) */
    float freq4;		/* Frequency 4.  (RADD) */
    float freq5;		/* Frequency 5.  (RADD) */
    float interpulse_per1;	/* Interpulse period 1.  (RADD) */
    float interpulse_per2;	/* Interpulse period 2.  (RADD) */
    float interpulse_per3;	/* Interpulse period 3.  (RADD) */
    float interpulse_per4;	/* Interpulse period 4.  (RADD) */
    float interpulse_per5;	/* Interpulse period 5.  (RADD) */
    unsigned n_cells;		/* Number of cells (CELV, also from PARM) */
    float *distP;		/* Distance from radar to cell, m,
				 * dimensions n_cells. (CELV)*/
    struct Dorade_ParmDesc *parmP;
    				/* Array of n_parms parameter descriptors
				 * (PARM) */
    struct Hash_Table parm_tbl;	/* This hash table maps parameter names
				 * to indeces in the parmP array */
    Angle azimuth_corr;		/* Correction added to azimuth (CFAC) */
    Angle elevation_corr;	/* Correction added to elevation (CFAC) */
    float range_delay_corr;	/* Correction used for range delay[m] (CFAC) */
    Angle longitude_corr;	/* Correction added to radar longitude (CFAC) */
    Angle latitude_corr;	/* Correction added to radar latitude (CFAC) */
    float pressure_alt_corr;	/* Correction added to pressure altitude
				 * (msl)[km] (CFAC) */
    float radar_alt_corr;	/* Correction added to radar altitude above
				 * ground level(agl) [km] (CFAC) */
    float ew_gndspd_corr;	/* Correction added to radar platform ground
				 * speed(E-W)[m/s] (CFAC) */
    float ns_gndspd_corr;	/* Correction added to radar platform ground
				 * speed(N-S)[m/s] (CFAC) */
    float vert_vel_corr;	/* Correction added to radar platform vertical
				 * velocity[m/s] (CFAC) */
    Angle heading_corr;		/* Correction added to radar platform heading
				 * (CFAC) */
    Angle roll_corr;		/* Correction added to radar platform roll
				 * (CFAC) */
    Angle pitch_corr;		/* Correction added to radar platform pitch
				 * (CFAC) */
    Angle drift_corr;		/* Correction added to radar platform drift
				 * (CFAC) */
    Angle rot_angle_corr;	/* Corrrection add to radar rotation angle
				 * (CFAC) */
    Angle tilt_corr;		/* Correction added to radar tilt angle
				 * (CFAC) */
    char swib_comment[9];	/* Comment (SWIB) */
    long  sweep_num;		/* Sweep number from the beginning of the volume
				 * (SWIB) */
    long  n_rays;		/* number of rays recorded in this sweep
				 * (SWIB) */
    long n_good_rays;		/* Number of rays flagged as "good"
				 * (See struct Dorade_RayHdr declaration) */
    Angle start_angle;		/* true start angle (SWIB) */
    Angle stop_angle;		/* true stop angle  (SWIB) */
    Angle fixed_angle;		/* (SWIB) */
    long  filter_flag;		/* (SWIB) */
    struct Dorade_RayHdr *rayHdrP;
				/* Array of n_rays ray headers */
    float ***dat;		/* Data, dimensioned
				 * [n_parms][n_rays][n_cells] */
    int interpolated;		/* If true, data is interpolated */
};

/*
 * Declarations for functions access and manage Dorade_Sweep struct's.
 */

void Dorade_InitSweep(struct Dorade_Sweep *);
float ***Dorade_AllocDat(long, long, unsigned);
void Dorade_FreeSweep(struct Dorade_Sweep *);
void Dorade_InitParm(struct Dorade_ParmDesc *);
void Dorade_InitRayHdr(struct Dorade_RayHdr *rayHdrP);
int Dorade_ReadSweepHeaders(FILE *, struct Dorade_Sweep *);
int Dorade_HashParms(struct Dorade_Sweep *);
struct Dorade_ParmDesc * Dorade_GetParm(struct Dorade_Sweep *, char *);
int Dorade_ReadSweep(FILE *, struct Dorade_Sweep *);
void Dorade_RmParm(struct Dorade_Sweep *, char *);
char *Dorade_SweepFName(char *, size_t *, struct Dorade_Sweep *);
int Dorade_WriteSweep(FILE *, struct Dorade_Sweep *);
int Dorade_ToASCII(FILE *, struct Dorade_Sweep *);
Angle Dorade_SweepElev(struct Dorade_Sweep *);
int Dorade_ShiftAz(struct Dorade_Sweep *, Angle);
int Dorade_IncrTime(struct Dorade_Sweep *, double);
int Dorade_CopyParm(struct Dorade_Sweep *, char *, char *, char *);
int Dorade_ParmPlus(struct Dorade_Sweep *, char *, char *);
int Dorade_ParmMinus(struct Dorade_Sweep *, char *, char *);
int Dorade_ParmMultiply(struct Dorade_Sweep *, char *, double);
int Dorade_ParmLog10(struct Dorade_Sweep *, char *);
int Dorade_SweepLatLonToBin(struct Dorade_Sweep *, GeoPt, int *, int *);
struct GeoLn *Dorade_SweepBinOutline(struct Dorade_Sweep *, unsigned, unsigned,
	unsigned);
int Dorade_SweepToGrnd(struct Dorade_Sweep *, unsigned,
	const struct Radar_PolarGrid *, struct Radar_PolarData *);
void Dorade_DestroyLib(void);
    
#ifdef __cplusplus
}
#endif

#endif
