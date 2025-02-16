/*
 * radar --
 *
 *	This file defines some utility functions and macros to use with radar
 *	data.
 * 
 * Copyright (c) 2007 Gordon D. Carrie
 *
 * Licensed under the Open Software License version 2.1
 *
 * Please send feedback to user0@tkgeomap.org
 *
 * @(#) $Id: radar.c,v 4.0 2013/02/08 18:31:45 ywang Exp $
 *
 **********************************************************************
 *
 */

#include <stdio.h>
#include <float.h>
#include <limits.h>
#include "errMsg.h"
#include "alloc.h"
#include "radar.h"

/*
 * This value indicates no data at a radar bin.
 */

#define NO_DATA -FLT_MAX

/*
 *------------------------------------------------------------------------
 *
 * Radar_NoData --
 *
 *	This function returns the value that indicates no data at a bin.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 * 	None.
 *
 *------------------------------------------------------------------------
 */

float Radar_NoData(void)
{
    return NO_DATA;
}

/*
 *------------------------------------------------------------------------
 *
 * Radar_ValIsNoData --
 *
 * 	This function determines whether a value is data or not.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 * 	None.
 *
 *------------------------------------------------------------------------
 */

int Radar_ValIsNoData(float v)
{
    return v == NO_DATA;
}

/*
 *------------------------------------------------------------------------
 *
 * Radar_ValIsData --
 *
 * 	This function determines whether a value is data or not.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 * 	None.
 * 	
 *------------------------------------------------------------------------
 */

int Radar_ValIsData(float v)
{
    return (v != NO_DATA);
}

/*
 *------------------------------------------------------------------------
 *
 * Radar_PrintVal --
 *
 * 	This convenience function prints a float value if possible.
 *
 * Arguments:
 *	FILE *out	- output stream
 *	float v		- value to print
 *	char *dat_fmt	- format to use in fprintf if v is data.
 *	char *nodat_fmt	- format to use in fprintf if v is not data.
 *
 * Results:
 * 	Return value is the return value of fprintf.
 *
 * Side effects:
 * 	The output stream is modified.
 * 	
 *------------------------------------------------------------------------
 */

int Radar_PrintVal(FILE *out, float v, char *dat_fmt, char *nodat_fmt)
{
    if (Radar_ValIsData(v)) {
	return fprintf(out, dat_fmt, v);
    } else {
	return fprintf(out, nodat_fmt);
    }
}

/*
 *------------------------------------------------------------------------
 *
 * Radar_GrndRng --
 *
 *	This function computes distance over ground given distance along a ray
 *	and ray tilt angle.
 *
 * Arguments:
 *	double d	- Distance along ray in meters
 *	double h	- Antenna height in meters
 *	Angle tilt	- Ray tilt
 *
 * Results:
 * 	See the user documentation.
 *
 * Side effects:
 * 	None.
 *
 *------------------------------------------------------------------------
 */

double Radar_GrndRng(double d, double h, Angle tilt)
{
    double t = AngleToRad(tilt);
    return REarth() * atan(d * cos(t) / (REarth() + h + d * sin(t)));
}

/*
 *------------------------------------------------------------------------
 *
 * Radar_RayRng --
 *
 *	This function computes distance along a ray given distance over ground 
 *	and ray tilt angle.
 *
 * Arguments:
 *	double d	- Distance along ground from radar to
 *			  point of interest, in meters 
 *	double h	- Antenna height in meters 
 *	Angle tilt	- Ray tilt 
 *
 * Results:
 * 	See the user documentation.
 *
 * Side effects:
 * 	None.
 *
 *------------------------------------------------------------------------
 */

double
Radar_RayRng(double d, double h, Angle tilt)
{
    double delta = d / REarth();
    return (REarth() + h) * sin(delta) / cos(AngleToRad(tilt) + delta);
}

/*
 *------------------------------------------------------------------------
 *
 * Radar_ArrSearch --
 *
 * 	This function searches a monotonic array.
 *	Ref: Press, et al., Numerical Recipes in C, 2nd Edition.  p. 117.
 *
 * Arguments:
 *	float *xx	- Array to search
 *	int n		- Number of elements in xx
 *	float x		- Value to search for
 *
 * Results:
 * 	See the user documentation.
 *
 * Side effects:
 * 	None.
 *
 *------------------------------------------------------------------------
 */

int Radar_ArrSearch(float *xx, int n, float x)
{
    int jl;		/* Index for lower bound */
    int jm;		/* Index for midpoint in bisection search */
    int ju;		/* Index for upper bound */

    if (Radar_ValIsNoData(x)) {
	return -1;
    }

    /*
     * Bisect to find the interval that contains the value.
     */

    jl = 0;
    ju = n - 1;
    if (xx[n - 1] > xx[0]) {
	/*
	 * xx is increasing
	 */

	if (x < xx[0] || x > xx[n - 1]) {
	    return -1;
	}
	while (ju - jl > 1) {
	    jm = (jl + ju) / 2;
	    if (x > xx[jm]) {
		jl = jm;
	    } else {
		ju = jm;
	    }
	}
    } else {
	/*
	 * xx is decreasing
	 */

	if (x > xx[0] || x < xx[n - 1]) {
	    return -1;
	}
	while (ju - jl > 1) {
	    jm = (jl + ju) / 2;
	    if (x > xx[jm]) {
		ju = jm;
	    } else {
		jl = jm;
	    }
	}
    }
    return jl;
}

/*
 *------------------------------------------------------------------------
 *
 * Radar_AzArrSearch --
 *
 * 	This function searches a monotonic array of angles.
 *	Ref: Press, et al., Numerical Recipes in C, 2nd Edition.  p. 117.
 *
 * Arguments:
 *	Angle *xx	- Array to search
 *	int n		- Number of elements in xx
 *	Angle x		- Value to search for
 *
 * Results:
 * 	See the user documentation.
 *
 * Side effects:
 * 	None.
 *
 *------------------------------------------------------------------------
 */

int Radar_AzArrSearch(Angle *xx, int n, Angle x)
{
    int jl;		/* Index for lower bound */
    int jm;		/* Index for midpoint in bisection search */
    int ju;		/* Index for upper bound */

    /*
     * Bisect to find the interval that contains the value.
     */

    jl = 0;
    ju = n - 1;
    if (LonCmp(xx[1], xx[0]) == East) {
	/*
	 * xx is increasing
	 */

	if (LonCmp(x, xx[0]) == West || LonCmp(x, xx[n - 1]) == East) {
	    return -1;
	}
	while (ju - jl > 1) {
	    jm = (jl + ju) / 2;
	    if (LonCmp(x, xx[jm]) == East) {
		jl = jm;
	    } else {
		ju = jm;
	    }
	}
    } else {
	/*
	 * xx is decreasing
	 */

	if (LonCmp(x, xx[0]) == East || LonCmp(x, xx[n - 1]) == West) {
	    return -1;
	}
	while (ju - jl > 1) {
	    jm = (jl + ju) / 2;
	    if (LonCmp(x, xx[jm]) == East) {
		ju = jm;
	    } else {
		jl = jm;
	    }
	}
    }
    return jl;
}

/*
 *------------------------------------------------------------------------
 *
 * Radar_BilinInterp --
 *
 * 	This is a bilinear interpolation function.  It handles azimuths
 * 	correctly.
 *
 * Arguments:
 *	Angle *az0	- Azimuths from 0.0 degrees clockwise to
 *			  points with known value, dimensioned 2.
 *	float *r0	- Ranges from origin to points with known value
 *			  dimensioned 2.
 *	float **z0	- Known values, dimensioned [2][2].
 *			  First index increments azimuth,
 *			  second index increments range.
 *	Angle az1	- Azimuth from 0.0 degrees clockwise to point for
 *			  which value is sought
 *	float r1	- Range from origin to point for which value is
 *			  sought.  Points are assumed to be arranged as
 *			  follows:
 *
 *			  	     |                |
 *			  	     |                |
 *			 r0[1]-- z0[0][1] ------ z0[1][1]
 *			  	     |                |
 *			 r1          |       *        |
 *			  	     |                |
 *			  	     |                |
 *			 r0[0]-- z0[0][0] ------ z0[1][0]
 *			  	     |                |
 *			  	     |                |
 *			   	  az0[0]    az1    az0[1]
 *
 *
 * Results:
 * 	See the user documentation.
 *
 * Side effects:
 * 	None.
 *
 *------------------------------------------------------------------------
 */

float Radar_BilinInterp(Angle *az0, float *r0, float **z0, Angle az1, float r1)
{
    double t, u;	/* Computational constants.
			 * See Press, et al.,
			 * Numerical Recipes in C, 2nd Edition.
			 * p. 123. */

    /*
     * Put angles within 180 degrees of az0[0].
     */

    az1 = DomainLon(az1, az0[0]);
    az0[1] = DomainLon(az0[1], az0[0]);

    if (Radar_ValIsNoData(z0[0][0])
	    || Radar_ValIsNoData(z0[0][1])
	    || Radar_ValIsNoData(z0[1][0])
	    || Radar_ValIsNoData(z0[1][1])
	    || fabs(az0[1] - az0[0]) < FLT_EPSILON
	    || fabs(r0[1] - r0[0]) < FLT_EPSILON) {
	return Radar_NoData();
    }
    t = (double)(az1 - az0[0]) / (double)(az0[1] - az0[0]);
    u = (r1 - r0[0]) / (r0[1] - r0[0]);
    return (1 - t) * (1 - u) * z0[0][0] + t * (1 - u) * z0[1][0]
	+ t * u * z0[1][1] + (1 - t) * u * z0[0][1];
}

/*
 *------------------------------------------------------------------------
 *
 * Radar_PolarDataInit --
 *
 * 	This constructor initializes a Radar_PolarData structure with bogus
 * 	values.
 *
 * Results:
 * 	None.
 *
 * Side effects:
 * 	Members in the structure are set to bogus values.
 *
 *------------------------------------------------------------------------
 *
 */

void
Radar_PolarDataInit(struct Radar_PolarData *rpdP)
{
    rpdP->ctr = GeoPtNowhere();
    rpdP->n_az = 0;
    rpdP->az = NULL;
    rpdP->n_rng = 0;
    rpdP->rng = NULL;
    rpdP->dat = NULL;
}

/*
 *------------------------------------------------------------------------
 *
 * Radar_PolarDataSetAlloc --
 *
 * 	This constructor sets values and allocates memory in a
 * 	Radar_PolarData structure.
 *
 * Arguments:
 *	struct Radar_PolarData *rpdP	- Structure to initialize
 *	GeoPt ctr			- Origin
 *	unsigned n_az			- Number of azimuths
 *	unsigned n_rng			- Number of ranges
 *
 * Results:
 * 	Members of a Radar_PolarData structure receive values.
 * 	Previous contents are assumed to be garbage.
 * 	Return value is true if function succeeds, otherwise false.
 *
 * Side effects:
 * 	Memory is allocated.  It should eventually be freed with a call to
 * 	Radar_PolarDataFree.
 *
 *------------------------------------------------------------------------
 */

int Radar_PolarDataSetAlloc(struct Radar_PolarData *rpdP, GeoPt ctr,
	unsigned n_az, unsigned n_rng)
{
    Angle *az = NULL;
    float *rng = NULL;
    float **dat = NULL;
    unsigned i;
    float *d, *de;

    az = (Angle *)MALLOC(n_az * sizeof(float));
    rng = (float *)MALLOC(n_rng * sizeof(float));
    dat = (float **)MALLOC(n_az * sizeof(float *));
    dat[0] = (float *)MALLOC(n_az * n_rng * sizeof(float));
    for (i = 1; i < n_az; i++) {
	dat[i] = dat[i - 1] + n_rng;
    }
    for (d = *dat, de = d + n_az * n_rng; d < de; d++) {
	*d = Radar_NoData();
    }
    rpdP->ctr = ctr;
    rpdP->n_az = n_az;
    rpdP->az = az;
    rpdP->n_rng = n_rng;
    rpdP->rng = rng;
    rpdP->dat = dat;
    return 1;
}

/*
 *------------------------------------------------------------------------
 *
 * Radar_PolarDataFree --
 *
 * 	This is the destructor for a Radar_PolarData structure.
 *
 *------------------------------------------------------------------------
 */

void Radar_PolarDataFree(struct Radar_PolarData *rpdP)
{
    FREE(rpdP->az);
    FREE(rpdP->rng);
    if (rpdP->dat) {
	if (rpdP->dat[0]) {
	    FREE(rpdP->dat[0]);
	}
	FREE(rpdP->dat);
    }
    Radar_PolarDataInit(rpdP);
}
