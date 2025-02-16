/*
 * radar --
 *
 * 	This file declares general purpose functions and macros to use with
 * 	radar data.
 * 
 * Copyright (c) 2007 Gordon D. Carrie
 *
 * Licensed under the Open Software License version 2.1
 *
 * Please send feedback to user0@tkgeomap.org
 *
 * @(#) $Id: radar.h,v 4.0 2013/02/08 18:31:45 ywang Exp $
 *
 **********************************************************************
 *
 */

#ifndef RADAR_H_
#define RADAR_H_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include "geoLines.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Speed of light in meters per second
 */

#define SPEED_OF_LIGHT 2.9979e8

/*
 * The following functions tell whether or not a value is usable.
 */

float Radar_NoData(void);
int Radar_ValIsNoData(float v);
int Radar_ValIsData(float v);
int Radar_PrintVal(FILE *out, float v, char *dat, char *nodat);

/*
 * These functions perform ray geometry calculations.
 */

double Radar_GrndRng(double d, double h, Angle tilt);
double Radar_RayRng(double d, double h, Angle tilt);

/*
 * Array search functions
 */

int Radar_ArrSearch(float *xx, int n, float x);
int Radar_AzArrSearch(Angle *xx, int n, Angle x);

/*
 * Interpolation functions
 */

typedef float (Radar_InterpProc)(Angle *, float *, float **, Angle, float);
Radar_InterpProc Radar_BilinInterp;

/*
 * This structure stores a polar grid.
 */

struct Radar_PolarGrid {
    unsigned n_az;	/* Number of azimuths */
    Angle *az;		/* Array of n_az azimuths */
    unsigned n_rng;	/* Number of range gates (or cells) */
    float *rng;		/* Ranges (meters) */
};

/*
 * This data structure stores polar gridded data.
 */

struct Radar_PolarData {
    unsigned n_az;	/* Number of azimuths */
    Angle *az;		/* Azimuth values */
    unsigned n_rng;	/* Number of ranges */
    float *rng;		/* Range values */
    GeoPt ctr;		/* Origin */
    float **dat;	/* Data, dimensioned [n_az][n_rng] */
};

void Radar_PolarDataInit(struct Radar_PolarData *rpdP);
int Radar_PolarDataSetAlloc(struct Radar_PolarData *rpdP, GeoPt ctr,
	unsigned n_az, unsigned n_rng);
void Radar_PolarDataFree(struct Radar_PolarData *rpdP);
    
#ifdef __cplusplus
}
#endif

#endif
