/*
 * geography.h --
 *
 *	This header file declares general structures, constants, and functions
 *	for basic geographic calculations.  See the Geography(3) man page
 *	for more information.
 * 
 * Copyright (c) 2007 Gordon D. Carrie.  All rights reserved.
 * 
 * Licensed under the Open Software License version 2.1
 * 
 * Please address questions and feedback to user0@tkgeomap.org
 *
 * @(#) $Id: geography.h,v 4.0 2013/02/08 18:31:45 ywang Exp $
 *
 ********************************************
 *
 */


/*
 * Note:
 * Default distance units (unless otherwise noted):
 *   great circle degrees on sphere
 *   meters on Earth's surface
 *   meters on plane
 */

#ifndef _GEOGRAPHY_H_
#define _GEOGRAPHY_H_

#include <stdio.h>
#include <float.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

/* 
 * Distance units
 */

#define NMI 1852.0
#define SMI 1609.344

/*
 * Angle conversion factors
 */

#define MPERDEG (REarth() * M_PI / 180.0)
#define DEGPERM (180.0 / (REarth() * M_PI))
#define NMIPERDEG (MPERDEG / NMI)
#define SMIPERDEG (MPERDEG / SMI)
#define KMPERDEG (MPERDEG / 1000.0)
#define RADPERDEG 0.017453292519943294892
#define DEGPERRAD 57.29577951308232088
#define DEG_INV 0.002777777777777778

/*
 * The following data type stores angle measurements.
 */

typedef int Angle;

/*
 * This structure stores a geopoint - point on the Earth's surface.
 */

typedef struct
{
  Angle			lat;		/* Latitude */
  Angle			lon;		/* Longitude */
} GeoPt;

/*
 * This structure stores a map-point - a point on a 2D map.
 */

typedef struct
{
  float			abs;		/* Abscissa */
  float			ord;		/* Ordinate */
} MapPt;

/*
 * A point in a 3D Cartesian coordinate system with origin at the Earth's,
 * center, in which Earth has a radius of 1.0.
 */

typedef struct
{
  double		x;
  double		y;
  double		z;
} CartPt;

/*
 * Global functions
 */

double	REarth (void);
void	SetREarth (double r);
Angle	AngleFmDeg (double deg);
double	AngleToDeg (Angle a);
Angle	AngleFmRad (double rad);
double	AngleToRad (Angle a);
Angle	BadAngle (void);
int	AngleIsOK (Angle a);
int	AngleIsBad (Angle a);
int	AnglePrintVal(FILE *out, Angle a, char *ok_fmt, char *bad_fmt);
double	ISin (Angle a);
double	ICos (Angle a);
GeoPt	GeoPtFmDeg (double dLat, double dLon);
GeoPt	GeoPtFmRad (double dLat, double dLon);
void	GeoPtGetDeg (GeoPt geoPt, double *dLatP, double *dLonP);
void	GeoPtGetRad (GeoPt geoPt, double *dLatP, double *dLonP);
CartPt	LatLonToCart (GeoPt geoPt);
int	GeoPtIsSomewhere (GeoPt geoPt);
int	GeoPtIsNowhere (GeoPt geoPt);
GeoPt	GeoPtNowhere (void);
int	MapPtIsSomewhere (MapPt mapPt);
int	MapPtIsNowhere (MapPt mapPt);
MapPt	MapPtNowhere (void);
MapPt	ScaleMapPt (MapPt mapPt, double scale);
GeoPt	GeoStep (GeoPt geoPt, Angle dir, Angle dist);
Angle	GeoDistance (GeoPt p1, GeoPt p2);
double	GeoQuickDistance (GeoPt p1, GeoPt p2);
Angle	Azimuth (GeoPt p1, GeoPt p2);
GeoPt	GCircleX (GeoPt ln1pt1, GeoPt ln1pt2, GeoPt ln2pt1, GeoPt ln2pt2);

/*
 * These constants specify relative Cardinal directions.
 */

enum LatSgn {North, Eq, South};
enum LonSgn {West, PrMd, East};

/*
 * Global functions for comparing angles and putting them into domains.
 */

Angle	DomainLat (Angle lat);
Angle	DomainLon (Angle lon, Angle refLon);
Angle	GwchLon (Angle lon);
GeoPt	DomainLonPt (GeoPt geoPt, Angle refLon);
GeoPt	GwchLonPt (GeoPt geoPt);
enum	LonSgn LonCmp (Angle lon0, Angle lon1);
enum	LatSgn LatCmp (Angle lat0, Angle lat1);
int	AngleCmp (Angle d0, Angle d1);
int	LonBtwn (Angle lon, Angle lon0, Angle lon1);
int	LonBtwn1 (Angle lon, Angle lon0, Angle lon1);

/*
 * This structure stores a time instant as a Julian day and seconds
 * since midnight.
 */

struct GeoTime_Jul {
    int   day;		/* Julian day */
    double second;	/* Seconds since midnight */
};

/*
 * This structure stores a time instant as a calendar date and clock time.
 */

struct GeoTime_Cal {
    int year;
    int month;
    int day;
    int hour;
    int minute;
    double second;
};

/*
 * The following functions convert between time representations.
 */

struct GeoTime_Cal GeoTime_CalSet(int year, int month, int day, int hour,
	int minute, double second);
struct GeoTime_Jul GeoTime_JulSet(int day, double second);
struct GeoTime_Jul GeoTime_CalToJul(struct GeoTime_Cal cal);
struct GeoTime_Cal GeoTime_JulToCal(struct GeoTime_Jul jul);
void GeoTime_Incr(struct GeoTime_Jul *jul, double ds);
int GeoTime_Cmp(struct GeoTime_Jul jul1, struct GeoTime_Jul jul2);
double GeoTime_Diff(struct GeoTime_Jul jul1, struct GeoTime_Jul jul2);
    
#ifdef __cplusplus
}
#endif

#endif
