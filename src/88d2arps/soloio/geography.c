/*
 * geography.c --
 *
 * 	This file provides general purpose geography functions.
 * 	See the user documentation for more information.
 * 
 * Copyright (c) 2007 Gordon D. Carrie.  All rights reserved.
 * 
 * Licensed under the Open Software License version 2.1
 * 
 * Please address questions and feedback to user0@tkgeomap.org
 *
 * @(#) $Id: geography.c,v 4.0 2013/02/08 18:31:45 ywang Exp $
 *
 ********************************************
 *
 */

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "geography.h"

/*
 * The following constants enable conversion between floating point degree
 * and radian measurements and the Angle data type.
 *
 * Angles are stored as integer microdegrees.  This helps reduce the effects
 * of round off.  A microdegree is about 0.11 meters on the spherical Earth,
 * so points closer than this are assumed to be at the same location.
 */

#define UNITPERDEG 1.0e+6
#define DEGPERUNIT 1.0e-6
#define RADPERUNIT 0.017453292519943294892e-6
#define UNITPERRAD 57.29577951308232088e+6
#define MAXDEG (INT_MAX * DEGPERUNIT)
#define MINDEG (INT_MIN * DEGPERUNIT)
#define MAXRAD (INT_MAX * RADPERUNIT)
#define MINRAD (INT_MIN * RADPERUNIT)

/*
 * The following constants give Angle representations of some frequently used
 * angles.
 */

#define D360 360000000
#define D270 270000000
#define D180 180000000
#define D90   90000000

static double rearth = 6366707.01949;

/*
 *----------------------------------------------------------------------
 *
 * REarth --
 *
 *	This procedure retrieves the radius of the Earth.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

double REarth(void)
{
    return rearth;
}

/*
 *----------------------------------------------------------------------
 *
 * SetREarth --
 *
 *	This procedure sets the assumed radius of the Earth.
 *
 * Results:
 * 	None.
 *
 * Side effects:
 *	A static variable is set.
 *
 *----------------------------------------------------------------------
 */

void SetREarth(double r)
    /*
     * r - New Earth radius
     */
{
    rearth = r;
}

/*
 *----------------------------------------------------------------------
 *
 * BadAngle --
 *
 *	This procedure returns a bogus angle.  It is used to indicate
 *	error conditions.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

Angle BadAngle(void)
{
    return -INT_MAX;
}

/*
 *----------------------------------------------------------------------
 *
 * AngleIsOK --
 *
 *	This procedure determines whether an angle value is valid.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

int AngleIsOK(Angle a)
    /* 
     * a - Angle to evaluate
     */
{
    return a != BadAngle();
}

/*
 *----------------------------------------------------------------------
 *
 * AngleIsBad --
 *
 *	This procedure determines whether an angle value is invalid.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

int AngleIsBad(Angle a)
    /*
     * a - Angle to evaluate
     */
{
    return a == BadAngle();
}

/*
 *------------------------------------------------------------------------
 *
 * AnglePrintVal --
 *
 * 	This convenience function prints an angle value, if possible.
 *
 * Arguments:
 *	FILE *out	- output stream
 *	Angle a		- angle value to print
 *	char *dat_fmt	- format to use if fprintf if a is an angle.
 *	char *nodat_fmt	- format to use if fprintf if a is not an angle.
 *
 * Results:
 * 	Return value is the return value of fprintf.  If a is angle, it's
 * 	value in degrees is printed to the output stream.
 *
 * Side effects:
 * 	The output stream is modified.
 * 	
 *------------------------------------------------------------------------
 */

int AnglePrintVal(FILE *out, Angle a, char *ok_fmt, char *bad_fmt)
{
    if (AngleIsOK(a)) {
	return fprintf(out, ok_fmt, AngleToDeg(a));
    } else {
	return fprintf(out, bad_fmt);
    }
}

/*
 *----------------------------------------------------------------------
 *
 * AngleFmDeg --
 *
 *	This procedure converts a degree value to an Angle value.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

Angle AngleFmDeg(double deg)
    /*
     * deg - Degree value to convert
     */
{
    return (deg > MAXDEG || deg < MINDEG)
	? BadAngle() : UNITPERDEG * deg + (deg > 0.0 ? 0.5 : -0.5);
}

/*
 *----------------------------------------------------------------------
 *
 * AngleToDeg --
 *
 * 	This procedure converts an Angle value to degrees.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

double AngleToDeg(Angle a)
    /*
     * a - Angle value to convert
     */
{
    return a * DEGPERUNIT;
}

/*
 *----------------------------------------------------------------------
 *
 * AngleFmRad --
 *
 *	This procedure converts a radian value to an Angle value.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

Angle AngleFmRad(double rad)
    /*
     * rad - Radian value to convert
     */
{
    return (rad > MAXRAD || rad < MINRAD)
	? BadAngle() : UNITPERRAD * rad + (rad > 0.0 ? 0.5 : -0.5);
}

/*
 *----------------------------------------------------------------------
 *
 * AngleToRad --
 *
 * 	This procedure converts an Angle value to radians.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

double AngleToRad(Angle a)
    /*
     * a - Angle value to convert
     */
{
    return a * RADPERUNIT;
}

/*
 *----------------------------------------------------------------------
 *
 * ISin --
 *
 *	This procedure computes the sine of an Angle value.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

double ISin(Angle a)
    /*
     * a - Angle value
     */
{
    return sin(a * RADPERUNIT);
}

/*
 *----------------------------------------------------------------------
 *
 * ICos --
 *
 *	This procedure computes the cosine of an Angle value.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

double ICos(Angle a)
    /*
     * a - Angle value
     */
{
    return cos(a * RADPERUNIT);
}

/*
 *----------------------------------------------------------------------
 *
 * GeoPtFmDeg --
 *
 *	This procedure sets a geopoint given latitude and longitude values
 *	in degrees.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

GeoPt GeoPtFmDeg(double dLat, double dLon)
    /*
     * dLat - Latitude, degrees
     * dLon - Longitude, degrees
     */
{
    GeoPt geoPt;			/* Return value */
    geoPt.lat = AngleFmDeg(dLat);
    geoPt.lon = AngleFmDeg(dLon);
    return (AngleIsBad(geoPt.lat) || AngleIsBad(geoPt.lon))
	? GeoPtNowhere() : geoPt;
}

/*
 *----------------------------------------------------------------------
 *
 * GeoPtFmRad --
 *
 *	This procedure sets a geopoint given latitude and longitude values
 *	in radians.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

GeoPt GeoPtFmRad(double dLat, double dLon)
    /*
     * dLat - Latitude, radians
     * dLon - Longitude, radians
     */
{
    GeoPt geoPt;
    geoPt.lat = AngleFmRad(dLat);
    geoPt.lon = AngleFmRad(dLon);
    return (AngleIsBad(geoPt.lat) || AngleIsBad(geoPt.lon))
	? GeoPtNowhere() : geoPt;
}

/*
 *----------------------------------------------------------------------
 *
 * GeoPtGetDeg --
 *
 * 	This procedure retrieves latitude and longitude values from a
 * 	geopoint in degrees.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

void GeoPtGetDeg(GeoPt geoPt, double *dLatP, double *dLonP)
    /*
     * geoPt	- Geographic point
     * dLatP	- Recipient of latitude
     * dLonP	- Recipient of longitude
     */
{
    *dLatP = AngleToDeg(geoPt.lat);
    *dLonP = AngleToDeg(geoPt.lon);
}

/*
 *----------------------------------------------------------------------
 *
 * GeoPtGetRad --
 *
 * 	This procedure retrieves latitude and longitude values from a
 * 	geopoint in degrees.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

void GeoPtGetRad(GeoPt geoPt, double *dLatP, double *dLonP)
    /*
     * geoPt	- Geographic point
     * dLatP	- Recipient of latitude
     * dLonP	- Recipient of longitude
     */
{
    *dLatP = AngleToRad(geoPt.lat);
    *dLonP = AngleToRad(geoPt.lon);
}

/*
 *----------------------------------------------------------------------
 *
 * GeoPtIsSomewhere --
 *
 * 	This procedure determines whether a geopoint is valid.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

int GeoPtIsSomewhere(GeoPt geoPt)
    /*
     * geoPt	- Geographic point to evaluate
     */
{
    return AngleIsOK(geoPt.lat) && AngleIsOK(geoPt.lon);
}
/*
 *----------------------------------------------------------------------
 *
 * GeoPtIsNowhere --
 *
 * 	This procedure determines whether a geopoint is invalid.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

int GeoPtIsNowhere(GeoPt geoPt)
    /*
     * geoPt - Geographic point to evaluate
     */
{
    return AngleIsBad(geoPt.lat) || AngleIsBad(geoPt.lon);
}

/*
 *----------------------------------------------------------------------
 *
 * GeoPtNowhere --
 *
 * 	This procedure returns an invalid geopoint.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

GeoPt GeoPtNowhere(void)
{
    GeoPt nowhere;		/* Return value */
    nowhere.lat = nowhere.lon = BadAngle();
    return nowhere;
}
/*
 *----------------------------------------------------------------------
 *
 * MapPtIsSomewhere --
 *
 * 	This procedure determines whether a map-point is valid.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */


int MapPtIsSomewhere(MapPt mapPt)
    /*
     * mapPt	- Map point to evaluate
     */
{
    return mapPt.abs != FLT_MAX && mapPt.ord != FLT_MAX;
}

/*
 *----------------------------------------------------------------------
 *
 * MapPtIsNowhere --
 *
 * 	This procedure determines whether a map-point is invalid.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

int MapPtIsNowhere(MapPt mapPt)
    /*
     * mapPt	- Map point to evaluate
     */
{
    return mapPt.abs == FLT_MAX || mapPt.ord == FLT_MAX;
}

/*
 *----------------------------------------------------------------------
 *
 * MapPtNowhere --
 *
 * 	This procedure returns an invalid map-point.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

MapPt MapPtNowhere(void)
{
    MapPt nowhere = {FLT_MAX, FLT_MAX};
    return nowhere;
}

/*
 *----------------------------------------------------------------------
 *
 * LatLonToCart --
 *
 *	This procedure convert a geopoint to 3D Geocentric Cartesian
 *	coordinates where Earth has unit radius
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

CartPt LatLonToCart(GeoPt geoPt)
    /*
     * geoPt - Geographic point value to convert
     */
{
    double lat, lon;			/* Latitude, longitude in radians */
    double coslat;
    CartPt cpt;				/* Return value */

    GeoPtGetRad(geoPt, &lat, &lon);
    coslat = cos(lat);
    cpt.x = coslat * cos(lon);
    cpt.y = coslat * sin(lon);
    cpt.z = sin(lat);
    return cpt;
}

/*
 *----------------------------------------------------------------------
 *
 * ScaleMapPt --
 *
 *	This procedure scales a map-point.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

MapPt ScaleMapPt(MapPt mapPt, double scale)
    /*
     * mapPt - Map point to scale
     * scale - Amount to scale it by
     */
{
    if (MapPtIsNowhere(mapPt) || scale <= 0.0) {
	return MapPtNowhere();
    }
    mapPt.abs *= scale;
    mapPt.ord *= scale;
    return mapPt;
}

/*
 *----------------------------------------------------------------------
 *
 * GeoDistance --
 *
 *	This procedure computes the distance between two geopoints.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

Angle GeoDistance(GeoPt p1, GeoPt p2)
    /*
     * p1 - First point
     * p2 - Second point
     */
{
    double lat1, lon1, lat2, lon2;
    double sin_dlon_2, sin_dlat_2;
    double a;

    /*
     * Reference:
     *
     *	R.W. Sinnott, "Virtues of the Haversine", Sky and Telescope,
     *  vol. 68, no. 2, 1984, p. 159
     *
     * cited in:
     *	http://www.census.gov/cgi-bin/geo/gisfaq?Q5.1
     */

    GeoPtGetRad(p1, &lat1, &lon1);
    GeoPtGetRad(p2, &lat2, &lon2);
    sin_dlon_2 = sin(0.5 * (lon2 - lon1));
    sin_dlat_2 = sin(0.5 * (lat2 - lat1));
    a = sqrt(sin_dlat_2 * sin_dlat_2
	+ cos(lat1) * cos(lat2) * sin_dlon_2 * sin_dlon_2);
    return AngleFmRad(a > 1.0 ? M_PI : 2.0 * asin(a));
}
/*
 *----------------------------------------------------------------------
 *
 * GeoQuickDistance --
 *
 *	This procedure computes a measure of the distance between two geopoints.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

double GeoQuickDistance(GeoPt p1, GeoPt p2)
    /*
     * p1 - First point
     * p2 - Second point
     */
{
    double lat1, lon1, lat2, lon2;
    double sin_dlon_2, sin_dlat_2;

    /*
     * Reference:
     *
     *	R.W. Sinnott, "Virtues of the Haversine", Sky and Telescope,
     *  vol. 68, no. 2, 1984, p. 159
     *
     * cited in:
     *	http://www.census.gov/cgi-bin/geo/gisfaq?Q5.1
     */

    GeoPtGetRad(p1, &lat1, &lon1);
    GeoPtGetRad(p2, &lat2, &lon2);
    sin_dlon_2 = sin(0.5 * (lon2 - lon1));
    sin_dlat_2 = sin(0.5 * (lat2 - lat1));
    return sin_dlat_2 * sin_dlat_2
	+ cos(lat1) * cos(lat2) * sin_dlon_2 * sin_dlon_2;
}

/*
 *----------------------------------------------------------------------
 *
 * Azimuth --
 *
 *	This procedure computes the azimuth (bearing) from one point to another.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

Angle Azimuth(GeoPt p1, GeoPt p2)
    /*
     * p1 - First point
     * p2 - Second point
     */
{
    double lat1, lon1, lat2, lon2;
    double cosDLon, sinDLon, sinDLat, sinMLat;

    GeoPtGetRad(p1, &lat1, &lon1);
    GeoPtGetRad(p2, &lat2, &lon2);
    cosDLon = cos(lon2 - lon1);
    sinDLon = sin(lon2 - lon1);
    sinDLat = sin(lat1 - lat2);
    sinMLat = sin(lat2 + lat1);
    return AngleFmRad(atan2(cos(lat2) * sinDLon,
		0.5 * (sinMLat - sinDLat - (sinMLat + sinDLat) * cosDLon)));
}

/*
 *----------------------------------------------------------------------
 *
 * GeoStep --
 *
 *	This procedure finds the point a given distance and direction from
 *	a another point.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

GeoPt GeoStep(GeoPt geoPt, Angle iDir, Angle iDist)
    /*
     * geoPt	- Starting point
     * iDir	- Which way to step 
     * iDist	- How far to step
     */
{
    double dist, dir;
    double lat, lon;
    double cos_dist, sin_dist, cos_dir, sin_dir;
    double cos_lat, cos_lon, sin_lon, sin_lat;
    double x, y, z, h_1, h_2, dh;

    dist = AngleToRad(iDist);
    cos_dist = cos(dist);
    sin_dist = sin(dist);

    dir = AngleToRad(iDir);
    cos_dir = cos(dir);
    sin_dir = sin(dir);

    GeoPtGetRad(geoPt, &lat, &lon);
    cos_lat = cos(lat);
    cos_lon = cos(lon);
    sin_lon = sin(lon);
    sin_lat = sin(lat);
    x = cos_dist * cos_lon * cos_lat - sin_dir * sin_dist * sin_lon 
	- cos_lon * cos_dir * sin_dist * sin_lat;
    y = sin_dir * cos_lon * sin_dist + cos_dist * cos_lat * sin_lon 
	- cos_dir * sin_dist * sin_lon * sin_lat;
    z = cos_lat * cos_dir * sin_dist + cos_dist * sin_lat;

    /*
     * Compute x^2  +  y^2
     */

    h_1 = cos_dist * cos_lat - cos_dir * sin_dist * sin_lat;
    h_2 = (sin_dir) * (sin_dist);
    dh = h_1 * h_1 + h_2 * h_2;

    lat = (dh == 0.0) ? (z > 0 ? M_PI_2 : -M_PI_2) : atan(z / sqrt(dh));
    lon = atan2(y, x);

    return GeoPtFmRad(lat, lon);
}

/*
 *----------------------------------------------------------------------
 *
 * GCircleX --
 *
 *	This procedure computes the intersection of two great circles.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

GeoPt GCircleX(GeoPt ln1pt1, GeoPt ln1pt2, GeoPt ln2pt1, GeoPt ln2pt2)
    /*
     * ln1pt1 - First point on first great circle
     * ln1pt2 - Second point on first great circle
     * ln2pt1 - First point on second great circle
     * ln2pt2 - Second point on second great circle
     */
{
    CartPt ln1pt1c, ln1pt2c, ln2pt1c, ln2pt2c; 
    				/* Points on the great circles in
				 * Cartesian Geocentric coordinates */
    CartPt p1;			/* Normal to plane containing ln1pt1 ln1pt2 */
    CartPt p2;			/* Normal to plane containing ln2pt2 ln2pt2 */
    CartPt c;			/* Vector perpendicular to p1 and p2 */
    double cs;			/* Normalizing factor */
    CartPt m;			/* Mean of the four input points */
    double dp, dm;		/* Distance from mean to +-(cx, cy, cz) */
    double lat, lon;		/* Lat and lon of return value */

    /*
     * Points in 3D Geocentric Cartesian Coordinates
     */

    ln1pt1c = LatLonToCart(ln1pt1);
    ln1pt2c = LatLonToCart(ln1pt2);
    ln2pt1c = LatLonToCart(ln2pt1);
    ln2pt2c = LatLonToCart(ln2pt2);

    /*
     * Normal to plane containing ln1pt1 and ln1pt2 (given by cross product)
     */

    p1.x = ln1pt1c.y * ln1pt2c.z - ln1pt1c.z * ln1pt2c.y;
    p1.y = ln1pt1c.z * ln1pt2c.x - ln1pt1c.x * ln1pt2c.z;
    p1.z = ln1pt1c.x * ln1pt2c.y - ln1pt1c.y * ln1pt2c.x;

    /*
     * Normal to plane containing ln2pt1 and ln2pt2 (given by cross product)
     */

    p2.x = ln2pt1c.y * ln2pt2c.z - ln2pt1c.z * ln2pt2c.y;
    p2.y = ln2pt1c.z * ln2pt2c.x - ln2pt1c.x * ln2pt2c.z;
    p2.z = ln2pt1c.x * ln2pt2c.y - ln2pt1c.y * ln2pt2c.x;

    /*
     * Perpendicular to p1 and p2 (another cross product)
     */

    c.x = p1.y * p2.z - p1.z * p2.y;
    c.y = p1.z * p2.x - p1.x * p2.z;
    c.z = p1.x * p2.y - p1.y * p2.x;

    /*
     * If all pts on one line, undefined intersection
     */

    if (c.x == 0.0 && c.y == 0.0 && c.z == 0.0) {
	return GeoPtNowhere();
    }

    /*
     * Normalize
     */

    cs = 1.0 / sqrt(c.x * c.x + c.y * c.y + c.z * c.z);
    c.x *= cs;
    c.y *= cs;
    c.z *= cs;

    /*
     * Intersection is at a pole
     */

    if (c.x == 0.0 && c.y == 0.0) {
	return GeoPtFmDeg((c.z > 0) ? 90.0 : -90.0, 0.0);
    }

    /*
     * Find which of the two intersections of the two great circles is
     * closer to the mean of the four input points		
     */

    m.x = 0.25 * (ln1pt1c.x + ln1pt2c.x + ln2pt2c.x + ln2pt2c.x);
    m.y = 0.25 * (ln1pt1c.y + ln1pt2c.y + ln2pt2c.y + ln2pt2c.y);
    m.z = 0.25 * (ln1pt1c.z + ln1pt2c.z + ln2pt2c.z + ln2pt2c.z);
    dp = (m.x-c.x) * (m.x-c.x) + (m.y-c.y) * (m.y-c.y) + (m.z-c.z) * (m.z-c.z);
    dm = (m.x+c.x) * (m.x+c.x) + (m.y+c.y) * (m.y+c.y) + (m.z+c.z) * (m.z+c.z);
    if (dm < dp) {
	/*
	 * Mean point is closer to (-c.x, -c.y, -c.z)
	 */

	c.x = -c.x;
	c.y = -c.y;
	c.z = -c.z;
    }

    lat = atan(c.z / sqrt(c.x * c.x + c.y * c.y));
    lon = atan2(c.y, c.x);
    return GeoPtFmRad(lat, lon);
}

/*
 *----------------------------------------------------------------------
 *
 * DomainLat --
 *
 *	This procedure puts an angle into the range -90.0 <= angle <= 90.0
 *	without changing its sine.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

Angle DomainLat(Angle lat)
    /*
     * lat - Latitude value
     */
{
    if (lat > D360) {
	lat = lat - D360 * (lat / D360);
    } else if (lat < 0) {
	lat = lat + D360 * (-lat / D360 + 1);
    }
    if (lat > D90 && lat < D270) {
        lat = D180 - lat;
    } else if (lat >= D270) {
        lat = lat - D360;
    }
    return lat;
}

/*
 *----------------------------------------------------------------------
 *
 * DomainLon --
 *
 *	This procedure puts an angle within 180 degrees of a reference angle.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

Angle DomainLon(Angle lon, Angle refLon)
    /*
     * lon	- Longitude value
     * refLon	- Reference longitude
     */
{
    if (lon == refLon) {
	return lon;
    }
    if (lon > refLon + D360) {
	lon -= D360 * ((lon - refLon) / D360);
    } else if (lon < refLon - D360) {
	lon += D360 * ((refLon - lon) / D360);
    }
    if (lon > refLon + D180) {
	return lon - D360;
    } else if (lon < refLon - D180) {
	return lon + D360;
    } else {
	return lon;
    }
}

/*
 *----------------------------------------------------------------------
 *
 * GwchLon --
 *
 *	This procedure is equivalent to DomainLon(lon, 0.0)
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

Angle GwchLon(Angle lon)
    /*
     * lon - Longitude value
     */
{
    if (lon == 0) {
	return lon;
    }
    if (lon > D360) {
	lon -= D360 * ((lon ) / D360);
    } else if (lon < -D360) {
	lon += D360 * ((-lon) / D360);
    }
    if (lon > D180) {
	return lon - D360;
    } else if (lon < -D180) {
	return lon + D360;
    } else {
	return lon;
    }
}

/*
 *----------------------------------------------------------------------
 *
 * DomainLonPt --
 *
 *	This procedure returns a geographic point whose longitude is within
 *	180 degrees of a reference longitude.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

GeoPt DomainLonPt(GeoPt geoPt, Angle refLon)
    /*
     * geoPt	- Geographic point
     * refLon	- Reference longitude
     */
{
    geoPt.lat = DomainLat(geoPt.lat);
    geoPt.lon = DomainLon(geoPt.lon, refLon);
    return geoPt;
}

/*
 *----------------------------------------------------------------------
 *
 * GwchLonPt --
 *
 *	This procedure is equivalent to DomainLonPt(geoPt, 0.0).
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

GeoPt GwchLonPt(GeoPt geoPt)
    /*
     * geoPt - Geographic point
     */
{
    return DomainLonPt(geoPt, 0);
}

/*
 *----------------------------------------------------------------------
 *
 * LonCmp --
 *
 *	This procedure compares two longitudes.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

enum LonSgn LonCmp(Angle lon0, Angle lon1)
    /*
     * Angle lon0 - First longitude
     * Angle lon1 - Second longitude
     */
{
    lon0 = DomainLon(lon0, lon1);
    return (lon0 < lon1) ? West : (lon0 > lon1) ? East : PrMd;
}

/*
 *----------------------------------------------------------------------
 *
 * LonBtwn --
 *
 *	This procedure determines whether a meridian lies between two others.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

int LonBtwn(Angle lon, Angle lon0, Angle lon1)
    /*
     * lon	- Longitude to evaluate
     * lon0	- Longitude at one end of interval
     * lon1	- Longitude at other end of interval
     */
{
    lon0 = DomainLon(lon0, lon);
    lon1 = DomainLon(lon1, lon);
    return ((lon0 > lon1) ? lon0 - lon1 : lon1 - lon0) < D180
	&& ((lon0 < lon && lon < lon1) || (lon1 < lon && lon < lon0));
}

/*
 *----------------------------------------------------------------------
 *
 * LonBtwn1 --
 *
 *	This procedure determines whether a meridian lies between or coincides
 *	with two others,
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

int LonBtwn1(Angle lon, Angle lon0, Angle lon1)
    /*
     * lon	- Longitude to evaluate
     * lon0	- Longitude at one end of interval
     * lon1	- Longitude at other end of interval
     */
{
    lon0 = DomainLon(lon0, lon);
    lon1 = DomainLon(lon1, lon);
    return    ((lon0 > lon1) ? lon0 - lon1 : lon1 - lon0) < D180
	&& ((lon0 < lon && lon <= lon1) || (lon1 < lon && lon <= lon0));
}

/*
 *----------------------------------------------------------------------
 *
 * LatCmp --
 *
 *	This procedure compares two latitudes.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

enum LatSgn LatCmp(Angle lat0, Angle lat1)
    /*
     * lat0	- First latitude
     * lat1	- Second latitude
     */
{
    return (lat0 < lat1) ? South : (lat0 > lat1) ? North : Eq;
}

/*
 *----------------------------------------------------------------------
 *
 * AngleCmp --
 *
 *	This procedure compares to angles.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

int AngleCmp(Angle a0, Angle a1)
    /*
     * a0 - First angle to compare 
     * a1 - Second angle to compare
     */
{
    return (a0 < a1) ? -1 : (a0 > a1) ? 1 : 0;
}
