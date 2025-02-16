/*
 * geoLines.c --
 *
 * 	This file defines functions that manage geolines and geolinearrays.
 * 	See the user documentation for more information.
 * 
 * Copyright (c) 2007 Gordon D. Carrie.  All rights reserved.
 * 
 * Licensed under the Open Software License version 2.1
 * 
 * Please address questions and feedback to user0@tkgeomap.org
 *
 * @(#) $Id: geoLines.c,v 4.0 2013/02/08 18:31:45 ywang Exp $
 *
 ********************************************
 *
 */

#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "alloc.h"
#include "errMsg.h"
#include "geoLines.h"

/*
 * Local function declaration.
 */

static void lnInit (struct GeoLn *lnP);


/*
 *----------------------------------------------------------------------
 *
 * GeoLnCreate --
 *
 *	This procedure creates and initializes a new geoline.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	See the user documentation.
 *
 *----------------------------------------------------------------------
 */

struct GeoLn *GeoLnCreate(unsigned nPtsMax)
{
    struct GeoLn *lnP = NULL;

    lnP = (struct GeoLn *)MALLOC(sizeof(*lnP));
    lnInit(lnP);
    if (nPtsMax == 0) {
	return lnP;
    }
    lnP->pts = (GeoPt *)MALLOC(nPtsMax * sizeof(GeoPt));
    lnP->nPtsMax = nPtsMax;
    return lnP;
}

/*
 *----------------------------------------------------------------------
 *
 * lnInit --
 *
 *	This procedure initializes a new geoline.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	The contents of a GeoLn structure are set.  The previous contents,
 *	which are assumed to be garbage, are ignored.
 *
 *----------------------------------------------------------------------
 */

static void
lnInit(struct GeoLn *lnP)
{
    /* Initialize a new Line struct. */

    lnP->nPts = lnP->nPtsMax = 0;
    lnP->lonMax = lnP->latMax = -INT_MAX;
    lnP->lonMin = lnP->latMin = INT_MAX;
    lnP->pts = NULL;
}

/*
 *----------------------------------------------------------------------
 *
 * GeoLnClear --
 *
 *	This procedure set the number of points in a geoline to zero without
 *	releasing its internal storage.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	See the user documentation.
 *
 *----------------------------------------------------------------------
 */

void
GeoLnClear(struct GeoLn *lnP)
{
    lnP->nPts = 0;
    lnP->lonMax = lnP->latMax = -INT_MAX;
    lnP->lonMin = lnP->latMin = INT_MAX;
}

/*
 *----------------------------------------------------------------------
 *
 * GeoLnDestroy --
 *
 *	This procedure frees internal storage in a GeoLn structure, and
 *	frees the structure.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	See the user documentation.
 *
 *----------------------------------------------------------------------
 */

void
GeoLnDestroy(struct GeoLn *lnP)
{
    if (!lnP) {
	return;
    }
    FREE((char *)lnP->pts);
    FREE((char *)lnP);
}

/*
 *----------------------------------------------------------------------
 *
 * GeoLnSetAlloc --
 *
 *	This procedure changes the number of points a geoline can store.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	This procedure reallocates a GeoLn structure's internal storage.
 *	If the number of points is reduced, the geoline's maxima and minima
 *	are recomputed.
 *
 *----------------------------------------------------------------------
 */

void GeoLnSetAlloc(struct GeoLn *lnP, unsigned nPtsMax)
{
    if (lnP->nPtsMax == nPtsMax) {
	return;
    }
    if (nPtsMax == 0) {
	FREE((char *)lnP->pts);
	lnInit(lnP);
	return;
    }
    lnP->pts = (GeoPt *)REALLOC(lnP->pts, nPtsMax * sizeof(GeoPt));
    lnP->nPtsMax = nPtsMax;
    if (lnP->nPts > lnP->nPtsMax) {
	/*
	 * If shrinking line, reset nPts and recompute maxima and minima.
	 */

	GeoPt *p, *pe;
	lnP->nPts = lnP->nPtsMax;
	for (p = lnP->pts, pe = lnP->pts + lnP->nPts; p < pe; p++) {
	    if (GeoPtIsSomewhere(*p)) {
		lnP->latMax
		    = (lnP->latMax > p->lat) ? lnP->latMax : p->lat;
		lnP->lonMax
		    = (lnP->lonMax > p->lon) ? lnP->lonMax : p->lon;
		lnP->latMin
		    = (lnP->latMin < p->lat) ? lnP->latMin : p->lat;
		lnP->lonMin
		    = (lnP->lonMin < p->lon) ? lnP->lonMin : p->lon;
	    }
	}
    }
}

/*
 *----------------------------------------------------------------------
 *
 * GeoLnAddPt --
 *
 *	This procedure adds a point to a geoline.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	The geoline's allocation might increase.  It's maxima and minima will
 *	be adjusted.
 *
 *----------------------------------------------------------------------
 */

void GeoLnAddPt(GeoPt p, struct GeoLn *lnP)
{
    if (lnP->nPts + 1 > lnP->nPtsMax) {
	GeoLnSetAlloc(lnP, ((lnP->nPtsMax + 4) * 5) / 4);
    }
    if (GeoPtIsSomewhere(p)) {
	lnP->latMax = (lnP->latMax > p.lat) ? lnP->latMax : p.lat;
	lnP->lonMax = (lnP->lonMax > p.lon) ? lnP->lonMax : p.lon;
	lnP->latMin = (lnP->latMin < p.lat) ? lnP->latMin : p.lat;
	lnP->lonMin = (lnP->lonMin < p.lon) ? lnP->lonMin : p.lon;
    }
    lnP->pts[lnP->nPts] = p;
    lnP->nPts++;
    return;
}

/*
 *----------------------------------------------------------------------
 *
 * GeoLnGetPt --
 *
 *	This procedure retrieve a geopoint from a geoline.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

GeoPt GeoLnGetPt(struct GeoLn *lnP, unsigned n)
{
    /* Return MultiPt n from lnP, or a bogus location if out of bounds. */
    return (n < lnP->nPts) ? *(lnP->pts + n) : GeoPtNowhere();
}

/*
 *----------------------------------------------------------------------
 *
 * GeoLnCtr --
 *
 *	This procedure computes the average position of the points in a
 *	geoline in 3-dimensional space.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

CartPt GeoLnCtr(struct GeoLn *lnP)
{
    CartPt ctr = {0.0, 0.0, 0.0};	/* Return value */
    CartPt cpt;
    GeoPt *p, *end;			/* Loop index */

    for (p = lnP->pts, end = p + lnP->nPts; p < end; p++) {
	cpt = LatLonToCart(*p);
	ctr.x += cpt.x;
	ctr.y += cpt.y;
	ctr.z += cpt.z;
    }
    ctr.x /= lnP->nPts;
    ctr.y /= lnP->nPts;
    ctr.z /= lnP->nPts;
    return ctr;
}

/*
 *----------------------------------------------------------------------
 *
 * GeoLnContainGeoPt --
 *
 *	This procedure determines if a geoline enclosed a geopoint.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

int GeoLnContainGeoPt(GeoPt geoPt, struct GeoLn *lnP)
{
    int mrdx;	/* Number of times line crosses meridian containing geoPt */
    int lnx;	/* Number of times line crosses line from geoPt to North pole */
    GeoPt *p0, *p1, *end;

    /*
     * Loop through segments in lnP, counting number of times lnP
     * crosses meridian, and number of crossings between geoPt and North pole
     */

    for (mrdx = 0, lnx = 0,
	    p0 = lnP->pts + lnP->nPts - 1,
	    p1 = lnP->pts,
	    end = lnP->pts + lnP->nPts;
	    p1 < end;
	    p0 = p1++) {
	/*
	 * Determine if segment defined by p0--p1 straddles meridian
	 * containing geoPt, or is on boundary.  Do not count segments on
	 * boundary as more than one crossing
	 */

	if (LonBtwn1(geoPt.lon, p1->lon, p0->lon)) {
	    double lat0 = AngleToDeg(p0->lat);
	    double lon0 = AngleToDeg(p0->lon);
	    double lat1 = AngleToDeg(p1->lat);
	    double lon1 = AngleToDeg(p1->lon);
	    double xlat;		/* Latitude of segment crossing */

	    mrdx++;
	    xlat = lat0 + (AngleToDeg(geoPt.lon) - lon0)
		* (lat1 - lat0) / (lon1 - lon0);
	    if (LatCmp(AngleFmDeg(xlat), geoPt.lat) == North) {
		lnx = !lnx;
	    }

	}
    }

    if (mrdx % 2 == 1) {
	/*
	 * Odd number of meridian crossings => region contains a pole */
	/* Assume pole is that of hemisphere containing lnP's mean
	 */

	CartPt ctr = GeoLnCtr(lnP);
	if (ctr.z > 0.0) {
	    /*
	     * Region contains North Pole
	     */

	    return !lnx;
	}

    }
    return lnx;
}

/*
 *----------------------------------------------------------------------
 *
 * GeoLnArrCreate --
 *
 *	This procedure creates and initializes an new geolinearray.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	See the user documentation.
 *
 *----------------------------------------------------------------------
 */

struct GeoLnArr *GeoLnArrCreate(unsigned nLinesMax)
{
    struct GeoLnArr *lnArrP;
    unsigned n;
    size_t sz;

    lnArrP = (struct GeoLnArr *)MALLOC(sizeof(*lnArrP));
    lnArrP->descr = NULL;
    lnArrP->lines = NULL;
    GeoLnArrSetDescr(lnArrP, "");
    lnArrP->nLines = lnArrP->nLinesMax = 0;
    lnArrP->nPts = lnArrP->nMax = 0;
    lnArrP->lonMax = lnArrP->latMax = -INT_MAX;
    lnArrP->lonMin = lnArrP->latMin = INT_MAX;
    lnArrP->lines = NULL;
    if (nLinesMax == 0) {
	return lnArrP;
    }
    sz = nLinesMax * sizeof(struct GeoLn *);
    lnArrP->lines = (struct GeoLn **)MALLOC(sz);
    lnArrP->nLinesMax = nLinesMax;
    for (n = 0; n < nLinesMax; n++) {
	lnArrP->lines[n] = NULL;
    }
    return lnArrP;
}

/*
 *----------------------------------------------------------------------
 *
 * GeoLnArrSetDescr --
 *
 *	This procedure sets the descriptor for a geolinearray.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	See the user documentation.
 *
 *----------------------------------------------------------------------
 */

void GeoLnArrSetDescr(struct GeoLnArr *lnArrP, char *descr)
{
    char *s, *s1, *d;

    lnArrP->descr = (char *)REALLOC(lnArrP->descr, strlen(descr) + 1);
    for (s = descr, s1 = s + strlen(descr), d = lnArrP->descr; s < s1; s++, d++) {
	*d = *s;
    }
    *d = '\0';
}

/*
 *----------------------------------------------------------------------
 *
 * GeoLnArrSetAlloc --
 *
 *	This procedure changes the number of geolines a geolinearray can store.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	See the user documentation.
 *
 *----------------------------------------------------------------------
 */

void GeoLnArrSetAlloc(struct GeoLnArr *lnArrP, unsigned nLinesMax)
{
    unsigned n;		/* Loop index */
    size_t sz;

    if (lnArrP->nLinesMax == nLinesMax) {
	return;
    }
    for (n = nLinesMax; n < lnArrP->nLinesMax; n++) {
	/*
	 * Free excess lines
	 */

	GeoLnDestroy(lnArrP->lines[n]);
    }
    sz = nLinesMax * sizeof(struct GeoLn *);
    lnArrP->lines = (struct GeoLn **)REALLOC(lnArrP->lines, sz);
    lnArrP->nLinesMax = nLinesMax;
    for (n = lnArrP->nLines; n < lnArrP->nLinesMax; n++) {
	/*
	 * Initialize additional lines.
	 */

	lnArrP->lines[n] = NULL;
    }

}

/*
 *----------------------------------------------------------------------
 *
 * GeoLnArrAddLine --
 *
 *	This procedure copies a geoline onto the end of a geolinearray.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	See the user documentation.
 *
 *----------------------------------------------------------------------
 */

int GeoLnArrAddLine(struct GeoLn *ln, struct GeoLnArr *lnArrP)
{
    unsigned nLines;		/* Initial number of lines in lnArrP */

    nLines = lnArrP->nLines;
    if (nLines + 1 > lnArrP->nLinesMax) {
	/*
	 * If lnArrP needs more space in its allocation, grow it by 25%
	 */

	GeoLnArrSetAlloc(lnArrP, ((lnArrP->nLinesMax + 4) * 5) / 4);
    }
    if ( !(lnArrP->lines[nLines] = GeoLnCreate(ln->nPts)) ) {
	return 0;
    }
    lnArrP->nPts += ln->nPts;
    lnArrP->nMax = (lnArrP->nMax > ln->nPts) ? lnArrP->nMax : ln->nPts;
    lnArrP->latMax
	= (lnArrP->latMax > ln->latMax) ? lnArrP->latMax : ln->latMax;
    lnArrP->lonMax
	= (lnArrP->lonMax > ln->lonMax) ? lnArrP->lonMax : ln->lonMax;
    lnArrP->latMin
	= (lnArrP->latMin < ln->latMin) ? lnArrP->latMin : ln->latMin;
    lnArrP->lonMin
	= (lnArrP->lonMin < ln->lonMin) ? lnArrP->lonMin : ln->lonMin;
    memcpy(lnArrP->lines[nLines]->pts, ln->pts, ln->nPts * sizeof(GeoPt));
    lnArrP->lines[nLines]->nPts = ln->nPts;
    lnArrP->lines[nLines]->lonMax = ln->lonMax;
    lnArrP->lines[nLines]->lonMin = ln->lonMin;
    lnArrP->lines[nLines]->latMax = ln->latMax;
    lnArrP->lines[nLines]->latMin = ln->latMin;
    lnArrP->nLines++;
    return 1;
}

/*
 *----------------------------------------------------------------------
 *
 * GeoLnArrPutLine --
 *
 *	This procedure gives a geoline to a geolinearray.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	See the user documentation.
 *
 *----------------------------------------------------------------------
 */

int GeoLnArrPutLine(struct GeoLn *ln, struct GeoLnArr *lnArrP)
{
    unsigned nLines;		/* Initial number of lines in lnArrP */

    nLines = lnArrP->nLines;
    if (nLines + 1 > lnArrP->nLinesMax) {
	/*
	 * If lnArrP needs more space in its allocation, grow it by 25%
	 */

	GeoLnArrSetAlloc(lnArrP, ((lnArrP->nLinesMax + 4) * 5) / 4);
    }
    lnArrP->nPts += ln->nPts;
    lnArrP->nMax = (lnArrP->nMax > ln->nPts) ? lnArrP->nMax : ln->nPts;
    lnArrP->latMax
	= (lnArrP->latMax > ln->latMax) ? lnArrP->latMax : ln->latMax;
    lnArrP->lonMax
	= (lnArrP->lonMax > ln->lonMax) ? lnArrP->lonMax : ln->lonMax;
    lnArrP->latMin
	= (lnArrP->latMin < ln->latMin) ? lnArrP->latMin : ln->latMin;
    lnArrP->lonMin
	= (lnArrP->lonMin < ln->lonMin) ? lnArrP->lonMin : ln->lonMin;
    lnArrP->lines[nLines] = ln;
    lnArrP->nLines++;
    return 1;
}

/*
 *----------------------------------------------------------------------
 *
 * GeoLnArrGetDescr --
 *
 *	This is a member access function.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

char * GeoLnArrGetDescr(struct GeoLnArr *lnArrP)
{
	return lnArrP->descr;
}

/*
 *----------------------------------------------------------------------
 *
 * GeoLnArrGetLine --
 *
 *	This procedure retrieves a geoline from a geolinearray.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

struct GeoLn *GeoLnArrGetLine(struct GeoLnArr *lnArrP, unsigned n)
{
    return (n < lnArrP->nLines) ? lnArrP->lines[n] : NULL;
}

/*
 *----------------------------------------------------------------------
 *
 * GeoLnArrFree --
 *
 *	This procedure frees internal storage in a geolinearray.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	See the user documentation.
 *
 *----------------------------------------------------------------------
 */

void GeoLnArrFree(struct GeoLnArr *lnArrP)
{
    unsigned n;

    if ( !lnArrP ) {
	return;
    }
    for (n = 0; n < lnArrP->nLines; n++) {
	GeoLnDestroy(lnArrP->lines[n]);
    }
    FREE((char *)lnArrP->lines);
    FREE((char *)lnArrP->descr);
}

/*
 *----------------------------------------------------------------------
 *
 * GeoLnArrDestroy --
 *
 *	This procedure frees internal storage in a geolinearray and frees the
 *	geolinearray.
 *
 * Results:
 *	None.
 *
 * Side effects:
 *	See the user documentation.
 *
 *----------------------------------------------------------------------
 */

void GeoLnArrDestroy(struct GeoLnArr *lnArrP)
{
    GeoLnArrFree(lnArrP);
    FREE((char *)lnArrP);
}

/*
 *----------------------------------------------------------------------
 *
 * GeoLnArrContainGeoPt --
 *
 *	This procedure determines whether any of the geolines in a geolinearray
 *	contain a geopoint.
 *
 * Results:
 *	See the user documentation.
 *
 * Side effects:
 *	None.
 *
 *----------------------------------------------------------------------
 */

int GeoLnArrContainGeoPt(GeoPt geoPt, struct GeoLnArr *lnArrP)
{
    unsigned n;

    for (n = 0; n < lnArrP->nLines; n++) {
	if ( GeoLnContainGeoPt(geoPt, lnArrP->lines[n]) ) {
	    return 1;
	}
    }
    return 0;
}
