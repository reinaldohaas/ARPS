/*
 * geoLines.h --
 *
 *	This header file declares structures and functions that
 *	manage GeoLn's (arrays of GeoPts) and GeoLnArr's (arrays of
 *	GeoLn's).  See the GeoLn(3) and GeoLnArr(3) man pages for details.
 * 
 * Copyright (c) 2007 Gordon D. Carrie.  All rights reserved.
 * 
 * Licensed under the Open Software License version 2.1
 * 
 * Please address questions and feedback to user0@tkgeomap.org
 *
 * @(#) $Id: geoLines.h,v 4.0 2013/02/08 18:31:45 ywang Exp $
 *
 ********************************************
 *
 */

#ifndef _GEOLINES_H_
#define _GEOLINES_H_

#include "geography.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Structure to store an array of geographic points
 */

struct GeoLn {
    unsigned	nPts;		/* Number of points */
    unsigned	nPtsMax;	/* Max number of points in current allocation */
    Angle	latMax;		/* Max longitude */
    Angle	lonMax;		/* Max latitude */
    Angle	latMin;		/* Min longitude */
    Angle	lonMin;		/* Min latitude */
    GeoPt	*pts;		/* Array of points in the line */
};

/*
 * Global functions that manage GeoLn structures
 */

struct GeoLn 	*GeoLnCreate (unsigned nptsMax);
void		GeoLnClear (struct GeoLn *geoLn);
void		GeoLnSetAlloc (struct GeoLn *geoLn, unsigned nptsMax);
void		GeoLnAddPt (GeoPt geoPt, struct GeoLn *geoLn);
GeoPt		GeoLnGetPt (struct GeoLn *geoLn, unsigned n);
void		GeoLnDestroy (struct GeoLn *geoLn);
CartPt		GeoLnCtr (struct GeoLn *ln);
int		GeoLnContainGeoPt (GeoPt geoPt, struct GeoLn *geoLn);

/*
 * Structure to store an array of GeoLn's
 */

struct GeoLnArr {
    char	*descr;		/* Descriptor */
    unsigned	nLines;		/* Number of lines */
    unsigned	nLinesMax;	/* Max number of lines in current allocation */
    unsigned	nPts;		/* Number of points for all lines */
    unsigned	nMax;		/* Number of points in longest line */
    Angle	latMax;		/* Max longitude for all lines */
    Angle	lonMax;		/* Max latitude for all lines */
    Angle	latMin;		/* Min longitude for all lines */
    Angle	lonMin;		/* Min latitude for all lines */
    struct GeoLn **lines;	/* Array of lines */
};

/*
 * Global functions to manage GeoLnArr's
 */

struct GeoLnArr *GeoLnArrCreate (unsigned nLinesMax);
void		GeoLnArrSetDescr (struct GeoLnArr *lnArr, char *descr);
void		GeoLnArrSetAlloc (struct GeoLnArr *lnArr, unsigned nLinesMax);
int		GeoLnArrAddLine (struct GeoLn *geoLn, struct GeoLnArr *lnArr);
int		GeoLnArrPutLine (struct GeoLn *geoLn, struct GeoLnArr *lnArr);
int		GeoLnArrContainGeoPt (GeoPt geoPt, struct GeoLnArr *lnArr);
char		*GeoLnArrGetDescr (struct GeoLnArr *lnArr);
struct GeoLn	*GeoLnArrGetLine (struct GeoLnArr *lnArr, unsigned n);
void		GeoLnArrFree (struct GeoLnArr *lnArr);
void		GeoLnArrDestroy (struct GeoLnArr *lnArr);
    
#ifdef __cplusplus
}
#endif

#endif
