/*
 * geoTime.c --
 *
 *	This file defines functions for manipulate date and times
 * 
 * Copyright (c) 2007 Gordon D. Carrie
 *
 * Licensed under the Open Software License version 2.1
 *
 * Please send feedback to user0@tkgeomap.org
 *
 * @(#) $Id: geoTime.c,v 4.0 2013/02/08 18:31:45 ywang Exp $
 *
 **********************************************************************
 */

#include <stdlib.h>
#include "geography.h"

/*
 *------------------------------------------------------------------------
 *
 * GeoTime_CalSet --
 *
 * 	This function assigns values to a GeoTime_Cal structure.
 *
 * Results:
 * 	Assigns values to a calendar time structure.
 *
 * Side effects:
 * 	None.
 *------------------------------------------------------------------------
 */

struct GeoTime_Cal
GeoTime_CalSet(int year, int month, int day, int hour, int minute, double second)
{
    struct GeoTime_Cal cal;
    int s;
    double fsec;

    s = floor(second);
    fsec = second - s;
    minute += s / 60;
    s %= 60;
    hour += minute / 60;
    minute %= 60;
    day += hour / 24;
    hour %= 24;

    cal.year = year;
    cal.month = month;
    cal.day = day;
    cal.hour = hour;
    cal.minute = minute;
    cal.second = s + fsec;
    return cal;
}

/*
 *------------------------------------------------------------------------
 *
 * GeoTime_JulSet --
 *
 * 	This function assigns values to a GeoTime_Jul structure.
 *
 * Results:
 * 	Assigns values to a Julian time structure.
 *
 * Side effects:
 * 	None.
 *
 *------------------------------------------------------------------------
 */

struct GeoTime_Jul GeoTime_JulSet(int day, double second)
{
    struct GeoTime_Jul jul;
    int s;
    double fsec;

    s = floor(second);
    fsec = second - s;
    day += s / 86400;
    s %= 86400;
    jul.day = day;
    jul.second = s + fsec;
    return jul;
}

/*
 *------------------------------------------------------------------------
 *
 * GeoTime_CalToJul --
 *
 *	This function converts a calendar time representation to Julian date.
 *
 * Arguments:
 *	cal	- callendar time to convert
 *
 * Results:
 * 	Return value gives the Julian representation of a calendar time.
 *
 * Side effects:
 *	None.
 *------------------------------------------------------------------------
 */

struct GeoTime_Jul GeoTime_CalToJul(struct GeoTime_Cal cal)
{
    int year = cal.year;
    int month = cal.month;
    int day = cal.day;
    int hour = cal.hour;
    int minute = cal.minute;
    double second = cal.second;
    struct GeoTime_Jul jul;	/* Return value */
    int d;			/* Number of days in second */

    /*
     * Compute Julian day.
     * Ref. Henry F. Fliegel and Thomas C. Van Flandern, from
     * http://serendipity.magnet.ch/hermetic/cal_stud/jdn.htm or
     * http://aa.usno.navy.mil/faq/docs/JD_Formula.html
     */

    jul.day =   (1461 * (year + 4800 + (month - 14) / 12)) / 4
	+ (367 * (month - 2 - 12 * ((month - 14) / 12))) / 12
	- (3 * ((year + 4900 + (month - 14) / 12) / 100)) / 4
	+ day - 32075;

    /*
     * Compute total seconds.  If greater than 1 day, increment days.
     */

    second = 3600.0 * hour + 60.0 * minute + second;
    d = floor(second / 86400.0);
    jul.day += d;
    jul.second = second - d * 86400.0;
    return jul;
}

/*
 *------------------------------------------------------------------------
 *
 * GeoTime_GetCal --
 *
 *	This function retrieves a calendar date and time from a GeoTime.
 *
 * Arguments:
 *	struct GeoTime_Jul jul - Time instant in calendar form
 *
 * Results:
 * 	Return value gives the calendar representation of a Julian time.
 *
 * Side effects:
 * 	None.
 *
 *------------------------------------------------------------------------
 */

struct GeoTime_Cal GeoTime_JulToCal(struct GeoTime_Jul jul)
{
    struct GeoTime_Cal cal;	/* Return value, calendar equivalent of jul */
    int l, n, i, j;		/* Intermediaries */
    int h, m;			/* Hour, minute */
    int s;			/* Second */
    double fsec;

    /*
     * Compute month, day, and year from Julian date
     * Ref. Henry F. Fliegel and Thomas C. Van Flandern,
     * http://serendipity.magnet.ch/hermetic/cal_stud/jdn.htm
     */

    l = jul.day + 68569;
    n = (4 * l) / 146097;
    l = l - (146097 * n + 3) / 4;
    i = (4000 * (l + 1)) / 1461001;
    l = l - (1461 * i) / 4 + 31;
    j = (80 * l) / 2447;
    cal.day = l - (2447 * j) / 80;
    l = j / 11;
    cal.month = j + 2 - (12 * l);
    cal.year = 100 * (n - 49) + i + l;

    /*
     * Convert seconds to clock time
     */

    s = floor(jul.second);
    fsec = jul.second - s;
    h = s / 3600;
    s %= 3600;
    m = s / 60;
    s %= 60;
    if (h > 23) {
	cal.day += h / 24;
	h %= 24;
    }
    cal.hour = h;
    cal.minute = m;
    cal.second = s + fsec;
    return cal;
}

/*
 *------------------------------------------------------------------------
 *
 * GeoTime_Incr --
 *
 *	This function increments a time instant.
 *
 * Arguments:
 *	struct GeoTime_Jul *jul	- GeoTime to increment
 *	double ds;		- Increment, in seconds
 *
 * Results:
 * 	This given time is incremented.
 *
 * Side effects:
 * 	None.
 *
 *------------------------------------------------------------------------
 */

void GeoTime_Incr(struct GeoTime_Jul *jul, double ds)
{
    int d;
    double s;

    s = jul->second + ds;
    d = floor(s / 86400.0);
    jul->day += d;
    jul->second = s - d * 86400.0;
}

/*
 *------------------------------------------------------------------------
 *
 * GeoTime_Cmp --
 *
 *	Compare two times.
 *
 * Arguments:
 *	struct GeoTime_Jul jul1;	GeoTimes to compare
 *	struct GeoTime_Jul jul2;
 *
 * Results:
 * 	Return value is a collating code.
 *	    -1 if jul1 is before jul2
 *	     0 if jul1 is same time as jul2
 *	     1 if jul1 is after jul2
 *
 * Side effects:
 * 	None.
 *
 *------------------------------------------------------------------------
 */

int GeoTime_Cmp(struct GeoTime_Jul jul1, struct GeoTime_Jul jul2)
{

    jul1 = GeoTime_JulSet(jul1.day, jul1.second);
    jul2 = GeoTime_JulSet(jul2.day, jul2.second);
    if (jul1.day > jul2.day) {
	return 1;
    }
    if (jul1.day < jul2.day) {
	return -1;
    }
    if (jul1.second > jul2.second) {
	return 1;
    }
    if (jul1.second < jul2.second) {
	return -1;
    }
    return 0;
}

/*
 *------------------------------------------------------------------------
 *
 * GeoTime_Diff --
 *
 *	Subtract one time from another.
 *
 * Arguments:
 *	struct GeoTime_Jul jul1
 *	struct GeoTime_JUL jul2		- GeoTimes to compare
 *
 * Results:
 * 	Return value is the difference in seconds between two times.
 *
 * Side effects:
 * 	None.
 *
 *------------------------------------------------------------------------
 */

double
GeoTime_Diff(struct GeoTime_Jul jul1, struct GeoTime_Jul jul2)
{

    jul1 = GeoTime_JulSet(jul1.day, jul1.second);
    jul2 = GeoTime_JulSet(jul2.day, jul2.second);
    return (jul1.day - jul2.day) * 86400.0 + jul1.second - jul2.second;
}
