#ifndef	_CONFIG_H
#define	_CONFIG_H

/*
#define	DEBUG
 */

/*
 *	Compile flags for the A2IO software.  We put these here so to make
 *	life a little easier.
 */

#pragma ident "@(#)config.h	5.3	03/10/04 CAPS"

/*
 *	BZIP2_OLD
 *		Support BZIP2 software prior to version 1.0.  Use
 *		"bzip2 --help" to find the version number.  If the version
 *		number is less than 1.0, define this, otherwise, don't!
 */

/* #define	BZIP2_OLD */

/*
 *	ALLOW_BZIP2
 *		Add support for processing BZIP2 input files.  This definition
 *		is required to read LDM 88D data files, which is in a modified
 *		BZIP2 format.
 *
 *		Libraries required:  -lbz2
 */

#define	ALLOW_BZIP2

/*
 *	ALLOW_COMPRESS
 *		Add support to read compressed (.Z) input files.
 *
 *		Libraries required:  None
 */

#define	ALLOW_COMPRESS

/*
 *	ALLOW_GZIP
 *		Add support to read gzipped input files.
 *
 *		Libraries required:  -lz
 */

/* #define	ALLOW_GZIP */

/*
 *	ALLOW_TDWR
 *		Add support to read TDWR format files and tapes.
 *		There is an old bug that may or may not exist that will cause
 *		the software to think that TDWR format data is 88D format data.
 *		When this happens, user software may get confused.
 *
 *		The TDWR code hasn't been tested recently, so it may not work!
 */

#define	ALLOW_TDWR

/*
 *     ALLOW_PKU
 *              Small modification to the header to read chinese netrad (CINRAD)
 *              radar data. Provided by Cunxi Zhang from Peking University.
 *              (07/22/2007).
 */

/* #define      ALLOW_PKU  */

/*
 *	NO_ATOLL
 *		Some systems do not support the atoll() call.  This disables
 *		the code, so the problem doesn't occur.  Individual data files
 *		can still be decoded, but directories of data files won't.
 */

/* #define	NO_ATOLL */

#endif	/* _CONFIG_H */
