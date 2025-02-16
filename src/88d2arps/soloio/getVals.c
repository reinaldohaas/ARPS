/*
 * getVals.c ---
 *
 *	This file defines functions that extracting values of various types
 *	stored in a char array.
 *
 * Copyright (c) 2007 Gordon D. Carrie
 *
 * Licensed under the Open Software License version 2.1
 *
 * Please send feedback to user0@tkgeomap.org
 *
 * @(#) $Id: getVals.c,v 4.0 2013/02/08 18:31:45 ywang Exp $
 *
 **********************************************************************
 *
 */

#include <stdlib.h>
#include <string.h>
#include "sizes.h"
#include "getVals.h"

/*
 * Flag - if true, we are byte swapping
 */

static int swapping;

/*
 * These macros swap bytes.  They assume 8 bits per byte.
 */

#define B0 0x000000FF
#define B1 0x0000FF00
#define B2 0x00FF0000
#define B3 0xFF000000
#define SWAP2(s) \
	((((uint16)s & B0) << 8) | (((uint16)s & B1) >> 8))
#define SWAP4(i) \
	((((uint32)i & B0) << 24) | (((uint32)i & B1) << 8) \
	 | (((uint32)i & B2) >> 8) | (((uint32)i & B3) >> 24))

/*
 *------------------------------------------------------------------------
 *
 * GetVals_SwapOff --
 *
 * 	This function forces i.o. to swap.
 *
 * Results:
 * 	None.
 *
 * Side effects:
 * 	The local swapping flag is set to true, enabling swapping.
 *
 *------------------------------------------------------------------------
 */

void GetVals_SwapOff(void)
{
    swapping = 0;
}

/*
 *------------------------------------------------------------------------
 *
 * GetVals_SwapOn --
 *
 * 	This function forces i.o. to swap.
 *
 * Results:
 * 	None.
 *
 * Side effects:
 * 	The local swapping flag is set to true, enabling swapping.
 *
 *------------------------------------------------------------------------
 */

void GetVals_SwapOn()
{
    swapping = 1;
}

/*
 *------------------------------------------------------------------------
 *
 * GetVals_ToggleSwap --
 *
 * 	This function toggles the swap flag.
 *
 * Results:
 * 	None.
 *
 * Side effects:
 * 	The local swapping flag is switched.  If swapping is enabled, this
 * 	function disables it, and vice versa.
 *
 *------------------------------------------------------------------------
 */

void GetVals_ToggleSwap()
{
    swapping = !swapping;
}

/*
 *------------------------------------------------------------------------
 *
 * GetVals_GetSINT16 --
 *
 * 	This function retrieves a signed, 16 bit integer from an address.
 *
 * Results:
 * 	Return value is a signed integer corresponding to the value
 * 	stored at the given address.
 *
 * Side effects:
 * 	None.
 *
 *------------------------------------------------------------------------
 */

int GetVals_GetSINT16(char *b)
{
    int16 s = *(int16 *)(b);
    return swapping ? SWAP2(s) : s;
}

/*
 *------------------------------------------------------------------------
 *
 * GetVals_GetUINT16 --
 *
 * 	This function retrieves an unsigned, 16 bit integer from an address.
 *
 * Arguments:
 *	char *b	- address of the 16 bit integer.
 *
 * Results:
 * 	Return value is an unsigned integer corresponding to the value
 * 	stored at the given address.
 *
 * Side effects:
 * 	None.
 *
 *------------------------------------------------------------------------
 */

unsigned GetVals_GetUINT16(char *b)
{
    uint16 s = *(uint16 *)(b);
    return swapping ? SWAP2(s) : s;
}

/*
 *------------------------------------------------------------------------
 *
 * GetVals_GetSINT32 --
 *
 * 	This function retrieves a signed, 32 bit integer from an address.
 *
 * Arguments:
 *	char *b	- address of the 32 bit integer.
 *
 * Results:
 * 	Return value is a signed integer corresponding to the value
 * 	stored at the given address.
 *
 * Side effects:
 * 	None.
 *
 *------------------------------------------------------------------------
 */

int GetVals_GetSINT32(char *b)
{
    int32 s = *(int32 *)(b);
    return swapping ? SWAP4(s) : s;
}

/*
 *------------------------------------------------------------------------
 *
 * GetVals_GetUINT32 --
 *
 * 	This function retrieves an unsigned, 32 bit integer from an address.
 *
 * Arguments:
 *	char *b	- address of the 32 bit integer.
 *
 * Results:
 * 	Return value is an unsigned integer corresponding to the value
 * 	stored at the given address.
 *
 * Side effects:
 * 	None.
 *
 *------------------------------------------------------------------------
 */

unsigned GetVals_GetUINT32(char *b)
{
    uint32 s = *(uint32 *)(b);
    return swapping ? SWAP4(s) : s;
}

/*
 *------------------------------------------------------------------------
 *
 * GetVals_GetFLOAT32 --
 *
 * 	This function retrieves a 32 bit floating point value from an address.
 *
 * Arguments:
 *	char *b	- address of the 32 bit float value.
 *
 * Results:
 * 	Return value is a float corresponding to the value
 * 	stored at the given address.
 *
 * Side effects:
 * 	None.
 *
 *------------------------------------------------------------------------
 */

double GetVals_GetFLOAT32(char *b)
{
    float f = *(float *)(b);
    return swapping ? SWAP4(f) : f;
}

/*
 *------------------------------------------------------------------------
 *
 * GetVals_SwapArr16 --
 *
 * 	This function applies swapping to an array of 16 byte integers.
 *
 * Arguments:
 *	unsigned short *r	- address of the array
 *	size_t nw		- number of integers to swap
 *
 * Results:
 * 	None.
 *
 * Side effects:
 * 	The input array is rearranged.
 *
 *------------------------------------------------------------------------
 */

void GetVals_SwapArr16(unsigned short *r, size_t nw)
{
    uint16 *p, *pe;
    if ( !swapping )
      return;
    else
      for (p = r, pe = r + nw; p < pe; p++)
        *p = ((*p & B0) << 8) | ((*p & B1) >> 8);
}

/*
 *------------------------------------------------------------------------
 *
 * GetVals_PutSINT16 --
 *
 * 	This function puts an integer at the address of a signed, 16 bit int.
 *
 * Arguments:
 *	int i	- integer to put at the address
 *	char *b	- the address
 *
 * Results:
 * 	A value is copied to the given address.
 *
 * Side effects:
 * 	None.
 *
 *------------------------------------------------------------------------
 */

void GetVals_PutSINT16(int i, char *b)
{
    int16 i2 = i;
    *(int16 *)b = swapping ? SWAP2(i2) : i2;
}

/*
 *------------------------------------------------------------------------
 *
 * GetVals_PutUINT16 --
 *
 * 	This function puts an integer at the address of an unsigned, 16 bit int.
 *
 * Arguments:
 *	unsigned u	- integer to put at the address
 *	char *b		- the address
 *
 * Results:
 * 	A value is copied to the given address.
 *
 * Side effects:
 * 	None.
 *
 *------------------------------------------------------------------------
 */

void GetVals_PutUINT16(unsigned u, char *b)
{
    uint16 u2 = u;
    *(uint16 *)b = swapping ? SWAP2(u2) : u2;
}

/*
 *------------------------------------------------------------------------
 *
 * GetVals_PutSINT32 --
 *
 * 	This function puts an integer at the address of a signed, 32 bit int.
 *
 * Arguments:
 *	int i	- integer to put at the address
 *	char *b	- the address
 *
 * Results:
 * 	A value is copied to the given address.
 *
 * Side effects:
 * 	None.
 *
 *------------------------------------------------------------------------
 */

void GetVals_PutSINT32(int i, char *b)
{
    int32 i4 = i;
    *(int32 *)b = swapping ? SWAP2(i4) : i4;
}

/*
 *------------------------------------------------------------------------
 *
 * GetVals_PutUINT32 --
 *
 * 	This function puts an integer at the address of an unsigned, 32 bit int.
 *
 * Arguments:
 *	unsigned u	- integer to put at the address
 *	char *b		- the address
 *
 * Results:
 * 	A value is copied to the given address.
 *
 * Side effects:
 * 	None.
 *
 *------------------------------------------------------------------------
 */

void GetVals_PutUINT32(unsigned u, char *b)
{
    uint32 u4 = u;
    *(uint32 *)b = swapping ? SWAP2(u4) : u4;
}

/*
 *------------------------------------------------------------------------
 *
 * GetVals_PutFLOAT32 --
 *
 * 	This function puts a floating point value at the address of a
 * 	32 bit float.
 *
 * Arguments:
 *	float f		- float value to put at the address
 *	char *b		- the address
 *
 * Results:
 * 	A value is copied to the given address.
 *
 * Side effects:
 * 	None.
 *
 *------------------------------------------------------------------------
 */

void GetVals_PutFLOAT32(float f, char *b)
{
    *(float *)b = swapping ? SWAP4(f) : f;
}
