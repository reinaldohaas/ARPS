/*
 * getVals.h --
 *
 * 	This file declares general purpose functions and macros to use with
 * 	read and write data.
 * 
 * Copyright (c) 2007 Gordon D. Carrie
 *
 * Licensed under the Open Software License version 2.1
 *
 * Please send feedback to user0@tkgeomap.org
 *
 * @(#) $Id: getVals.h,v 4.0 2013/02/08 18:31:45 ywang Exp $
 *
 **********************************************************************
 *
 */

#ifndef GETVALS_H_
#define GETVALS_H_

#ifdef __cplusplus
extern "C" {
#endif

/*
 * The following functions extract values of various types stored in a
 * char buffer.
 */

void GetVals_SwapOn(void);
void GetVals_SwapOff(void);
void GetVals_ToggleSwap(void);
int GetVals_GetSINT16(char *b);
unsigned GetVals_GetUINT16(char *b);
int GetVals_GetSINT32(char *b);
unsigned GetVals_GetUINT32(char *b);
double GetVals_GetFLOAT32(char *b);
void GetVals_SwapArr16(unsigned short *r, size_t nw);

/*
 * The following functions put values of various types into a char buffer.
 */

void GetVals_PutSINT16(int i, char *b);
void GetVals_PutUINT16(unsigned u, char *b);
void GetVals_PutSINT32(int i, char *b);
void GetVals_PutUINT32(unsigned u, char *b);
void GetVals_PutFLOAT32(float f, char *b);
    
#ifdef __cplusplus
}
#endif

#endif
