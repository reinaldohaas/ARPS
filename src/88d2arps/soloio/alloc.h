/*
 * alloc.c --
 *
 * 	This file declares memory allocators.
 * 
 * Copyright (c) 2007, 2008 Gordon D. Carrie
 *
 * Licensed under the Open Software License version 2.1
 *
 * Please send feedback to user0@tkgeomap.org
 *
 * @(#) $Id: alloc.h,v 4.0 2013/02/08 18:31:44 ywang Exp $
 *
 **********************************************************************
 *
 */

#ifndef ALLOC_H_
#define ALLOC_H_

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef MEM_DEBUG
#define MALLOC(x) DB_Malloc(x, __FILE__, __LINE__)
#define REALLOC(x,s) DB_Realloc((x),(s), __FILE__, __LINE__)
#define FREE(x) DB_Free((x), __FILE__, __LINE__)
#else
#define MALLOC(x) malloc(x)
#define REALLOC(x,s) realloc((x),(s))
#define FREE(x) free((x))
#endif

void *Malloc(size_t sz);
void *Realloc(void *m, size_t sz);
void *DB_Malloc(size_t sz, char *fl_nm, int ln);
void *DB_Realloc(void *m, size_t sz, char *fl_nm, int ln);
void DB_Free(void *m, char *fl_nm, int ln);
float ***Alloc_Arr3(size_t i, size_t j, size_t k);
void Alloc_Free3(float ***d);
float ****Alloc_Arr4(size_t i, size_t j, size_t k, size_t l);
void Alloc_Free4(float ****d);
    
#ifdef __cplusplus
}
#endif

#endif
