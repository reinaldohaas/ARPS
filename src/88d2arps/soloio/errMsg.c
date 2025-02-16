/*
 * errMsg --
 *
 * 	This file defines general purpose functions and macros to
 * 	generate globally visible error messages.
 * 
 * Copyright (c) 2007 Gordon D. Carrie
 *
 * Licensed under the Open Software License version 2.1
 *
 * Please send feedback to user0@tkgeomap.org
 *
 * @(#) $Id: errMsg.c,v 4.0 2013/02/08 18:31:45 ywang Exp $
 *
 **********************************************************************
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "alloc.h"
#include "errMsg.h"

static char *errMsg;	/* Current error message */
static size_t alloc;	/* Allocation at errMsg */
static size_t len;	/* strlen(errMsg) */

/*
 *------------------------------------------------------------------------
 *
 * ErrMsg_Append --
 *
 * 	This function appends a string to the current error message.
 *
 * Arguments:
 *	const char *msg	- string to append to the error message.
 *
 * Results:
 * 	None.
 *
 * Side effects:
 * 	A string is appended to the string at errMsg, if any.  Memory is
 * 	reallocated if necessary.  If allocation fails, the application exits.
 *
 *------------------------------------------------------------------------
 */

void ErrMsg_Append(const char *msg)
{
    size_t new_len, new_alloc;
    char *e, *e1;
    const char *m;

    if ( !msg || strlen(msg) == 0 ) {
	return;
    }
    new_len = len + strlen(msg);
    new_alloc = new_len + 1;
    if (new_alloc > alloc) {
	errMsg = REALLOC(errMsg, new_alloc);
	alloc = new_alloc;
    }
    for (e = errMsg + len, m = msg, e1 = e + strlen(msg); e < e1; e++, m++) {
	*e  = *m;
    }
    *e = '\0';
    len = new_len;
}

/*
 *------------------------------------------------------------------------
 *
 * ErrMsg_Get --
 *
 * 	This function returns the current error message.
 *
 *
 * Results:
 *	Return value is the address of the error message.  User should not
 *	modify it.

 * Side effects:
 *	Error message is set to zero length, although it remains allocated.
 *
 *------------------------------------------------------------------------
 */

char *ErrMsg_Get(void)
{
    if (errMsg) {
	len = 0;
	return errMsg;
    } else {
	return "";
    }
}

/*
 *------------------------------------------------------------------------
 *
 * ErrMsg_Destroy --
 *
 * 	This clean up function frees memory allocated in this file.
 *	It should be used when the application exits.
 *
 * Side effects:
 *	The error message is freed.
 *
 *------------------------------------------------------------------------
 */

void ErrMsg_Destroy(void)
{
    if (errMsg) {
	FREE(errMsg);
    }
    errMsg = NULL;
    len = 0;
}
