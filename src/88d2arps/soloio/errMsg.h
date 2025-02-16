/*
 * errMsg --
 *
 * 	This file declares general purpose functions and macros to
 * 	generate globally visible error messages.
 * 
 * Copyright (c) 2007 Gordon D. Carrie
 *
 * Licensed under the Open Software License version 2.1
 *
 * Please send feedback to user0@tkgeomap.org
 *
 * @(#) $Id: errMsg.h,v 4.0 2013/02/08 18:31:45 ywang Exp $
 *
 **********************************************************************
 *
 */

#ifndef ERRMSG_H_
#define ERRMSG_H_

#ifdef __cplusplus
extern "C" {
#endif

/*
 * The following functions access the global error message.
 */

void ErrMsg_Append(const char *msg);
char *ErrMsg_Get(void);
void ErrMsg_Destroy(void);
    
#ifdef __cplusplus
}
#endif

#endif
