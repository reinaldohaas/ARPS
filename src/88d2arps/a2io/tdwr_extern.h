/*
 *	Tdwr_extern.h
 *
 *	TDWR externals...
 */

#pragma ident	"@(#)tdwr_extern.h	5.1	10/10/01	CAPS"

#ifndef	_TDWR_EXTERN_H
#define	_TDWR_EXTERN_H

#include "config.h"

#ifdef	ALLOW_TDWR

extern struct tdwr_msg_hdr *tdwr_pmsg;
extern struct tdwr_basedata_hdr *tdwr_phdr;
extern struct tdwr_prod_2b00 *tdwr_p2b00;
extern struct tdwr_prod_2b01 *tdwr_p2b01;
extern struct tdwr_prod_2b02 *tdwr_p2b02;

#endif	/* ALLOW_TDWR */

#endif	/* _TDWR_EXTERN_H */
