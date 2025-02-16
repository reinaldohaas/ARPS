/*
 *	Tdwr_struct.h
 *
 *	TDWR structures...
 */

#pragma ident	"@(#)tdwr_struct.h	5.1	10/10/01	CAPS"

#ifndef	_TDWR_STRUCT_H
#define	_TDWR_STRUCT_H

#include "config.h"

#ifdef	ALLOW_TDWR

struct tdwr_msg_hdr *tdwr_pmsg;
struct tdwr_basedata_hdr *tdwr_phdr;
struct tdwr_prod_2b00 *tdwr_p2b00;
struct tdwr_prod_2b01 *tdwr_p2b01;
struct tdwr_prod_2b02 *tdwr_p2b02;

#endif	/* ALLOW_TDWR */

#endif	/* _TDWR_STRUCT_H */
