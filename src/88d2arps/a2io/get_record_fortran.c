/*
 *	Get_record_fortran.c
 *
 *	This is a Fortran interface version.
 *
 *	Do *NOT* include the "len" argument in any Fortran calls!
 *	It is supplied automatically.
 */
 
#pragma ident	"@(#)get_record_fortran.c	5.1	03/01/05	CAPS"

#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <strings.h>

get_record_f( mybuf, irc, len )
char *mybuf;
int *irc;
int len;
{
	*irc = get_record( mybuf );
	return;
}

get_record_f_( mybuf, irc, len )
char *mybuf;
int *irc;
int len;
{
	get_record_f( mybuf, irc, len );
	return;
}

set_record_format_f( n, irc )
int *n;
int *irc;
{
	*irc = set_record_format( *n );
	return;
}

set_record_format_f_( n, irc )
int *n;
int *irc;
{
	set_record_format_f( n, irc );
	return;
}

get_record_format_f( irc )
int *irc;
{
	*irc = get_record_format();
	return;
}

get_record_format_f_( irc )
int *irc;
{
	get_record_format_f( irc );
	return;
}
